#include "SimpleMRFSegmenter.h"

#include <itkExpectationMaximizationMixtureModelEstimator.h>
#include <itkGaussianMixtureModelComponent.h>
#include <itkImageToListSampleFilter.h>

#include <mrf.h>
#include <MaxProdBP.h>

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

namespace vt
{
// ------------------------------------------------------------------------
SimpleMRFSegmenter::SimpleMRFSegmenter()
{
	m_InitialPositiveMean = 45;
	m_InitialPositiveVariance = 1000;
	m_InitialNegativeMean = 245;
	m_InitialNegativeVariance = 1000;
	m_SmoothnessCost = 1.0;
	m_SmallestCheck = 0.00000001;
	m_OutputValue  = 1;
}

// ------------------------------------------------------------------------
void SimpleMRFSegmenter::Segment()
{
	// compute the gaussians for the classifiers
	CreateMembershipFunction();


	// create the data costs
	unsigned int numberOfPixels = m_Image->GetLargestPossibleRegion().GetNumberOfPixels();
	unsigned int width = m_Image->GetLargestPossibleRegion().GetSize()[0];
	unsigned int height = m_Image->GetLargestPossibleRegion().GetSize()[1];

	MRF::EnergyVal * dCosts = new MRF::EnergyVal[numberOfPixels*2];

	// iterate through the image
	typedef itk::ImageRegionConstIterator<ImageType> ImageIteratorType;
	typedef itk::ImageRegionConstIterator<MaskType> MaskIteratorType;
	ImageIteratorType imageIt(m_Image, m_Image->GetLargestPossibleRegion());

	unsigned int pcount = 0;

	while(!imageIt.IsAtEnd())
	{
		unsigned short pixelValue = imageIt.Get();
		MeanType measurement;
		measurement[0] = pixelValue;

		double positiveProbability = m_PositiveFunction->Evaluate(measurement);
		double negativeProbability = m_NegativeFunction->Evaluate(measurement);

		dCosts[pcount*2] = -log(std::max(m_SmallestCheck, positiveProbability));
		dCosts[pcount*2+1] = -log(std::max(m_SmallestCheck, negativeProbability));

		++imageIt; ++pcount;
	}

	DataCost * data = new DataCost(dCosts);
	SmoothnessCost * smoothness = new SmoothnessCost(2, 2000, m_SmoothnessCost);
	EnergyFunction *eng = new EnergyFunction(data, smoothness);
	
	MRF * mrf = new MaxProdBP(width, height, 2, eng);
	mrf->initialize();
	mrf->clearAnswer();
	float t;
	mrf->optimize(5,t);


	//std::cout << "Data Cost: " << mrf->dataEnergy() << std::endl;
	//std::cout << "Smoothness Cost: " << mrf->smoothnessEnergy() << std::endl;


	// now we create the output
	m_Output = MaskType::New();
	m_Output->SetSpacing(m_Mask->GetSpacing());
	m_Output->SetDirection(m_Mask->GetDirection());
	m_Output->SetOrigin(m_Mask->GetOrigin());
	m_Output->SetRegions(m_Mask->GetLargestPossibleRegion());
	m_Output->Allocate();
	m_Output->FillBuffer(0);

	MaskIteratorType maskIt(m_Mask, m_Image->GetLargestPossibleRegion());
	typedef itk::ImageRegionIterator<MaskType> OutputIteratorType;
	OutputIteratorType outputIt(m_Output, m_Output->GetLargestPossibleRegion());

	pcount = 0;
	while(!outputIt.IsAtEnd())
	{
		if(mrf->getLabel(pcount) == 0)
			outputIt.Set(m_OutputValue);
		else
			outputIt.Set(0);
		++outputIt; ++pcount;
	}

	delete mrf;

}


// ------------------------------------------------------------------------
void SimpleMRFSegmenter::CreateMembershipFunction()
{
	// get the image samples
	typedef itk::Statistics::ImageToListSampleFilter<ImageType, MaskType> ListSampleFilterType;
	ListSampleFilterType::Pointer listSampleFilter = ListSampleFilterType::New();
	listSampleFilter->SetInput(m_Image);
	listSampleFilter->SetMaskImage(m_Mask);
	listSampleFilter->SetMaskValue(1);

	try 
	{
		listSampleFilter->Update();
	}
	catch(itk::ExceptionObject &e)
	{
		std::cout << "Error in list sample filter: " << e << std::endl;
		exit(1);
	}


	// create the parameters
	typedef itk::Array<double> ParametersType;
	ParametersType params(2);

	std::vector<ParametersType> initialParameters(2);
	params[0] = m_InitialPositiveMean;
	params[1] = m_InitialPositiveVariance;
	initialParameters[0] = params;

	params[0] = m_InitialNegativeMean;
	params[1] = m_InitialNegativeVariance;
	initialParameters[1] = params;


	// create the components
	typedef itk::Statistics::GaussianMixtureModelComponent<SampleType> ComponentType;
	std::vector<ComponentType::Pointer> components;
	components.push_back(ComponentType::New());
	components[0]->SetSample(listSampleFilter->GetOutput());
	components[0]->SetParameters(initialParameters[0]);

	components.push_back(ComponentType::New());
	components[1]->SetSample(listSampleFilter->GetOutput());
	components[1]->SetParameters(initialParameters[1]);

	//std::cout << listSampleFilter->GetOutput()->GetTotalFrequency() << std::endl;


	// create the estimator
	typedef itk::Statistics::ExpectationMaximizationMixtureModelEstimator<SampleType> EstimatorType;
	EstimatorType::Pointer estimator = EstimatorType::New();

	estimator->SetSample(listSampleFilter->GetOutput());
	estimator->SetMaximumIteration(200);

	itk::Array< double > initialProportions(2);
	initialProportions[0] = 0.5;
	initialProportions[1] = 0.5;
	
	estimator->SetInitialProportions(initialProportions);
	estimator->AddComponent((ComponentType::Superclass*) (components[0]).GetPointer());
	estimator->AddComponent((ComponentType::Superclass*) (components[1]).GetPointer());

	try 
	{
		estimator->Update();
	}
	catch(itk::ExceptionObject &e)
	{
		std::cout << "Error in EM algorithm: " << e << std::endl;
		exit(1);
	}

	//std::cout << "Positive: " << components[0]->GetFullParameters() << std::endl;
	//std::cout << "Negative: " << components[1]->GetFullParameters() << std::endl;

	// build the membership functions
	MeanType positiveMean, negativeMean;
	positiveMean[0] = components[0]->GetFullParameters()[0];
	negativeMean[0] = components[1]->GetFullParameters()[0];

	VarianceType positiveVariance(1,1), negativeVariance(1,1);
	positiveVariance(0,0) = components[0]->GetFullParameters()[1];
	negativeVariance(0,0) = components[1]->GetFullParameters()[1];

	m_PositiveFunction = MembershipFunctionType::New();
	m_PositiveFunction->SetMean(positiveMean);
	m_PositiveFunction->SetCovariance(positiveVariance);

	m_NegativeFunction = MembershipFunctionType::New();
	m_NegativeFunction->SetMean(negativeMean);
	m_NegativeFunction->SetCovariance(negativeVariance);
}


}
