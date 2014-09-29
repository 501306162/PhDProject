#include <iostream>
#include <FilenamesReader.h>
#include <CommonDefinitions.h>


#include <itkListSample.h>
#include <itkVector.h>
#include <itkImageToListSampleFilter.h>
#include <itkGaussianMixtureModelComponent.h>
#include <itkExpectationMaximizationMixtureModelEstimator.h>
#include <MatrixWriter.h>


using namespace utils;
int main(int , char *argv[])
{
	// read the input filenames
	FilenamesReader::FilenamesType labelFilenames = FilenamesReader::Read(argv[1]);
	FilenamesReader::FilenamesType imageFilenames = FilenamesReader::Read(argv[2]);
	

	typedef itk::Vector<unsigned short, 1> MeasurementVectorType;
	typedef itk::Statistics::ListSample<MeasurementVectorType> SampleType;
	SampleType::Pointer positiveSamples = SampleType::New();
	SampleType::Pointer negativeSamples = SampleType::New();
	positiveSamples->SetMeasurementVectorSize(1);
	negativeSamples->SetMeasurementVectorSize(1);

	// build the positive and negative sample list
	for(unsigned int i = 0; i < labelFilenames.size(); i++)
	{
		LabelVolume::Pointer label = LabelVolumeIO::Read(labelFilenames[i]);
		ImageVolume::Pointer image = ImageVolumeIO::Read(imageFilenames[i]);


		// get the samples from the positive samples
		typedef itk::Statistics::ImageToListSampleFilter<ImageVolume, LabelVolume> ListSampleGenerator;
		ListSampleGenerator::Pointer generator = ListSampleGenerator::New();
		generator->SetInput(image);
		generator->SetMaskImage(label);
		generator->SetMaskValue(1);
		generator->Update();

		const ListSampleGenerator::ListSampleType * partList = generator->GetOutput();

		// add to the set of positive values
		for(unsigned int j = 0; j < partList->Size(); j++)
		{
			MeasurementVectorType mv;
			mv[0] = partList->GetMeasurementVector(j)[0];
			positiveSamples->PushBack(mv);			
		}


		// do it for the negative samples
		generator->SetMaskValue(0);
		generator->Modified();
		generator->Update();


		const ListSampleGenerator::ListSampleType * partList2 = generator->GetOutput();

		// add to the set of positive values
		for(unsigned int j = 0; j < partList2->Size(); j++)
		{
			MeasurementVectorType mv;
			mv[0] = partList2->GetMeasurementVector(j)[0];
			negativeSamples->PushBack(mv);			
		}

	}



	// now we compute the gmm for the positive class
	const unsigned int positiveModes = 1;
	typedef itk::Array<double> ParametersType;
	ParametersType positiveParams(2);
	std::vector<ParametersType> initialParameters(positiveModes);
	positiveParams[0] = 160;
	positiveParams[1] = 1500;

	initialParameters[0] = positiveParams;


	// create the positive component
	typedef itk::Statistics::GaussianMixtureModelComponent<SampleType> ComponentType;
	std::vector<ComponentType::Pointer> components;
	ComponentType::Pointer positiveComponent = ComponentType::New();
	components.push_back(positiveComponent);
	positiveComponent->SetSample(positiveSamples);
	positiveComponent->SetParameters(initialParameters[0]);


	// do the EM
	typedef itk::Statistics::ExpectationMaximizationMixtureModelEstimator<SampleType> EstimatorType;
	EstimatorType::Pointer positiveEstimator = EstimatorType::New();
	positiveEstimator->SetSample(positiveSamples);
	positiveEstimator->SetMaximumIteration(200);

	itk::Array<double> initialProportions(1);
	initialProportions[0] = 1.0;

	positiveEstimator->SetInitialProportions(initialProportions);
	positiveEstimator->AddComponent((ComponentType::Superclass*) positiveComponent.GetPointer());

	std::cout << positiveComponent->GetFullParameters() << std::endl;

	positiveEstimator->Update();

	std::cout << positiveComponent->GetFullParameters() << std::endl;



	// do the negative em
	const unsigned int negativeModes = 2;
	ParametersType negativeParameters(2);
	std::vector<ParametersType> initialNegativeParams(negativeModes);

	negativeParameters[0] = 57.0;
	negativeParameters[1] = 2450.0;
	initialNegativeParams[0] = negativeParameters;
	negativeParameters[0] = 237.0;
	negativeParameters[1] = 20510.0;
	initialNegativeParams[1] = negativeParameters;


	std::vector<ComponentType::Pointer> negativeComponents;
	for ( unsigned int i = 0; i < negativeModes; i++ )
	{
		negativeComponents.push_back( ComponentType::New() );
		(negativeComponents[i])->SetSample( negativeSamples );
		(negativeComponents[i])->SetParameters( initialNegativeParams[i] );
	}

	EstimatorType::Pointer negativeEstimator = EstimatorType::New();
	negativeEstimator->SetSample(negativeSamples);
	negativeEstimator->SetMaximumIteration(200);
	
	itk::Array<double> negativeProportions(negativeModes);
	negativeProportions[0] = 0.5;
	negativeProportions[1] = 0.5;
	negativeEstimator->SetInitialProportions(negativeProportions);

	for ( unsigned int i = 0; i < negativeModes; i++)
	{
		negativeEstimator->AddComponent( (ComponentType::Superclass*)
				(negativeComponents[i]).GetPointer() );
	}

	negativeEstimator->Update();

	for ( unsigned int i = 0; i < 2; i++ )
	{
		std::cout << "Cluster[" << i << "]" << std::endl;
		std::cout << "    Parameters:" << std::endl;
		std::cout << "         " << (negativeComponents[i])->GetFullParameters()
			<< std::endl;
		std::cout << "    Proportion: ";
		std::cout << "         " << negativeEstimator->GetProportions()[i] << std::endl;
	}


	typedef utils::MatrixWriter MatrixWriter;
	typedef utils::MatrixDataSet DataSetType;
	utils::DoubleMatrixType negOutput = utils::DoubleMatrixType::Zero(negativeModes, 3);
	negOutput(0,0) = negativeComponents[0]->GetFullParameters()[0];
	negOutput(0,1) = negativeComponents[0]->GetFullParameters()[1];
	negOutput(0,2) = negativeEstimator->GetProportions()[0];
	
	negOutput(1,0) = negativeComponents[1]->GetFullParameters()[0];
	negOutput(1,1) = negativeComponents[1]->GetFullParameters()[1];
	negOutput(1,2) = negativeEstimator->GetProportions()[1];


	utils::DoubleMatrixType posOutput = utils::DoubleMatrixType::Zero(1, 3);
	posOutput(0,0) = positiveComponent->GetFullParameters()[0];
	posOutput(0,1) = positiveComponent->GetFullParameters()[1];
	posOutput(0,2) = positiveEstimator->GetProportions()[0];


	utils::MatrixDataSet::Pointer outputData = utils::MatrixDataSet::New();
	outputData->AddData("NegativeParams", negOutput);
	outputData->AddData("PositiveParams", posOutput);

	MatrixWriter::Pointer matWriter = MatrixWriter::New();
	matWriter->SetInput(outputData);
	matWriter->SetFilename(argv[3]);
	matWriter->Update();


	return 0;
}
