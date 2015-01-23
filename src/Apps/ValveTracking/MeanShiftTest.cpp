#include <iostream>
#include <itkImageMomentsCalculator.h>
#include <itkDiscreteGaussianImageFilter.h>

#include "FittingFunctions.h"

#include <opencv2/video/tracking.hpp>

#include <ValveIO.h>
#include <ParameterHelpers.h>
#include <PatchTrainingData.h>
#include <FlipChecker.h>
#include <ValveNormaliser.h>
#include <OpenCVValve.h>
#include <itkImageRegionIterator.h>
#include <CommonDefinitions.h>
#include <itkGaussianKernelFunction.h>
#include <itkGaussianBlurImageFunction.h>

#include <PatchExtractor.h>

using namespace vt;


typedef itk::Image<double, 3> RealImageType;
typedef ValveLine<3> ValveType;
typedef ImageType::RegionType RegionType;
typedef ImageType::PointType PointType;

void alignSequence(const ValveSequence<3>::Pointer &input, ValveSequence<3>::Pointer &output, 
		bool flipImage, bool flipLine);

void extractMask(const ValveType::Pointer &input, const unsigned int pIndex, ImageType::RegionType &region);

void computeProbImage(const PatchParams &params, const unsigned int pIndex, const ClassifierList &cls, 
		const ValveLine<3>::Pointer &valve, RealImageType::Pointer &prob);
void kdeEstimate(const RealImageType::Pointer &prob, const double h, const RegionType &region, RealImageType::Pointer &output);


void track(const ValveType::Pointer &valve, const PointType &start, const RealImageType::Pointer &prob, PointType &end);


int main(int argc, char ** argv)
{
	// load the params 
	PatchParams params(argv[1]);
	const std::string type = argv[4];
	unsigned int num = atoi(argv[3]);


	// load the classifiers 
	PatchTrainingData::Pointer patchData = PatchTrainingData::New();



	// load the valve line 
	ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
	reader->SetFileName(argv[2]);
	ValveSequence<3>::Pointer sequence = reader->GetOutput();

	FlipChecker::Pointer checker = FlipChecker::New();
	bool flipImage = checker->FlipImage(type, num);
	bool flipPoints = checker->FlipImage(type, num);


	ValveLine<3>::Pointer startLine = sequence->GetValveLine(0);
	std::vector<PointType> points;
	points.push_back(startLine->GetP1());
	points.push_back(startLine->GetP2());


	for(unsigned int ts = 1; ts < sequence->GetNumberOfLines(); ts++)
	{
		// set up the test valve 
		ValveType::Pointer testValve = sequence->GetValveLine(ts);
		testValve->SetP1(points[0]);
		testValve->SetP2(points[1]);
		testValve->UpdateIndexs();

		ValveNormaliser::Pointer normaliser = ValveNormaliser::New();
		normaliser->SetInput(testValve);
		normaliser->SetFlip(flipImage);
		normaliser->SetFlipPoints(flipPoints);
		normaliser->Normalise();

		ValveType::Pointer aligned = normaliser->GetOutput();

		params.timeStep = ts;
		ClassifierMap classifiers;
		loadOrTrainClassifiers(params, 1000, patchData, classifiers);
		
			
		// loop through the points
		for(unsigned int p = 0; p < 2; p++)
		{
			// compute the probability image 
			RealImageType::Pointer prob = RealImageType::New();
			computeProbImage(params, p, classifiers[type], aligned, prob);

			PointType end;
			track(aligned, aligned->GetPoint(p+1), prob, end);

			points[p] = end;
		}

		aligned->SetP1(points[0]);
		aligned->SetP2(points[1]);


		normaliser->SetInput(aligned);
		normaliser->UnNormalise();
		ValveType::Pointer output = normaliser->GetOutput();

		points[0] = output->GetP1();
		points[1] = output->GetP2();

		std::stringstream ss;
		ss << "o" << ts;
		ValveWriter<3>::Pointer writer = ValveWriter<3>::New();
		writer->SetInput(output);
		writer->SetFileName(ss.str());
		writer->Write();




	}

	return 0;
}

// ------------------------------------------------------------------------
void track(const ValveType::Pointer &valve, const PointType &start, const RealImageType::Pointer &prob, PointType &end)
{
	const double eps = 1.0;
	const unsigned int ts = 100;
	end = start;
	for(unsigned int i = 0; i < ts; i++)
	{
		typedef itk::ImageMomentsCalculator<RealImageType> MomentsCalculatorType;
		MomentsCalculatorType::Pointer calc = MomentsCalculatorType::New();


		// get the patch region 
		ContIndexType contInd;
		prob->TransformPhysicalPointToContinuousIndex(end, contInd);
		
		RealPatchExtractor::Pointer extractor =RealPatchExtractor::New();
		extractor->SetImage(prob);
		
		RealPatchExtractor::SizeType size;
		size.Fill(10);

		extractor->SetPatchSize(size);
		extractor->SetPatchCenter(contInd);

		calc->SetImage(extractor->ExtractPatch());
		calc->Compute();

		std::cout << calc->GetCenterOfGravity() << std::endl;
		
		PointType centroid = calc->GetCenterOfGravity();
		double diff = centroid.SquaredEuclideanDistanceTo(end);

		//for(unsigned int j = 0; j < 2; j++)
		//{
			//end[j] = centroid[j];
		//}
		end = centroid;

		if(diff < eps)
			break;


	}



	
}

// ------------------------------------------------------------------------
void extractMask(const ValveType::Pointer &input, const unsigned int pIndex, ImageType::RegionType &region)
{
	ImagePatchExtractor::Pointer extractor = ImagePatchExtractor::New();
	extractor->SetImage(input->GetImage());
	extractor->SetPatchCenter(input->GetIndex(pIndex));	

	MaskPatchExtractor::SizeType size;
	size.Fill(50);
	extractor->SetPatchSize(size);
	extractor->SetMaskValue(255);
	
	region = extractor->ExtractionRegion();
}




// ------------------------------------------------------------------------
void computeProbImage(const PatchParams &params, const unsigned int pIndex, const ClassifierList &cls, 
		const ValveLine<3>::Pointer &valve, RealImageType::Pointer &prob)
{
	// create the output image 
	ImageType::Pointer image = valve->GetImage();
	prob->SetOrigin(image->GetOrigin());
	prob->SetRegions(image->GetLargestPossibleRegion());
	prob->SetSpacing(image->GetSpacing());
	prob->SetDirection(image->GetDirection());
	prob->Allocate();
	prob->FillBuffer(0.0);


	// get the extraction region
	ImageType::RegionType region;
	extractMask(valve, pIndex+1, region);


	// get all the indices
	itk::ImageRegionIterator<RealImageType> imIt(prob, region);
	while(!imIt.IsAtEnd())
	{
		ImageType::IndexType index = imIt.GetIndex();
		ContIndexType contInd;
		for(unsigned int i = 0; i < 3; i++)
		{
			contInd[i] = index[i];
		}

		MatrixType feature;
		//extractLBPFeature(params, image, contInd, feature);

		MatrixType probs;
		IntMatrixType classes;
		cls[pIndex]->PredictProbability(feature, classes, probs);
		imIt.Set(probs(0,1));

		++imIt;
	}


	typedef itk::DiscreteGaussianImageFilter<RealImageType, RealImageType> BlurType;
	BlurType::Pointer blur = BlurType::New();
	blur->SetInput(prob);
	BlurType::ArrayType vars;
	vars.Fill(1.0);
	blur->SetVariance(vars);
	blur->Update();

	prob = blur->GetOutput();
}

// ------------------------------------------------------------------------
void kdeEstimate(const RealImageType::Pointer &prob, const double h, 
		const RegionType &region, RealImageType::Pointer &output)
{
	output->SetDirection(prob->GetDirection());
	output->SetSpacing(prob->GetSpacing());
	output->SetOrigin(prob->GetOrigin());
	output->SetRegions(prob->GetLargestPossibleRegion());
	output->Allocate();
	output->FillBuffer(0);

	unsigned int numPixels = region.GetNumberOfPixels();

	itk::ImageRegionIterator<RealImageType> probIt(prob, region);
	itk::ImageRegionIterator<RealImageType> outIt(output, region);

	while(!outIt.IsAtEnd())
	{
		ImageType::IndexType ind1 = outIt.GetIndex();
		ImageType::PointType p1;
		output->TransformIndexToPhysicalPoint(ind1, p1);

		double val = 0;
		probIt.GoToBegin();
		while(!probIt.IsAtEnd())
		{
			ImageType::IndexType ind2 = probIt.GetIndex();
			ImageType::PointType p2;
			prob->TransformIndexToPhysicalPoint(ind2, p2);

			itk::GaussianKernelFunction<double>::Pointer kernal = itk::GaussianKernelFunction<double>::New();

			ImageType::PointType pdiff = p1-p2;
			double psum = 1.0;
			for(unsigned int i = 0; i < 2; i++)
			{
				psum *= kernal->Evaluate(pdiff[i]/h);
				psum *= (1.0/h);
			}

			

			val += psum;
			++probIt;
		}

		outIt.Set(val/ static_cast<double>(numPixels));
		
		++outIt;
	}

}

// ------------------------------------------------------------------------
void alignSequence(const ValveSequence<3>::Pointer &input, ValveSequence<3>::Pointer &output, 
		bool flipImage, bool flipLine)
{
	for(unsigned int i = 0; i < input->GetNumberOfLines(); i++)
	{
		ValveNormaliser::Pointer normaliser = ValveNormaliser::New();
		normaliser->SetInput(input->GetValveLine(i));
		normaliser->SetFlip(flipImage);
		normaliser->SetFlipPoints(flipLine);
		normaliser->Normalise();

		output->AddValveLine(normaliser->GetOutput());
	}

}
