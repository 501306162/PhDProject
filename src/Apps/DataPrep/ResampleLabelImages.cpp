#include <iostream>

#include <FilenamesReader.h>
#include <CommonDefinitions.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkIdentityTransform.h>
#include <itkConstantPadImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageSliceConstIteratorWithIndex.h>
#include <itkSliceBySliceImageFilter.h>
#include <itkRecursiveGaussianImageFilter.h>

using namespace utils;

void CreateIntermediateImage(LabelVolume::Pointer &input, LabelVolume::Pointer &output);
void Resample(RealVolume::Pointer input, RealVolume::Pointer &output);
int main(int argc, char *argv[])
{

	// read the input filenames
	FilenamesReader::FilenamesType inputFilenames = FilenamesReader::Read(argv[1]);
	FilenamesReader::FilenamesType outputFilenames = FilenamesReader::Read(argv[2]);
	for(unsigned int i = 0; i < inputFilenames.size(); i++)
	{
		// load the image
		LabelVolume::Pointer label = LabelVolumeIO::Read(inputFilenames[i]);

		typedef itk::SliceBySliceImageFilter<LabelVolume, RealVolume> SliceBySliceImageFilterType;

		typedef itk::SignedMaurerDistanceMapImageFilter<LabelSlice, RealSlice> DistanceMapType;
		DistanceMapType::Pointer distanceMap = DistanceMapType::New();
		distanceMap->SetBackgroundValue(0);
		distanceMap->SetUseImageSpacing(true);
		distanceMap->SetSquaredDistance(false);

		SliceBySliceImageFilterType::Pointer sliceBySliceFilter = SliceBySliceImageFilterType::New();
		sliceBySliceFilter->SetInput(label);
		sliceBySliceFilter->SetFilter(distanceMap);
		sliceBySliceFilter->Update();

		RealVolume::Pointer map = sliceBySliceFilter->GetOutput();
		RealVolume::Pointer map2 = RealVolume::New();
		Resample(map, map2);



		typedef itk::BinaryThresholdImageFilter<RealVolume, LabelVolume> ThresholderType;
		ThresholderType::Pointer thresholder = ThresholderType::New();
		thresholder->SetInput(map2);
		thresholder->SetUpperThreshold(0.0);
		thresholder->SetOutsideValue(0.0);
		thresholder->SetInsideValue(1.0);

		thresholder->Update();


		LabelVolumeIO::Write(outputFilenames[i], thresholder->GetOutput());
		
		
		


	

		
	}
	
	return 0;
}

void Resample(RealVolume::Pointer input, RealVolume::Pointer &output)
{
	// first we smooth along the x and y dimension
	typedef itk::RecursiveGaussianImageFilter<RealVolume, RealVolume> SmootherType;
	SmootherType::Pointer smootherX = SmootherType::New();
	SmootherType::Pointer smootherY = SmootherType::New();

	smootherX->SetInput(input);
	smootherY->SetInput(smootherX->GetOutput());
	
	// set up the resampler 
	RealVolume::SpacingType inputSpacing = input->GetSpacing();
	double isoSpacing = std::sqrt(inputSpacing[0] * inputSpacing[2]);
	isoSpacing = inputSpacing[0];

	smootherX->SetDirection(0);
	smootherY->SetDirection(1);

	smootherX->SetSigma(isoSpacing);
	smootherY->SetSigma(isoSpacing);


	typedef itk::ResampleImageFilter<RealVolume, RealVolume> ResamplerType;
	ResamplerType::Pointer resampler = ResamplerType::New();
	

	// set the identity transform
	typedef itk::IdentityTransform< double, 3 >  TransformType;
	TransformType::Pointer transform = TransformType::New();
	transform->SetIdentity();
	resampler->SetTransform( transform );

	typedef itk::LinearInterpolateImageFunction<RealVolume, double>  InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	resampler->SetInterpolator(interpolator);
	
	resampler->SetDefaultPixelValue(100);

	RealVolume::SpacingType spacing;
	spacing.Fill(isoSpacing);

	resampler->SetOutputSpacing(spacing);
	resampler->SetOutputOrigin(input->GetOrigin());
	resampler->SetOutputDirection(input->GetDirection());


	// get the new size
	RealVolume::SizeType inputSize = input->GetLargestPossibleRegion().GetSize();
	const double dx = inputSize[0] * inputSpacing[0] / isoSpacing;
	const double dy = inputSize[1] * inputSpacing[1] / isoSpacing;
	const double dz = (inputSize[2] - 1 ) * inputSpacing[2] / isoSpacing;


	typedef RealVolume::SizeValueType SizeValueType;
	RealVolume::SizeType   size;
	size[0] = static_cast<SizeValueType>( dx );
	size[1] = static_cast<SizeValueType>( dy );
	size[2] = static_cast<SizeValueType>( dz );
	resampler->SetSize( size );	


	resampler->SetInput(input);
	resampler->Update();

	output = resampler->GetOutput();




}

void CreateIntermediateImage(LabelVolume::Pointer &input, LabelVolume::Pointer &output)
{
	// compute the output spacing etc
	double inputZSpacing = input->GetSpacing()[2];
	double outputZSpacing = input->GetSpacing()[1];
	int inputSize = input->GetLargestPossibleRegion().GetSize()[2];
	double inputDistance  =inputZSpacing*inputSize;

	int outputSize = (int) inputDistance / outputZSpacing;

	LabelVolume::SizeType size;
	size = input->GetLargestPossibleRegion().GetSize();
	size[2] = outputSize;

	LabelVolume::SpacingType spacing = input->GetSpacing();
	spacing[2] = outputZSpacing;

	LabelVolume::RegionType region = input->GetLargestPossibleRegion();
	region.SetSize(size);


	// get the z positions of the current slices
	typedef itk::ImageSliceConstIteratorWithIndex<LabelVolume> IteratorType;
	IteratorType it(input, input->GetLargestPossibleRegion());
	it.SetFirstDirection(0);
	it.SetSecondDirection(1);
	it.GoToBegin();

	std::vector<double> sliceLocations;
	while(!it.IsAtEnd())
	{
		int pxSum = 0;
		while(!it.IsAtEndOfSlice())
		{
			while(!it.IsAtEndOfLine())
			{
				LabelVolume::ValueType val = it.Get();
				if(val == 1)
				{
					pxSum++;
				}
				++it;
			}
			it.NextLine();
		}

		if(pxSum > 0)
		{
			LabelVolume::IndexType index = it.GetIndex();
			LabelVolume::PointType point;
			input->TransformIndexToPhysicalPoint(index, point);
			sliceLocations.push_back(point[2]);
		}
		it.NextSlice();
	}




}
