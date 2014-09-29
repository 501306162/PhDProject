#include <iostream>

#include <FilenamesReader.h>
#include <CommonDefinitions.h>

#include <RegistrationHelpers.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

using namespace utils;
int main(int argc, char *argv[])
{
	// read the filenames
	typedef utils::FilenamesReader FilenamesReader;
	FilenamesReader::FilenamesType labelFilenames = FilenamesReader::Read(argv[1]);
	FilenamesReader::FilenamesType imageFilenames = FilenamesReader::Read(argv[2]);
	FilenamesReader::FilenamesType outputImageFilenames = FilenamesReader::Read(argv[3]);
	FilenamesReader::FilenamesType outputLabelFilenames = FilenamesReader::Read(argv[4]);
	FilenamesReader::FilenamesType outputDistanceFilenames = FilenamesReader::Read(argv[5]);


	// load all of the images
	utils::LabelVolumeList labels;
	utils::ImageVolumeList images;
	for(unsigned int i = 0; i < labelFilenames.size(); i++)
	{
		utils::LabelVolume::Pointer label = utils::LabelVolumeIO::Read(labelFilenames[i]);
		labels.push_back(label);

		utils::ImageVolume::Pointer image = utils::ImageVolumeIO::Read(imageFilenames[i]);
		images.push_back(image);
	}


	// select the fixed image for the first set of alignments
	utils::LabelVolume::Pointer target = labels[0];
	for(unsigned int i = 0; i < labels.size(); i++)
	{
		std::cout << "Processing: " << i+1 << " / " << labels.size() << std::endl;
		utils::Rigid3DTransformType::Pointer transform = utils::Rigid3DTransformType::New();
		utils::RegisterLabelVolumeImages(target, labels[i], transform);

		typedef itk::ResampleImageFilter<utils::ImageVolume, utils::ImageVolume> ResamplerType;
		typedef itk::LinearInterpolateImageFunction<utils::ImageVolume, double> InterpolatorType;

		// resample the output image
		ResamplerType::Pointer resampler = ResamplerType::New();
		InterpolatorType::Pointer interpolator = InterpolatorType::New();

		resampler->SetTransform(transform);
		resampler->SetInput(images[i]);
		resampler->SetInterpolator(interpolator);
		resampler->SetOutputParametersFromImage(images[0]);
		resampler->SetUseReferenceImage(true);
		resampler->SetDefaultPixelValue(0);
		resampler->Update();

		utils::ImageVolumeIO::Write(outputImageFilenames[i], resampler->GetOutput());


		// resample the output label
		typedef itk::ResampleImageFilter<utils::LabelVolume, utils::LabelVolume> LabelResampler;
		LabelResampler::Pointer labelResampler = LabelResampler::New();
		typedef itk::NearestNeighborInterpolateImageFunction<utils::LabelVolume, double> LabelInterpolatorType;
		LabelInterpolatorType::Pointer labelInterpolator = LabelInterpolatorType::New();
		labelResampler->SetInput(labels[i]);
		labelResampler->SetInterpolator(labelInterpolator);
		labelResampler->SetTransform(transform);
		labelResampler->SetOutputParametersFromImage(images[0]);
		labelResampler->SetDefaultPixelValue(0);
		labelResampler->Update();


		utils::LabelVolumeIO::Write(outputLabelFilenames[i], labelResampler->GetOutput());



		// compute the distance maps
		typedef itk::SignedMaurerDistanceMapImageFilter<utils::LabelVolume, utils::RealVolume> DistanceMapFilterType;
		DistanceMapFilterType::Pointer distanceFilter = DistanceMapFilterType::New();
		distanceFilter->SetInput(labelResampler->GetOutput());
		distanceFilter->SetSquaredDistance(false);
		distanceFilter->SetUseImageSpacing(true);
		distanceFilter->SetBackgroundValue(0);
		distanceFilter->Update();


		utils::RealVolumeIO::Write(outputDistanceFilenames[i], distanceFilter->GetOutput());
	}


	
	return 0;
}
