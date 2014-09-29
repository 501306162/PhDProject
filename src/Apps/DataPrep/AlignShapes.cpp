#include <iostream>

#include <FilenamesReader.h>
#include <CommonDefinitions.h>

#include <RegistrationHelpers.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>

using namespace utils;
int main(int argc, char *argv[])
{
	// read the filenames
	typedef utils::FilenamesReader FilenamesReader;
	FilenamesReader::FilenamesType inputFilenames = FilenamesReader::Read(argv[1]);
	FilenamesReader::FilenamesType outputFilenames = FilenamesReader::Read(argv[2]);
	std::string meanFilename = argv[3];


	// load all of the images
	utils::LabelVolumeList labels;
	for(unsigned int i = 0; i < inputFilenames.size(); i++)
	{
		utils::LabelVolume::Pointer label = utils::LabelVolumeIO::Read(inputFilenames[i]);
		labels.push_back(label);
	}

	RealVolume::Pointer sumImage = RealVolume::New();


	// select the fixed image for the first set of alignments
	utils::LabelVolume::Pointer target = labels[0];
	for(unsigned int i = 0; i < labels.size(); i++)
	{
		std::cout << "Processing: " << i+1 << " / " << labels.size() << std::endl;
		utils::Rigid3DTransformType::Pointer transform = utils::Rigid3DTransformType::New();
		utils::RegisterLabelVolumeImages(target, labels[i], transform);

		utils::LabelVolume::Pointer output = utils::LabelVolume::New();
		utils::Apply3DRigidTransformTOLabelVolume(labels[i], target, transform, output);


		// compute the signed distance function of the shape
		typedef itk::SignedMaurerDistanceMapImageFilter<LabelVolume, RealVolume> DistanceMapType;
		DistanceMapType::Pointer distanceMapFilter = DistanceMapType::New();
		distanceMapFilter->SetUseImageSpacing(true);
		distanceMapFilter->SetSquaredDistance(false);
		distanceMapFilter->SetBackgroundValue(0);
		distanceMapFilter->SetInput(output);

		distanceMapFilter->Update();

		RealVolume::Pointer distance = distanceMapFilter->GetOutput();

		if(i == 0)
		{
			sumImage = distance;
		}
		else
		{
			// add the images
			typedef itk::AddImageFilter<RealVolume, RealVolume, RealVolume> AddImageType;
			AddImageType::Pointer adder = AddImageType::New();
			adder->SetInput1(sumImage);
			adder->SetInput2(distance);
			adder->Update();
			sumImage = adder->GetOutput();
		}

		LabelVolumeIO::Write(outputFilenames[i], output);
	}

	// do the image division to get the mean
	typedef itk::DivideImageFilter<RealVolume, RealVolume, RealVolume> DividerType;
	DividerType::Pointer divider = DividerType::New();
	divider->SetInput(sumImage);
	divider->SetConstant(static_cast<double>(labels.size()));
	divider->Update();


	typedef itk::BinaryThresholdImageFilter<RealVolume, LabelVolume> ThresholderType;
	ThresholderType::Pointer thresholder = ThresholderType::New();
	thresholder->SetInput(divider->GetOutput());
	thresholder->SetUpperThreshold(0.0);
	thresholder->SetOutsideValue(0.0);
	thresholder->SetInsideValue(1.0);

	thresholder->Update();


	LabelVolumeIO::Write(meanFilename, thresholder->GetOutput());

	
	return 0;
}
