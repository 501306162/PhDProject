#include <iostream>

#include <FilenamesReader.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <CommonDefinitions.h>


using namespace utils;
int main(int argc, char *argv[])
{
	// load up the images
	FilenamesReader::FilenamesType inputFilenames = FilenamesReader::Read(argv[1]);
	FilenamesReader::FilenamesType outputFilenames = FilenamesReader::Read(argv[2]);

	// apply to all of the images
	for(unsigned int i = 0; i < inputFilenames.size(); i++)
	{
		LabelVolume::Pointer label = LabelVolumeIO::Read(inputFilenames[i]);

		typedef itk::SignedMaurerDistanceMapImageFilter<LabelVolume, RealVolume> DistanceMapType;
		DistanceMapType::Pointer distanceMap = DistanceMapType::New();
		distanceMap->SetInput(label);
		distanceMap->SetUseImageSpacing(true);
		distanceMap->SetSquaredDistance(false);
		distanceMap->SetBackgroundValue(0);
		distanceMap->Update();

		RealVolumeIO::Write(outputFilenames[i], distanceMap->GetOutput());
	}
	
	
	return 0;
}
