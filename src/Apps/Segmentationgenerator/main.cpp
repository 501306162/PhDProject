#include <iostream>

#include "common.h"
#include "functions.h"

#include <itkSimilarity3DTransform.h>


int main(int, char ** argv)
{
	std::string optionsFilename = argv[1];
	std::string dicomFolder = argv[2];


	// first we need to compute the transform to the origin
	typedef itk::Similarity3DTransform<double> TransformType;

	


	// load the input file
	SeriesTransform::Map transforms;
	readOptionsFile(optionsFilename, transforms);

	// scan the dicom directory
	gdcm::Directory::FilenamesType filenames;
	scanDicomFiles(dicomFolder, filenames);	


	// group the dicom files by the series
	groupDicomFiles(filenames, transforms);

	// now iterate through the images 
	SeriesTransform::Map::iterator mapIt = transforms.begin();
	unsigned int count = 1;
	while(mapIt != transforms.end())
	{
		SeriesTransform &trans = mapIt->second;

		std::cout << "Loading series: " << count  << std::endl;

		// load the images
		ImageList images;
		loadImageSeries(trans.imageFilenames, images);
		
		
		// apply the transforms
		
		
		
		++mapIt; ++count;
	}



	return 0;

}
