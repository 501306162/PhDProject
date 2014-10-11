#include <iostream>

#include "common.h"
#include "functions.h"
#include "transforms.h"
#include "read_series_transforms.h"
#include "InitialTransformExtractor.h"

#include <itkSimilarity3DTransform.h>
#include <itkImageRegionIterator.h>
#include <itkLinearInterpolateImageFunction.h>

#include <itkScalableAffineTransform.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkNrrdImageIO.h>
#include <itkNiftiImageIO.h>


int main(int, char ** argv)
{
	// get all the input strings
	std::string segmentationDirectory = argv[1];
	std::string xmlFilename = segmentationDirectory + "/simple_example.xml";
	std::string registrationFilename = segmentationDirectory + "/registration.txt";
	std::string seriesLookupFilename = segmentationDirectory + "/series_lookup.txt";
	std::string levelSetFilename = segmentationDirectory + "/phiContour0_t=0_final.nii";


	// load the level set
	typedef itk::ImageFileReader<LevelSetType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(levelSetFilename);
	reader->SetImageIO(itk::NiftiImageIO::New());
	reader->Update();

	LevelSetType::Pointer levelSet = reader->GetOutput();


	// read the xml file
	std::cout << " --> Reading XML file" << std::endl;
	OptionsData xmlOptions;
	readXMLValues(xmlFilename, xmlOptions);
	
	// get the number of instances by looking at the number of output level sets
	xmlOptions.instanceNumber = getNumberOfInstances(segmentationDirectory);



	// set up the series lookup
	std::cout << " --> Reading the transformations" << std::endl;
	SeriesTransform::Map transforms;
	readSeriesTransforms(registrationFilename, seriesLookupFilename, transforms);


	// find the reference image
	std::cout << " --> Finding initial transform" << std::endl;
	InitialTransformExtractor extractor;
	extractor.SetOptions(xmlOptions);
	extractor.SetDicomDir(xmlOptions.dataDirectory);
	extractor.Compute();
	ImageType::Pointer reference = extractor.GetReference();


	// scan the dicom directory
	std::cout << " --> Scanning and grouping dicom images" << std::endl;
	gdcm::Directory::FilenamesType filenames;
	scanDicomFiles(xmlOptions.dataDirectory, filenames);	
	groupDicomFiles(filenames, transforms);


	// iterate through the series
	SeriesTransform::Map::iterator mapIt = transforms.begin();
	unsigned int count = 1;

	// st up the interpolator
	typedef itk::LinearInterpolateImageFunction<LevelSetType, double> InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	interpolator->SetInputImage(levelSet);


	while(mapIt != transforms.end())
	{
		SeriesTransform &trans = mapIt->second;
		std::cout << "Loading series: " << count  << std::endl;

		// load the images
		loadImageSeries(trans, xmlOptions.instanceNumber);

		for(unsigned int i = 0; i < xmlOptions.instanceNumber; i++)
		{
			ImageType::Pointer label = ImageType::New();
			createLabelImage(trans, reference, levelSet, xmlOptions.roiOffset, label);
			trans.labelImages.push_back(label);
		}

		++mapIt; ++count;
	}


	// group the series
	std::vector<std::vector<SeriesTransform> > groupedTransforms;
	groupImageSeries(transforms, groupedTransforms);

	for(unsigned int i = 0; i < groupedTransforms.size(); i++)
	{
		SeriesTransform::List series = groupedTransforms[i];
		ImageType::Pointer label = ImageType::New();
		ImageType::Pointer image = ImageType::New();

		buildOutput(series, image, label, 0);

		std::stringstream ss;
		ss << "image_" << i << ".nrrd";
		std::stringstream ss2;
		ss2 << "label_" << i << ".nrrd";

		typedef itk::ImageFileWriter<ImageType> WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetInput(label);
		writer->SetFileName(ss2.str());
		writer->SetImageIO(itk::NrrdImageIO::New());
		writer->Update();

		writer->SetInput(image);
		writer->SetFileName(ss.str());
		writer->Modified();
		writer->Update();

	}




	return 0;

}


