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
#include <itkJoinSeriesImageFilter.h>


int main(int, char ** argv)
{
	// get all the input strings
	std::string segmentationDirectory = argv[1];
	std::string xmlFilename = segmentationDirectory + "/simple_example.xml";
	std::string registrationFilename = segmentationDirectory + "/registration.txt";
	std::string seriesLookupFilename = segmentationDirectory + "/series_lookup.txt";
	std::string outputDirectory = argv[2];




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


	while(mapIt != transforms.end())
	{
		SeriesTransform &trans = mapIt->second;
		std::cout << "Loading series: " << count  << std::endl;

		// load the images
		loadImageSeries(trans, xmlOptions.instanceNumber);

		for(unsigned int i = 0; i < xmlOptions.instanceNumber; i++)
		{
			std::stringstream ssLevelSet;
			ssLevelSet << segmentationDirectory << "/phiContour0_t=" << i << "_final.nii";
			std::string levelSetFilename = ssLevelSet.str();

			// load the level set for this instance
			typedef itk::ImageFileReader<LevelSetType> ReaderType;
			ReaderType::Pointer reader = ReaderType::New();
			reader->SetFileName(levelSetFilename);
			reader->SetImageIO(itk::NiftiImageIO::New());
			reader->Update();

			LevelSetType::Pointer levelSet = reader->GetOutput();

			// st up the interpolator
			typedef itk::LinearInterpolateImageFunction<LevelSetType, double> InterpolatorType;
			InterpolatorType::Pointer interpolator = InterpolatorType::New();
			interpolator->SetInputImage(levelSet);



			ImageType::Pointer label = ImageType::New();
			createLabelImage(trans, reference, levelSet, xmlOptions.roiOffset, label, i);
			trans.labelImages.push_back(label);
		}

		++mapIt; ++count;
	}


	// group the series
	std::vector<std::vector<SeriesTransform> > groupedTransforms;
	groupImageSeries(transforms, groupedTransforms);

	// create the lookup file
	std::string lookupFilename = outputDirectory + "/lookup.txt";
	std::ofstream lookupFile;
	lookupFile.open(lookupFilename.c_str());

	if(!lookupFile.is_open())
	{
		std::cout << "Couldn't open the output file" << std::endl;
		exit(1);
	}





	for(unsigned int i = 0; i < groupedTransforms.size(); i++)
	{

		// create the 4d output image
		typedef itk::Image<unsigned short, 4> OutputImageType;
		typedef itk::JoinSeriesImageFilter<ImageType, OutputImageType> JoinerType;
		JoinerType::Pointer imageJoiner = JoinerType::New();
		JoinerType::Pointer labelJoiner = JoinerType::New();

		SeriesTransform::List series = groupedTransforms[i];

		for(unsigned int time = 0; time < xmlOptions.instanceNumber; time++)
		{
			ImageType::Pointer label = ImageType::New();
			ImageType::Pointer image = ImageType::New();

			buildOutput(series, image, label, time);


			imageJoiner->SetInput(time, image);
			labelJoiner->SetInput(time, label);
		}



		std::stringstream ss, ss2;
		ss << outputDirectory << "/label_" << series.front().description << "_" << series.front().dcmSeries << ".nrrd"; 
		ss2 << outputDirectory << "/image_" << series.front().description << "_" << series.front().dcmSeries << ".nrrd"; 


		for(unsigned int j = 0; j < series.size(); j++)
		{
			lookupFile << "label_" << series.front().description << "_" << series.front().dcmSeries <<  ".nrrd" << ":" << j << ":" << series[j].imageFilenames[0] << "\n";			
		}

		labelJoiner->Update();
		imageJoiner->Update();


		typedef itk::ImageFileWriter<OutputImageType> WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetInput(labelJoiner->GetOutput());
		writer->SetFileName(ss.str());
		writer->SetImageIO(itk::NrrdImageIO::New());
		writer->Update();

		writer->SetInput(imageJoiner->GetOutput());
		writer->SetFileName(ss2.str());
		writer->Modified();
		writer->Update();

	}

	lookupFile.close();



	return 0;

}


