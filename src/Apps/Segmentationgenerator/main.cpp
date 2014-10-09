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
	
	// TODO update this
	xmlOptions.instanceNumbers.push_back(0);


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
		loadImageSeries(trans, xmlOptions.instanceNumbers);


		// create the output image 
		ImageType::Pointer output = ImageType::New();
		output->SetDirection(trans.images[0]->GetDirection());
		output->SetSpacing(trans.images[0]->GetSpacing());
		output->SetOrigin(trans.images[0]->GetOrigin());
		output->SetRegions(trans.images[0]->GetLargestPossibleRegion());
		output->Allocate();
		output->FillBuffer(0);
		
		
		typedef itk::ImageRegionIterator<ImageType> ItType;
		ItType it(output, output->GetLargestPossibleRegion());
		while(!it.IsAtEnd())
		{
			ImageType::IndexType index = it.GetIndex();
			ImageType::PointType point;
			output->TransformIndexToPhysicalPoint(index, point);

			VectorType p1 = point.GetVnlVector();
			VectorType p2;

			transformPointToPatient(reference, trans, xmlOptions, p1, p2);

			ImageType::PointType point2;
			point2[0] = p2[0];
			point2[1] = p2[1];
			point2[2] = p2[2];
			
			if(interpolator->IsInsideBuffer(point2))
			{
				float val = interpolator->Evaluate(point2);
				if(val > 0)
				{
					it.Set(1);
				}
			}

			++it;
		}


		std::stringstream ss;
		ss << "label_" << trans.description << "_" << trans.dcmSeries << ".nrrd";

		typedef itk::ImageFileWriter<ImageType> WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName(ss.str());
		writer->SetInput(output);
		writer->SetImageIO(itk::NrrdImageIO::New());
		writer->Update();

		std::stringstream ss2;
		ss2 << "image_" << trans.description << "_" << trans.dcmSeries << ".nrrd";

		WriterType::Pointer writer2 = WriterType::New();
		writer2->SetFileName(ss2.str());
		writer2->SetInput(trans.images[0]);
		writer2->SetImageIO(itk::NrrdImageIO::New());
		writer2->Update();


		

		++mapIt; ++count;

		/*

		// normalise the images
		ImageType::Pointer image = trans.images[0];
		ImageType::Pointer normalised = ImageType::New();
		NormaliseImage(image, reference, normalised, trans);

		trans.normalisedImages.push_back(normalised);

		// compute the bounds of the normalised image

		std::cout << trans.description << " " << trans.dcmSeries << " " << normalised->GetOrigin() << std::endl;
		*/
	}
	return 0;
	



	// compute the translation using the min maxs and the rio offset
	std::cout << "Applying Transforms" << std::endl;
	mapIt = transforms.begin();
	while(mapIt != transforms.end())
	{
		SeriesTransform trans = mapIt->second;
		std::stringstream ss;
		ss << "output_" << trans.description << "_" << trans.dcmSeries << ".nrrd";

		ImageType::Pointer test = trans.normalisedImages[0];
		ImageType::PointType orig = test->GetOrigin();
		std::cout << orig << std::endl;
		orig[0] -= xmlOptions.roiOffset[0];
		orig[1] -= xmlOptions.roiOffset[1];
		orig[2] -= xmlOptions.roiOffset[2];


		std::cout << trans.toString() << std::endl;
		orig[0] -= (trans.translation[0]);// * (test->GetSpacing()[0] / reference->GetSpacing()[0]));
		orig[1] -= (trans.translation[1]);// * (test->GetSpacing()[1] / reference->GetSpacing()[1]));
		orig[2] -= (trans.translation[2]);// * (test->GetSpacing()[2] / reference->GetSpacing()[2]));



		test->SetOrigin(orig);
		std::cout << test->GetOrigin() << std::endl;
		typedef itk::ImageFileWriter<ImageType> WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetInput(test);
		writer->SetImageIO(itk::NrrdImageIO::New());
		writer->SetFileName(ss.str());
		writer->Update();

		++mapIt;
	}


	return 0;

}


