#include <iostream>

#include "common.h"
#include "functions.h"
#include "InitialTransformExtractor.h"

#include <itkSimilarity3DTransform.h>

#include <itkScalableAffineTransform.h>
#include <itkPermuteAxesImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkNrrdImageIO.h>

void NormaliseImage(const ImageType::Pointer &input, 
		const ImageType::Pointer &reference,	
		ImageType::Pointer &output);


typedef struct _bounds
{
	ImageType::PointType corners[4];
} BoundsType;

void getImageBounds(const ImageType::Pointer &image, 
		BoundsType &bounds);

void getMinMaxBounds(const std::vector<BoundsType> &bounds, double * minMaxs);

int main(int, char ** argv)
{
	std::string xmlFilename = argv[3];
	std::string optionsFilename = argv[1];
	std::string dicomFolder = argv[2];


	// read the xml file
	std::cout << "Reading XML file" << std::endl;
	OptionsData xmlOptions;
	readXMLValues(xmlFilename, xmlOptions);


	// first we need to compute the transform to the origin
	typedef itk::Similarity3DTransform<double> TransformType;

	std::cout << "Finding initial transform" << std::endl;
	InitialTransformExtractor extractor;
	extractor.SetOptions(xmlOptions);
	extractor.SetDicomDir(dicomFolder);
	extractor.Compute();


	ImageType::Pointer reference = extractor.GetReference();
	



	

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

	std::vector<ImageType::Pointer> normalisedImages;
	std::vector<BoundsType> boundsList;

	while(mapIt != transforms.end())
	{
		SeriesTransform &trans = mapIt->second;

		std::cout << "Loading series: " << count  << std::endl;

		// load the images
		ImageList images;
		loadImageSeries(trans.imageFilenames, images);
		
		// apply the transforms
		ImageType::Pointer image = images[0];
		ImageType::Pointer normalised = ImageType::New();
		NormaliseImage(image, reference, normalised);
		
		std::stringstream ss;
		ss << "output" << trans.dcmSeries << ".nrrd";

		typedef itk::ImageFileWriter<ImageType> WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetInput(normalised);
		writer->SetImageIO(itk::NrrdImageIO::New());
		writer->SetFileName(ss.str());
		writer->Update();

		normalisedImages.push_back(normalised);

		BoundsType bounds;
		getImageBounds(normalised, bounds);
		boundsList.push_back(bounds);

		++mapIt; ++count;
	}


	// get the min max bounds
	double minMaxs[6];
	getMinMaxBounds(boundsList, minMaxs);

	


	// compute the translation using the min maxs and the rio offset
	ImageType::PointType offset;
	offset[0] = minMaxs[0] - xmlOptions.roiOffset[0];
	offset[1] = minMaxs[2] - xmlOptions.roiOffset[1];
	offset[2] = minMaxs[4] - xmlOptions.roiOffset[2];



	

	for(unsigned int i = 0; i < normalisedImages.size(); i++)
	{
		std::stringstream ss;
		ss << "output_better" << i << ".nrrd";

		ImageType::Pointer test = normalisedImages[i];
		ImageType::PointType orig = test->GetOrigin();
		orig[0] -= minMaxs[0] -  xmlOptions.roiOffset[0];
		orig[1] -= minMaxs[2] - xmlOptions.roiOffset[1];
		orig[2] -= minMaxs[4] - xmlOptions.roiOffset[2];


		std::cout << orig << std::endl;

		test->SetOrigin(orig);
		typedef itk::ImageFileWriter<ImageType> WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetInput(test);
		writer->SetImageIO(itk::NrrdImageIO::New());
		writer->SetFileName(ss.str());
		writer->Update();
	}


	for(unsigned int i = 0; i < 6; i++)
	{
		std::cout << minMaxs[i] << std::endl;
	}

	
	for(unsigned int i = 0; i < 3; i++)
	{
		std::cout << xmlOptions.roiOffset[i] << std::endl;
	}





	return 0;

}


// ------------------------------------------------------------------------
void getMinMaxBounds(const std::vector<BoundsType> &bounds, double * minMaxs)
{
	// initialise the min maxs 
	minMaxs[0] = 1000;
	minMaxs[1] = -1000;
	minMaxs[2] = 1000;
	minMaxs[3] = -1000;
	minMaxs[4] = 1000;
	minMaxs[5] = -1000;

	for(unsigned int i = 0; i < bounds.size(); i++)
	{
		for(unsigned int j = 0; j < 4; j++)
		{
			ImageType::PointType corner = bounds[i].corners[j];
			double x = corner[0];
			double y = corner[1];
			double z = corner[2];
			
			if(x < minMaxs[0]) minMaxs[0] = x;
			if(x > minMaxs[1]) minMaxs[1] = x;
			if(y < minMaxs[2]) minMaxs[2] = y;
			if(y > minMaxs[3]) minMaxs[3] = y;
			if(z < minMaxs[4]) minMaxs[4] = z;
			if(z > minMaxs[5]) minMaxs[5] = z;
		}	
	}
}

// ------------------------------------------------------------------------
void getImageBounds(const ImageType::Pointer &image, BoundsType &bounds)
{
	// get the image size
	ImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
	
	typedef ImageType::IndexType IndexType;
	IndexType ind1, ind2, ind3, ind4;
	ind1.Fill(0); ind2.Fill(0); ind3.Fill(0); ind4.Fill(0);

	ind2[0] = size[0]-1;
	ind3[1] = size[1]-1;
	ind4[0] = size[0]-1;
	ind4[1] = size[1]-1;

	image->TransformIndexToPhysicalPoint(ind1, bounds.corners[0]);
	image->TransformIndexToPhysicalPoint(ind2, bounds.corners[1]);
	image->TransformIndexToPhysicalPoint(ind3, bounds.corners[2]);
	image->TransformIndexToPhysicalPoint(ind4, bounds.corners[3]);

}

void NormaliseImage(const ImageType::Pointer &input, 
		const ImageType::Pointer &reference,	
		ImageType::Pointer &output)
{
	// flip the x and y axis 
	typedef itk::PermuteAxesImageFilter<ImageType> FlipperType;
	FlipperType::Pointer flipper = FlipperType::New();
	itk::FixedArray<unsigned int, 3> order;
	order[0] = 1;
	order[1] = 0;
	order[2] = 2;

	flipper->SetOrder(order);
	flipper->SetInput(input);
	flipper->Update();

	output = flipper->GetOutput();

	ImageType::DirectionType refDirection = reference->GetDirection();
	ImageType::DirectionType direction = output->GetDirection();
	ImageType::PointType refOrigin = reference->GetOrigin();
	ImageType::PointType origin = output->GetOrigin();

	//origin = origin - refOrigin;

	vnl_vector<double> org;
	org.set_size(3);
	org[0] = origin[0] - refOrigin[0];
	org[1] = origin[1] - refOrigin[1];
	org[2] = origin[2] - refOrigin[2];

	ImageType::SpacingType inSpacing = output->GetSpacing();
	org /= inSpacing[0];

	org =  refDirection.GetTranspose() * org;


	//origin[0] = org[0] * (1.0 / inSpacing[0]);
	//origin[1] = org[1] * (1.0 / inSpacing[0]);
	//origin[2] = org[2] * (1.0 / inSpacing[0]);
	origin[0] = org[0];
	origin[1] = org[1];
	origin[2] = org[2];
	std::cout << org << std::endl;



	ImageType::DirectionType inv;
	inv = refDirection.GetTranspose();

	direction =  inv * direction;

	ImageType::SpacingType spacing;
	spacing.Fill(1.0);

	//std::cout << direction << std::endl;

	output->SetSpacing(spacing);
	output->SetDirection(direction);
	output->SetOrigin(origin);
}
