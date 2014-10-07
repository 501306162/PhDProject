#include <iostream>


#include <gdcmDirectory.h>
#include <gdcmScanner.h>
#include <gdcmTag.h>
#include <gdcmSorter.h>
#include <gdcmAttribute.h>
#include <gdcmDataSet.h>
#include <gdcmIPPSorter.h>
#include <QString>

#include <itkGDCMImageIO.h>
#include <itkNrrdImageIO.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImage.h>
#include <itkFlipImageFilter.h>
#include <itkPermuteAxesImageFilter.h>

int main(int, char ** argv)
{

	// parse the dicom directory
	gdcm::Directory dir;
	dir.Load(argv[1], true);
	gdcm::Directory::FilenamesType filenames = dir.GetFilenames();

	
	gdcm::Tag seriesDescription = gdcm::Tag(0x0008,0x103e);
	gdcm::Tag seriesNumber = gdcm::Tag(0x0020,0x0011);
	gdcm::Tag instanceNumber = gdcm::Tag(0x0020,0x0013);
	gdcm::Tag sliceThickness = gdcm::Tag(0x0018,0x0050);
	gdcm::Tag triggerTime = gdcm::Tag(0x0018,0x1060);
	gdcm::Tag numberOfImages = gdcm::Tag(0x0018,0x1090);
	gdcm::Tag slicePosition = gdcm::Tag(0x0019,0x1015);
	gdcm::Tag imagePosition = gdcm::Tag(0x0020,0x0032);
	gdcm::Tag imageOrientation = gdcm::Tag(0x0020,0x0037);
	gdcm::Tag sliceLocation = gdcm::Tag(0x0020,0x1041);
	gdcm::Tag rows = gdcm::Tag(0x0028,0x0010);
	gdcm::Tag cols = gdcm::Tag(0x0028,0x0011);
	gdcm::Tag pixelSpacing = gdcm::Tag(0x0028,0x0030);

	gdcm::Scanner scanner;
	scanner.AddTag(seriesDescription);
	scanner.AddTag(seriesNumber);
	scanner.AddTag(instanceNumber);
	scanner.AddTag(sliceThickness);
	scanner.AddTag(triggerTime);
	scanner.AddTag(numberOfImages);
	scanner.AddTag(slicePosition);
	scanner.AddTag(imagePosition);
	scanner.AddTag(imageOrientation);
	scanner.AddTag(sliceLocation);
	scanner.AddTag(rows);
	scanner.AddTag(cols);
	scanner.AddTag(pixelSpacing);

	scanner.Scan(filenames);
	gdcm::Scanner::MappingType mapping = scanner.GetMappings();


	// extract all the images that are Short axis and are the first instance
	gdcm::Directory::FilenamesType targetFilenames;
	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		const char * fname = filenames[i].c_str();
		if(mapping.count(fname) == 0) continue;


		// extract the image information
		std::string descriptionStr = mapping[fname][seriesDescription];
		std::string instanceNumberStr = mapping[fname][instanceNumber];
		QString description = QString::fromStdString(descriptionStr);
		unsigned int instanceNumber = QString::fromStdString(instanceNumberStr).toInt();

		// check that the description is a short axis one and is instance 1
		if(instanceNumber == 1 && description.contains("sa", Qt::CaseInsensitive))
		{
			targetFilenames.push_back(filenames[i]);
		}
	}


	// sort the images based on their slice position
	/*
	gdcm::Sorter sorter;
	sorter.SetSortFunction(position_sort);
	sorter.StableSort(targetFilenames);

	gdcm::Directory::FilenamesType sorted = sorter.GetFilenames();
	for(unsigned int i = 0; i < sorted.size(); i++)
	{
		const char * fname = sorted[i].c_str();
		std::string position = mapping[fname][imagePosition];
		std::cout << position << std::endl;
		std::cout << mapping[fname][sliceLocation] << std::endl;
	}
	*/

	std::cout << targetFilenames.size() << std::endl;


	// find the slice with the smallest slice position
	double smallest = 1000000.0;
	int smallestIndex;
	for(unsigned int i = 0; i < targetFilenames.size(); i++)
	{
		const char * fname = targetFilenames[i].c_str();
		std::string slicePosition = mapping[fname][sliceLocation];
		double pos = QString::fromStdString(slicePosition).toDouble();

		std::cout << pos << std::endl;
		if(pos < smallest)
		{	
			smallest = pos;
			smallestIndex = i;
		}
	}

	// load the image
	typedef itk::Image<unsigned short, 3> ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(targetFilenames[smallestIndex]);
	reader->SetImageIO(itk::GDCMImageIO::New());
	reader->Update();

	
	// flip the x and y axis
	typedef itk::PermuteAxesImageFilter<ImageType> FlipperType;
	FlipperType::Pointer flipper = FlipperType::New();
	itk::FixedArray<unsigned int, 3> order;
	order[0] = 1;
	order[1] = 0;
	order[2] = 2;
	flipper->SetOrder(order);
	flipper->SetInput(reader->GetOutput());
	flipper->Update();


	ImageType::Pointer referenceImage = flipper->GetOutput();
	ImageType::DirectionType direction = referenceImage->GetDirection();
	direction.SetIdentity();
	referenceImage->SetDirection(direction);
	
	ImageType::PointType origin = referenceImage->GetOrigin();
	origin.Fill(20.0);
	referenceImage->SetOrigin(origin);

	ImageType::SpacingType spacing;
	spacing.Fill(1.0);

	referenceImage->SetSpacing(spacing);


	
	flipper->GetOutput()->Print(std::cout);

	typedef itk::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(referenceImage);
	writer->SetImageIO(itk::NrrdImageIO::New());
	writer->SetFileName("test.nrrd");
	writer->Update();

	std::cout << targetFilenames[smallestIndex] << std::endl;

	return 0;
}
