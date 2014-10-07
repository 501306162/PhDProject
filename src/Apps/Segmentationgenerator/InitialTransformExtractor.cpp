#include "InitialTransformExtractor.h"

#include <QString>
#include <QStringList>

#include <itkImageFileReader.h>
#include <itkGDCMImageIO.h>
#include <itkNrrdImageIO.h>
#include <itkImageFileWriter.h>

#include <itkPermuteAxesImageFilter.h>

gdcm::Tag InitialTransformExtractor::seriesDescription = gdcm::Tag(0x0008,0x103e);
gdcm::Tag InitialTransformExtractor::seriesNumber = gdcm::Tag(0x0020,0x0011);
gdcm::Tag InitialTransformExtractor::instanceNumber = gdcm::Tag(0x0020,0x0013);
gdcm::Tag InitialTransformExtractor::sliceThickness = gdcm::Tag(0x0018,0x0050);
gdcm::Tag InitialTransformExtractor::triggerTime = gdcm::Tag(0x0018,0x1060);
gdcm::Tag InitialTransformExtractor::numberOfImages = gdcm::Tag(0x0018,0x1090);
gdcm::Tag InitialTransformExtractor::slicePosition = gdcm::Tag(0x0019,0x1015);
gdcm::Tag InitialTransformExtractor::studyId = gdcm::Tag(0x0020,0x0010);
gdcm::Tag InitialTransformExtractor::imagePosition = gdcm::Tag(0x0020,0x0032);
gdcm::Tag InitialTransformExtractor::imageOrientation = gdcm::Tag(0x0020,0x0037);
gdcm::Tag InitialTransformExtractor::sliceLocation = gdcm::Tag(0x0020,0x1041);
gdcm::Tag InitialTransformExtractor::rows = gdcm::Tag(0x0028,0x0010);
gdcm::Tag InitialTransformExtractor::cols = gdcm::Tag(0x0028,0x0011);
gdcm::Tag InitialTransformExtractor::pixelSpacing = gdcm::Tag(0x0028,0x0030);


// ------------------------------------------------------------------------
InitialTransformExtractor::InitialTransformExtractor()
{
	this->m_DicomDir = "";
}

// ------------------------------------------------------------------------
void InitialTransformExtractor::SetOptions(const OptionsData &options)
{
	this->m_Options = options;
}

// ------------------------------------------------------------------------
void InitialTransformExtractor::SetDicomDir(const std::string &dir)
{
	this->m_DicomDir = dir;
}


// ------------------------------------------------------------------------
void InitialTransformExtractor::Compute()
{
	GetFilenames(m_Filenames);
	BuildMapping(m_Mapping);
	std::string referenceFilename = GetReferenceImageFilename();


	// load the reference file
	typedef itk::ImageFileReader<ImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(referenceFilename);
	reader->SetImageIO(itk::GDCMImageIO::New());
	reader->Update();
	ImageType::Pointer refImage = reader->GetOutput();


	// flip the x and y axis 
	typedef itk::PermuteAxesImageFilter<ImageType> FlipperType;
	FlipperType::Pointer flipper = FlipperType::New();
	itk::FixedArray<unsigned int, 3> order;
	order[0] = 1;
	order[1] = 0;
	order[2] = 2;

	flipper->SetOrder(order);
	flipper->SetInput(refImage);
	flipper->Update();

	refImage = flipper->GetOutput();
	
	// transformation is the negative origin of the reference
	ImageType::PointType origin = refImage->GetOrigin();
	m_Translation[0] = -origin[0];
	m_Translation[1] = -origin[1];
	m_Translation[2] = -origin[2];

	m_Rotation = refImage->GetDirection().GetInverse();
	m_Reference = refImage;
}

// ------------------------------------------------------------------------
InitialTransformExtractor::TranslationType InitialTransformExtractor::ComputeTranslation(
		const ImageType::Pointer &image,
	   	const OptionsData &options)
{
	// set the offsets
	TranslationType translation;
	translation.Fill(0);
	translation[0] -= options.roiOffset[0];
	translation[1] -= options.roiOffset[1];
	translation[2] -= options.roiOffset[2];


	ImageType::PointType origin = image->GetOrigin();
	translation[0] -= origin[0];
	translation[1] -= origin[1];
	translation[2] -= origin[2];

	return translation;
}


// ------------------------------------------------------------------------
std::string InitialTransformExtractor::GetReferenceImageFilename()
{
	FilenamesType targetFilenames;
	for(unsigned int i = 0; i < m_Filenames.size(); i++)
	{
		const char * fname = m_Filenames[i].c_str();
		if(m_Mapping.count(fname) == 0) continue;

		// extract the image information
		std::string descriptionStr = m_Mapping[fname][seriesDescription];
		std::string instanceNumberStr = m_Mapping[fname][instanceNumber];
		QString description = QString::fromStdString(descriptionStr);
		unsigned int instanceNumber = QString::fromStdString(instanceNumberStr).toInt();

		// check that the description is a short axis one and is instance 1
		if(instanceNumber == 1 && description.contains("sa", Qt::CaseInsensitive))
		{
			targetFilenames.push_back(m_Filenames[i]);
		}
	}

	// find the slice with the smallest slice position
	double smallest = 1000000.0;
	int smallestIndex;
	for(unsigned int i = 0; i < targetFilenames.size(); i++)
	{
		const char * fname = targetFilenames[i].c_str();
		std::string slicePosition = m_Mapping[fname][sliceLocation];
		double pos = QString::fromStdString(slicePosition).toDouble();

		if(pos < smallest)
		{	
			smallest = pos;
			smallestIndex = i;
		}
	}

	return targetFilenames[smallestIndex];
}

// ------------------------------------------------------------------------
void InitialTransformExtractor::GetFilenames(FilenamesType &filenames)
{
	m_Dir.Load(m_DicomDir, true);
	filenames = m_Dir.GetFilenames();
}

// ------------------------------------------------------------------------
void InitialTransformExtractor::BuildMapping(gdcm::Scanner::MappingType &mapping)
{
	// set up the scanner
	m_Scanner.AddTag(seriesDescription);
	m_Scanner.AddTag(instanceNumber);
	m_Scanner.AddTag(numberOfImages);
	m_Scanner.AddTag(sliceLocation);
	m_Scanner.Scan(m_Filenames);



	mapping = m_Scanner.GetMappings();
}
