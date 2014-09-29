#include "SeriesExtractor.h"

#include <gdcmDirectory.h>
#include <QString>


gdcm::Tag SeriesExtractor::seriesDescription = gdcm::Tag(0x0008,0x103e);
gdcm::Tag SeriesExtractor::seriesNumber = gdcm::Tag(0x0020,0x0011);
gdcm::Tag SeriesExtractor::instanceNumber = gdcm::Tag(0x0020,0x0013);
gdcm::Tag SeriesExtractor::sliceThickness = gdcm::Tag(0x0018,0x0050);
gdcm::Tag SeriesExtractor::triggerTime = gdcm::Tag(0x0018,0x1060);
gdcm::Tag SeriesExtractor::numberOfImages = gdcm::Tag(0x0018,0x1090);
gdcm::Tag SeriesExtractor::slicePosition = gdcm::Tag(0x0019,0x1015);
gdcm::Tag SeriesExtractor::studyId = gdcm::Tag(0x0020,0x0010);
gdcm::Tag SeriesExtractor::imagePosition = gdcm::Tag(0x0020,0x0032);
gdcm::Tag SeriesExtractor::imageOrientation = gdcm::Tag(0x0020,0x0037);
gdcm::Tag SeriesExtractor::sliceLocation = gdcm::Tag(0x0020,0x1041);
gdcm::Tag SeriesExtractor::rows = gdcm::Tag(0x0028,0x0010);
gdcm::Tag SeriesExtractor::cols = gdcm::Tag(0x0028,0x0011);
gdcm::Tag SeriesExtractor::pixelSpacing = gdcm::Tag(0x0028,0x0030);

// ------------------------------------------------------------------------
SeriesExtractor::SeriesExtractor()
{
	m_Directory = "";
}


// ------------------------------------------------------------------------
SeriesExtractor::SeriesExtractor(const std::string &directory) : 
	m_Directory(directory)
{
}

// ------------------------------------------------------------------------
void SeriesExtractor::SetDirectory(const std::string &directory)
{
	this->m_Directory = directory;
}

// ------------------------------------------------------------------------
void SeriesExtractor::SetDescriptionFilters(const std::vector<std::string> &filters)
{
	m_DescriptionFilters = filters;
}


// ------------------------------------------------------------------------
bool SeriesExtractor::ExtractSeries()
{
	if(!InitialCheck())
	{
		std::cerr << "The directory has not been set for the scanning" << std::endl;
		return false;
	}

	// get the set of filenames from the directory
	std::vector<std::string> filenames = GetAllFilenames();
	
	// extract all the relevant image tags
	std::vector<DicomImage> images = ScanFiles(filenames);
	
	


	return true;
}



// ------------------------------------------------------------------------
std::vector<DicomImage> SeriesExtractor::ScanFiles(const std::vector<std::string> &filenames)
{
	// set up the scanner
	gdcm::Scanner scanner;
	scanner.AddTag(seriesDescription);
	scanner.AddTag(seriesNumber);
	scanner.AddTag(instanceNumber);
	scanner.AddTag(sliceThickness);
	scanner.AddTag(triggerTime);
	scanner.AddTag(numberOfImages);
	scanner.AddTag(slicePosition);
	scanner.AddTag(studyId);
	scanner.AddTag(imagePosition);
	scanner.AddTag(imageOrientation);
	scanner.AddTag(imagePosition);
	scanner.AddTag(sliceLocation);
	scanner.AddTag(rows);
	scanner.AddTag(cols);
	scanner.AddTag(pixelSpacing);

	scanner.Scan(filenames);

	gdcm::Scanner::MappingType mapping = scanner.GetMappings();

	// assign the scanned values to the filenames to create the dicom images
	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		// is this file a key 
		if(mapping.count(filenames[i].c_str()) == 0) continue;

		// is this file an image
		if(mapping[filenames[i].c_str()].count(instanceNumber) == 0) continue;

		// build the dicom image
		DicomImage image = BuildDicomImage(filenames[i], mapping[filenames[i].c_str()]);
		
		
	}

	std::vector<DicomImage> images;

	return images;
}

// ------------------------------------------------------------------------
DicomImage SeriesExtractor::BuildDicomImage(const std::string &filename, const gdcm::Scanner::TagToValue &mappings)
{
	DicomImage image;
	image.seriesNumber = ConvertInt(mappings.at(seriesNumber));
	image.imageNumber = ConvertInt(mappings.at(instanceNumber));
	image.sliceThickness = ConvertDouble(mappings.at(sliceThickness));

	


	std::cout << image.sliceThickness << std::endl;
	


	return image;

}



/*
typedef struct _dicom_image
{
	unsigned int seriesNumber;
	unsigned int imageNumber;
	unsigned int sliceThickness;
	unsigned int rows;
	unsigned int cols;
	double triggerTime;
	double sliceLocation;
	double slicePosition[3];
	double sliceOrientation[6];

	std::string filename;
} DicomImage;
*/


// ------------------------------------------------------------------------
double SeriesExtractor::ConvertDouble(const std::string &value)
{
	QString qvalue = QString::fromStdString(value);
	return qvalue.toDouble();
}


// ------------------------------------------------------------------------
int SeriesExtractor::ConvertInt(const std::string &value)
{
	QString qvalue = QString::fromStdString(value);
	return qvalue.toInt();
}



// ------------------------------------------------------------------------
bool SeriesExtractor::HasTag(const std::string &filename, const gdcm::Tag &tag, const gdcm::Scanner &scanner)
{
	if(!scanner.IsKey(filename.c_str()))
		return false;

	return true;
}



// ------------------------------------------------------------------------
std::vector<std::string> SeriesExtractor::GetAllFilenames()
{
	gdcm::Directory gdcmDir;
	gdcmDir.Load(m_Directory, true);	
	
	return gdcmDir.GetFilenames();
}



// ------------------------------------------------------------------------
bool SeriesExtractor::InitialCheck()
{
	if(m_Directory.empty())
		return false;

	

	return true;
}
