#include "SeriesExtractor.h"

#include <gdcmDirectory.h>
#include <QString>
#include <QStringList>


gdcm::Tag SeriesExtractor::studyUID = gdcm::Tag(0x0020, 0x000D);
gdcm::Tag SeriesExtractor::seriesUID = gdcm::Tag(0x0020,0x000E);
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
DicomSeriesList SeriesExtractor::GetOutput() const
{
	return m_Series;
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
	
	// group the images into series
	m_Series = GroupSeries(images);

	return true;
}

// ------------------------------------------------------------------------
std::vector<DicomSeries> SeriesExtractor::GroupSeries(const std::vector<DicomImage> &images)
{
	// group the images on the series number
	typedef std::map<std::string, std::vector<DicomImage> > ImageMapType;
	ImageMapType imageMap;
	for(unsigned int i = 0; i < images.size(); i++)
	{
		std::stringstream ss;
		ss << images[i].seriesNumber << ":" << images[i].seriesUID;
		std::string key = ss.str();
		if(imageMap.count(key) == 0)
		{
			std::vector<DicomImage> imageSet;
			imageSet.push_back(images[i]);
			std::pair<std::string, std::vector<DicomImage> > newSet(key, imageSet);
			imageMap.insert(newSet);
		}
		else
		{
			imageMap[key].push_back(images[i]);
		}
	}

	// iterate through the image map and build the series data structure
	std::vector<DicomSeries> seriesList;
	
	ImageMapType::iterator mapIt = imageMap.begin();
	while(mapIt != imageMap.end())
	{
		DicomSeries series;
		series.numImages = mapIt->second.size();
		series.images = mapIt->second;
		series.description = series.images.front().seriesDescription;
		series.number = series.images.front().seriesNumber;

		seriesList.push_back(series);
		
		++mapIt;
	}
	
	return seriesList;
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
	scanner.AddTag(seriesUID);
	scanner.AddTag(studyUID);

	scanner.Scan(filenames);

	gdcm::Scanner::MappingType mapping = scanner.GetMappings();
	std::vector<DicomImage> images;


	

	// assign the scanned values to the filenames to create the dicom images
	bool studyFound = false;
	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		// is this file a key 
		if(mapping.count(filenames[i].c_str()) == 0) continue;

		// is this file an image
		if(mapping[filenames[i].c_str()].count(instanceNumber) == 0) continue;
		if(mapping[filenames[i].c_str()].count(imagePosition) == 0) continue;

		std::string uid = mapping[filenames[i].c_str()][studyUID];
		if(uidMap.count(uid) == 0)
		{
			uidMap[uid] = uidMap.size();
		}

		if(!studyFound)
		{
			m_FirstStudy = mapping[filenames[i].c_str()][studyUID];
			studyFound = true;
		}

		//if(mapping[filenames[i].c_str()][studyUID] != m_FirstStudy) continue;

		// check for the filtering based on series name
		if(FilterOnSeriesDescription(mapping[filenames[i].c_str()][seriesDescription])) continue;
		


		// build the dicom image
		DicomImage image = BuildDicomImage(filenames[i], mapping[filenames[i].c_str()]);

		images.push_back(image);
	}


	return images;
}

// ------------------------------------------------------------------------
bool SeriesExtractor::FilterOnSeriesDescription(const std::string &description)
{
	// if filters are empty then ignore
	if(m_DescriptionFilters.empty()) return false;	

	for(unsigned int i = 0; i < m_DescriptionFilters.size(); i++)
	{
		if(description.find(m_DescriptionFilters[i]) != std::string::npos)
		{
			return false;
		}
	}

	return true;
}


// ------------------------------------------------------------------------
DicomImage SeriesExtractor::BuildDicomImage(const std::string &filename, const gdcm::Scanner::TagToValue &mappings)
{
	DicomImage image;
	image.seriesDescription = mappings.at(seriesDescription);
	image.seriesNumber = ConvertInt(mappings.at(seriesNumber));
	image.imageNumber = ConvertInt(mappings.at(instanceNumber));
	image.seriesUID = mappings.at(seriesUID);
	image.studyUID = QString::number(uidMap[mappings.at(studyUID)]).toStdString();


	if(mappings.count(imagePosition) > 0)
		image.imagePosition = ConvertDoubleArray(mappings.at(imagePosition));



	//image.sliceThickness = ConvertDouble(mappings.at(sliceThickness));
	image.rows = ConvertInt(mappings.at(rows));
	image.cols = ConvertInt(mappings.at(cols));

	

	if(mappings.count(triggerTime) == 0)
		image.triggerTime = 0.0;
	else
		image.triggerTime = ConvertDouble(mappings.at(triggerTime));
	//image.sliceLocation = ConvertDouble(mappings.at(sliceLocation));
	//image.imagePosition = ConvertDoubleArray(mappings.at(imagePosition));
	//image.imageOrientation = ConvertDoubleArray(mappings.at(imageOrientation));
	image.filename = filename;

	return image;
}


// ------------------------------------------------------------------------
std::vector<double> SeriesExtractor::ConvertDoubleArray(const std::string &value)
{
	QString qvalue = QString::fromStdString(value);
	QStringList tokens = qvalue.split("\\");

	std::vector<double> output;
	for(int i = 0; i < tokens.size(); i++)
	{
		output.push_back(tokens[i].toDouble());
	}

	return output;
}

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
