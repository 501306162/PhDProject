#ifndef SERIES_EXTRACTOR_H
#define SERIES_EXTRACTOR_H


#include <gdcmScanner.h>
#include <gdcmTag.h>


typedef struct _dicom_image
{
	std::string seriesDescription;
	std::string studyUID;
	unsigned int seriesNumber;
	unsigned int imageNumber;
	double sliceThickness;
	unsigned int rows;
	unsigned int cols;
	double triggerTime;
	double sliceLocation;
	std::vector<double> imagePosition;
	std::vector<double> imageOrientation;
	std::string filename;
	std::string seriesUID;

	std::string toString()
	{
		std::stringstream ss;
		ss << "Series Number: " << seriesNumber << "\n";
		ss << "Image Number: " << imageNumber << "\n";
		ss << "Slice Thickness: " << sliceThickness << "\n";
		ss << "Dimensions: " << rows << "-" << cols << "\n";
		ss << "Trigger Time: " << triggerTime << "\n";
		ss << "Slice location: " << sliceLocation << "\n";
		ss << "Filename: " << filename << "\n";

		return ss.str();
	}


	bool operator < (const _dicom_image &im) const
	{
		return triggerTime < im.triggerTime;
	}

} DicomImage;


/**
 * Definition of a dicom series
 */
typedef struct _dicom_series
{
	std::string description;
	std::string number;
	unsigned int numImages;
	std::vector<DicomImage> images;

} DicomSeries;
typedef std::vector<DicomSeries> DicomSeriesList;





/**
 * Class that scans the dicom directory and extracts the set of s
 * series from it
 */
class SeriesExtractor
{
public:

	/**
	 * Constructor, takes the directory as an argument
	 */
	SeriesExtractor(const std::string &directory);
	SeriesExtractor();

	/**
	 * Set the directory that will be scanned
	 */
	void SetDirectory(const std::string &directory);

	/**
	 * Set a filter list, this will restrict the series descriptions
	 * to having one of these search strings within the name
	 */
	void SetDescriptionFilters(const std::vector<std::string> &filters);

	/**
	 * Get all the output of the series extraction
	 */
	DicomSeriesList GetOutput() const;

	/**
	 * Get a single output given an index
	 */
	DicomSeries GetOutput(unsigned int index) const;

	/**
	 * Function that extracts the series
	 */
	bool ExtractSeries();

private:
	
	/**
	 * Function to check if things are ok
	 */
	bool InitialCheck();


	/**
	 * Retrieves all filenames from the scanning
	 */
	std::vector<std::string> GetAllFilenames();


	/**
	 * Function to scan all the filenames to extract the relevant information
	 */
	std::vector<DicomImage> ScanFiles(const std::vector<std::string> &filenames);


	/**
	 * Function to check if an image has a particular mapping
	 */
	bool HasTag(const std::string &filename, const gdcm::Tag &tag, const gdcm::Scanner &scanner);

	/**
	 * Function to build the dicom image given the dicom tags
	 */
	DicomImage BuildDicomImage(const std::string &filename, const gdcm::Scanner::TagToValue &mappings);

	/**
	 * Function to group the extracted images by series
	 */
	std::vector<DicomSeries> GroupSeries(const std::vector<DicomImage> &images);


	/**
	 * Function to check if the image satisfies the filter values
	 */
	bool FilterOnSeriesDescription(const std::string &description);


	/**
	 * Conversion functions to turn the strings to their appropriate values
	 */
	int ConvertInt(const std::string &value);
	double ConvertDouble(const std::string &value);
	std::vector<double> ConvertDoubleArray(const std::string &value);


	std::string m_Directory;
	std::vector<std::string> m_DescriptionFilters;

	DicomSeriesList m_Series;
	std::map<std::string, int> uidMap;

	static gdcm::Tag seriesDescription;
	static gdcm::Tag seriesNumber;
	static gdcm::Tag instanceNumber;
	static gdcm::Tag sliceThickness;
	static gdcm::Tag triggerTime;
	static gdcm::Tag numberOfImages;
	static gdcm::Tag slicePosition;
	static gdcm::Tag studyId;
	static gdcm::Tag imagePosition;
	static gdcm::Tag imageOrientation;
	static gdcm::Tag sliceLocation;
	static gdcm::Tag rows;
	static gdcm::Tag cols;
	static gdcm::Tag pixelSpacing;
	static gdcm::Tag studyUID;
	static gdcm::Tag seriesUID;

	std::string m_FirstStudy;
};

#endif
