#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <vector>
#include <sstream>

typedef struct _dicom_image
{
	std::string seriesDescription;
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




#endif
