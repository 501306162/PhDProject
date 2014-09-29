#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <vector>


typedef struct _dicom_image
{
	unsigned int seriesNumber;
	unsigned int imageNumber;
	double sliceThickness;
	unsigned int rows;
	unsigned int cols;
	double triggerTime;
	double sliceLocation;
	double slicePosition[3];
	double sliceOrientation[6];

	std::string filename;
} DicomImage;


/**
 * Definition of a dicom series
 */
typedef struct _dicom_series
{
	std::string description;
	std::string number;
	unsigned int numImages;
	std::vector<DicomImage> imageFilenames;

} DicomSeries;
typedef std::vector<DicomSeries> DicomSeriesList;




#endif
