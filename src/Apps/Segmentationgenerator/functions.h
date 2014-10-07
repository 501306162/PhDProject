#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <fstream>
#include <gdcmDirectory.h>


#include "common.h"


/**
 * Function to read the xml values from the iresisd options file
 */
void readXMLValues(const std::string &input, OptionsData &options);


/**
 * Function to load the input options file
 */
void readOptionsFile(const std::string &filename, SeriesTransform::Map &transforms);

/**
 * Function to scan the dicom files and extract them
 */
void scanDicomFiles(const std::string &folder, gdcm::Directory::FilenamesType &filenames);

/**
 * Funciton to extract the filenames and put them into the transform objects. It will also, sort 
 * the series images by their trigger time
 */
void groupDicomFiles(const gdcm::Directory::FilenamesType &filenames,
	   	SeriesTransform::Map &transforms);


/**
 * Load the set of images for each of tge registration outputs
 */
void loadImageSeries(const gdcm::Directory::FilenamesType &filenames,
		ImageList &images);

/**
 * Function to apply the transform to the images
 */
void applyTransform(const ImageList &inputImages,
	   	const SeriesTransform &transform,
	   	ImageList &transformedImages);

#endif
