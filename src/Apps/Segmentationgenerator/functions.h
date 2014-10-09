#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <fstream>
#include <gdcmDirectory.h>
#include <itkSimilarity3DTransform.h>
#include "common.h"

/**
 * Normalises an image given the reference image
 */
void NormaliseImage(const ImageType::Pointer &input, 
		const ImageType::Pointer &reference,	
		ImageType::Pointer &output, SeriesTransform &trans);

/**
 * Function to find the bounds defined by the set of images
 */
void getImageBounds(const ImageType::Pointer &image, 
		BoundsType &bounds);

/**
 * Get the min max bounds 
 */
void getMinMaxBounds(const std::vector<BoundsType> &bounds, double * minMaxs);

/**
 * Compute transform matrix from the reference images
 */
void computeTransform(const ImageType::Pointer &image, vnl_matrix<double> &transform);



/**
 * Function to read the xml values from the iresisd options file
 */
void readXMLValues(const std::string &input, OptionsData &options);

void buildShortAxisVolume(const SeriesTransform::Map &transforms, 
		const unsigned int seriesNumber,
		ImageType::Pointer &saVolume);

/**
 * Function to build the output volumes from the inputs
 */
void buildOutputVolumes(SeriesTransform::Map &transforms);

/**
 * Function to get the series number for the short axis
 */
int getShortAxisSeriesNumber(SeriesTransform::Map &transforms);

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
void loadImageSeries(SeriesTransform &series,
		const std::vector<int> &instanceNumbers);

/**
 * Function to apply the transform to the images
 */
void applyTransform(const ImageList &inputImages,
	   	const SeriesTransform &transform,
	   	ImageList &transformedImages);

#endif
