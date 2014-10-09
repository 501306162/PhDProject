#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <sstream>
#include <vector>
#include <map>


#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include <itkImage.h>

typedef itk::Image<unsigned short, 3> ImageType;
typedef std::vector<ImageType::Pointer> ImageList;

typedef vnl_matrix<double> MatrixType;
typedef vnl_vector<double> VectorType;

typedef itk::Image<float, 3> LevelSetType;


/**
 * Structure defining image bounds
 */
typedef struct _bounds
{
	ImageType::PointType corners[4];
} BoundsType;

/**
 * Structure to hold the options values from the 
 * xml file
 */
typedef struct _options_data
{
	std::vector<int> volumeSize;
	std::vector<int> roiOffset;
	std::vector<int> border;
	std::string dataDirectory;
	std::vector<int> instanceNumbers;

	_options_data()
	{
		volumeSize.resize(3);
		roiOffset.resize(3);
	}

} OptionsData;

// define the structure to hold the input information 
typedef struct _trans
{
	unsigned int dcmSeries;
	unsigned int series;
	unsigned int slice;
	double sliceThickness;
	std::string description;

	std::vector<double> translation;
	std::vector<double> rotation;


	std::vector<std::string> imageFilenames;
	std::vector<ImageType::Pointer> images;
	std::vector<ImageType::Pointer> normalisedImages;

	typedef std::vector<_trans> List;
	typedef std::map<int, _trans> Map;


	std::string toString()
	{
		std::stringstream ss;
		ss << "Description: " << description << std::endl;
		ss << "DCM Series: " << dcmSeries << std::endl;
		ss << "Series: " << series << std::endl;
		ss << "Slice: " << slice << std::endl;
		ss << "Translation: ";
		for(unsigned int i = 0; i < translation.size(); i++)
		{
			ss << translation[i] << " ";
		}
		ss << std::endl;

		ss << "Rotation: ";
		for(unsigned int i = 0; i < rotation.size(); i++)
		{
			ss << rotation[i] << " ";
		}
		ss << std::endl;


		return ss.str();
	}

} SeriesTransform;

#endif
