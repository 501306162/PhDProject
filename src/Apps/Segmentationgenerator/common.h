#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <sstream>
#include <vector>
#include <map>

//#include <CommonDefinitions.h>


#include <itkImage.h>

typedef itk::Image<unsigned short, 3> ImageType;
typedef std::vector<ImageType::Pointer> ImageList;


// define the structure to hold the input information 
typedef struct _trans
{
	unsigned int dcmSeries;
	unsigned int series;
	unsigned int slice;

	std::vector<double> translation;
	std::vector<double> rotation;


	std::vector<std::string> imageFilenames;

	typedef std::vector<_trans> List;
	typedef std::map<int, _trans> Map;


	std::string toString()
	{
		std::stringstream ss;
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
