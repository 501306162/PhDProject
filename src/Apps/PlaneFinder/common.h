#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <vector>

#include <vtkImageData.h>
#include <vtkSmartPointer.h>


typedef vtkSmartPointer<vtkImageData> ImagePointer;
typedef std::vector<ImagePointer> ImageVolume;
typedef std::vector<ImageVolume> ImageSequence;


typedef struct image_data_
{
	ImageSequence images;
	std::string filename;
} ImageData;

typedef std::vector<ImageData> ImageDataList;


#endif

