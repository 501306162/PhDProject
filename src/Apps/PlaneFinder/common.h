#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <vector>

#include <vtkImageData.h>
#include <vtkSmartPointer.h>


typedef vtkSmartPointer<vtkImageData> ImagePointer;
typedef std::vector<ImagePointer> ImageList;


typedef struct image_data_
{
	ImageList images;
	std::string filename;
} ImageData;

typedef std::vector<ImageData> ImageDataList;


#endif

