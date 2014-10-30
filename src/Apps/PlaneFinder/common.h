#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <vector>

#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <itkImage.h>

#include "line.h"
#include <vtkPolyData.h>
#include <vtkImageSlice.h>
#include <vtkImageSliceMapper.h>



typedef vtkSmartPointer<vtkImageData> ImagePointer;
typedef std::vector<Line::List> LinesPerVolume;
typedef std::vector<LinesPerVolume> LinesPerSequence;

typedef itk::Image<unsigned short, 3> ImageType;



typedef struct data_holder_
{
	ImagePointer vtkImage;
	ImageType::Pointer itkImage;
	vtkSmartPointer<vtkImageSlice> actor;
	vtkSmartPointer<vtkImageSliceMapper> mapper;
	
	unsigned int slice;
	unsigned int timestep;
	Line::Map lines;

} DataHolder;

typedef std::vector<DataHolder> VolumeData;
typedef std::vector<VolumeData> SequenceData;

typedef struct data_instance_
{
	std::string filename;
	SequenceData images;
} DataInstance;

typedef std::vector<DataInstance> DataList;





#endif

