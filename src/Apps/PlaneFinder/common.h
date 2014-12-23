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


#include <QStringList>


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

	enum ImageType { UNKNOWN, T2C, T3C, T4C, TR2C, TR3C, TRVOT, TLVOT };

	static QStringList getTypes()
	{
		QStringList types;
		types << "2C";
		types << "3C";
		types << "4C";
		types << "R2C";
		types << "R3C";
		types << "RVOT";
		types << "LVOT";

		return types;
	}

	static ImageType getType(const std::string &name)
	{
		if(name == "2C")
			return T2C;
		else if (name == "3C")
			return T3C;
		else if (name  == "4C")
			return T4C;
		else if (name == "R3C")
			return TR3C;
		else if (name == "R2C")
			return TR2C;
		else if (name == "RVOT")
			return TRVOT;
		else if (name == "LVOT")
			return TLVOT;
		else
			return UNKNOWN;
	}

	ImageType imageType;

	std::string getImageTypeString()
	{
		switch(this->imageType)
		{
			case T2C:
				return "2C";
				break;
			case T3C:
				return "3C";
				break;
			case T4C:
				return "4C";
				break;
			case TR2C:
				return "R2C";
				break;
			case TR3C:
				return "R3C";
				break;
			case TLVOT:
				return "LVOT";
				break;
			case TRVOT:
				return "RVOT";
				break;
			default:
				return "Unknown";
				break;
		}
	}


	std::string filename;
	SequenceData images;

	data_instance_() : imageType(UNKNOWN) {}
} DataInstance;

typedef std::vector<DataInstance> DataList;





#endif

