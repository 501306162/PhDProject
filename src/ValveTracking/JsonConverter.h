#ifndef JSON_CONVERTER_H
#define JSON_CONVERTER_H

#include <iostream>
#include <itkObject.h>
#include <itkObjectFactory.h>

#include "JsonReader.h"
#include "ValveLine.h"


#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <itkImage.h>

namespace vt
{
class JsonConverter : public itk::Object
{
public:
	typedef JsonConverter Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(JsonConverter, Object);
	itkNewMacro(Self);

	typedef itk::Image<unsigned short, 4> ReferenceType;
	typedef itk::Image<unsigned short, 3> ImageType;
	typedef ImageType::PointType PointType;

	typedef JsonReader::LineSequence LineSequence;
	typedef JsonReader::Line Line;


	void SetInput(const LineSequence & seq) { m_Input = seq; }
	itkGetMacro(Input, LineSequence);

	ValveSequence<3>::Pointer GetOutput() const { return m_Output; }

	void Convert();



protected:

	JsonConverter() {};
	virtual ~JsonConverter() {}

	void ConvertLine(const Line & line, ImageType::Pointer &image, ValveLine<3>::Pointer & output);
	void LoadImages(const std::string &filename, ReferenceType::Pointer &ref, vtkImageData* image);
	void ConvertPoint(ImageType::Pointer &image, const double * input, double * output);
	void FinishUp(const unsigned int &timeStep, ValveLine<3>::Pointer &valve);
private:	

	LineSequence m_Input;

	ReferenceType::Pointer m_Reference;
	vtkSmartPointer<vtkImageData> m_VtkImage;

	ValveSequence<3>::Pointer m_Output;

	JsonConverter(const Self&);
	void operator=(const Self&);

	

};



}


#endif
