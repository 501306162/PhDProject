#ifndef LINE_H
#define LINE_H


#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>

class Line 
{
public:
	typedef std::vector<Line*> List;

	enum Type { MV, TP, AV };
	
	Line();

	 static Line * NewLine(vtkImageData * image, Type type);

	 Type getType() { return this->type; }



	 void getColour(double * col);
	 vtkSmartPointer<vtkPolyData> getPoly() { return this->poly; }

	 static Type getTypeEnum(unsigned int btnIndex);
	 static std::string getTypeString(Type type);
	 vtkActor * getActor() { return actor; }

private:

	Type type;
	double x;
   	double y; 
	double z;
	double x2;
   	double y2; 
	double z2;

	vtkSmartPointer<vtkPolyData> poly;
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkPolyDataMapper> mapper;


};





#endif
