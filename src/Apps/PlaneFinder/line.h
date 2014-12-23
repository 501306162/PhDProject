#ifndef LINE_H
#define LINE_H

#include <map>

#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>

class Line 
{
public:
	typedef std::vector<Line*> List;
	

	enum Type { MV, TP, AV, PV };
	
	typedef std::map<Type, Line*> Map;

	Line();

	 static Line * NewLine(vtkImageData * image, Type type);
	 static Line * NewLine(Type type, double * p1, double * p2);


	 Type getType() { return this->type; }

	 Line * copy();

	 bool isLocked() { return this->locked; }
	 void setLocked(bool val) { this->locked = val; }
	 void setPoint1(double x, double y, double z);
	 void setPoint2(double x, double y, double z);



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

	bool locked;

};





#endif
