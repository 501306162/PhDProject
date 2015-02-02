#ifndef POINT_LOCATIONS_H
#define POINT_LOCATIONS_H

#include <itkObject.h>
#include <itkObjectFactory.h>

#include <MatrixCommon.h>

namespace vt
{
class PointLocations : public itk::Object 
{
public:
	typedef PointLocations Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	typedef utils::DoubleMatrixType MatrixType;
	typedef std::map<std::string, MatrixType> DataMapType;


	itkTypeMacro(PointLocations, Object);
	itkNewMacro(Self);

	static PointLocations::Pointer Load(const std::string &data);
	void LoadData(const std::string &data);
	void GetPoints(const int exclude, const std::string type, unsigned int timeStep, MatrixType &p1, MatrixType &p2);

protected:
	PointLocations() {}
	virtual ~PointLocations() {}

private:
	PointLocations(const Self&);
	void operator=(const Self&);

	DataMapType m_Data;
};



}

#endif
