#ifndef VALVE_PLANE_H
#define VALVE_PLANE_H

#include <itkObject.h>
#include <itkPoint.h>
#include <itkVector.h>
#include <itkCovariantVector.h>
#include <itkObjectFactory.h>
#include <QVariant>

namespace vt 
{
// ------------------------------------------------------------------------
class ValvePlane : public itk::Object 
{
public:
	typedef ValvePlane Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(ValvePlane, Object);
	itkNewMacro(Self);

	typedef itk::Point<double, 3> PointType;
	typedef std::vector<PointType> PointListType;
	typedef std::map<std::string, PointListType> PointMapType;
	typedef itk::CovariantVector<double, 3> VectorType;

	void Initialise(const QVariantMap &planeData);
	itkGetMacro(Points, PointMapType);
	itkGetMacro(Normal, VectorType);
	itkGetMacro(Center, PointType);

	PointListType GetPoints(const std::string &key) const { return m_Points.at(key); }
	PointType GetPoint(const std::string &key, const unsigned int index) const { return m_Points.at(key)[index]; }


protected:
	ValvePlane() {}
	virtual ~ValvePlane() {}

private:
	ValvePlane(const Self&);
	void operator=(const Self&);

	PointMapType m_Points;
	VectorType m_Normal;
	PointType m_Center;


};


// ------------------------------------------------------------------------
class ValvePlaneSequence : public itk::Object 
{
public:
	typedef ValvePlaneSequence Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(ValvePlaneSequence, Object);
	itkNewMacro(Self);

	typedef std::vector<ValvePlane::Pointer> PlaneList;

	static Pointer Load(const std::string &filename);

	void LoadFromFile(const std::string &filename);

	ValvePlane::Pointer GetValvePlane(const unsigned int index) { return m_Planes[index]; }

	itkGetMacro(Name, std::string);

	unsigned int GetNumberOfPlanes() const  { return m_Planes.size(); }
protected:
	ValvePlaneSequence() {}
	virtual ~ValvePlaneSequence() {}

private:
	ValvePlaneSequence(const Self&);
	void operator=(const Self&);

	PlaneList m_Planes;
	std::string m_Name;
	

};

}

#endif
