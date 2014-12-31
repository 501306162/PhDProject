#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include <itkObject.h>
#include <itkObjectFactory.h>

#include <itkPoint.h>
#include <itkBoundingBox.h>
#include <itkVectorContainer.h>

#include <itkTransform.h>


namespace vt
{
class BoundingBox : public itk::Object 
{
public:
	typedef BoundingBox Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkNewMacro(Self);
	itkTypeMacro(BoundingBox, Object);

	typedef itk::Point<double, 3> PointType;
	typedef itk::VectorContainer<int, PointType> PointContainerType;
	typedef itk::BoundingBox<int, 3, double, PointContainerType> BoundingBoxType;

	void Load(const std::string &filename, const int exclude=-1);

	typedef itk::Transform<double, 3, 3> TransformType;

	void TransformBoundingBox(const TransformType::Pointer &transform);

protected:
	BoundingBox() {}
	virtual ~BoundingBox() {}

private:
	BoundingBox(const Self&);
	void operator=(const Self&);

	void SetBoundingBox(const PointContainerType::Pointer &points, BoundingBoxType::Pointer &box);

	BoundingBoxType::Pointer m_BoundingBoxP1;
	BoundingBoxType::Pointer m_BoundingBoxP2;


};


}

#endif
