#ifndef VALVE_ORIGIN_FINDER_H
#define VALVE_ORIGIN_FINDER_H

#include <itkImage.h>

namespace vt
{
class ValveOriginFinder : public itk::Object 
{
public:
	typedef ValveOriginFinder Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(ValveOriginFinder, Object);
	itkNewMacro(Self);

	typedef itk::Vector<double, 3> VectorType;
	typedef struct plane_
	{
		VectorType normal;		
		double d;
	} Plane;

	typedef unsigned short PixelType;
	typedef itk::Image<PixelType, 3> ImageType;
	typedef ImageType::PointType PointType;
	typedef std::vector<PointType> PointListType;
	typedef itk::Matrix<double, 3, 3> RotationType;

	void Set2CImage(const ImageType::Pointer &image) { m_2CImage = image; }
	void Set3CImage(const ImageType::Pointer &image) { m_3CImage = image; }
	void SetImageStack(const ImageType::Pointer &image) { m_Stack = image; }

	void Compute();
	
	void ComputePlane(const ImageType::Pointer &image, Plane &plane);
	void ExtractPlanePoints(const ImageType::Pointer &image, PointListType &points);

	void ExtractSASlice(const ImageType::Pointer &image, ImageType::Pointer &slice);
	
	PointType GetOrigin() const { return m_Origin; }
	VectorType GetAxis() const { return m_YAxis; }
	RotationType GetRotation() const { return m_Rotation; }

protected:
	ValveOriginFinder() {}
	virtual ~ValveOriginFinder() {}

private:
	ValveOriginFinder(const Self&);
	void operator=(const Self&);

	ImageType::Pointer m_2CImage;
	ImageType::Pointer m_3CImage;
	ImageType::Pointer m_Stack;

	PointType m_Origin;
	VectorType m_YAxis;

	RotationType m_Rotation;

};
}

#endif
