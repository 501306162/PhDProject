#ifndef LINE_INTERSECTION_FINDER_H
#define LINE_INTERSECTION_FINDER_H

#include <itkImage.h>
#include <vtkBoundingBox.h>
#include <Eigen/Dense>
#include <Eigen/StdVector>

namespace vt
{
class LineIntersectionFinder : public itk::Object
{
public:
	typedef LineIntersectionFinder Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	typedef itk::Image<unsigned short, 3> ImageType;

	typedef Eigen::Vector3d VectorType;
	typedef std::vector<VectorType, Eigen::aligned_allocator<VectorType> > VectorListType;
	typedef Eigen::Hyperplane<double, 3> PlaneType;
	typedef Eigen::ParametrizedLine<double, 3> LineType;

	itkTypeMacro(ValveIntersectionFinder, Object);
	itkNewMacro(Self);

	typedef struct line_
	{
		ImageType::PointType p1, p2;
		VectorType direction;
	} OutputLineType;

	void SetImage(const ImageType::Pointer &image) { m_Image = image; }
	void SetBoundingBox(const vtkBoundingBox &box) { m_BoundingBox = box; }
	void SetPlane(const VectorType &normal, const VectorType &point) { m_PlaneNormal = normal; m_PlanePoint = point; }
	void Compute();
	OutputLineType GetOutput() const { return m_Output; }

protected:
	LineIntersectionFinder() {}
	virtual ~LineIntersectionFinder() {}

private:

	void FindIntersectingLine(const PlaneType &p1, const PlaneType &p2, LineType &line);
	void GetBoundingPlaneIntersections(const LineType &line, vtkBoundingBox &box, double &minT, double &maxT);

	LineIntersectionFinder(const Self&);
	void operator=(const Self&);

	ImageType::Pointer m_Image;
	vtkBoundingBox m_BoundingBox;
	VectorType m_PlaneNormal;
	VectorType m_PlanePoint;

	OutputLineType m_Output;

};



}

#endif
