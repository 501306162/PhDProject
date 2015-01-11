#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include <itkObject.h>
#include <itkObjectFactory.h>

#include <itkPoint.h>
#include <itkBoundingBox.h>
#include <itkVectorContainer.h>

#include <itkSimilarity3DTransform.h>
#include <itkImage.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkBoundingBox.h>
#include <vtkPoints.h>


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
	
	typedef itk::Image<unsigned short, 3> ImageType;
	typedef ImageType::IndexType IndexType;
	typedef itk::Image<unsigned char, 3> MaskType;

	void Load(const std::string &filename, const int exclude=-1);

	typedef itk::Similarity3DTransform<double> TransformType;

	void TransformBoundingBox(const TransformType::Pointer &transform);
	void ApplyInflation(BoundingBoxType::Pointer &box, double &inflation);
	void ComputeImageMask(const ImageType::Pointer &image, unsigned int point, MaskType::Pointer &mask);
	void SetInfation(const double &inflation) { m_Inflation = inflation; }
	bool IsInside(const PointType &point, vtkSmartPointer<vtkPolyData> &box);


protected:
	BoundingBox() : m_Inflation(0.5) {}
	virtual ~BoundingBox() {}

private:
	BoundingBox(const Self&);
	void operator=(const Self&);

	void ExtractMaskMinMaxIndex(const ImageType::Pointer &image, 
			const BoundingBoxType::Pointer &box,
			ImageType::RegionType &region);

	void SetBoundingBox(vtkSmartPointer<vtkPoints> &points, vtkSmartPointer<vtkPolyData> &box);
	void ApplyTransform(vtkSmartPointer<vtkPolyData> &input, 
			const TransformType::Pointer &transform, 
			vtkSmartPointer<vtkPolyData> &output);

	vtkSmartPointer<vtkPolyData> m_BoundingBoxP1;
	vtkSmartPointer<vtkPolyData> m_BoundingBoxP2;

	double m_Inflation;

};


}

#endif
