#ifndef RESULT_VIEWER_H
#define RESULT_VIEWER_H

#include <itkObject.h>
#include <itkObjectFactory.h>
#include <itkImage.h>
#include <Eigen/Dense>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkVolume.h>
#include <ValveLine.h>
#include <vtkBoundingBox.h>
#include <TestData.h>

namespace vt
{
class ResultViewer : public itk::Object 
{
public:
	typedef ResultViewer Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkTypeMacro(ResultViewer, Object);
	itkNewMacro(Self);

	typedef itk::Image<unsigned short, 3> ImageType;
	typedef Eigen::Vector3d VectorType;
	typedef ImageType::PointType PointType;

	void SetImages(const std::vector<ImageType::Pointer> &images) { m_Images = images; }
	void SetPlane(VectorType &point, VectorType &normal) { m_Normal = normal; m_Point = point; }
	void SetStartPlane(VectorType &point, VectorType &normal) { m_StartNormal = normal; m_StartPoint = point; }
	void SetBoundingBox(vtkBoundingBox &box) { m_BoundingBox = box; }
	void SetTransform(const TestData::TransformType::Pointer &trans) { m_Transform = trans; }
	void AddPoints(const std::vector<std::vector<std::pair<PointType, PointType > > > &points) { m_Points = points; }
	void SetUp();
	void View();
	void SetWeights(const std::vector<std::vector<double> > &weights) { m_Weights = weights; }
	void Save(const std::string &filename);
	void SetGT(const std::vector<std::pair<PointType, PointType> > &gt) { m_GT = gt; }
	void SetRes(const std::vector<std::pair<PointType, PointType> > &res) { m_Res = res; }
protected:
	ResultViewer() {}
	virtual ~ResultViewer() {}

private:
	ResultViewer(const Self&);
	void operator=(const Self&);

	void GetPlane(VectorType &point, VectorType &normal, ImageType::Pointer &image, vtkSmartPointer<vtkPolyData> &poly);
	void GetImage(const ImageType::Pointer &in, vtkSmartPointer<vtkImageData> &image);
	void GetLine(const ImageType::Pointer &image, VectorType &point, VectorType &normal, 
			vtkSmartPointer<vtkImageData> &vtkImage,
		   	vtkSmartPointer<vtkPolyData> &line);

	TestData::TransformType::Pointer m_Transform;
	std::vector<ImageType::Pointer> m_Images;
	VectorType m_Point;
	VectorType m_Normal;
	VectorType m_StartPoint;
	VectorType m_StartNormal;
	std::vector<std::vector<std::pair<PointType, PointType > > > m_Points;
	std::vector<std::vector<double> > m_Weights;
	std::vector<std::pair<PointType, PointType> > m_GT;
	std::vector<std::pair<PointType, PointType> > m_Res;

	vtkBoundingBox m_BoundingBox;
	vtkSmartPointer<vtkRenderWindow> m_RenderWindow;
	vtkSmartPointer<vtkRenderWindowInteractor> m_Interactor;

};


}

#endif
