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

	void SetImages(const std::vector<ImageType::Pointer> &images) { m_Images = images; }
	void SetPlane(VectorType &point, VectorType &normal) { m_Normal = normal; m_Point = point; }
	void SetStartPlane(VectorType &point, VectorType &normal) { m_StartNormal = normal; m_StartPoint = point; }
	void SetBoundingBox(vtkBoundingBox &box) { m_BoundingBox = box; }
	void SetUp();
	void View();
	void Save(const std::string &filename);
protected:
	ResultViewer() {}
	virtual ~ResultViewer() {}

private:
	ResultViewer(const Self&);
	void operator=(const Self&);

	void GetPlane(vtkSmartPointer<vtkPolyData> &poly);
	void GetImage(const ImageType::Pointer &in, vtkSmartPointer<vtkImageData> &image);
	void GetLine(const ImageType::Pointer &image, VectorType &point, VectorType &normal, 
			vtkSmartPointer<vtkImageData> &vtkImage,
		   	vtkSmartPointer<vtkPolyData> &line);

	std::vector<ImageType::Pointer> m_Images;
	VectorType m_Point;
	VectorType m_Normal;
	VectorType m_StartPoint;
	VectorType m_StartNormal;

	vtkBoundingBox m_BoundingBox;
	vtkSmartPointer<vtkRenderWindow> m_RenderWindow;
	vtkSmartPointer<vtkRenderWindowInteractor> m_Interactor;

};


}

#endif
