#ifndef TEST_DATA_H
#define TEST_DATA_H

#include <itkImage.h>
#include <CMRFileExtractor.h>
#include <itkSimilarity3DTransform.h>
#include <vtkBoundingBox.h>
#include <Eigen/Dense>
#include <LineIntersectionFinder.h>

namespace vt 
{
class TestData : public itk::Object 
{
public:
	typedef TestData Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;


	typedef itk::Image<unsigned short, 3> ImageType;
	typedef ImageType::PointType PointType;
	typedef itk::Image<unsigned char, 3> LabelType;
	typedef itk::Similarity3DTransform<double> TransformType;

	typedef Eigen::Vector3d VectorType;

	typedef struct image_group_
	{
		ImageType::Pointer image;
		LabelType::Pointer seg;
		LabelType::Pointer contours;
		LabelType::Pointer mask;
	} ImageGroup;

	typedef std::map<std::string, ImageGroup> TestGroup;
	typedef LineIntersectionFinder::OutputLineType LineType;


	typedef std::vector<LineType> LineTypeList;
	typedef std::map<std::string, LineTypeList> LineGroup;


	itkTypeMacro(TestData, Object);
	itkNewMacro(Self);

	TransformType::Pointer GetTransform() { return m_InitialTransform; }

	static TestData::Pointer Initialise(const std::string &directory);
	void Prepare(unsigned int timeStep);
	void SetBoundingBox(vtkBoundingBox &box) { m_BoundingBox = box; };

	void SaveImageGroup(std::string type, std::string name);
	void GetImageGroup(const std::string type, ImageGroup &group);
	void ExtractTestLines(const VectorType &normal, const VectorType &point, LineGroup &lines);
	void GetLinesFromPlane(const VectorType &normal, const VectorType &point, LineGroup &lines);

	void GetImages(std::vector<ImageType::Pointer> &images);
	void Load(const std::string &directory);
	void ProcessImages();
	itkGetMacro(Id, unsigned int);
protected:
	TestData() {}
	virtual ~TestData() {}

private:
	TestData(const Self&);
	void operator=(const Self&);
	void SaveGroup(const ImageGroup &group, std::string name);

	void SaveLabel(const LabelType::Pointer &label, std::string filename);
	void SaveImage(const ImageType::Pointer &label, std::string filename);
	void GetTestLines(const ImageGroup &images, const VectorType &normal, const VectorType &point, LineTypeList &lines);


	void ComputeInitialTransform();
	void TransformInitialGroup();
	void TransformImage(const ImageType::Pointer &input, TransformType::Pointer &transform, ImageType::Pointer &output);
	void ComputeImageMask(ImageGroup &group);
	void ComputeSegmentation(ImageGroup &group);
	void ComputeSegmentation2(ImageGroup &group);
	void ComputeContours(ImageGroup &group);

	unsigned int m_Id;

	CMRFileExtractor::Pointer m_Extractor;

	ImageType::Pointer m_StackImage;
	TestGroup m_InitialGroup;
	TestGroup m_AlignedGroup;
	vtkBoundingBox m_BoundingBox;
	TransformType::Pointer m_InitialTransform;
	vtkBoundingBox m_NewBox;
};



}

#endif
