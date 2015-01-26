#ifndef COMPUTE_PROB_IMAGE_H
#define COMPUTE_PROB_IMAGE_H

#include <itkImageToImageFilter.h>
#include <SVMClassifier.h>
#include <PatchExtractor2.h>

using namespace itk;
namespace vt
{
class ComputeProbImage : public itk::ImageToImageFilter<
						 itk::Image<unsigned short, 3>, 
						 itk::Image<double, 3> >
{
public: 
	typedef ComputeProbImage Self;
	typedef itk::ImageToImageFilter<
		itk::Image<unsigned short, 3>, 
		itk::Image<double, 3> > Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkNewMacro(Self);

	typedef itk::Image<unsigned short, 3> ImageType;
	typedef itk::Image<unsigned char, 3> MaskType;
	typedef itk::Image<double, 3> OutputType;
	typedef OutputType::RegionType RegionType;
	typedef OutputType::IndexType IndexType;
	
	typedef PatchExtractor2<ImageType> ExtractorType;
	typedef ExtractorType::PointType PointType;
	typedef ExtractorType::DistanceType DistanceType;
	typedef ExtractorType::VectorType VectorType;
	typedef ExtractorType::SizeType SizeType;

	typedef SVMClassifier::MatrixType MatrixType;
	typedef SVMClassifier::IntMatrixType IntMatrixType;

	void SetClassifier(const SVMClassifier::Pointer &cls) { m_Classifier = cls; }
	void SetMask(const MaskType::Pointer &mask) { m_Mask = mask; m_HasMask = true; }

	itkSetMacro(PatchDistance, DistanceType);
	itkSetMacro(PatchSize, SizeType);
	itkSetMacro(Line, VectorType);
	itkSetMacro(Radius, double);
	itkSetMacro(GridSize, unsigned int);
	itkSetMacro(Neighbors, unsigned int);
	
	typedef std::vector<IndexType> IndexList;
	typedef std::vector<IndexList> ThreadIndexList;


protected:
	ComputeProbImage();
	virtual ~ComputeProbImage() {}

	virtual void ThreadedGenerateData(const RegionType & outputRegionForThread, ThreadIdType threadId);

private:
	ComputeProbImage(const Self&);
	void operator=(const Self&);

	void ExtractLBPFeature(const PointType &point, MatrixType &feature);

	MaskType::Pointer m_Mask;
	SVMClassifier::Pointer m_Classifier;
	DistanceType m_PatchDistance;
	VectorType m_Line;
	SizeType m_PatchSize;

	MatrixType m_Features;
	bool m_HasMask;
	double m_Radius;
	unsigned int m_GridSize;
	unsigned int m_Neighbors;
	ThreadIndexList m_Indices;
	std::vector<std::vector<double> > m_Values;

};




}

#endif
