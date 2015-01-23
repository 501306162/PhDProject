#ifndef SIMPLE_MRF_SEGMENTER_H
#define SIMPLE_MRF_SEGMENTER_H

#include <itkImage.h>
#include <itkGaussianMembershipFunction.h>
#include <itkListSample.h>

namespace vt
{
class SimpleMRFSegmenter : public itk::Object
{
public:
	typedef SimpleMRFSegmenter Self;
   	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;	

	itkNewMacro(Self);
	itkTypeMacro(SimpleMRFSegmenter, Object);

	typedef itk::Image<unsigned short, 3> ImageType;
	typedef itk::Image<unsigned char, 3> MaskType;

	typedef itk::FixedArray<unsigned short,1> MeasurementVectorType;
	typedef itk::Statistics::GaussianMembershipFunction<MeasurementVectorType> MembershipFunctionType;
	typedef itk::Statistics::ListSample<MeasurementVectorType> SampleType;
	typedef itk::Vector<double,1> MeanType;
	typedef itk::VariableSizeMatrix<double> VarianceType;

	void SetImage(const ImageType::Pointer &image) { m_Image = image; }
	void SetMask(const MaskType::Pointer &mask) { m_Mask = mask; }

	void Segment();

	MaskType::Pointer GetOutput() const { return m_Output; }
	itkSetMacro(SmoothnessCost, double);
	itkSetMacro(OutputValue, MaskType::ValueType); 

	itkSetMacro(MaskValue, unsigned char);


protected:
	SimpleMRFSegmenter();
	virtual ~SimpleMRFSegmenter() {}


private:
	void CreateMembershipFunction();


	SimpleMRFSegmenter(const Self&);
	void operator=(const Self&);

	ImageType::Pointer m_Image;
	MaskType::Pointer m_Mask;
	MaskType::Pointer m_Output;
	MaskType::ValueType m_OutputValue;
	unsigned char m_MaskValue;

	MembershipFunctionType::Pointer m_PositiveFunction;
	MembershipFunctionType::Pointer m_NegativeFunction;
	
	double m_InitialPositiveMean;
	double m_InitialNegativeMean;
	double m_InitialPositiveVariance;
	double m_InitialNegativeVariance;
	double m_SmoothnessCost;
	double m_SmallestCheck;
	



};

}

#endif
