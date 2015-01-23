#ifndef PATCH_EXTRACTOR2_H
#define PATCH_EXTRACTOR2_H

#include <itkObject.h>
#include <itkImageToImageFilter.h>
#include <itkInterpolateImageFunction.h>

using namespace itk;

namespace vt
{
template<typename TImageType>
class PatchExtractor2 : public itk::ImageToImageFilter<TImageType, TImageType> 
{
public:
	typedef PatchExtractor2 Self;
	typedef itk::ImageToImageFilter<TImageType, TImageType> Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	typedef TImageType ImageType;
	typedef typename ImageType::Pointer ImagePointer;
	typedef typename ImageType::RegionType ImageRegionType;
	typedef typename ImageType::PointType PointType;
	typedef typename ImageType::DirectionType DirectionType;
	typedef typename ImageType::IndexType IndexType;
	typedef typename ImageType::SizeType SizeType;
	typedef typename ImageType::SpacingType SpacingType;
	typedef itk::Vector<double, ImageType::ImageDimension> DistanceType;
	typedef itk::Vector<double, ImageType::ImageDimension> VectorType;

	typedef itk::InterpolateImageFunction<ImageType, double> InterpolatorType;
	typedef typename InterpolatorType::Pointer InterpolatorPointer;
	
	itkTypeMacro(PatchExtractor2, itk::ImageToImageFilter);
	itkNewMacro(Self);

	itkSetMacro(Distance, DistanceType);
	itkSetMacro(Size, SizeType);
	itkSetMacro(Center, PointType);
	itkSetMacro(Line, VectorType);
	itkSetObjectMacro(Interpolator, InterpolatorType);


protected:
	PatchExtractor2();
	virtual ~PatchExtractor2() {}

	virtual void ThreadedGenerateData(const ImageRegionType &,
                                    itk::ThreadIdType);
	virtual void GenerateOutputInformation();
	virtual void BeforeThreadedGenerateData();
	virtual void VerifyInputInformation() {}
	virtual void GenerateInputRequestedRegion();


private:
	PatchExtractor2(const Self&);
	void operator=(const Self&);

	DistanceType m_Distance;
	SizeType m_Size;
	PointType m_Center;
	VectorType m_Line;
	InterpolatorPointer m_Interpolator;




};



}

#include "PatchExtractor2.hpp"
#endif
