#ifndef IMAGE_PATCH_EXTRACTOR_H
#define IMAGE_PATCH_EXTRACTOR_H

#include <itkImageToImageFilter.h>
#include <itkTransform.h>
#include <itkInterpolateImageFunction.h>

namespace filter
{
template<typename TInputType>
class ImagePatchExtractor : public itk::ImageToImageFilter<TInputType, TInputType>
{
public:
	typedef ImagePatchExtractor Self;
	typedef itk::ImageToImageFilter<TInputType, TInputType> Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkTypeMacro(ImagePatchExtractor, ImageToImageFilter);
	itkNewMacro(Self);

	typedef TInputType InputType;
	typedef typename InputType::SizeType SizeType;
	typedef typename InputType::PointType PointType;
	typedef typename InputType::IndexType IndexType;
	typedef typename InputType::SpacingType SpacingType;
	typedef typename InputType::PointType RadiusType;
	typedef typename InputType::PixelType PixelType;

	typedef std::vector<PointType> PointSet;

	typedef itk::Transform<double, InputType::ImageDimension, InputType::ImageDimension> TransformType;
	typedef typename TransformType::Pointer TransformPointer;

	typedef itk::InterpolateImageFunction<InputType, double> InterpolatorType;
	typedef typename InterpolatorType::Pointer InterpolatorPointer;


	void SetInterpolator(InterpolatorType * interpolator) { m_Interpolator = interpolator; }
	itkSetObjectMacro(Transform, TransformType);

	itkSetMacro(PatchSize, SizeType);
	itkSetMacro(InputPoint, PointType);
	itkSetMacro(Scale, double);


	bool PatchInsideBuffer(const PointType &point);

protected:
	ImagePatchExtractor();
	virtual ~ImagePatchExtractor() {}


	void GenerateData();

	void GetPointsToExtract(PointSet &points);
	

private:
	
	PointType m_InputPoint;
	double m_Scale;

	InterpolatorPointer m_Interpolator;
	TransformPointer m_Transform;
	SizeType m_PatchSize;


	ImagePatchExtractor(const Self&);
	void operator=(const Self&);

};
} /* filter */ 

#include "ImagePatchExtractor.hpp"


#endif
