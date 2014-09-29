#ifndef MRF_MEAN_SQUARE_UNARY_POTENTIALS_METRIC_H
#define MRF_MEAN_SQUARE_UNARY_POTENTIALS_METRIC_H

#include "MRFUnaryPotentialsMetricBase.h"

namespace myitk
{
using namespace itk;
template< typename TFixedImage, typename TMovingImage >
class MRFMeanSquaredUnaryPotentialsMetric : public MRFUnaryPotentialsMetricBase<TFixedImage, TMovingImage>
{
public:
	typedef MRFMeanSquaredUnaryPotentialsMetric                     Self;
	typedef MRFUnaryPotentialsMetricBase<TFixedImage, TMovingImage> Superclass;
	typedef SmartPointer< Self >                                    Pointer;
	typedef SmartPointer< const Self >                              ConstPointer;


	itkNewMacro(Self);
	itkTypeMacro(MRFMeanSquaredUnaryPotentialsMetric, MRFUnaryPotentialsMetricBase);
	itkStaticConstMacro(ImageDimension, unsigned int, TMovingImage::ImageDimension);

	/** Typedefs for the moving image */
	typedef TMovingImage                            MovingImageType;
	typedef typename MovingImageType::Pointer       MovingImagePointer;
	typedef typename MovingImageType::ConstPointer  MovingImageConstPointer;
	typedef typename MovingImageType::IndexType     MovingImageIndexType;
	typedef typename MovingImageType::PointType     MovingImagePointType;
	typedef typename MovingImageType::ValueType     MovingImageValueType;

	typedef TFixedImage                             FixedImageType;
	typedef typename FixedImageType::Pointer        FixedImagePointer;
	typedef typename FixedImageType::ConstPointer   FixedImageConstPointer;
	typedef typename FixedImageType::IndexType      FixedImageIndexType;
	typedef typename FixedImageType::PointType      FixedImagePointType;
	typedef typename FixedImageType::ValueType      FixedImageValueType;

	/** Typedefs for the samples list */
	typedef typename Superclass::FixedImageSampleValue       FixedImageSampleValue;
	typedef typename Superclass::FixedImageSamples           FixedImageSamples;
	typedef typename Superclass::FixedImageSamplesList       FixedImageSamplesList;

	typedef typename Superclass::TranslationTransformType    TranslationTransformType;
	typedef typename TranslationTransformType::Pointer       TranslationTransformPointer;
	typedef typename TranslationTransformType::ParametersType ParametersType;

	typedef typename Superclass::MeasureType  MeasureType;
	typedef typename Superclass::DerivativeType DerivativeType;


	typedef typename Superclass::SampleRegion SampleRegion;



protected:
	MRFMeanSquaredUnaryPotentialsMetric() {}
	virtual ~MRFMeanSquaredUnaryPotentialsMetric() {}

	using Superclass::m_OutputValue;

	virtual void ComputeUnaryPotentials( const SampleRegion &region, int threadId );

private:
};


}

#ifndef ITK_MANUAL_INSTANTIATION
#include "MRFMeanSquaredUnaryPotentialsMetric.hpp"
#endif


#endif
