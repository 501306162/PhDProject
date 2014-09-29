#ifndef IMAGE_PATCH_EXTRACTOR_HPP
#define IMAGE_PATCH_EXTRACTOR_HPP

#include "ImagePatchExtractor.h"
#include <itkIdentityTransform.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>
#include <itkNumericTraits.h>

namespace filter
{
// ------------------------------------------------------------------------
template<typename TInputType>
ImagePatchExtractor<TInputType>::
ImagePatchExtractor()
{

	m_Interpolator = 0;
	m_Transform = 0;
	m_Scale = 1.0;
}

// ------------------------------------------------------------------------
template<typename TInputType>
void
ImagePatchExtractor<TInputType>::
GenerateData()
{
	if(!m_Transform)
	{
		// no transform set so just use identity
		m_Transform = itk::IdentityTransform<double, InputType::ImageDimension>::New();
	}	

	if(!m_Interpolator)
	{
		m_Interpolator = itk::LinearInterpolateImageFunction<InputType, double>::New();
	}


	

	typedef itk::ResampleImageFilter<InputType,InputType> ResamplerType;
	typename ResamplerType::Pointer resampler = ResamplerType::New();
	resampler->SetInput(this->GetInput());
	resampler->SetTransform(m_Transform);
	resampler->SetInterpolator(m_Interpolator);
	resampler->SetSize(m_PatchSize);
	resampler->SetOutputDirection(this->GetInput()->GetDirection());


	PixelType zeroValue = itk::NumericTraits<PixelType>::ZeroValue();
	resampler->SetDefaultPixelValue(zeroValue);


	// compute the output spacing and origin
	SpacingType outputSpacing;
	PointType outputOrigin;

	for(unsigned int i = 0; i < InputType::ImageDimension; i++)
	{
		double outputSize = m_Scale * 2.0;
		double spacing = outputSize / static_cast<double>(m_PatchSize[i]);

		unsigned int outputRadius = m_PatchSize[i] / 2;
		outputSpacing[i] = spacing;
		outputOrigin[i] = m_InputPoint[i] - (outputRadius*outputSpacing[i]) + (0.5 * outputSpacing[i]);
	}

	resampler->SetOutputSpacing(outputSpacing);
	resampler->SetOutputOrigin(outputOrigin);
	resampler->Update();


	this->GraftOutput(resampler->GetOutput());

}



// ------------------------------------------------------------------------
template<typename TInputType>
bool
ImagePatchExtractor<TInputType>::
PatchInsideBuffer(const PointType &point)
{
	
}

} /* filter */ 

#endif
