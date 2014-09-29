#ifndef IMAGE_TO_DO_G_IMAGE_FILTER_HPP
#define IMAGE_TO_DO_G_IMAGE_FILTER_HPP

#include "ImageToDoGImageFilter.h"
#include <itkDiscreteGaussianImageFilter.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkSubtractImageFilter.h>


namespace filter
{
// ------------------------------------------------------------------------
template<typename TInputType, typename TOutputType>
ImageToDoGImageFilter<TInputType, TOutputType>::
ImageToDoGImageFilter()
{
	m_Sigma1 = 1.0;
	m_Sigma2 = 2.0;
}

// ------------------------------------------------------------------------
template<typename TInputType, typename TOutputType>
void
ImageToDoGImageFilter<TInputType, TOutputType>::
GenerateData()
{
	const InputType * input = this->GetInput(0);

	typedef itk::RecursiveGaussianImageFilter<InputType, RealImageType> GaussianFilterType;
	typename GaussianFilterType::Pointer smoother1 = GaussianFilterType::New();
	smoother1->SetInput(input);
	smoother1->SetSigma(m_Sigma1);
	//smoother1->SetUseImageSpacing(m_UseImageSpacing);

	typename GaussianFilterType::Pointer smoother2 = GaussianFilterType::New();
	smoother2->SetInput(input);
	smoother2->SetSigma(m_Sigma2);
	//smoother2->SetUseImageSpacing(m_UseImageSpacing);


	typedef itk::SubtractImageFilter<RealImageType, RealImageType, OutputType> SubtractorType;
	typename SubtractorType::Pointer subtractor = SubtractorType::New();
	subtractor->SetInput1(smoother2->GetOutput());
	subtractor->SetInput2(smoother1->GetOutput());

	subtractor->Update();

	this->GraftOutput(subtractor->GetOutput());
}

} /* filter */ 



#endif
