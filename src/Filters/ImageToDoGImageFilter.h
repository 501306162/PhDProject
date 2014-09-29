#ifndef IMAGE_TO_DO_G_IMAGE_FILTER_H
#define IMAGE_TO_DO_G_IMAGE_FILTER_H

#include <itkImageToImageFilter.h>

namespace filter
{
template<typename TInputType, typename TOutputType>
class ImageToDoGImageFilter : public itk::ImageToImageFilter<TInputType, TOutputType>
{
public:
	typedef ImageToDoGImageFilter Self;
	typedef itk::ImageToImageFilter<TInputType, TOutputType> Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkTypeMacro(ImageToDoGImageFilter, ImageToImageFilter);
	itkNewMacro(Self);

	typedef TInputType InputType;
	typedef TOutputType OutputType;
	typedef itk::Image<double, InputType::ImageDimension> RealImageType;


	itkSetMacro(Sigma1, double);
	itkSetMacro(Sigma2, double);
	itkSetMacro(UseImageSpacing, bool);

protected:
	ImageToDoGImageFilter();
	virtual ~ImageToDoGImageFilter() {}

	void GenerateData();

private:
	ImageToDoGImageFilter(const Self&);
	void operator=(const Self&);


	double m_Sigma1;
	double m_Sigma2;
	bool m_UseImageSpacing;

};
} /* filter */ 

#include "ImageToDoGImageFilter.hpp"

#endif
