#ifndef IMAGE_TO_IMAGE_DISTANCE_MEASURE_HPP
#define IMAGE_TO_IMAGE_DISTANCE_MEASURE_HPP

#include "ImageToImageDistanceMeasure.h"


namespace manifold
{
// ------------------------------------------------------------------------
template<typename TInputImage>
ImageToImageDistanceMeasure<TInputImage>::ImageToImageDistanceMeasure()
{
	this->SetNumberOfRequiredInputs(2);	
}

// ------------------------------------------------------------------------
template<typename TInputImage>
void ImageToImageDistanceMeasure<TInputImage>::SetInput1(const ImageType * image)
{
	this->ProcessObject::SetNthInput(0, 
			const_cast<ImageType *>(image));
}

// ------------------------------------------------------------------------
template<typename TInputImage>
void ImageToImageDistanceMeasure<TInputImage>::SetInput2(const ImageType * image)
{
	this->ProcessObject::SetNthInput(1, 
			const_cast<ImageType *>(image));
}


// ------------------------------------------------------------------------
template<typename TInputImage>
const typename ImageToImageDistanceMeasure<TInputImage>::ImageType *
ImageToImageDistanceMeasure<TInputImage>::GetInput1() const
{
	return itkDynamicCastInDebugMode<const ImageType *>(
			this->ProcessObject::GetInput(0));
}


// ------------------------------------------------------------------------
template<typename TInputImage>
const typename ImageToImageDistanceMeasure<TInputImage>::ImageType *
ImageToImageDistanceMeasure<TInputImage>::GetInput2() const
{
	return itkDynamicCastInDebugMode<const ImageType *>(
			this->ProcessObject::GetInput(1));
}

// ------------------------------------------------------------------------
template<typename TInputImage>
void ImageToImageDistanceMeasure<TInputImage>::PreCheck()
throw(itk::ExceptionObject)
{
	if(this->ProcessObject::GetNumberOfInputs() != 2)
	{
		itkExceptionMacro(<< "Inputs Have Not Been Set");	
	}

	// check that the images are the same size
	ImageConstPointer input1 = this->GetInput1();
	ImageConstPointer input2 = this->GetInput2();

	if(input1->GetLargestPossibleRegion().GetSize() != input2->GetLargestPossibleRegion().GetSize())
	{
		itkExceptionMacro(<< "Inputs are Not of the Same Size");
	}
}

// ------------------------------------------------------------------------
template<typename TInputImage>
ImageToImageDistanceMeasure<TInputImage>::~ImageToImageDistanceMeasure()
{
}

} /* manifold */ 


#endif
