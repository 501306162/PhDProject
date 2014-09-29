#ifndef IMAGE_TO_IMAGE_DISTANCE_MEASURE_H
#define IMAGE_TO_IMAGE_DISTANCE_MEASURE_H


#include <itkProcessObject.h>
#include <itkImage.h>

namespace manifold
{

template<typename TInputImage>
class ImageToImageDistanceMeasure : public itk::ProcessObject
{
public:

	typedef ImageToImageDistanceMeasure<TInputImage>        Self;
	typedef itk::ProcessObject                              Superclass;
	typedef itk::SmartPointer<Self>  						Pointer;
	typedef itk::SmartPointer<const Self>                   ConstPointer;
	
	itkTypeMacro(ImageToImageDistanceMeasure, itk::ProcessObject);


	/** Typedefs for the image inputs */
	typedef TInputImage ImageType;
	typedef typename ImageType::Pointer InputImage1Pointer;
	typedef typename ImageType::ConstPointer ImageConstPointer;

	virtual ~ImageToImageDistanceMeasure ();

	virtual double GetValue() = 0;

	void SetInput1(const ImageType * image);
	void SetInput2(const ImageType * image);
	const ImageType * GetInput1() const;
	const ImageType * GetInput2() const;


protected:
	ImageToImageDistanceMeasure ();

	virtual void PreCheck() throw (itk::ExceptionObject);

private:
	/* data */
};	

} /* Manifold */ 



#ifndef ITK_MANUAL_INSTANTIATION
#include "ImageToImageDistanceMeasure.hpp"
#endif


#endif
