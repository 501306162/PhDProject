#ifndef MRF_LABEL_ALIGNMENT_H
#define MRF_LABEL_ALIGNMENT_H

#include <itkImageToImageFilter.h>
#include <itkTransform.h>


namespace segmentation
{
template<typename TFixedType, typename TMovingType, typename TTransformType>
class MRFLabelAlignment : public itk::ImageToImageFilter<TMovingType, TMovingType>
{
public:
	typedef MRFLabelAlignment Self;
	typedef itk::ImageToImageFilter<TMovingType, TMovingType> Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkTypeMacro(MRFLabelAlignment, ImageToImageFilter);
	itkNewMacro(Self);


	/** typedefs for the images */
	typedef TFixedType FixedType;
	void SetFixedImage(const FixedType * fixed);
	typedef TMovingType MovingType;
	typedef typename MovingType::Pointer MovingImagePointer;
	void SetMovingImage(const MovingType * moving);


	typedef TTransformType TransformType;
	typedef typename TransformType::ParametersType ParametersType;
	typedef typename TransformType::Pointer TransformPointer;

	void SetTransform(const TransformType * transform);
	ParametersType GetFinalParameters() const;
	TransformPointer GetTransform() const;

protected:
	MRFLabelAlignment();
	virtual ~MRFLabelAlignment() {};

	virtual void GenerateData();
	void ComputeTransform();
	void ResampleImage();

private:
	MRFLabelAlignment(const Self&);
	void operator=(const Self&);


	MovingImagePointer m_MovingImage;
	TransformPointer m_Transform;

};

} /* segmentation */ 

#include "MRFLabelAlignment.hpp"


#endif
