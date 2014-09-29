#ifndef STRUCTURE_TENSOR_KEY_POINT_VALIDATOR_H
#define STRUCTURE_TENSOR_KEY_POINT_VALIDATOR_H

#include <itkProcessObject.h>
#include "FeatureCommon.h"


namespace filter
{
template<typename TImageType>
class StructureTensorKeyPointValidator : public itk::ProcessObject
{
public:
	typedef StructureTensorKeyPointValidator Self;
	typedef itk::ProcessObject Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkTypeMacro(StructureTensorKeyPointValidator, ProcessObject);
	itkNewMacro(Self);


	typedef filter::OrientatedKeyPoint<TImageType::ImageDimension> KeyPointType;

	typedef TImageType ImageType;
	typedef typename ImageType::PixelType PixelType;

	void SetInput(const ImageType * input);
	void SetKeyPoint(const KeyPointType &keyPoint);
	const ImageType * GetInput() const;

	bool IsValid();

	void Update();

	itkSetMacro(Ratio, double);
	 

protected:
	StructureTensorKeyPointValidator();
	virtual ~StructureTensorKeyPointValidator() {}


private:
	StructureTensorKeyPointValidator(const Self&);
	void operator=(const Self&);
	KeyPointType m_KeyPoint;
	double m_Ratio;
	bool m_Valid;
};
} /* filter */ 

#include "StructureTensorKeyPointValidator.hpp"


#endif
