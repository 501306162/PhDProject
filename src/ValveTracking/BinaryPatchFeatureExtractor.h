#ifndef BINARY_PATCH_FEATURE_EXTRACTOR_H
#define BINARY_PATCH_FEATURE_EXTRACTOR_H

#include "FeatureExtractor.h"

namespace vt
{
class BinaryPatchFeatureExtractor : public FeatureExtractor<itk::Image<unsigned char, 3> >
{
public:
	typedef itk::Image<unsigned char, 3> ImageType;
	typedef BinaryPatchFeatureExtractor Self;
	typedef FeatureExtractor<ImageType> Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkTypeMacro(BinaryPatchFeatureExtractor, FeatureExtractor);
	itkNewMacro(Self);


	virtual void Extract(FeatureType &feature);


protected:
	BinaryPatchFeatureExtractor() {}
	virtual ~BinaryPatchFeatureExtractor() {}



private:
	BinaryPatchFeatureExtractor(const Self&);
	void operator=(const Self&);


};


}

#endif
