#ifndef LBP_FEATURE_EXTRACTOR_H
#define LBP_FEATURE_EXTRACTOR_H

#include <FeatureExtractor.h>


namespace vt
{
template<typename TImageType>
class LBPFeatureExtractor : public FeatureExtractor<TImageType>
{
public:
	typedef LBPFeatureExtractor Self;
	typedef FeatureExtractor<TImageType> Superclass;	
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	typedef itk::Image<unsigned int, 3> ImageType;
	

	itkTypeMacro(LBPFeatureExtractor, FeatureExtractor);
	itkNewMacro(Self);

protected:
	LBPFeatureExtractor() {}
	virtual ~LBPFeatureExtractor() {}


private:
	LBPFeatureExtractor(const Self&);
	void operator=(const Self&);

};


}

#endif
