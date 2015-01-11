#ifndef FEATURE_EXTRACTOR_H
#define FEATURE_EXTRACTOR_H

#include <itkImage.h>
#include <Eigen/Dense>
#include <MatrixCommon.h>

namespace vt 
{
template<typename TInputType>
class FeatureExtractor : public itk::Object 
{
public:
	typedef FeatureExtractor Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(FeatureExtractor, Object);


	typedef TInputType ImageType;
	typedef typename ImageType::Pointer ImagePointer;
	typedef utils::DoubleMatrixType FeatureType;
	
	void SetInput(const ImagePointer &input) { m_Input = input; }
	virtual void Extract(FeatureType &feature) = 0;



protected:
	FeatureExtractor() {}
	virtual ~FeatureExtractor() {}

	ImagePointer m_Input;

private:
	FeatureExtractor(const Self&);
	void operator=(const Self&);

};

}

#endif
