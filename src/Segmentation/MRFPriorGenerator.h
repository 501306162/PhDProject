#ifndef MRF_PRIOR_GENERATOR_H
#define MRF_PRIOR_GENERATOR_H

#include <itkImageToImageFilter.h>

namespace segmentation
{
template<typename TInputType, typename TOutputType>
class MRFPriorGenerator : public itk::ImageToImageFilter<TInputType, TOutputType>
{
public:
	typedef MRFPriorGenerator Self;
	typedef itk::ImageToImageFilter<TInputType, TOutputType> Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;


	typedef TInputType InputType;
	

	typedef TOutputType OutputType;
	typedef typename OutputType::IndexType OutputIndexType;
	typedef typename OutputType::RegionType OutputRegionType;



	typedef itk::Image<double, InputType::ImageDimension> RealImageType;
	typedef typename RealImageType::Pointer RealImagePointer;

	itkTypeMacro(MRFPriorGenerator, ImageToImageFilter);
	itkNewMacro(Self);

	itkSetMacro(Spread, double);
	void SetDistanceMap(const RealImageType * map);


protected:
	MRFPriorGenerator();
	virtual ~MRFPriorGenerator() {}

	virtual void BeforeThreadedGenerateData();
	virtual void ThreadedGenerateData(const OutputRegionType &outputRegion,
			itk::ThreadIdType threadId);

private:
	MRFPriorGenerator(const Self&);
	void operator=(const Self&);

	bool m_DistanceMapInput;
	double m_Spread;

	RealImagePointer m_DistanceImage;

};
} /* segmenation */ 


#include "MRFPriorGenerator.hpp"

#endif
