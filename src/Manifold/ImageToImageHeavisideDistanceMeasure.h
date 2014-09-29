#ifndef IMAGE_TO_IMAGE_HEAVISIDE_DISTANCE_MEASURE_H
#define IMAGE_TO_IMAGE_HEAVISIDE_DISTANCE_MEASURE_H


#include "ImageToImageDistanceMeasure.h"
#include <itkDomainThreader.h>
#include <itkThreadedImageRegionPartitioner.h>


namespace manifold
{
template<typename TAssociate, int TDimension>
class ImageToImageHeavisideDistanceMeasureThreader : 
	public itk::DomainThreader<itk::ThreadedImageRegionPartitioner<TDimension>, TAssociate>
{
public:
	typedef ImageToImageHeavisideDistanceMeasureThreader Self;
	typedef itk::DomainThreader<
		itk::ThreadedImageRegionPartitioner<TDimension>, TAssociate> Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	
	typedef typename Superclass::DomainType DomainType;
	itkNewMacro(Self);

protected:
	virtual void BeforeThreadedExecution();
	virtual void AfterThreadedExecution();
	virtual void ThreadedExecution(const DomainType &domain, const itk::ThreadIdType id);
	ImageToImageHeavisideDistanceMeasureThreader() {}	
};




template<typename TInputImage>
class ImageToImageHeavisideDistanceMeasure : public ImageToImageDistanceMeasure<TInputImage>
{
public:
	typedef ImageToImageHeavisideDistanceMeasure<TInputImage> 	Self;
	typedef ImageToImageDistanceMeasure<TInputImage>          	Superclass;
	typedef itk::SmartPointer<Self>								Pointer;
	typedef itk::SmartPointer<const Self>						ConstPointer;
	itkNewMacro(Self);
	itkTypeMacro(ImageToImageHeavisideDistanceMeasure, ImageToImageDistanceMeasure);


	/** Typedefs for the image inputs */
	typedef TInputImage ImageType;
	typedef typename ImageType::Pointer InputImage1Pointer;
	typedef typename ImageType::ConstPointer ImageConstPointer;
	typedef typename ImageType::ValueType ValueType;

	virtual double GetValue();

	itkSetMacro(BoundryValue, ValueType);
	itkSetMacro(NumberOfThreads, unsigned int);
	itkGetMacro(NumberOfThreads, unsigned int);
	

protected:
	ImageToImageHeavisideDistanceMeasure();
	virtual ~ImageToImageHeavisideDistanceMeasure();

	ValueType m_BoundryValue;




	struct ThreadData
	{
		double m_Value;
	};

	itkAlignedTypedef(64, ThreadData, AlignedThreadData);
	AlignedThreadData * m_ThreadData;

	unsigned int m_NumberOfThreads;

	void Initialise();

	typedef ImageToImageHeavisideDistanceMeasureThreader<Self, ImageType::ImageDimension> ThreaderType;
	typedef typename ThreaderType::DomainType DomainType;
	friend class ImageToImageHeavisideDistanceMeasureThreader<Self, ImageType::ImageDimension>;
	typename ThreaderType::Pointer m_Threader;
	
	// threaded compute function
	virtual void ComputeThreadedRegion(const DomainType &domain, const itk::ThreadIdType id);

};


} /* manifold */ 


#include "ImageToImageHeavisideDistanceMeasure.hpp"

#endif
