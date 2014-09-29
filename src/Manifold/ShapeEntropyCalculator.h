#ifndef SHAPE_ENTROPY_CALCULATOR_H
#define SHAPE_ENTROPY_CALCULATOR_H

#include <itkProcessObject.h>
#include <itkDomainThreader.h>
#include <itkThreadedImageRegionPartitioner.h>

#include <itkMaskedImageToHistogramFilter.h>

namespace manifold
{
template<typename TAssociate, int TDimension>
class ShapeEntropyCalculatorThreader : 
	public itk::DomainThreader<itk::ThreadedImageRegionPartitioner<TDimension>, TAssociate>
{
public:
	typedef ShapeEntropyCalculatorThreader Self;
	typedef itk::DomainThreader<
		itk::ThreadedImageRegionPartitioner<TDimension>, TAssociate> Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkNewMacro(Self);
	typedef typename Superclass::DomainType DomainType;

protected:
	virtual void BeforeThreadedExecution();
	virtual void ThreadedExecution(const DomainType &domain, 
			itk::ThreadIdType id);


	ShapeEntropyCalculatorThreader() {}
	~ShapeEntropyCalculatorThreader() {}

private:
	ShapeEntropyCalculatorThreader(const Self&);
	void operator=(const Self&);	
};	

// ------------------------------------------------------------------------
template<int TDimension>
class ShapeEntropyCalculator : public itk::ProcessObject
{
public:
	typedef ShapeEntropyCalculator Self;
	typedef itk::ProcessObject Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkNewMacro(Self);


	typedef itk::Image<unsigned char, TDimension> LabelType;
	typedef itk::Image<unsigned short, TDimension> ImageType;

	// define the histograms
	typedef itk::Statistics::MaskedImageToHistogramFilter<ImageType, LabelType> HistogramFilterType;
	typedef typename HistogramFilterType::HistogramType HistogramType;
	typedef typename HistogramFilterType::HistogramSizeType HistogramSizeType;
	typedef typename HistogramFilterType::HistogramMeasurementVectorType HistogramMeasurementVectorType;


	void SetInputImage(const ImageType * image);
	const ImageType * GetInputImage() const; 
	void SetInputShape(const LabelType * label);
	const LabelType * GetInputShape() const;


	itkSetMacro(NumberOfThreads, unsigned int);
	itkGetMacro(NumberOfThreads, unsigned int);


	virtual double GetValue();


protected:
	ShapeEntropyCalculator();
	~ShapeEntropyCalculator() {}

	void ComputeHistograms();

	typename HistogramFilterType::Pointer fgHistogramFilter;
	typename HistogramFilterType::Pointer bgHistogramFilter;



private:
	ShapeEntropyCalculator(const Self&);
	void operator=(const Self&);


	typedef ShapeEntropyCalculatorThreader<Self, TDimension> ThreaderType;
	typedef typename ThreaderType::DomainType DomainType;
	friend class ShapeEntropyCalculatorThreader<Self, TDimension>;

	typename ThreaderType::Pointer m_Threader;
	unsigned int m_NumberOfThreads;

	virtual void ComputeRegionEntropy(const DomainType &domain, 
			const itk::ThreadIdType id);


	struct ThreadData
	{
		double m_Value;
	};

	itkAlignedTypedef(64, ThreadData, AlignedThreadData);
	AlignedThreadData * m_ThreadData;

};



} /* manifold */ 


#include "ShapeEntropyCalculator.hpp"

#endif
