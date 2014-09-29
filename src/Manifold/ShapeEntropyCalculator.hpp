#ifndef SHAPE_ENTROPY_CALCULATOR_HPP
#define SHAPE_ENTROPY_CALCULATOR_HPP

#include "ShapeEntropyCalculator.h"
#include <itkResampleImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkIdentityTransform.h>
#include <itkImageRegionConstIterator.h>


namespace manifold
{

// ------------------------------------------------------------------------
template<typename TAssociate, int TDimension>
void
ShapeEntropyCalculatorThreader<TAssociate, TDimension>::
BeforeThreadedExecution()
{
	// clear the threaded values
	if(this->m_Associate->m_ThreadData != NULL)
		delete[] this->m_Associate->m_ThreadData;

	const unsigned int numThreads = this->GetNumberOfThreadsUsed();
	this->m_Associate->m_ThreadData = new typename TAssociate::AlignedThreadData[numThreads];

	for(unsigned int i = 0; i < numThreads; i++)
	{
		this->m_Associate->m_ThreadData[i].m_Value = 0.0;
	}
}

// ------------------------------------------------------------------------
template<typename TAssociate, int TDimension>
void
ShapeEntropyCalculatorThreader<TAssociate, TDimension>::
ThreadedExecution(const DomainType &domain, const itk::ThreadIdType id)
{
	this->m_Associate->ComputeRegionEntropy(domain, id);
}



// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
template<int TDimension>
ShapeEntropyCalculator<TDimension>::
ShapeEntropyCalculator()
{
	this->SetNumberOfRequiredInputs(2);
	m_Threader = ThreaderType::New();
	m_NumberOfThreads = m_Threader->GetMaximumNumberOfThreads();
	m_ThreadData = NULL;
	
}


// ------------------------------------------------------------------------
template<int TDimension>
void
ShapeEntropyCalculator<TDimension>::
SetInputShape(const LabelType * shape)
{
	this->SetNthInput(1, const_cast<LabelType*>(shape));
}

// ------------------------------------------------------------------------
template<int TDimension>
void
ShapeEntropyCalculator<TDimension>::
SetInputImage(const ImageType * image)
{
	this->SetNthInput(0, const_cast<ImageType*>(image));
}

// ------------------------------------------------------------------------
template<int TDimension>
const typename ShapeEntropyCalculator<TDimension>::ImageType *
ShapeEntropyCalculator<TDimension>::
GetInputImage() const
{
	return itkDynamicCastInDebugMode<const ImageType*>(this->GetInput(0));
}

// ------------------------------------------------------------------------
template<int TDimension>
const typename ShapeEntropyCalculator<TDimension>::LabelType *
ShapeEntropyCalculator<TDimension>::
GetInputShape() const
{
	return itkDynamicCastInDebugMode<const LabelType*>(this->GetInput(1));
}

// ------------------------------------------------------------------------
template<int TDimension>
void
ShapeEntropyCalculator<TDimension>::
ComputeHistograms()
{
	// resample the mask image
	const LabelType * mask = this->GetInputShape();
	const ImageType * image = this->GetInputImage();

	typedef itk::ResampleImageFilter<LabelType, LabelType> ResamplerType;
	typedef itk::NearestNeighborInterpolateImageFunction<LabelType, double> InterpolatorType;
	typedef itk::IdentityTransform<double, TDimension> IdentityType;

	typename ResamplerType::Pointer resampler = ResamplerType::New();
	typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
	typename IdentityType::Pointer transform = IdentityType::New();
	transform->SetIdentity();

	resampler->SetInput(mask);
	resampler->SetOutputSpacing(image->GetSpacing());
	resampler->SetOutputDirection(image->GetDirection());
	resampler->SetOutputOrigin(image->GetOrigin());
	resampler->SetSize(image->GetLargestPossibleRegion().GetSize());
	resampler->SetInterpolator(interpolator);
	resampler->SetTransform(transform);
	resampler->Update();






	HistogramMeasurementVectorType binMin(1), binMax(1);
	binMin[0] = -0.5;
	binMax[0] = 500.5;

	HistogramSizeType histSize(1);
	histSize[0] = 50;
	

	fgHistogramFilter = HistogramFilterType::New();
	fgHistogramFilter->SetInput(image);
	fgHistogramFilter->SetHistogramSize(histSize);
	fgHistogramFilter->SetHistogramBinMaximum(binMax);
	fgHistogramFilter->SetHistogramBinMinimum(binMin);
	fgHistogramFilter->SetMaskImage(resampler->GetOutput());
	fgHistogramFilter->SetMaskValue(1);
	fgHistogramFilter->Update();

	bgHistogramFilter = HistogramFilterType::New();
	bgHistogramFilter->SetInput(image);
	bgHistogramFilter->SetHistogramSize(histSize);
	bgHistogramFilter->SetHistogramBinMaximum(binMax);
	bgHistogramFilter->SetHistogramBinMinimum(binMin);
	bgHistogramFilter->SetMaskImage(resampler->GetOutput());
	bgHistogramFilter->SetMaskValue(0);
	bgHistogramFilter->Update();


}

// ------------------------------------------------------------------------
template<int TDimension>
double
ShapeEntropyCalculator<TDimension>::
GetValue()
{
	// first we compute the histograms
	ComputeHistograms();
	

	// now we set up the threader
	m_Threader->SetMaximumNumberOfThreads(m_NumberOfThreads);
	m_Threader->Execute(this, this->GetInputImage()->GetLargestPossibleRegion());


	double entropy = 0.0;
	for(unsigned int i = 0; i < m_Threader->GetNumberOfThreadsUsed(); i++)
	{
		entropy += m_ThreadData[i].m_Value;
	}

	delete[] m_ThreadData;
	m_ThreadData = NULL;

	return entropy;
}

// ------------------------------------------------------------------------
template<int TDimension>
void
ShapeEntropyCalculator<TDimension>::
ComputeRegionEntropy(const DomainType &domain, 
			const itk::ThreadIdType id)
{
	const ImageType * image = this->GetInputImage();
	const HistogramType * fgHist = fgHistogramFilter->GetOutput();
	const HistogramType * bgHist = bgHistogramFilter->GetOutput();

	itk::ImageRegionConstIterator<ImageType> it(image, domain);
	
	while(!it.IsAtEnd())
	{
		HistogramMeasurementVectorType mv(1);
		mv[0] = it.Get();

		typedef typename HistogramType::IndexType HistIndexType;
		HistIndexType histIndex;
		fgHist->GetIndex(mv, histIndex);


		double bgProb = (double) bgHist->GetFrequency(histIndex) / (double) bgHist->GetTotalFrequency();
		double fgProb = (double) fgHist->GetFrequency(histIndex) / (double) fgHist->GetTotalFrequency();

		fgProb = std::max(0.000000001, fgProb);
		bgProb = std::max(0.000000001, bgProb);

		double entropy = -1 * ((fgProb * log(fgProb)) + (bgProb*log(bgProb)));


		m_ThreadData[id].m_Value += entropy; 

		++it;
	}
}

} /* manifold */ 

#endif
