#ifndef MRF_PRIOR_GENERATOR_HPP
#define MRF_PRIOR_GENERATOR_HPP

#include "MRFPriorGenerator.h"
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>

namespace segmentation
{
// ------------------------------------------------------------------------
template<typename TInputType, typename TOutputType>
MRFPriorGenerator<TInputType, TOutputType>::
MRFPriorGenerator()
{
	m_Spread = 1.0;
	m_DistanceMapInput = false;
}

// ------------------------------------------------------------------------
template<typename TInputType, typename TOutputType>
void
MRFPriorGenerator<TInputType, TOutputType>::
BeforeThreadedGenerateData()
{
	if(m_DistanceMapInput)
	{
		return;
	}



	// compute the distance image 
	typedef itk::SignedMaurerDistanceMapImageFilter<InputType, RealImageType> DistanceFilterType;
	typename DistanceFilterType::Pointer distanceFilter = DistanceFilterType::New();
	distanceFilter->SetInput(this->GetInput(0));
	distanceFilter->SetBackgroundValue(0);
	distanceFilter->SetUseImageSpacing(true);
	distanceFilter->SquaredDistanceOff();
	distanceFilter->Update();

	m_DistanceImage = distanceFilter->GetOutput();

}

// ------------------------------------------------------------------------
template<typename TInputType, typename TOutputType>
void
MRFPriorGenerator<TInputType, TOutputType>::
SetDistanceMap(const RealImageType * map)
{
	m_DistanceImage = map;
	m_DistanceMapInput = true;
}

// ------------------------------------------------------------------------
template<typename TInputType, typename TOutputType>
void
MRFPriorGenerator<TInputType, TOutputType>::
ThreadedGenerateData(const OutputRegionType &outputRegion,
		itk::ThreadIdType threadId)
{
	typedef itk::ImageRegionConstIterator<RealImageType> RealIt;
	RealIt realIt(m_DistanceImage, outputRegion);
	typedef itk::ImageRegionIterator<OutputType> OutIt;
	OutIt outIt(this->GetOutput(), outputRegion);


	while(!outIt.IsAtEnd())
	{
		double dist = realIt.Get();
		if(dist <= 0.0)
		{
			outIt.Set(1.0);
		}
		else
		{
			double outVal = exp(-dist / m_Spread);
			outIt.Set(outVal);
		}
		
		++realIt; ++outIt;
	}

}

} /* segmentation */ 

#endif
