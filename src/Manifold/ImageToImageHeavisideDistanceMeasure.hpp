#ifndef IMAGE_TO_IMAGE_HEAVISIDE_DISTANCE_MEASURE_HPP
#define IMAGE_TO_IMAGE_HEAVISIDE_DISTANCE_MEASURE_HPP

#include "ImageToImageHeavisideDistanceMeasure.h"
#include <itkImageRegionConstIterator.h>

#include <itkDomainThreader.h>


namespace manifold
{
// ------------------------------------------------------------------------
template<typename TAssociate, int TDimension>
void
ImageToImageHeavisideDistanceMeasureThreader<TAssociate, TDimension>::
BeforeThreadedExecution()
{
	// clear all of the input data
	delete[] this->m_Associate->m_ThreadData;


	// initialise new count data
	const itk::ThreadIdType numThreads = this->GetNumberOfThreadsUsed();
	this->m_Associate->m_ThreadData = new typename TAssociate::AlignedThreadData[numThreads];
	for(unsigned int i = 0; i < numThreads; i++)
	{
		this->m_Associate->m_ThreadData[i].m_Value = 0.0;
	}

}


// ------------------------------------------------------------------------
template<typename TAssociate, int TDimension>
void
ImageToImageHeavisideDistanceMeasureThreader<TAssociate, TDimension>::
AfterThreadedExecution()
{
	// clear all of the input data
	//delete[] this->m_Associate->m_ThreadData;
}



// ------------------------------------------------------------------------
template<typename TAssociate, int TDimension>
void
ImageToImageHeavisideDistanceMeasureThreader<TAssociate, TDimension>::
ThreadedExecution(const DomainType &domain, const itk::ThreadIdType id)
{
	// clear all of the input data
	this->m_Associate->ComputeThreadedRegion(domain, id);
}





// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
template<typename TInputImage>
void
ImageToImageHeavisideDistanceMeasure<TInputImage>::
ComputeThreadedRegion(const DomainType &domain, const itk::ThreadIdType id)
{
	typedef typename itk::ImageRegionConstIterator<ImageType> IteratorType;
	ImageConstPointer input1 = this->GetInput1();
	ImageConstPointer input2 = this->GetInput2();

	IteratorType it1(input1, domain);
	IteratorType it2(input2, domain);
	it1.GoToBegin(); it2.GoToBegin();
	
	while(!it1.IsAtEnd())
	{
		ValueType v1 = it1.Get();
		ValueType v2 = it2.Get();
		unsigned int h1 = 0, h2 = 0;


		// get the heaviside step function value
		if(v1 >= m_BoundryValue) h1 = 1;
		if(v2 >= m_BoundryValue) h2 = 1;


		// add the difference to the output
		m_ThreadData[id].m_Value += (h1-h2)*(h1-h2);


		++it1;++it2;
	}
}



// ------------------------------------------------------------------------
template<typename TImageType>
ImageToImageHeavisideDistanceMeasure<TImageType>::ImageToImageHeavisideDistanceMeasure()
{
	m_BoundryValue = 1;
	m_Threader = ThreaderType::New();
	m_NumberOfThreads = m_Threader->GetMaximumNumberOfThreads();
}


// ------------------------------------------------------------------------
template<typename TInputImage>
double ImageToImageHeavisideDistanceMeasure<TInputImage>::GetValue()
{
	Superclass::PreCheck();
	
	m_Threader->SetMaximumNumberOfThreads(m_NumberOfThreads);
	m_Threader->Execute(this, this->GetInput1()->GetLargestPossibleRegion());

	double output = 0.0;
	for(unsigned int i = 0; i < m_Threader->GetNumberOfThreadsUsed(); i++)
	{
		output+= m_ThreadData[i].m_Value;
	}

	
	return output;
}


// ------------------------------------------------------------------------
template<typename TImageType>
ImageToImageHeavisideDistanceMeasure<TImageType>::
~ImageToImageHeavisideDistanceMeasure()
{
	delete[] m_ThreadData;
   	m_ThreadData = NULL;	
}

} /* manifold */ 

#endif
