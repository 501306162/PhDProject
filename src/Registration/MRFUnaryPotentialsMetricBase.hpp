#ifndef MRF_UNARY_POTENTIALS_METRIC_BASE_HPP
#define MRF_UNARY_POTENTIALS_METRIC_BASE_HPP

#include "MRFUnaryPotentialsMetricBase.h"

namespace myitk
{
using namespace itk;
template< typename TFixedImage, typename TMovingImage >
MRFUnaryPotentialsMetricBase<TFixedImage, TMovingImage>::
MRFUnaryPotentialsMetricBase()
{
	// initailise the moving image
	this->m_MovingImage = 0;
	this->m_FixedImage  = 0;
	this->m_FixedImageSamplesList.clear();

	// initalise the transform
	this->m_Transform    = TranslationTransformType::New();

	this->m_Threader = MultiThreader::New();
	this->m_NumberOfThreads = this->m_Threader->GetNumberOfThreads();
}


template< typename TFixedImage, typename TMovingImage >
void
MRFUnaryPotentialsMetricBase<TFixedImage, TMovingImage>::
Initialise()
throw(ExceptionObject &)
{
	if(!m_MovingImage)
	{
		itkExceptionMacro( << " Moving Image is not Assigned");
	}

	if(!m_FixedImage)
	{
		itkExceptionMacro( << " Fixed Image is not Assigned");
	}

	if(m_FixedImageSamplesList.empty())
	{
		itkExceptionMacro( << " Fixed Image Samples are not set");
	}

	if( !m_Interpolator )
	{
		itkExceptionMacro( << " Interpolator is not set");
	}

}

template < typename TFixedImage, typename TMovingImage >
unsigned int
MRFUnaryPotentialsMetricBase<TFixedImage, TMovingImage>::
GetNumberOfValues() const
{
	if(!m_FixedImageSamplesList.empty())
	{
		return m_FixedImageSamplesList.size();
	}
	else
	{
		return 0;
	}
}


template < typename TFixedImage, typename TMovingImage >
typename MRFUnaryPotentialsMetricBase<TFixedImage, TMovingImage>::MeasureType
MRFUnaryPotentialsMetricBase<TFixedImage, TMovingImage>::
GetNonConstValue( const ParametersType & parameters )
{
	// make sure that the values have been set
	try
	{
		this->Initialise();
	}
	catch (ExceptionObject &e)
	{
		throw e;
	}

	// clear the output value and set it up again
	this->m_OutputValue.clear();
	this->m_OutputValue.SetSize(this->m_FixedImageSamplesList.size());


	// set up the transform
	this->m_Transform->SetParameters(parameters);
	this->m_Transform->Modified();

	// Set the interpolator
	this->m_Interpolator->SetInputImage(m_MovingImage);


	// set up the multi threading
	//
	//
	//
	ThreadStruct str;
	str.Metric = this;

	this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
	this->GetMultiThreader()->SetSingleMethod(this->ThreaderCallback, &str);


	// execute the threader
	//
	//
	//
	this->GetMultiThreader()->SingleMethodExecute();


	// return
	return m_OutputValue;
}



template < typename TFixedImage, typename TMovingImage >
typename MRFUnaryPotentialsMetricBase<TFixedImage, TMovingImage>::MeasureType
MRFUnaryPotentialsMetricBase<TFixedImage, TMovingImage>::
GetValue( const ParametersType & parameters ) const
{
	return const_cast<Self*>(this)->GetNonConstValue(parameters);
}



template<typename TFixedImage, typename TMovingImage>
ITK_THREAD_RETURN_TYPE
MRFUnaryPotentialsMetricBase<TFixedImage, TMovingImage>::
ThreaderCallback(void * arg)
{
	ThreadStruct *str;
	int total, threadId, threadCount;

	// Cast out the information about the threader
	//
	//
	//
	threadId    = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
	threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;
	str = (ThreadStruct *)(((MultiThreader::ThreadInfoStruct *)(arg))->UserData);


	// Now we need to split the input samples into the
	// Correct amnount, given the number of threads
	//
	//
	//
	SampleRegion sampleRegion;
	total = str->Metric->SplitFixedSampleList(threadId, threadCount, sampleRegion);


	// Now we can actually calculate the values
	//
	//
	//
	if(threadId < total)
	{
		str->Metric->ComputeUnaryPotentials(sampleRegion, threadId);
	}



	return ITK_THREAD_RETURN_VALUE;


}

template<typename TFixedImage, typename TMovingImage>
int
MRFUnaryPotentialsMetricBase<TFixedImage, TMovingImage>::
SplitFixedSampleList(int i, int num, SampleRegion &sampleRegion)
{
	// Find out the number of peices that the list will be split into
	//
	//
	//
	unsigned int range  = this->m_FixedImageSamplesList.size();
	int valuesPerThread = (int) vcl_ceil( range / (double) num );
	int maxThreadIdUsed  = (int) vcl_ceil( range / (double) valuesPerThread ) - 1;




	// Next we split the data. The last thread has to take on the "rest"
	//
	//
	//
	int threadStartIndex = 0;
	int threadRange = 0;;
	if(i < maxThreadIdUsed)
	{
		threadStartIndex = i*valuesPerThread;
		threadRange      = valuesPerThread;
	}

	if(i == maxThreadIdUsed)
	{
		threadStartIndex = i*valuesPerThread;
		threadRange      = this->GetNumberOfValues() - threadStartIndex;
	}



	// Now we set the values of the sample region
	//
	//
	//
	sampleRegion.threadStartIndex = threadStartIndex;
	sampleRegion.threadRange      = threadRange;


	return maxThreadIdUsed + 1;

}

}



#endif
