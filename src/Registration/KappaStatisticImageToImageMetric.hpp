/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkKappaStatisticImageToImageMetric_hxx
#define __itkKappaStatisticImageToImageMetric_hxx

#include "KappaStatisticImageToImageMetric.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace registration
{
using namespace itk;
/**
 * Constructor
 */
template <class TFixedImage, class TMovingImage>
KappaStatisticImageToImageMetric<TFixedImage, TMovingImage>
::KappaStatisticImageToImageMetric()
{
	itkDebugMacro("Constructor");

	this->SetComputeGradient(true);
	this->m_WithinThreadPreProcess = false;
	this->m_WithinThreadPostProcess = false;
	m_PerThread = NULL;

	m_ForegroundValue = 255;
	m_Complement = false;
}

// ---------------------------------------------------------------------------------
template <class TFixedImage, class TMovingImage>
void
KappaStatisticImageToImageMetric<TFixedImage, TMovingImage>
::Initialize()
throw( ExceptionObject )
{
	this->Superclass::Initialize();
	this->Superclass::MultiThreadingInitialize();

	delete[] m_PerThread;

	this->m_PerThread = new AlignedPerThreadType[this->GetNumberOfThreads()];


}

/**
 * Get the match Measure
 */
template <class TFixedImage, class TMovingImage>
typename KappaStatisticImageToImageMetric<TFixedImage, TMovingImage>::MeasureType
KappaStatisticImageToImageMetric<TFixedImage, TMovingImage>
::GetValue(const TransformParametersType & parameters) const
{
	itkDebugMacro("GetValue( " << parameters << " ) ");

	this->SetTransformParameters(parameters);

	// Get the fixed image
	//
	//
	FixedImageConstPointer fixedImage = this->m_FixedImage;
	if( !fixedImage )
	{
		itkExceptionMacro(<< "Fixed image has not been assigned");
	}

	// Get the moving image
	//
	//
	MovingImageConstPointer movingImage = this->m_MovingImage;
	if( !movingImage )
	{
		itkExceptionMacro(<< "Moving image has not been assigned");
	}

	

	this->m_Transform->SetParameters( parameters );
	

	// initailise the thread holder
	for( ThreadIdType id = 0; id < this->GetNumberOfThreads(); id++ )
	{
		m_PerThread[id].m_MovingForeground = 0;
		m_PerThread[id].m_FixedForeground = 0;
		m_PerThread[id].m_Intersection = 0;
		m_PerThread[id].m_Sum1 = ArrayType( this->GetNumberOfParameters() );
		m_PerThread[id].m_Sum1.Fill(0);
		m_PerThread[id].m_Sum2 = ArrayType( this->GetNumberOfParameters() );
		m_PerThread[id].m_Sum2.Fill(0);
	}

	this->GetValueMultiThreadedInitiate();


	MeasureType measure;

	int fixedArea = 0;
	int movingArea = 0;
	int intersection = 0;
	for( unsigned int i = 0; i < this->GetNumberOfThreads(); i++ )
	{
		fixedArea += m_PerThread[i].m_FixedForeground;
		movingArea += m_PerThread[i].m_MovingForeground;
		intersection += m_PerThread[i].m_Intersection;
	}

	// Compute the final metric value
	//
	//
	if( !m_Complement )
	{
		measure = 2.0 * ( intersection ) / ( fixedArea + movingArea );
	}
	else
	{
		measure = 1.0 - 2.0 * ( intersection ) / ( fixedArea + movingArea );
	}

	return measure;
}

// ---------------------------------------------------------------------------------
template <class TFixedImage, class TMovingImage>
bool
KappaStatisticImageToImageMetric<TFixedImage, TMovingImage>
::GetValueThreadProcessSample(ThreadIdType threadID,
			SizeValueType fixedImageSample,
			const MovingImagePointType & mappedPoint,
			double movingImageValue) const
{

	FixedImageSample fixedSample = this->m_FixedImageSamples[fixedImageSample];
	
	// get the fixed image point
	InputPointType fixedInputPoint = fixedSample.point;
	


	// check if it is valid ( if there is a mask )
	if( this->m_FixedImageMask && !this->m_FixedImageMask->IsInside(fixedInputPoint) )
	{
		return true;
	}

	// get the fixed value
	const RealType fixedValue = fixedSample.value;

	// Increment 'fixedForegroundArea'
	//
	//
	if( fixedValue == m_ForegroundValue )
	{
		m_PerThread[threadID].m_FixedForeground++;
	}

	// now we check the moving point
	if( this->m_MovingImageMask && !this->m_MovingImageMask->IsInside(mappedPoint) )
	{
		return true;
	}

	// Compute movingForegroundArea and intersection
	//
	//
	if( static_cast<MovingImagePixelType>(movingImageValue) == m_ForegroundValue )
	{
		m_PerThread[threadID].m_MovingForeground++;
		if( fixedValue == m_ForegroundValue )
		{
			m_PerThread[threadID].m_Intersection++;
		}
	}

	return true;
}


// ---------------------------------------------------------------------------------
template <class TFixedImage, class TMovingImage>
bool
KappaStatisticImageToImageMetric<TFixedImage, TMovingImage>
::GetValueAndDerivativeThreadProcessSample(ThreadIdType threadID,
			SizeValueType fixedImageSample,
			const MovingImagePointType & mappedPoint,
			double movingImageValue,
			const ImageDerivativesType &
			movingImageGradientValue) const
{

	FixedImageSample fixedSample = this->m_FixedImageSamples[fixedImageSample];
	
	// get the fixed image point
	InputPointType fixedInputPoint = fixedSample.point;
	


	// check if it is valid ( if there is a mask )
	if( this->m_FixedImageMask && !this->m_FixedImageMask->IsInside(fixedInputPoint) )
	{
		return true;
	}

	// get the fixed value
	const RealType fixedValue = fixedSample.value;

	// Increment 'fixedForegroundArea'
	//
	//
	if( fixedValue == m_ForegroundValue )
	{
		m_PerThread[threadID].m_FixedForeground++;
	}

	// now we check the moving point
	if( this->m_MovingImageMask && !this->m_MovingImageMask->IsInside(mappedPoint) )
	{
		return true;
	}

	// Compute movingForegroundArea and intersection
	//
	//
	if( static_cast<MovingImagePixelType>(movingImageValue) == m_ForegroundValue )
	{
		m_PerThread[threadID].m_MovingForeground++;
		if( fixedValue == m_ForegroundValue )
		{
			m_PerThread[threadID].m_Intersection++;
		}
	}


	TransformJacobianType jacobian;
	this->m_Transform->ComputeJacobianWithRespectToParameters(
			fixedInputPoint, jacobian);

	this->m_NumberOfPixelsCounted++;

	// Get the gradient by NearestNeighboorInterpolation:
	// which is equivalent to round up the point components.
	typedef typename OutputPointType::CoordRepType CoordRepType;
	typedef ContinuousIndex<CoordRepType, MovingImageType::ImageDimension>
		MovingImageContinuousIndexType;

	MovingImageContinuousIndexType tempIndex;
	this->m_MovingImage->TransformPhysicalPointToContinuousIndex(mappedPoint, tempIndex);

	typename MovingImageType::IndexType mappedIndex;
	mappedIndex.CopyWithRound(tempIndex);

	unsigned int ParametersDimension = this->GetNumberOfParameters();
	unsigned int ImageDimension = FixedImageType::ImageDimension;
	for( unsigned int par = 0; par < ParametersDimension; par++ )
	{
		for( unsigned int dim = 0; dim < ImageDimension; dim++ )
		{
			m_PerThread[threadID].m_Sum2[par] += jacobian(dim, par) * movingImageGradientValue[dim];
			if( fixedValue == m_ForegroundValue )
			{
				m_PerThread[threadID].m_Sum1[par] += 2.0 * jacobian(dim, par) * movingImageGradientValue[dim];
			}
		}
	}


	return true;
}







/**
 * Get the Derivative Measure
 */
template <class TFixedImage, class TMovingImage>
void
KappaStatisticImageToImageMetric<TFixedImage, TMovingImage>
::GetDerivative(const TransformParametersType & parameters,
		DerivativeType & derivative) const
{
	MeasureType m;
	this->GetValueAndDerivative( parameters, m, derivative );
}

/*
 * Compute the image gradient and assign to m_GradientImage.
 */
template <class TFixedImage, class TMovingImage>
void
KappaStatisticImageToImageMetric<TFixedImage, TMovingImage>
::ComputeGradient()
{
	const unsigned int dim = MovingImageType::ImageDimension;

	typename GradientImageType::Pointer tempGradientImage = GradientImageType::New();
	tempGradientImage->SetRegions( this->m_MovingImage->GetBufferedRegion().GetSize() );
	tempGradientImage->Allocate();
	tempGradientImage->Update();

	typedef  ImageRegionIteratorWithIndex<GradientImageType>    GradientIteratorType;
	typedef  ImageRegionConstIteratorWithIndex<MovingImageType> MovingIteratorType;

	GradientIteratorType git( tempGradientImage, tempGradientImage->GetBufferedRegion() );
	MovingIteratorType   mit( this->m_MovingImage, this->m_MovingImage->GetBufferedRegion() );

	git.GoToBegin();
	mit.GoToBegin();

	typename MovingImageType::IndexType minusIndex;
	typename MovingImageType::IndexType plusIndex;
	typename MovingImageType::IndexType currIndex;
	typename GradientImageType::PixelType tempGradPixel;
	typename MovingImageType::SizeType movingSize = this->m_MovingImage->GetBufferedRegion().GetSize();
	while( !mit.IsAtEnd() )
	{
		currIndex = mit.GetIndex();
		minusIndex = mit.GetIndex();
		plusIndex = mit.GetIndex();
		for( unsigned int i = 0; i < dim; i++ )
		{
			if( ( currIndex[i] == 0 )
					|| ( static_cast<typename MovingImageType::SizeType::SizeValueType>( currIndex[i] ) == ( movingSize[i] - 1 ) ) )
			{
				tempGradPixel[i] = 0;
			}
			else
			{
				minusIndex[i] = currIndex[i] - 1;
				plusIndex[i] = currIndex[i] + 1;
				double minusVal = double( this->m_MovingImage->GetPixel(minusIndex) );
				double plusVal  = double( this->m_MovingImage->GetPixel(plusIndex) );
				if( ( minusVal != m_ForegroundValue ) && ( plusVal == m_ForegroundValue ) )
				{
					tempGradPixel[i] = 1;
				}
				else if( ( minusVal == m_ForegroundValue ) && ( plusVal != m_ForegroundValue ) )
				{
					tempGradPixel[i] = -1;
				}
				else
				{
					tempGradPixel[i] = 0;
				}
			}
			minusIndex = currIndex;
			plusIndex  = currIndex;
		}
		git.Set(tempGradPixel);
		++git;
		++mit;
	}

	this->m_GradientImage = tempGradientImage;
}

/**
 * Get both the match Measure and theDerivative Measure
 */
template <class TFixedImage, class TMovingImage>
void
KappaStatisticImageToImageMetric<TFixedImage, TMovingImage>
::GetValueAndDerivative(const TransformParametersType & parameters,
		MeasureType & Value, DerivativeType  & Derivative) const
{

	itkDebugMacro("GetValue( " << parameters << " ) ");

	this->SetTransformParameters(parameters);

	// Get the fixed image
	//
	//
	FixedImageConstPointer fixedImage = this->m_FixedImage;
	if( !fixedImage )
	{
		itkExceptionMacro(<< "Fixed image has not been assigned");
	}

	// Get the moving image
	//
	//
	MovingImageConstPointer movingImage = this->m_MovingImage;
	if( !movingImage )
	{
		itkExceptionMacro(<< "Moving image has not been assigned");
	}

	
	Derivative = DerivativeType( this->m_NumberOfParameters );
	Derivative.Fill(0.0);

	this->m_Transform->SetParameters( parameters );
	

	// initailise the thread holder
	for( ThreadIdType id = 0; id < this->GetNumberOfThreads(); id++ )
	{
		m_PerThread[id].m_MovingForeground = 0;
		m_PerThread[id].m_FixedForeground = 0;
		m_PerThread[id].m_Intersection = 0;
		m_PerThread[id].m_Sum1 = ArrayType( this->GetNumberOfParameters() );
		m_PerThread[id].m_Sum1.Fill(0);
		m_PerThread[id].m_Sum2 = ArrayType( this->GetNumberOfParameters() );
		m_PerThread[id].m_Sum2.Fill(0);
	}


	this->GetValueAndDerivativeMultiThreadedInitiate();

	MeasureType measure;

	int fixedArea = 0;
	int movingArea = 0;
	int intersection = 0;


	ArrayType sum1(this->GetNumberOfParameters());
	sum1.Fill(0);
	ArrayType sum2(this->GetNumberOfParameters());
	sum2.Fill(0);

	for( unsigned int i = 0; i < this->GetNumberOfThreads(); i++ )
	{
		fixedArea += m_PerThread[i].m_FixedForeground;
		movingArea += m_PerThread[i].m_MovingForeground;
		intersection += m_PerThread[i].m_Intersection;
		sum1 += m_PerThread[i].m_Sum1;
		sum2 += m_PerThread[i].m_Sum2;
	}

	// Compute the final metric value
	//
	//
	if( !m_Complement )
	{
		measure = 2.0 * ( intersection ) / ( fixedArea + movingArea );
	}
	else
	{
		measure = 1.0 - 2.0 * ( intersection ) / ( fixedArea + movingArea );
	}

	Value = measure;

	double areaSum = double(fixedArea) + double(movingArea);
	unsigned int ParametersDimension = this->GetNumberOfParameters();
	for( unsigned int par = 0; par < ParametersDimension; par++ )
	{
		Derivative[par] = -( areaSum * sum1[par] - 2.0 * intersection * sum2[par] ) / ( areaSum * areaSum );
	}
}

/**
 * PrintSelf
 */
template <class TFixedImage, class TMovingImage>
void
KappaStatisticImageToImageMetric<TFixedImage, TMovingImage>
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
	Superclass::PrintSelf(os, indent);
	os << indent << "Complement: "         << ( m_Complement ? "On" : "Off" )  << std::endl;
	os << indent << "ForegroundValue: "    << m_ForegroundValue << std::endl;
}

} // end namespace itk

#endif
