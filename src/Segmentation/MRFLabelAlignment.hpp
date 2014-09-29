#ifndef MRF_LABEL_ALIGNMENT_HPP
#define MRF_LABEL_ALIGNMENT_HPP

#include "MRFLabelAlignment.h"
#include <itkCenteredTransformInitializer.h>
#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkImageRegistrationMethod.h>
#include "KappaStatisticImageToImageMetric.h"
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>

namespace segmentation
{
// ------------------------------------------------------------------------
template<typename TFixedType, typename TMovingType, typename TTransformType>
MRFLabelAlignment<TFixedType, TMovingType, TTransformType>::
MRFLabelAlignment()
{
	this->SetNumberOfRequiredInputs(1);
}

// ------------------------------------------------------------------------
template<typename TFixedType, typename TMovingType, typename TTransformType>
void
MRFLabelAlignment<TFixedType, TMovingType, TTransformType>::
SetMovingImage(const MovingType * moving)
{
	m_MovingImage = const_cast<MovingType*>(moving);
}


// ------------------------------------------------------------------------
template<typename TFixedType, typename TMovingType, typename TTransformType>
void
MRFLabelAlignment<TFixedType, TMovingType, TTransformType>::
SetFixedImage(const FixedType * fixed)
{
	this->SetNthInput(0, const_cast<MovingType*>(fixed));
}

// ------------------------------------------------------------------------
template<typename TFixedType, typename TMovingType, typename TTransformType>
void
MRFLabelAlignment<TFixedType, TMovingType, TTransformType>::
GenerateData()
{
	// compute the transform
	m_Transform = TransformType::New();
	ComputeTransform();
	ResampleImage();

}

// ------------------------------------------------------------------------
template<typename TFixedType, typename TMovingType, typename TTransformType>
void
MRFLabelAlignment<TFixedType, TMovingType, TTransformType>::
ResampleImage()
{
	typedef itk::ResampleImageFilter<TMovingType, TMovingType> ResamplerType;
	typedef itk::NearestNeighborInterpolateImageFunction<TMovingType, double> InterpolatorType;

	typename ResamplerType::Pointer resampler = ResamplerType::New();
	typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

	resampler->SetTransform( m_Transform );
	resampler->SetInterpolator( interpolator );
	resampler->SetInput( m_MovingImage );
	resampler->SetOutputSpacing( this->GetInput(0)->GetSpacing() );
	resampler->SetOutputDirection( this->GetInput(0)->GetDirection() );
	resampler->SetOutputOrigin( this->GetInput(0)->GetOrigin() );
	resampler->SetSize( this->GetInput(0)->GetLargestPossibleRegion().GetSize() );
	resampler->SetDefaultPixelValue(0);
	resampler->Update();

	this->GraftOutput(resampler->GetOutput());	
}

// ------------------------------------------------------------------------
template<typename TFixedType, typename TMovingType, typename TTransformType>
typename MRFLabelAlignment<TFixedType, TMovingType, TTransformType>::TransformPointer
MRFLabelAlignment<TFixedType, TMovingType, TTransformType>::
GetTransform() const
{
	return m_Transform;

}


// ------------------------------------------------------------------------
template<typename TFixedType, typename TMovingType, typename TTransformType>
void
MRFLabelAlignment<TFixedType, TMovingType, TTransformType>::
ComputeTransform()
{
	// intialise the transform
	typedef itk::CenteredTransformInitializer<TransformType, FixedType, MovingType> InitialiserType;
	typename InitialiserType::Pointer initialiser = InitialiserType::New();
	initialiser->SetMovingImage(m_MovingImage);
	initialiser->SetFixedImage(this->GetInput(0));
	initialiser->SetTransform(m_Transform);
	initialiser->MomentsOn();
	initialiser->InitializeTransform();
	

	// set up the registration
	// 
	// Now we set up the registration
	//  
	typedef itk::RegularStepGradientDescentOptimizer OptimiserType;
	typedef itk::ImageRegistrationMethod< FixedType, MovingType > RegistrationType;
	typedef itk::NearestNeighborInterpolateImageFunction< FixedType, double > InterpolatorType;
	typedef registration::KappaStatisticImageToImageMetric< FixedType, MovingType > MetricType;

	OptimiserType::Pointer optimiser = OptimiserType::New();
	typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
	typename RegistrationType::Pointer registration = RegistrationType::New();
	typename MetricType::Pointer metric = MetricType::New();

	registration->SetOptimizer( optimiser );
	registration->SetMetric( metric );
	registration->SetTransform( m_Transform );
	registration->SetInterpolator( interpolator );
	registration->SetFixedImage( this->GetInput(0) );
	registration->SetMovingImage(m_MovingImage );
	registration->SetFixedImageRegion( this->GetInput(0)->GetLargestPossibleRegion() );
	registration->SetInitialTransformParameters( m_Transform->GetParameters() );


	// set up the metric
	metric->SetForegroundValue(1);
	metric->SetComplement( true );


	// set up the optimiser
	typedef OptimiserType::ScalesType OptimiserScaleType;
	OptimiserScaleType optimiserScales ( m_Transform->GetNumberOfParameters() );
	const double translationScale = 1.0 / 500.0;
	const double scaleScale = 0.1;
	const double rotationScale = 0.1;

	optimiserScales[0] = rotationScale;
	optimiserScales[1] = rotationScale;
	optimiserScales[2] = rotationScale;
	optimiserScales[3] = translationScale;
	optimiserScales[4] = translationScale;
	optimiserScales[5] = translationScale;
	optimiserScales[6] = scaleScale;

	optimiser->SetScales( optimiserScales );
	optimiser->SetMaximumStepLength( 0.05 );
	optimiser->SetMinimumStepLength( 0.0001 );
	optimiser->SetNumberOfIterations( 300 );


	try
	{
		registration->Update();
	}
	catch ( itk::ExceptionObject &e )
	{
		std::cout << "Error in registration: " << std::endl;
		std::cout << e << std::endl;
		exit(1);
	}

	m_Transform->SetParameters( registration->GetLastTransformParameters() );

}


// ------------------------------------------------------------------------
template<typename TFixedType, typename TMovingType, typename TTransformType>
typename MRFLabelAlignment<TFixedType, TMovingType, TTransformType>::ParametersType
MRFLabelAlignment<TFixedType, TMovingType, TTransformType>::
GetFinalParameters() const
{
	return m_Transform->GetParameters();
}
} /* segmentation */ 

#endif
