#include "RegistrationHelpers.h"

#include "KappaStatisticImageToImageMetric.h"
#include <itkCenteredTransformInitializer.h>
#include <itkResampleImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkImageRegistrationMethod.h>

namespace utils
{

// ------------------------------------------------------------------------
void RegisterLabelVolumeImages(
		LabelVolume::Pointer &fixed, 
		LabelVolume::Pointer &moving,
		Rigid3DTransformType::Pointer &transform)
{
	// first we initialise the process using the centered transform initialiser
	typedef itk::CenteredTransformInitializer<Rigid3DTransformType, LabelVolume, LabelVolume> InitialiserType;
	InitialiserType::Pointer initialiser = InitialiserType::New();
	initialiser->SetFixedImage(fixed);
	initialiser->SetMovingImage(moving);
	initialiser->SetTransform(transform);
	initialiser->MomentsOn();

	initialiser->InitializeTransform();

		
	// set up the registration
// 
	// Now we set up the registration
	//  
	typedef itk::RegularStepGradientDescentOptimizer OptimiserType;
	typedef itk::ImageRegistrationMethod< LabelVolume, LabelVolume > RegistrationType;
	typedef itk::NearestNeighborInterpolateImageFunction< LabelVolume, double > InterpolatorType;
	typedef registration::KappaStatisticImageToImageMetric< LabelVolume, LabelVolume > MetricType;

	OptimiserType::Pointer optimiser = OptimiserType::New();
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	RegistrationType::Pointer registration = RegistrationType::New();
	MetricType::Pointer metric = MetricType::New();

	registration->SetOptimizer( optimiser );
	registration->SetMetric( metric );
	registration->SetTransform( transform );
	registration->SetInterpolator( interpolator );
	registration->SetFixedImage( fixed );
	registration->SetMovingImage( moving );
	registration->SetFixedImageRegion( fixed->GetLargestPossibleRegion() );
	registration->SetInitialTransformParameters( transform->GetParameters() );


	// set up the metric
	metric->SetForegroundValue(1);
	metric->SetComplement( true );


	// set up the optimiser
	typedef OptimiserType::ScalesType OptimiserScaleType;
	OptimiserScaleType optimiserScales ( transform->GetNumberOfParameters() );
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

	transform->SetParameters( registration->GetLastTransformParameters() );
	
}

// ------------------------------------------------------------------------
void RegisterLabelSliceImages(
		LabelSlice::Pointer &fixed, 
		LabelSlice::Pointer &moving,
		Rigid2DTransformType::Pointer &transform)
{
	// first we initialise the process using the centered transform initialiser
	typedef itk::CenteredTransformInitializer<Rigid2DTransformType, LabelSlice, LabelSlice> InitialiserType;
	InitialiserType::Pointer initialiser = InitialiserType::New();
	initialiser->SetFixedImage(fixed);
	initialiser->SetMovingImage(moving);
	initialiser->SetTransform(transform);
	initialiser->MomentsOn();

	initialiser->InitializeTransform();

		
	// set up the registration
// 
	// Now we set up the registration
	//  
	typedef itk::RegularStepGradientDescentOptimizer OptimiserType;
	typedef itk::ImageRegistrationMethod< LabelSlice, LabelSlice > RegistrationType;
	typedef itk::NearestNeighborInterpolateImageFunction< LabelSlice, double > InterpolatorType;
	typedef registration::KappaStatisticImageToImageMetric< LabelSlice, LabelSlice > MetricType;

	OptimiserType::Pointer optimiser = OptimiserType::New();
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	RegistrationType::Pointer registration = RegistrationType::New();
	MetricType::Pointer metric = MetricType::New();

	registration->SetOptimizer( optimiser );
	registration->SetMetric( metric );
	registration->SetTransform( transform );
	registration->SetInterpolator( interpolator );
	registration->SetFixedImage( fixed );
	registration->SetMovingImage( moving );
	registration->SetFixedImageRegion( fixed->GetLargestPossibleRegion() );
	registration->SetInitialTransformParameters( transform->GetParameters() );


	// set up the metric
	metric->SetForegroundValue(1);
	metric->SetComplement( true );


	// set up the optimiser
	typedef OptimiserType::ScalesType OptimiserScaleType;
	OptimiserScaleType optimiserScales ( transform->GetNumberOfParameters() );
	const double translationScale = 1.0 / 500.0;
	const double scaleScale = 0.1;
	const double rotationScale = 0.1;

	optimiserScales[0] = scaleScale;
	optimiserScales[1] = rotationScale;
	optimiserScales[2] = translationScale;
	optimiserScales[3] = translationScale;
	
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

	transform->SetParameters( registration->GetLastTransformParameters() );
	
}

// ------------------------------------------------------------------------
void Apply2DRigidTransformToLabelSlice(
		LabelSlice::Pointer &input,
		LabelSlice::Pointer &reference,
		Rigid2DTransformType::Pointer &transform,
		LabelSlice::Pointer &output)
{

	typedef itk::ResampleImageFilter<LabelSlice, LabelSlice> ResamplerType;
	typedef itk::NearestNeighborInterpolateImageFunction<LabelSlice, double> InterpolatorType;

	ResamplerType::Pointer resampler = ResamplerType::New();
	InterpolatorType::Pointer interpolator = InterpolatorType::New();

	resampler->SetTransform( transform );
	resampler->SetInterpolator( interpolator );
	resampler->SetInput( input );
	resampler->SetOutputSpacing( reference->GetSpacing() );
	resampler->SetOutputDirection( reference->GetDirection() );
	resampler->SetOutputOrigin( reference->GetOrigin() );
	resampler->SetSize( reference->GetLargestPossibleRegion().GetSize() );
	resampler->SetDefaultPixelValue(0);
	resampler->Update();

	output = resampler->GetOutput();


}


// ------------------------------------------------------------------------
void Apply3DRigidTransformTOLabelVolume(
		LabelVolume::Pointer &input,
		LabelVolume::Pointer &reference,
		Rigid3DTransformType::Pointer &transform,
		LabelVolume::Pointer &output)
{
	

	typedef itk::ResampleImageFilter<LabelVolume, LabelVolume> ResamplerType;
	typedef itk::NearestNeighborInterpolateImageFunction<LabelVolume, double> InterpolatorType;

	ResamplerType::Pointer resampler = ResamplerType::New();
	InterpolatorType::Pointer interpolator = InterpolatorType::New();

	resampler->SetTransform( transform );
	resampler->SetInterpolator( interpolator );
	resampler->SetInput( input );
	resampler->SetOutputSpacing( reference->GetSpacing() );
	resampler->SetOutputDirection( reference->GetDirection() );
	resampler->SetOutputOrigin( reference->GetOrigin() );
	resampler->SetSize( reference->GetLargestPossibleRegion().GetSize() );
	resampler->SetDefaultPixelValue(0);
	resampler->Update();

	output = resampler->GetOutput();
}


// ------------------------------------------------------------------------
void Apply3DRigidTransformToImageVolume(
		ImageVolume::Pointer &input,
		ImageVolume::Pointer &reference,
		Rigid3DTransformType::Pointer &transform,
		ImageVolume::Pointer &output)
{
	

	typedef itk::ResampleImageFilter<ImageVolume, ImageVolume> ResamplerType;
	typedef itk::LinearInterpolateImageFunction<ImageVolume, double> InterpolatorType;

	ResamplerType::Pointer resampler = ResamplerType::New();
	InterpolatorType::Pointer interpolator = InterpolatorType::New();

	resampler->SetTransform( transform );
	resampler->SetInterpolator( interpolator );
	resampler->SetInput( input );
	resampler->SetOutputSpacing( reference->GetSpacing() );
	resampler->SetOutputDirection( reference->GetDirection() );
	resampler->SetOutputOrigin( reference->GetOrigin() );
	resampler->SetSize( reference->GetLargestPossibleRegion().GetSize() );
	resampler->SetDefaultPixelValue(0);
	resampler->Update();

	output = resampler->GetOutput();
}

}
