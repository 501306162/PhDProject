#include <iostream>

#include <CommonDefinitions.h>
#include <RegistrationHelpers.h>
#include <itkImageRegistrationMethod.h>
#include <itkBSplineTransform.h>
#include <itkRegularStepGradientDescentOptimizer.h>
#include <KappaStatisticImageToImageMetric.h>
#include <itkBSplineTransformInitializer.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>
#include <itkImageRegionIterator.h>

#include <itkMultiResolutionPyramidImageFilter.h>
#include <itkComposeDisplacementFieldsImageFilter.h>
#include <itkWarpImageFilter.h>

#include <MRFImageRegistrationMethod.h>
#include <MRFMeanSquaredUnaryPotentialsMetric.h>

using namespace utils;


int main(int, char ** argv)
{

	typedef unsigned char PixelType;
	const unsigned int Dimensions = 2;
	typedef itk::Image<PixelType, Dimensions> LabelType;
	typedef itk::Image<unsigned short, Dimensions> ImageType;

	std::string fixedImageFilename  = argv[1];
	std::string fixedLabelFilename  = argv[2];
	std::string movingImageFilename = argv[3];
	std::string movingLabelFilename = argv[4];
	std::string outputImageFilename = argv[5];


	// load the images 
	ImageType::Pointer fixedImage  = ImageSliceIO::Read(fixedImageFilename);
	LabelType::Pointer fixedLabel  = LabelSliceIO::Read(fixedLabelFilename);
	ImageType::Pointer movingImage = ImageSliceIO::Read(movingImageFilename);
	LabelType::Pointer movingLabel = LabelSliceIO::Read(movingLabelFilename);


	// do the rigid registration between the label volumes
	Rigid2DTransformType::Pointer rigidTransform = Rigid2DTransformType::New();
	RegisterLabelSliceImages(fixedLabel, movingLabel, rigidTransform);


	// apply the transform to the image and label
	ImageVolume::Pointer rigidlyTransformedImage = ImageVolume::New();
	LabelType::Pointer rigidlyTransformedLabel = LabelType::New();
	//Apply3DRigidTransformToImageVolume(movingImage, fixedImage, rigidTransform, rigidlyTransformedImage);
	Apply2DRigidTransformToLabelSlice(movingLabel, fixedLabel, rigidTransform, rigidlyTransformedLabel);



	const unsigned int numberOfLevels = 3;


	std::vector<int> sampleNumbers;
	sampleNumbers.push_back(1);
	sampleNumbers.push_back(1);
	sampleNumbers.push_back(1);


	std::vector<int> gridSizes;
	gridSizes.push_back(5);
	gridSizes.push_back(10);
	gridSizes.push_back(25);
	

	// resample the images using an image pyramid
	typedef itk::MultiResolutionPyramidImageFilter<LabelType, LabelType> PyramidType;
	PyramidType::Pointer fixedPyramid = PyramidType::New();
	fixedPyramid->SetInput(fixedLabel);
	fixedPyramid->SetNumberOfLevels(numberOfLevels);
	fixedPyramid->Update();


	PyramidType::Pointer movingPyramid = PyramidType::New();
	movingPyramid->SetInput(rigidlyTransformedLabel);
	movingPyramid->SetNumberOfLevels(numberOfLevels);
	movingPyramid->Update();

	typedef itk::BSplineTransform<double, 2, 3> NNRTransformType;
	NNRTransformType::ParametersType finalParameters;



	// create the displacement field 
	typedef itk::Vector<double, 2> DisplacementType;
	typedef itk::Image<DisplacementType, 2> DisplacementImageType;
	DisplacementImageType::Pointer displacement = DisplacementImageType::New();
	displacement->SetOrigin(fixedLabel->GetOrigin());
	displacement->SetSpacing(fixedLabel->GetSpacing());
	displacement->SetDirection(fixedLabel->GetDirection());
	displacement->SetRegions(fixedLabel->GetLargestPossibleRegion());
	displacement->Allocate();

	DisplacementType zeroDisplacement;
	zeroDisplacement.Fill(0.0);

	displacement->FillBuffer(zeroDisplacement);



	std::vector<NNRTransformType::Pointer> transforms;


	for(unsigned int i = 0; i < numberOfLevels; i++)
	{
		LabelType::Pointer fixed  = fixedPyramid->GetOutput(i);
		LabelType::Pointer moving = movingPyramid->GetOutput(i);
	
		// set up the non rigid registration
		NNRTransformType::Pointer nnrTransform = NNRTransformType::New();

		typedef itk::BSplineTransformInitializer< NNRTransformType, LabelType > BSplineInitialiserType;
		BSplineInitialiserType::Pointer initialiser = BSplineInitialiserType::New();
		NNRTransformType::MeshSizeType  meshSize;
		meshSize.Fill(gridSizes[i]);

		initialiser->SetTransform(nnrTransform);
		initialiser->SetImage(fixed);
		initialiser->SetTransformDomainMeshSize(meshSize);
		initialiser->InitializeTransform();

		LabelSliceIO::Write("starting.nrrd", moving);


		for(unsigned int j = 0; j < transforms.size(); j++)
		{
			typedef itk::ResampleImageFilter<LabelType,LabelType>    ResampleFilterType;
			ResampleFilterType::Pointer resample = ResampleFilterType::New();
			resample->SetTransform( transforms[j] );
			resample->SetInput( moving );
			resample->SetSize(    fixed->GetLargestPossibleRegion().GetSize() );
			resample->SetOutputOrigin(  fixed->GetOrigin() );
			resample->SetOutputSpacing( fixed->GetSpacing() );
			resample->SetOutputDirection( fixed->GetDirection() );
			resample->SetDefaultPixelValue( 0 );
			resample->Update();


			moving = resample->GetOutput();
		}



		
	

		std::cout << nnrTransform->GetParameters() << std::endl;



		typedef myitk::MRFImageRegistrationMethod<LabelType, LabelType, NNRTransformType> RegistrationType;
		RegistrationType::Pointer registration = RegistrationType::New();

		typedef myitk::MRFMeanSquaredUnaryPotentialsMetric<LabelType, LabelType> MetricType;
		MetricType::Pointer metric = MetricType::New();

		typedef itk::NearestNeighborInterpolateImageFunction<LabelType, double> InterpolatorType;
		InterpolatorType::Pointer interpolator = InterpolatorType::New();

		registration->SetBSplineTransform(nnrTransform);
		registration->SetFixedImage(fixed);
		registration->SetMovingImage(moving);
		registration->SetInitialTransformParameters(nnrTransform->GetParameters());
		registration->SetLambda(1);
		registration->SetMetric(metric);
		registration->SetInterpolator(interpolator);
		registration->SetLabelSteps(5);
		registration->SetDenseSampling(true);
		registration->SetSamplingRate(sampleNumbers[i]);
		registration->SetOptimiserLevels(5);
		registration->SetMaxLabelDisplacement(0.8);
		registration->SetMaxDisplacementScaleFactor(0.6);
		registration->SetVerbose(true);

		registration->Update();


		transforms.push_back(nnrTransform);



		// add to the displacement
		DisplacementImageType::Pointer newDisplacement = DisplacementImageType::New();
		newDisplacement->SetOrigin(fixedLabel->GetOrigin());
		newDisplacement->SetSpacing(fixedLabel->GetSpacing());
		newDisplacement->SetDirection(fixedLabel->GetDirection());
		newDisplacement->SetRegions(fixedLabel->GetLargestPossibleRegion());
		newDisplacement->Allocate();


		typedef itk::ImageRegionIterator<DisplacementImageType> DispIteratorType;
		DispIteratorType it(newDisplacement, newDisplacement->GetLargestPossibleRegion());

		while(!it.IsAtEnd())
		{
			DisplacementImageType::IndexType index = it.GetIndex();
			DisplacementImageType::PointType point, tpoint;

			newDisplacement->TransformIndexToPhysicalPoint(index, point);
			tpoint = nnrTransform->TransformPoint(point);

			DisplacementType disp;
			disp[0] = tpoint[0] - point[0];
			disp[1] = tpoint[1] - point[1];

			it.Set(disp);

			++it;
		}

		typedef itk::ComposeDisplacementFieldsImageFilter<DisplacementImageType, DisplacementImageType> ComposerType;
		ComposerType::Pointer composer = ComposerType::New();
		composer->SetDisplacementField(newDisplacement);
		composer->SetWarpingField(displacement);
		composer->Update();

		displacement = composer->GetOutput();







		typedef itk::ResampleImageFilter<LabelType,LabelType>    ResampleFilterType;
		ResampleFilterType::Pointer resample = ResampleFilterType::New();
		resample->SetTransform( nnrTransform );
		resample->SetInput( moving );
		resample->SetSize(    fixed->GetLargestPossibleRegion().GetSize() );
		resample->SetOutputOrigin(  fixed->GetOrigin() );
		resample->SetOutputSpacing( fixed->GetSpacing() );
		resample->SetOutputDirection( fixed->GetDirection() );
		resample->SetDefaultPixelValue( 0 );
		resample->Update();

		LabelSliceIO::Write("warped.nrrd", resample->GetOutput());

		std::stringstream ss;
		ss << i;

		ImageSliceIO::Write("fixed" + ss.str() + ".nrrd", fixedImage);
		ImageSliceIO::Write("moving" + ss.str() + ".nrrd", movingImage);
		//ImageSliceIO::Write("rigidTransformImage" + ss.str() + ".nrrd", rigidlyTransformedImage);
		LabelSliceIO::Write("fixedLabel" + ss.str() + ".nrrd", fixed);
		LabelSliceIO::Write("movingLabel" + ss.str() + ".nrrd", moving);
		LabelSliceIO::Write("rigidTransformLabel" + ss.str() + ".nrrd", rigidlyTransformedLabel);



	}

	typedef itk::NearestNeighborInterpolateImageFunction<LabelSlice, double> InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();

	// apply the displacement 
	typedef itk::WarpImageFilter<LabelSlice, LabelSlice, DisplacementImageType> WarperType;
	WarperType::Pointer warper = WarperType::New();
	warper->SetInput(rigidlyTransformedLabel);
	warper->SetDisplacementField(displacement);
	warper->SetInterpolator(interpolator);
	warper->SetOutputParametersFromImage(fixedLabel);
	warper->SetEdgePaddingValue(0);
	warper->Update();

	LabelSliceIO::Write("final.nrrd", warper->GetOutput());





	/*
	
	typedef itk::ImageRegistrationMethod<LabelVolume, LabelVolume> RegistrationType;
	RegistrationType::Pointer registration = RegistrationType::New();

	typedef registration::KappaStatisticImageToImageMetric<LabelVolume, LabelVolume> MetricType;
	MetricType::Pointer metric = MetricType::New();
	metric->SetForegroundValue(1);
	metric->SetComplement( true );

	typedef itk::RegularStepGradientDescentOptimizer OptimiserType;
	OptimiserType::Pointer optimiser = OptimiserType::New();
	optimiser->SetMaximumStepLength( 10 );
	optimiser->SetMinimumStepLength( 0.001 );
	optimiser->SetNumberOfIterations( 300 );

	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  	optimiser->AddObserver( itk::IterationEvent(), observer );




	registration->SetFixedImage(fixedLabel);
	registration->SetMovingImage(rigidlyTransformedLabel);
	registration->SetTransform(nnrTransform);
	registration->SetMetric(metric);
	registration->SetInterpolator(interpolator);
	registration->SetOptimizer(optimiser);
	registration->SetFixedImageRegion(fixedLabel->GetLargestPossibleRegion());

	registration->SetInitialTransformParameters(nnrTransform->GetParameters());

	registration->Update();


	std::cout << nnrTransform->GetParameters() << std::endl;

	
	*/

	


	return 0;
}
