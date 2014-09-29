#ifndef MYITKMRFIMAGEREGISTRATIONMETHOD_HPP
#define MYITKMRFIMAGEREGISTRATIONMETHOD_HPP

#include "MRFImageRegistrationMethod.h"
#include "MRFImageRegistrationMethodLabelGenerator.h"

#include <itkResampleImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkNrrdImageIO.h>

//#include <Fast_PD.h>

#include <mrf.h>
#include <GCoptimization.h>

namespace myitk
{
// ---------------------------------------------------------------------------------
template <typename TFixedImage, typename TMovingImage, typename TTransform >
MRFImageRegistrationMethod< TFixedImage, TMovingImage, TTransform >
::MRFImageRegistrationMethod()
{
	// set some of the initial values
	this->SetNumberOfRequiredOutputs(1);

	m_FixedImage = 0;
	m_MovingImage = 0;
	m_BSplineTransform = 0;
	m_Interpolator = 0;
	m_Metric = 0;
	m_Verbose = false;
	m_FixedRegionSet = false;
	m_SamplingRate = 1;

	m_MaxLabelDisplacement = 0.4;
	m_OptimiserLevels = 5;
	m_MaxDisplacementScaleFactor = 0.4;
	m_LabelSteps = 5;
	m_DenseSampling = false;
	m_Lambda = 1.0;

	m_InitialTransformParameters = ParametersType(1);
	m_InitialTransformParameters.Fill(0.0);
	m_LastTransformParameters    = ParametersType(1);
	m_LastTransformParameters.Fill(0.0);

	TransformOutputPointer transformDecorator = 
		itkDynamicCastInDebugMode< TransformOutputType * >( this->MakeOutput(0).GetPointer() );
	this->ProcessObject::SetNthOutput(0, transformDecorator.GetPointer() );
	this->SetNumberOfThreads( this->GetMultiThreader()->GetNumberOfThreads() );
}

// ---------------------------------------------------------------------------------
template< typename TFixedImage, typename TMovingImage, typename TTransform >
void
MRFImageRegistrationMethod< TFixedImage, TMovingImage, TTransform >
::SetFixedImageRegion(const FixedImageRegionType & region)
{
	m_FixedImageRegion = region;
	m_FixedRegionSet = true;
}

// ---------------------------------------------------------------------------------
template< typename TFixedImage, typename TMovingImage, typename TTransform >
void
MRFImageRegistrationMethod< TFixedImage, TMovingImage, TTransform >
::SetInitialTransformParameters(const ParametersType & param)
{
	m_InitialTransformParameters = param;
	this->Modified();
}

// ---------------------------------------------------------------------------------
template< typename TFixedImage, typename TMovingImage, typename TTransform>
void
MRFImageRegistrationMethod< TFixedImage, TMovingImage, TTransform >
::Initialize()
throw ( ExceptionObject )
{
	if ( !m_FixedImage )
	{
		itkExceptionMacro(<< "FixedImage is not present");
	}

	if ( !m_MovingImage )
	{
		itkExceptionMacro(<< "MovingImage is not present");
	}

	if ( !m_Metric )
	{
		itkExceptionMacro(<< "Metric is not present");
	}

	if ( !m_BSplineTransform )
	{
		itkExceptionMacro(<< "Transform is not present");
	}

	//
	// Connect the transform to the Decorator.
	//
	TransformOutputType *transformOutput =
	static_cast< TransformOutputType * >( this->ProcessObject::GetOutput(0) );

	transformOutput->Set( m_BSplineTransform.GetPointer() );

	if ( !m_Interpolator )
	{
		itkExceptionMacro(<< "Interpolator is not present");
	}

	// 
	// Compute the connections
	// 
	if( m_Verbose )
	{
		std::cout << "Computing Connections" << std::endl;
	}

	GridImagePointer gridImage = m_BSplineTransform->GetCoefficientImages()[0];
	typename ConnectivityFilterType::Pointer connectionFilter = ConnectivityFilterType::New();
	connectionFilter->SetGridImage( gridImage );
	m_Connections = connectionFilter->ComputeConnections();




	// 
	// compute the fixed image samples
	// 
	if( !m_FixedRegionSet )
	{
		m_FixedImageRegion = m_FixedImage->GetLargestPossibleRegion();
	}


	if( m_Verbose )
	{
		std::cout << "Computing Grid Influences" << std::endl;
	}
	this->ComputeFixedImageSampleList( gridImage );

	//
	// Set up the metric
	//
	this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
	this->m_Metric->SetNumberOfThreads( this->GetNumberOfThreads() );
	m_Metric->SetMovingImage(m_MovingImage);
	m_Metric->SetFixedImage(m_FixedImage);
	m_Metric->SetInterpolator(m_Interpolator);
	m_Metric->SetFixedImageSamplesList( m_SampleList );



	// Validate initial transform parameters
	if ( m_InitialTransformParameters.Size() != m_BSplineTransform->GetNumberOfParameters() )
	{
		itkExceptionMacro(<< "Size mismatch between initial parameters and transform."
				<< "Expected " << m_BSplineTransform->GetNumberOfParameters() << " parameters and received "
				<<  m_InitialTransformParameters.Size() << " parameters");
	}

}

// ---------------------------------------------------------------------------------
template< typename TFixedImage, typename TMovingImage, typename TTransform >
void
MRFImageRegistrationMethod< TFixedImage, TMovingImage, TTransform>
::ComputeFixedImageSampleList( const  GridImagePointer & gridImage )
{
	// initailise the sample list variable
	unsigned int numGridNodes = gridImage->GetLargestPossibleRegion().GetNumberOfPixels();
	m_SampleList.resize( numGridNodes, FixedImageSamples() );

	unsigned int numWeights = m_BSplineTransform->GetNumberOfWeights();

	// first we iterate through all the points in the fixed image
	typedef itk::ImageRegionConstIterator< FixedImageType > IteratorType;
	IteratorType it( m_FixedImage, m_FixedImageRegion );
	it.GoToBegin();
	int pcount = 0;
	while( !it.IsAtEnd() )
	{
		
		if( pcount % m_SamplingRate == 0 )
		{

		// find the location of the iterator
		typename FixedImageType::IndexType index = it.GetIndex();		
		typename FixedImageType::PointType point;
		m_FixedImage->TransformIndexToPhysicalPoint( index, point );

		// try transforming the point using the bspline transform
		typename FixedImageType::PointType outputPoint;
		typename BSplineTransformType::WeightsType weights( numWeights );
		typename BSplineTransformType::ParameterIndexArrayType indices( numWeights );

		bool isInside;
		m_BSplineTransform->TransformPoint( point, outputPoint, weights, indices, isInside );

		if( isInside )
		{
			for( unsigned int i = 0; i < numWeights; i++ )
			{
				int gridIndex = indices[i];
				double gridWeight = weights[i];

				FixedImageSampleValue value;
				value.first = index;
				value.second = gridWeight;

				m_SampleList[gridIndex].push_back( value );

			}
		}
		}

		pcount++;
		++it;		
	}
}

// ---------------------------------------------------------------------------------
template< typename TFixedImage, typename TMovingImage, typename TTransform >
void
MRFImageRegistrationMethod< TFixedImage, TMovingImage, TTransform>
::StartOptimization(void)
{
	// 
	// Set some useful variables
	//
	double currentMaxDisplacement = m_MaxLabelDisplacement;
	GridImagePointer gridImage = m_BSplineTransform->GetCoefficientImages()[0];
	int numPoints = gridImage->GetLargestPossibleRegion().GetNumberOfPixels();


	// 
	// Set the working parameters
	//
	ParametersType workingParameters;
	workingParameters = m_InitialTransformParameters;
	

	// 
	// Start iterating through the optimiser levels
	// 
	for( unsigned int level = 0; level < m_OptimiserLevels; level++ )
	{



		if( m_Verbose )
		{
			std::cout << "Optimising Level: " << level+1 << std::endl;
		}


		// 
		// Set up the working image by resampling the moving image by the 
		// initial transform parameters
		//
		typedef itk::ResampleImageFilter< MovingImageType, MovingImageType > ResamplerType;
		typename ResamplerType::Pointer resampler = ResamplerType::New();
		resampler->SetInput( this->GetMovingImage() );
		m_BSplineTransform->SetParametersByValue( workingParameters );
		resampler->SetTransform( m_BSplineTransform );
		resampler->SetReferenceImage( m_MovingImage );
		resampler->SetUseReferenceImage( true );
		resampler->SetInterpolator( m_Interpolator );
		resampler->SetDefaultPixelValue( 0 );

		try
		{
			resampler->Update();
		}
		catch ( itk::ExceptionObject &e )
		{
			std::cout << "Error in resampler: Couldn't resampler the initial working image" << std::endl;
			std::cout << e << std::endl;
			return;
		}

		m_Metric->SetMovingImage( resampler->GetOutput() );

		
		typedef itk::ImageFileWriter< MovingImageType > WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetInput( resampler->GetOutput() );
		writer->SetFileName("test.nrrd");
		writer->SetImageIO( itk::NrrdImageIO::New() );
		writer->Update();




		



		// 
		// Generate the labels for this optimiser level
		//
		typedef myitk::MRFImageRegistrationMethodLabelGenerator<ImageDimension> LabelGeneratorType;
		typedef typename LabelGeneratorType::LabelType LabelType;
		typedef typename LabelGeneratorType::LabelListType LabelListType;
		
		typename LabelGeneratorType::Pointer labelGenerator = LabelGeneratorType::New();
		labelGenerator->SetDenseSampling( m_DenseSampling );
		labelGenerator->SetMaxDisplacement( currentMaxDisplacement );
		labelGenerator->SetLabelSteps( m_LabelSteps );
		labelGenerator->SetSpacing( gridImage->GetSpacing() );

		LabelListType labels = labelGenerator->GenerateLabels();



		// 
		// Set up the variables for the optimisation
		//
		int numLabels = labels.size();
		MRF::CostVal * dCost = new MRF::CostVal[numLabels*numPoints];


		// 
		// compute the unary potentials
		//
		if( m_Verbose )
		{
			std::cout << "Computing Unary Potentials" << std::endl;
		}

		for( unsigned int labelNum = 0; labelNum < labels.size(); labelNum ++ )
		{
			typedef typename MetricType::MeasureType MeasureType;
			MeasureType value = m_Metric->GetValue( labels[labelNum] );

			//std::cout << labels[labelNum] << std::endl;
			for( unsigned int nodeNum = 0; nodeNum < value.size(); nodeNum++ )
			{
				dCost[nodeNum*numLabels+labelNum] = static_cast<MRF::CostVal> (value[nodeNum]);
				//std::cout << dCost[nodeNum*numLabels+labelNum] << std::endl;
			}
		}

		DataCost * data = new DataCost(dCost);


		// 
		// Compute the smoothness costs
		//
		MRF::CostVal *sCost = new MRF::CostVal[numLabels*numLabels];
		
		for( unsigned int l = 0; l < numLabels; l++ )
		{
			for( unsigned int k = 0; k < numLabels; k++ )
			{
				// compute the vector difference between the two labels
				LabelType l1 = labels[l];
				LabelType l2 = labels[k];

				double sum = 0.0;
				for( unsigned int dim = 0; dim < ImageDimension; dim++ )
				{
					sum += pow( (l2(dim) - l1(dim) ), 2.0 );				
				}

				sCost[l+numLabels*k] = static_cast<int>(m_Lambda * sqrt( sum ));
				//std::cout << sCost[l+numLabels*k] << std::endl;
			}
		}

		SmoothnessCost * smoothness = new SmoothnessCost(sCost);
		EnergyFunction * eng = new EnergyFunction( data, smoothness );

		float t;
		MRF * mrf = new Expansion( numPoints, numLabels, eng );
		for( unsigned int i = 0; i < m_Connections.size(); i++ )
		{
			mrf->setNeighbors(m_Connections[i].first, m_Connections[i].second, 1.0);
		}

		if( m_Verbose )
		{
			std::cout << "MRF Optimisation" << std::endl;
		}
		mrf->initialize();  
		mrf->clearAnswer();
		mrf->optimize(5,t);  // run for 5 iterations, store time t it took 

		MRF::EnergyVal E_smooth = mrf->smoothnessEnergy();
		MRF::EnergyVal E_data   = mrf->dataEnergy(); 

		for (int pix =0; pix < numPoints; pix++ ) 
		{
			int labelIndex = mrf->getLabel(pix);
			LabelType bestLabel = labels[labelIndex];
			for( unsigned int dim = 0; dim < ImageDimension; dim++ )
			{
				workingParameters(dim*numPoints+pix) += bestLabel[dim];
			}

			//std::cout << pix << "  " <<  labelIndex << std::endl;

		}
		std::cout << " --> Level " << level << " ( smoothness: " << E_smooth << " / data: " << E_data << " )"  << std::endl;
		
		delete mrf;


		/*
		// 
		// Set up the MRF optimisation data
		//
		int numLabels = labels.size();
		Graph::Real * lcosts = new Graph::Real[numLabels*numPoints];
		int numPairs = m_Connections.size();
		int *  pairs = new int[2*numPairs];
		Graph::Real * dist  = new Graph::Real[numLabels*numLabels];
		int max_iters = 10;
		Graph::Real * wcosts = new Graph::Real[numPairs];


		// 
		// Set up the pairs 
		//
		for( unsigned int p = 0; p < m_Connections.size(); p++ )
		{
			ConnectionType con = m_Connections[p];
			pairs[2*p] = con.first;
			pairs[2*p+1] = con.second;
			wcosts[p] = 1;
		}

		// 
		// Set up the weights
		//
		for( unsigned int l = 0; l < numLabels; l++ )
		{
			for( unsigned int k = 0; k < numLabels; k++ )
			{
				// compute the vector difference between the two labels
				LabelType l1 = labels[l];
				LabelType l2 = labels[k];

				double sum = 0.0;
				for( unsigned int dim = 0; dim < ImageDimension; dim++ )
				{
					sum += pow( (l2(dim) - l1(dim) ), 2.0 );				
				}

				dist[k*numLabels+l] = static_cast<int>(m_Lambda * sqrt( sum ));
				std::cout << dist[k*numLabels+l] << std::endl;
			}
		}


		// 
		// compute the unary potentials
		//
		if( m_Verbose )
		{
			std::cout << "Computing Unary Potentials" << std::endl;
		}

		for( unsigned int labelNum = 0; labelNum < labels.size(); labelNum ++ )
		{
			typedef typename MetricType::MeasureType MeasureType;
			MeasureType value = m_Metric->GetValue( labels[labelNum] );

			for( unsigned int nodeNum = 0; nodeNum < value.size(); nodeNum++ )
			{
				lcosts[labelNum*numPoints+nodeNum] = static_cast<int> (value[nodeNum]);
				//std::cout << lcosts[labelNum*numPoints+nodeNum] << std::endl;
			}
		}


		if( m_Verbose )
		{
			std::cout << "MRF Optimisation" << std::endl;
		}

		// 
		// Now we can actually do the optimisation, fingers crossed...
		//
		CV_Fast_PD pd( numPoints, numLabels, 
			   lcosts, numPairs, pairs, dist, max_iters, wcosts );

		pd.run();



		// 
		// Now we update the parameters
		//
		for( unsigned int i = 0; i < numPoints; i++ )
		{
			LabelType newLabel = labels[(int) pd._pinfo[i].label];
			for( unsigned int dim = 0; dim < ImageDimension; dim ++ )
			{
				workingParameters(dim*numPoints+i) += newLabel(dim);			
			}
		}
		*/


		// 
		// Update the max displacement for the next round of optimisation
		//
		currentMaxDisplacement *= m_MaxDisplacementScaleFactor;
	}


	m_LastTransformParameters = workingParameters;

}

// ---------------------------------------------------------------------------------
template< typename TFixedImage, typename TMovingImage, typename TTransform >
void
MRFImageRegistrationMethod< TFixedImage, TMovingImage, TTransform >
::GenerateData()
{
	ParametersType empty(1);
	empty.Fill(0.0);
	try
	{
		// initialize the interconnects between components
		this->Initialize();
	}
	catch ( ExceptionObject & err )
	{
		m_LastTransformParameters = empty;

		// pass exception to caller
		throw err;
	}

	this->StartOptimization();
}

// ---------------------------------------------------------------------------------
template< typename TFixedImage, typename TMovingImage, typename TTransform>
const typename MRFImageRegistrationMethod< TFixedImage, TMovingImage, TTransform>::TransformOutputType *
MRFImageRegistrationMethod< TFixedImage, TMovingImage, TTransform >
::GetOutput() const
{
	return static_cast< const TransformOutputType * >( this->ProcessObject::GetOutput(0) );
}

// ---------------------------------------------------------------------------------
template< typename TFixedImage, typename TMovingImage, typename TTransform >
DataObject::Pointer
MRFImageRegistrationMethod< TFixedImage, TMovingImage, TTransform >
::MakeOutput(DataObjectPointerArraySizeType output)
{
	switch ( output )
	{
		case 0:
			return TransformOutputType::New().GetPointer();
			break;
		default:
			itkExceptionMacro("MakeOutput request for an output number larger than the expected number of outputs");
			return 0;
	}
}

// ---------------------------------------------------------------------------------
template< typename TFixedImage, typename TMovingImage, typename TTransform >
void
MRFImageRegistrationMethod< TFixedImage, TMovingImage, TTransform >
::SetFixedImage(const FixedImageType *fixedImage)
{
	itkDebugMacro("setting Fixed Image to " << fixedImage);

	if ( this->m_FixedImage.GetPointer() != fixedImage )
	{
		this->m_FixedImage = fixedImage;

		// Process object is not const-correct so the const_cast is required here
		this->ProcessObject::SetNthInput( 0,
				const_cast< FixedImageType * >( fixedImage ) );

		this->Modified();
	}
}

template< typename TFixedImage, typename TMovingImage, typename TTransform >
void
MRFImageRegistrationMethod< TFixedImage, TMovingImage, TTransform>
::SetMovingImage(const MovingImageType *movingImage)
{
	itkDebugMacro("setting Moving Image to " << movingImage);

	if ( this->m_MovingImage.GetPointer() != movingImage )
	{
		this->m_MovingImage = movingImage;

		// Process object is not const-correct so the const_cast is required here
		this->ProcessObject::SetNthInput( 1,
				const_cast< MovingImageType * >( movingImage ) );

		this->Modified();
	}
}



}

#endif /* End of MYITKMRFIMAGEREGISTRATIONMETHOD_HPP */
