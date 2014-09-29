#ifndef MYITKMRFIMAGEREGISTRATIONMETHOD_H
#define MYITKMRFIMAGEREGISTRATIONMETHOD_H

#include <itkProcessObject.h>
#include <itkDataObjectDecorator.h>

#include "MRFUnaryPotentialsMetricBase.h"
#include "MRFImageRegistrationMethodConnectivityFilter.h"

namespace myitk
{
/**
 * @brief This is a class that performes MRF based registration
 */
template< typename TFixedImage, typename TMovingImage, typename TTransformType >
class MRFImageRegistrationMethod : public itk::ProcessObject 
{
public:
	/** Standard typedefs */
	typedef MRFImageRegistrationMethod		Self;
	typedef itk::ProcessObject				Superclass;
	typedef itk::SmartPointer< Self >		Pointer;
	typedef itk::SmartPointer< const Self > ConstPointer;

	/** factory macros */
	itkNewMacro( Self );
	itkTypeMacro( MRFImageRegistrationMethod, itk::ProcessObject );

	/** typedefs for the images */
	typedef	TFixedImage FixedImageType;
	typedef typename FixedImageType::ConstPointer FixedImageConstPointer;
	typedef typename FixedImageType::RegionType FixedImageRegionType;
	typedef TMovingImage MovingImageType;
	typedef typename MovingImageType::ConstPointer MovingImageConstPointer;

	
	itkStaticConstMacro(ImageDimension,
			unsigned int,
			TMovingImage::ImageDimension);
	
	/** Typedef for the metric */
	typedef myitk::MRFUnaryPotentialsMetricBase< FixedImageType, MovingImageType > MetricType;
	typedef typename MetricType::Pointer MetricPointer;
	typedef typename MetricType::FixedImageSampleValue FixedImageSampleValue;
	typedef typename MetricType::FixedImageSamples FixedImageSamples;
	typedef typename MetricType::FixedImageSamplesList FixedImageSamplesList;

	/** typedef the transform */
	typedef TTransformType BSplineTransformType;
	typedef typename BSplineTransformType::Pointer BSplineTransformPointer;
	typedef typename BSplineTransformType::ParametersType ParametersType;

	typedef itk::DataObjectDecorator< BSplineTransformType > TransformOutputType;
	typedef typename TransformOutputType::Pointer TransformOutputPointer;
	typedef typename TransformOutputType::ConstPointer TransformOutputConstPointer;
		
	/** typedef the grid image */
	typedef typename BSplineTransformType::ImageType GridImageType;
	typedef typename GridImageType::Pointer GridImagePointer;

	/** typedef for the interpolator */
	typedef typename MetricType::InterpolatorType InterpolatorType;
	typedef typename InterpolatorType::Pointer InterpolatorPointer;


	/** typedefs for the data object pointer */
	typedef typename itk::DataObject::Pointer DataObjectPointer;


	/** typedefs for the connectivity */
	typedef myitk::MRFImageRegistrationMethodConnectivityFilter< GridImageType > ConnectivityFilterType;
	typedef typename ConnectivityFilterType::Pointer	ConnectivityFilterPointer;
	typedef typename ConnectivityFilterType::ConnectionType ConnectionType;
	typedef typename ConnectivityFilterType::ConnectionListType ConnectionListType;

	/** Set and Get for the images */
	void SetFixedImage( const FixedImageType * image );
	itkGetConstObjectMacro( FixedImage, FixedImageType );
	void SetMovingImage( const MovingImageType * image );
	itkGetConstObjectMacro( MovingImage, MovingImageType );

	/** Set and Get the fixed image region */
	void SetFixedImageRegion( const FixedImageRegionType &region );

	/** Set and get for the object objects */
	itkSetObjectMacro( Metric, MetricType );
	itkGetConstObjectMacro( Metric, MetricType );
	itkSetObjectMacro( BSplineTransform, BSplineTransformType );
	itkGetConstObjectMacro( BSplineTransform, BSplineTransformType );
	itkSetObjectMacro( Interpolator, InterpolatorType );
	itkGetConstObjectMacro( Interpolator, InterpolatorType );

	/** Set/Get the initial transformation parameters. */
	virtual void SetInitialTransformParameters(const ParametersType & param);
	itkGetConstReferenceMacro(InitialTransformParameters, ParametersType);

	/** Get the last transformation parameters visited by
	 * the optimizer. */
	itkGetConstReferenceMacro(LastTransformParameters, ParametersType);

	/** Initialize by setting the interconnects between the components. */
	virtual void Initialize()
		throw ( ExceptionObject );

	/** Returns the transform resulting from the registration process  */
	const TransformOutputType * GetOutput() const;

	/** Set / Get for Verbose */
	itkSetMacro( Verbose, bool );
	itkGetConstMacro( Verbose, bool );

	/** Make a DataObject of the correct type to be used as the specified
	 * output. */
	typedef ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;
	using Superclass::MakeOutput;
	virtual DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx);


	/** Set / Get for OptimiserLevels */
	itkSetMacro( OptimiserLevels, int );
	itkGetConstMacro( OptimiserLevels, int );
	/** Set / Get for MaxLabelDisplacement */
	itkSetMacro( MaxLabelDisplacement, double );
	itkGetConstMacro( MaxLabelDisplacement, double );
	/** Set / Get for MaxDisplacementScaleFactor */
	itkSetMacro( MaxDisplacementScaleFactor, double );
	itkGetConstMacro( MaxDisplacementScaleFactor, double );
	/** Set / Get for LabelSteps */
	itkSetMacro( LabelSteps, int );
	itkGetConstMacro( LabelSteps, int );
	/** Set / Get for DenseSampling */
	itkSetMacro( DenseSampling, bool );
	itkGetConstMacro( DenseSampling, bool );
	/** Set / Get for SamplingRate */
	itkSetMacro( SamplingRate, int );
	itkGetConstMacro( SamplingRate, int );

	/** Set / Get for Lambda */
	itkSetMacro( Lambda, double );
	itkGetConstMacro( Lambda, double );


protected:
	MRFImageRegistrationMethod();
	virtual ~MRFImageRegistrationMethod() {};

	void GenerateData();
	void StartOptimization();

private:
	MRFImageRegistrationMethod( const Self & );
	void operator=( const Self &);

	/** Function to compute the grid point influences */
	void ComputeFixedImageSampleList( const GridImagePointer &gridImage );
	FixedImageSamplesList m_SampleList;

	/** objects for the registration */
	MovingImageConstPointer m_MovingImage;
	FixedImageConstPointer m_FixedImage;
	MetricPointer m_Metric;
	BSplineTransformPointer m_BSplineTransform;
	InterpolatorPointer m_Interpolator;

	ParametersType m_InitialTransformParameters;
	ParametersType m_LastTransformParameters;

	FixedImageRegionType m_FixedImageRegion;
	bool m_FixedRegionSet;

	ConnectionListType m_Connections;
	bool m_Verbose;

	int m_OptimiserLevels;
	int m_LabelSteps;
	bool m_DenseSampling;
	double m_MaxLabelDisplacement;
	double m_MaxDisplacementScaleFactor;
	double m_Lambda;

	int m_SamplingRate;

};


} // end namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "MRFImageRegistrationMethod.hpp"
#endif

#endif /* End of MYITKMRFIMAGEREGISTRATIONMETHOD_H */
