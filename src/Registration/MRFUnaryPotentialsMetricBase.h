#ifndef MRF_POTENTIALS_METRIC_BASE_H
#define MRF_POTENTIALS_METRIC_BASE_H

#include <itkImage.h>
#include <itkMultipleValuedCostFunction.h>
#include <itkTranslationTransform.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkSmartPointer.h>
#include <itkMultiThreader.h>

namespace myitk
{
using namespace itk;
template< typename TFixedImage, typename TMovingImage >
class  MRFUnaryPotentialsMetricBase : public itk::MultipleValuedCostFunction
{
public:
	typedef MRFUnaryPotentialsMetricBase Self;
	typedef itk::MultipleValuedCostFunction   Superclass;
	typedef SmartPointer< Self >         Pointer;
	typedef SmartPointer< const Self >   ConstPointer;



	itkTypeMacro(MRFUnaryPotentialsMetricBase, MultipleValuedCostFunction);
	itkStaticConstMacro(ImageDimension, unsigned int, TMovingImage::ImageDimension);

	/** Typedefs for the moving image */
	typedef TMovingImage                            MovingImageType;
	typedef typename MovingImageType::Pointer       MovingImagePointer;
	typedef typename MovingImageType::ConstPointer  MovingImageConstPointer;
	typedef typename MovingImageType::IndexType     MovingImageIndexType;
	typedef typename MovingImageType::PointType     MovingImagePointType;

	typedef TFixedImage                             FixedImageType;
	typedef typename FixedImageType::Pointer        FixedImagePointer;
	typedef typename FixedImageType::ConstPointer   FixedImageConstPointer;
	typedef typename FixedImageType::IndexType      FixedImageIndexType;
	typedef typename FixedImageType::PointType      FixedImagePointType;


	/** Typedefs for the samples list */
	typedef std::pair<FixedImageIndexType, double>  FixedImageSampleValue;
	typedef std::vector<FixedImageSampleValue>      FixedImageSamples;
	typedef std::vector<FixedImageSamples>          FixedImageSamplesList;


	/** typedefs for the transformation */
	typedef TranslationTransform<double, ImageDimension>      TranslationTransformType;
	typedef typename TranslationTransformType::Pointer        TranslationTransformPointer;
	typedef typename TranslationTransformType::ParametersType ParametersType;


	/** typedef for the interpolator */
	typedef InterpolateImageFunction< MovingImageType, double> InterpolatorType;
	typedef typename InterpolatorType::Pointer                 InterpolatorPointer;

	typedef Superclass::MeasureType  MeasureType;
	typedef Superclass::DerivativeType DerivativeType;

	virtual void GetDerivative(const ParametersType &params, DerivativeType &derivative) const {}

	virtual MeasureType GetValue( const ParametersType & parameters ) const;


	/** Get/Set the number of threads to create when executing. */
	itkSetClampMacro( NumberOfThreads, int, 1, ITK_MAX_THREADS );
	itkGetConstReferenceMacro( NumberOfThreads, int );

	/** Getters and setters */
	itkSetConstObjectMacro(MovingImage, MovingImageType);
	itkGetConstObjectMacro(MovingImage, MovingImageType);

	itkSetConstObjectMacro(FixedImage, FixedImageType);
	itkGetConstObjectMacro(FixedImage, FixedImageType);

	itkGetConstObjectMacro(Transform, TranslationTransformType);
	itkSetObjectMacro(Interpolator, InterpolatorType);
	itkGetConstObjectMacro(Interpolator, InterpolatorType);

	void SetFixedImageSamplesList(FixedImageSamplesList &list) { this->m_FixedImageSamplesList = list; }
	FixedImageSamplesList & GetFixedImageSamplesList() { return this->m_FixedImageSamplesList; }

	itkSetMacro(UseDistanceWeighting, bool);
	itkGetMacro(UseDistanceWeighting, bool);

	MultiThreader * GetMultiThreader() { return m_Threader; }

	unsigned int GetNumberOfParameters() const { return m_Transform->GetNumberOfParameters(); }

	virtual unsigned int GetNumberOfValues() const;


protected:

	struct SampleRegion
	{
		int threadStartIndex;
		int threadRange;
	};


	MRFUnaryPotentialsMetricBase();
	virtual ~MRFUnaryPotentialsMetricBase() {}

	virtual void Initialise() throw(ExceptionObject &);

	virtual MeasureType GetNonConstValue(const ParametersType & params);

	virtual void ComputeUnaryPotentials( const SampleRegion &region,  int threadId ) = 0;

	MovingImageConstPointer     m_MovingImage;
	FixedImageConstPointer      m_FixedImage;
	TranslationTransformPointer m_Transform;
	FixedImageSamplesList       m_FixedImageSamplesList;
	InterpolatorPointer         m_Interpolator;
	bool                        m_UseDistanceWeighting;
	MeasureType                 m_OutputValue;


	struct ThreadStruct
	{
		Pointer Metric;
	};

	static ITK_THREAD_RETURN_TYPE ThreaderCallback( void *arg );

	virtual int SplitFixedSampleList(int i, int num, SampleRegion &region);

private:
	MRFUnaryPotentialsMetricBase(const Self&);
	void operator=(const Self&);



	/*************************************************
	 * Multi threading stuff.... eak
	 */
	typedef MultiThreader MultithreaderType;

	MultithreaderType::Pointer m_Threader;
	int m_NumberOfThreads;



};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "MRFUnaryPotentialsMetricBase.hpp"
#endif

#endif
