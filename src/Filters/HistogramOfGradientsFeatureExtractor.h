#ifndef HISTOGRAM_OF_GRADIENTS_FEATURE_EXTRACTOR_H
#define HISTOGRAM_OF_GRADIENTS_FEATURE_EXTRACTOR_H

#include <itkProcessObject.h>
#include <itkImage.h>
#include <itkVariableLengthVector.h>
#include <itkHistogram.h>
#include "FeatureCommon.h"

namespace filter
{
template<typename TGradientImageType>
class HistogramOfGradeintsFeatureExtractor : public itk::ProcessObject
{
public:
	typedef HistogramOfGradeintsFeatureExtractor Self;
	typedef itk::ProcessObject Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkTypeMacro(HistogramOfGradeintsFeatureExtractor, ProcessObject);
	itkNewMacro(Self);



	typedef TGradientImageType GradientImageType;
	typedef typename GradientImageType::SizeType SizeType;
	typedef typename GradientImageType::OffsetValueType OffsetValueType;
	typedef typename GradientImageType::IndexType IndexType;
	typedef typename GradientImageType::PixelType PixelType;

	typedef std::vector<double> OrientationListType;
	typedef std::vector<OrientationListType> HistogramType;





	typedef OrientatedKeyPoint<GradientImageType::ImageDimension> KeyPointType;
	typedef HoGFeature<GradientImageType::ImageDimension> FeatureType;
	typedef FeatureType OutputType;
	




	itkSetMacro(SpatialBinSize, SizeType);
	itkSetMacro(OrientationBins, unsigned int);

	void SetInput(const GradientImageType *image);
	void SetKeyPoint(const KeyPointType &keyPoint);
	const GradientImageType * GetInput() const;

	
	OutputType GetOutput() const;


	void Update();

protected:
	HistogramOfGradeintsFeatureExtractor();
	~HistogramOfGradeintsFeatureExtractor() {}
	

	void GetHistogramIndex(const int count, const double angle, 
			int &spatialBin, int &orientationBin);


	void GetSpatialIndex(const int count, int &spatialIndex);


	void InitialiseHistogram(HistogramType &histogram);
	//void RankOrder(const HistogramType &histogram, 
	//		OutputType &output);
private:
	HistogramOfGradeintsFeatureExtractor(const Self&);
	void operator=(const Self&);


	SizeType m_SpatialBinSize;
	unsigned int m_OrientationBins;
	unsigned int m_NumPixels;
	unsigned int m_SpatialBins;


	KeyPointType m_KeyPoint;

	OffsetValueType m_OffsetTable[GradientImageType::ImageDimension+1];
	OutputType m_Output;


};
} /* filter */ 

#include "HistogramOfGradientsFeatureExtractor.hpp"


#endif
