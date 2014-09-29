#ifndef HISTOGRAM_OF_GRADIENTS_FEATURE_EXTRACTOR_HPP
#define HISTOGRAM_OF_GRADIENTS_FEATURE_EXTRACTOR_HPP

#include "HistogramOfGradientsFeatureExtractor.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageHelper.h>

namespace filter
{
bool hist_sort(const std::pair<int, double> &a, const std::pair<int, double> &b)
{
	return a.second < b.second;
}


// ------------------------------------------------------------------------
template<typename TGradientImageType>
HistogramOfGradeintsFeatureExtractor<TGradientImageType>::
HistogramOfGradeintsFeatureExtractor()
{
	this->SetNumberOfRequiredInputs(1);
	m_SpatialBinSize.Fill(2);
	m_OrientationBins = 8;		
}


// ------------------------------------------------------------------------
template<typename TGradientImageType>
const typename HistogramOfGradeintsFeatureExtractor<TGradientImageType>::GradientImageType *
HistogramOfGradeintsFeatureExtractor<TGradientImageType>::
GetInput() const
{
	return itkDynamicCastInDebugMode<const GradientImageType*>(this->ProcessObject::GetInput(0));
}

// ------------------------------------------------------------------------
template<typename TGradientImageType>
void 
HistogramOfGradeintsFeatureExtractor<TGradientImageType>::
SetKeyPoint(const KeyPointType &keyPoint)
{
	m_KeyPoint = keyPoint;	
}


// ------------------------------------------------------------------------
template<typename TGradientImageType>
void 
HistogramOfGradeintsFeatureExtractor<TGradientImageType>::
SetInput(const GradientImageType * image)
{
	this->SetNthInput(0, const_cast<GradientImageType*>(image));
}


// ------------------------------------------------------------------------
template<typename TGradientImageType>
void 
HistogramOfGradeintsFeatureExtractor<TGradientImageType>::
Update()
{
	// create an offset table
	unsigned int start = 1;
	m_OffsetTable[0] = start;
	for(unsigned int i = 0; i < GradientImageType::ImageDimension; i++)
	{
		start *= m_SpatialBinSize[i];
		m_OffsetTable[i+1] = start;		
	}
	


	HistogramType histogram;
	InitialiseHistogram(histogram);


	// start iterating through the input patch
	typedef itk::ImageRegionConstIterator<GradientImageType> IteratorType;
	IteratorType it(GetInput(), GetInput()->GetLargestPossibleRegion());
	it.GoToBegin();

	unsigned int count = 0;
	while(!it.IsAtEnd())
	{
		PixelType grad = it.Get();
		double angle = std::atan2(grad[1], grad[0]);

		int sb, ob;
		GetHistogramIndex(count, angle, sb, ob);

		

		histogram[sb][ob] += grad.GetNorm();

		++count;
		++it;
	}

	// prepare the output
	m_Output.keyPoint = m_KeyPoint;
	m_Output.histogram = histogram;
}


// ------------------------------------------------------------------------
template<typename TGradientImageType>
typename  HistogramOfGradeintsFeatureExtractor<TGradientImageType>::OutputType
HistogramOfGradeintsFeatureExtractor<TGradientImageType>::
GetOutput() const
{
	return m_Output;
}


/*
// ------------------------------------------------------------------------
template<typename TGradientImageType>
void 
HistogramOfGradeintsFeatureExtractor<TGradientImageType>::
RankOrder(const HistogramType &histogram, 
			OutputType &output)
{
	unsigned int spatialBins = histogram.size();
	unsigned int orientationBins = histogram.front().size();
	unsigned int totalBins = spatialBins*orientationBins;
	
	std::vector<std::pair<int, double> > sortingList;

	for(unsigned int i = 0; i < totalBins; i++)
	{
		unsigned int sb = i / orientationBins;		
		unsigned int ob = i % orientationBins;

		std::pair<int, double> p(i, histogram[sb][ob]);
		sortingList.push_back(p);
	}

	// sort the list
	
	for(unsigned int i = 0; i < sortingList.size(); i++)
	{
		std::cout << sortingList[i].first << " " << sortingList[i].second << std::endl;
	}

	std::sort(sortingList.begin(), sortingList.end(), hist_sort);
	
	for(unsigned int i = 0; i < sortingList.size(); i++)
	{
		std::cout << sortingList[i].first << " " << sortingList[i].second << std::endl;
	}
}
*/

// ------------------------------------------------------------------------
template<typename TGradientImageType>
void 
HistogramOfGradeintsFeatureExtractor<TGradientImageType>::
GetSpatialIndex(const int count, int &spatialIndex)
{
	// get the image index
	IndexType index = GetInput()->ComputeIndex(count);
	SizeType imageSize = GetInput()->GetLargestPossibleRegion().GetSize();
	IndexType binIndex;
	for(unsigned int dim = 0; dim < GradientImageType::ImageDimension; dim++)
	{
		unsigned int pixelPerBinDimension = imageSize[dim] / m_SpatialBinSize[dim];  
		binIndex[dim] = index[dim] / pixelPerBinDimension;	
	}

	IndexType tmpIndex;
	tmpIndex.Fill(0);
	OffsetValueType offset = 0;
	itk::ImageHelper<GradientImageType::ImageDimension,
		GradientImageType::ImageDimension>::ComputeOffset(tmpIndex, binIndex, m_OffsetTable, offset);


	spatialIndex = offset;
}

// ------------------------------------------------------------------------
template<typename TGradientImageType>
void 
HistogramOfGradeintsFeatureExtractor<TGradientImageType>::
GetHistogramIndex(const int count, const double angle, 
		int &spatialBin, int &orientationBin)
{
	GetSpatialIndex(count, spatialBin);


	if(angle <= -M_PI)
	{
		orientationBin = 0;
	}
	else if(angle >= M_PI)
	{
		orientationBin = m_OrientationBins-1;
	}
	else
	{
		// align to positive
		double posAngle = angle + M_PI;	
		double orientationSplit = (2*M_PI) / static_cast<double>(m_OrientationBins);
		double position = posAngle / orientationSplit;
		orientationBin = static_cast<int>(std::floor(position));
	}
}

// ------------------------------------------------------------------------
template<typename TGradientImageType>
void 
HistogramOfGradeintsFeatureExtractor<TGradientImageType>::
InitialiseHistogram(HistogramType &histogram)
{
	m_NumPixels = GetInput()->GetLargestPossibleRegion().GetNumberOfPixels();
	// compute the spatial bin num
	m_SpatialBins = 1;
	for(unsigned int i = 0; i < GradientImageType::ImageDimension; i++)
	{
		m_SpatialBins *= m_SpatialBinSize[i];		
	}
	
	for(unsigned int i = 0; i < m_SpatialBins; i++)
	{
		OrientationListType list(m_OrientationBins,0.0);
		histogram.push_back(list);
	}
}


} /* filter */ 

#endif
