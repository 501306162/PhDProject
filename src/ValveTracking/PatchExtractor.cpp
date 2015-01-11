#include "PatchExtractor.h"
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkImageRegionIterator.h>

namespace vt
{
// ------------------------------------------------------------------------
template<typename PValue>
typename PatchExtractor<PValue>::ImageType::Pointer 
PatchExtractor<PValue>::ExtractPatch()
{
	ContIndexType startIndex;
	for(unsigned int i = 0; i < 2; i++)
	{
		startIndex[i] = m_Center[i] - (m_Size[i] / 2.0);
	}
	startIndex[2] = 0;

	PointType newOrigin;
	m_Image->TransformContinuousIndexToPhysicalPoint(startIndex, newOrigin);


	typename ImageType::Pointer patch = ImageType::New();
	patch->SetDirection(m_Image->GetDirection());
	patch->SetOrigin(newOrigin);
	patch->SetSpacing(m_Image->GetSpacing());
	
	typename ImageType::RegionType region;
	m_Size[2] = 1;
	region.SetSize(m_Size);
	
	patch->SetRegions(region);
	patch->Allocate();

	
	typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType;
	typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
	interpolator->SetInputImage(m_Image);


	itk::ImageRegionIterator<ImageType> imageIt(patch, region);
	while(!imageIt.IsAtEnd())
	{
		typename ImageType::IndexType index = imageIt.GetIndex();
		PointType point;
		patch->TransformIndexToPhysicalPoint(index, point);

		if(interpolator->IsInsideBuffer(point))
			imageIt.Set(interpolator->Evaluate(point));
		else
			imageIt.Set(0);
		++imageIt;
	}

	return patch;
	
}


// ------------------------------------------------------------------------
template<typename PValue>
typename PatchExtractor<PValue>::RegionType PatchExtractor<PValue>::ExtractionRegion()
{
	ContIndexType startIndex;
	for(unsigned int i = 0; i < 2; i++)
	{
		startIndex[i] = m_Center[i] - (m_Size[i] / 2.0);
	}
	startIndex[2] = 0;


	typename ImageType::IndexType properIndex;
	for(unsigned int i = 0; i < 3; i++)
	{
		properIndex[i] = round(startIndex[i]);
	}

	for(unsigned int i = 0; i < 2; i++)
	{
		if(properIndex[i] < 0) properIndex[i] = 0;
		if((m_Size[i] + properIndex[i]) > (m_Image->GetLargestPossibleRegion().GetSize()[i]-1))
		{
			m_Size[i] = m_Image->GetLargestPossibleRegion().GetSize()[i]-1-properIndex[i];
		}
	}
	
	MaskType::RegionType region;
	region.SetIndex(properIndex);
	m_Size[2] = 1;
	region.SetSize(m_Size);

	return region;
}

// ------------------------------------------------------------------------
template<typename PValue>
typename PatchExtractor<PValue>::MaskType::Pointer PatchExtractor<PValue>::ExtractMask()
{
	ContIndexType startIndex;
	for(unsigned int i = 0; i < 2; i++)
	{
		startIndex[i] = m_Center[i] - (m_Size[i] / 2.0);
	}
	startIndex[2] = 0;


	typename ImageType::IndexType properIndex;
	for(unsigned int i = 0; i < 3; i++)
	{
		properIndex[i] = round(startIndex[i]);
	}

	for(unsigned int i = 0; i < 2; i++)
	{
		if(properIndex[i] < 0) properIndex[i] = 0;
		if((m_Size[i] + properIndex[i]) > (m_Image->GetLargestPossibleRegion().GetSize()[i]-1))
		{
			m_Size[i] = m_Image->GetLargestPossibleRegion().GetSize()[i]-1-properIndex[i];
		}
	}
	
	MaskType::RegionType region;
	region.SetIndex(properIndex);
	m_Size[2] = 1;
	region.SetSize(m_Size);

	MaskType::Pointer mask = MaskType::New();
	mask->SetDirection(m_Image->GetDirection());
	mask->SetSpacing(m_Image->GetSpacing());
	mask->SetOrigin(m_Image->GetOrigin());
	mask->SetRegions(m_Image->GetLargestPossibleRegion());
	mask->Allocate();
	mask->FillBuffer(0);

	itk::ImageRegionIterator<MaskType> imageIt(mask, region);
	while(!imageIt.IsAtEnd())
	{
		imageIt.Set(m_MaskValue);
		++imageIt;
	}

	return mask;
}


template class PatchExtractor<unsigned short>;
template class PatchExtractor<float>;
template class PatchExtractor<double>;
template class PatchExtractor<unsigned char>;

}
