#ifndef PATCH_EXTRACTOR2_HPP
#define PATCH_EXTRACTOR2_HPP


#include "PatchExtractor2.h"
#include <itkLinearInterpolateImageFunction.h> 
#include <itkImageRegionIterator.h>


using namespace itk;

namespace vt
{
// ------------------------------------------------------------------------
template<typename TImageType>
PatchExtractor2<TImageType>::
PatchExtractor2()
{
}


// ------------------------------------------------------------------------
template<typename TImageType>
void
PatchExtractor2<TImageType>::
ThreadedGenerateData(const ImageRegionType & region,
                               itk::ThreadIdType id)
{
	itk::ImageRegionIterator<ImageType> it(this->GetOutput(), region);
	
	while(!it.IsAtEnd())
	{
		IndexType index = it.GetIndex();
		PointType point;
		this->GetOutput()->TransformIndexToPhysicalPoint(index, point);

		if(m_Interpolator->IsInsideBuffer(point))
		{
			it.Set(m_Interpolator->Evaluate(point));
		}
		else 
		{
			it.Set(0);
		}

		++it;
	}
}


// ------------------------------------------------------------------------
template<typename TImageType>
void
PatchExtractor2<TImageType>::
BeforeThreadedGenerateData()
{
	if(!m_Interpolator)
	{
		m_Interpolator = itk::LinearInterpolateImageFunction<ImageType, double>::New();
	}
	
	m_Interpolator->SetInputImage(this->GetInput(0));

}

// ------------------------------------------------------------------------
template<typename TImageType>
void
PatchExtractor2<TImageType>::
GenerateInputRequestedRegion()
{
	Superclass::GenerateInputRequestedRegion();
	ImageType * input = const_cast<ImageType*>(this->GetInput());
	input->SetRequestedRegionToLargestPossibleRegion();
}


// ------------------------------------------------------------------------
template<typename TImageType>
void
PatchExtractor2<TImageType>::
GenerateOutputInformation()
{
	const ImageType * input = this->GetInput();
	ImageType * output = this->GetOutput();

	// compute the output direction 
	DirectionType imDirection = input->GetDirection();
	VectorType xDir = m_Line;
	xDir.Normalize();
	VectorType zDir;

	for(unsigned int i = 0; i < ImageType::ImageDimension; i++)
	{
		zDir[i] = imDirection(i,2);		
	}

	
	VectorType yDir = itk::CrossProduct(xDir, zDir);
	double trace = 0;
	for(unsigned int i = 0; i < 3; i++)
	{
		trace += imDirection(i,i);
	}

	if(trace < 0.0)
	{
		yDir = -yDir;
	}

	
	DirectionType outDirection;
	PointType origin;
	SpacingType spacing;
	for(unsigned int i = 0; i < ImageType::ImageDimension; i++)
	{
		outDirection(i,0) = xDir[i];		
		outDirection(i,1) = yDir[i];		
		outDirection(i,2) = zDir[i];		

		origin[i] = m_Center[i];

	}
	

	for(unsigned int i = 0; i < 3; i++)
	{
		for(unsigned int j = 0; j < 3; j++)
		{
			origin[i] -= m_Distance[j]*0.5*outDirection(i,j);
		}

		spacing[i] = std::max(m_Distance[i], input->GetSpacing()[i]) / static_cast<double>(m_Size[i]);
	}


	ImageRegionType region;
	IndexType index;
	index.Fill(0);
	region.SetIndex(index);
	region.SetSize(m_Size);

	output->SetDirection(outDirection);
	output->SetOrigin(origin);
	output->SetSpacing(spacing);
	output->SetRegions(region);


}


}


#endif


