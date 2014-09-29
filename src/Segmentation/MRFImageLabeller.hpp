#ifndef MRF_IMAGE_LABELLER_HPP
#define MRF_IMAGE_LABELLER_HPP

#include "MRFImageLabeller.h"
#include <itkImageRegionIterator.h>
#include <itkProgressReporter.h>

namespace segmentation
{
template<typename TOutputImage>
MRFImageLabeller<TOutputImage>::
MRFImageLabeller()
{
	
}


// ------------------------------------------------------------------------
template<typename TOutputImage>
void 
MRFImageLabeller<TOutputImage>::
UseImageParameters(const ImageBaseType * image)
{
	SetSpacing(image->GetSpacing());
	SetSize(image->GetLargestPossibleRegion().GetSize());
	SetOrigin(image->GetOrigin());
	SetDirection(image->GetDirection());
}

// ------------------------------------------------------------------------
template<typename TOutputImage>
void
MRFImageLabeller<TOutputImage>::
SetLabelSet(const LabelSetType & labels)
{
	m_LabelSet = labels;
}


//----------------------------------------------------------------------------
template< typename TOutputImage >
void
MRFImageLabeller<TOutputImage>
::GenerateOutputInformation()
{
	TOutputImage *output;
	IndexType     index;

	index.Fill(0);

	output = this->GetOutput(0);

	typename TOutputImage::RegionType largestPossibleRegion;
	largestPossibleRegion.SetSize(this->m_Size);
	largestPossibleRegion.SetIndex(index);
	output->SetLargestPossibleRegion(largestPossibleRegion);

	output->SetSpacing(m_Spacing);
	output->SetOrigin(m_Origin);
	output->SetDirection(m_Direction);
}

//----------------------------------------------------------------------------
template< typename TOutputImage >
void
MRFImageLabeller< TOutputImage >
::ThreadedGenerateData(const ImageRegionType & outputRegionForThread,
                       itk::ThreadIdType threadId)
{

	// Support progress methods/callbacks
	itk::ProgressReporter progress( this, threadId, outputRegionForThread.GetNumberOfPixels() );

	typedef typename TOutputImage::PixelType scalarType;
	typename TOutputImage::Pointer image = this->GetOutput(0);

	itk::ImageRegionIterator< TOutputImage > it(image, outputRegionForThread);

	// Random number seed
	while(!it.IsAtEnd())
	{
		// get the index
		IndexType index = it.GetIndex();
		
		// get the linear index
		int idx = 0;
		if(ImageType::ImageDimension == 2)
		{
			idx = (index[1]*m_Size[0])+index[0];
		}
		else if(ImageType::ImageDimension == 3)
		{
			
			idx = (index[2]*m_Size[1]*m_Size[0])+(index[1]*m_Size[0])+index[0];
		}


		LabelType label = m_LabelSet[idx];

		it.Set(label);

		progress.CompletedPixel();

		++it;
	}
}


} /* segmentation */ 

#endif
