#ifndef HELPERS_HPP
#define HELPERS_HPP

#include "Helpers.h"
#include <itkExtractImageFilter.h>


namespace utils
{
template<typename TInputType, typename TOutputType>
void
ExtractSlice(const TInputType * input, const unsigned int slice,
	   typename TOutputType::Pointer &output)
{
	typedef itk::ExtractImageFilter<TInputType, TOutputType> ExtractorType;
	typename ExtractorType::Pointer extractor = ExtractorType::New();
	extractor->SetInput(input);

	// get the output region
	typedef TInputType InputType;
	typedef typename InputType::RegionType RegionType;
	typedef typename InputType::IndexType IndexType;
	typedef typename InputType::SizeType SizeType;

	RegionType inputRegion = input->GetLargestPossibleRegion();
	IndexType index = inputRegion.GetIndex();
	SizeType size = inputRegion.GetSize();

	RegionType outputRegion;
	index[2] = slice;
	size[2] = 0;

	outputRegion.SetSize(size);
	outputRegion.SetIndex(index);

	extractor->SetExtractionRegion(outputRegion);
	extractor->SetDirectionCollapseToIdentity();
	extractor->Update();

	output = extractor->GetOutput();
	
}
} /* utils */ 

#endif
