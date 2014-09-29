#ifndef SHAPE_GENERATOR_HPP
#define SHAPE_GENERATOR_HPP

#include "ShapeGenerator.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkNumericTraits.h>

namespace manifold
{
// ------------------------------------------------------------------------
template<typename TInputType, typename TOutputType>
ShapeGenerator<TInputType, TOutputType>::
ShapeGenerator()
{

}

// ------------------------------------------------------------------------
template<typename TInputType, typename TOutputType>
void
ShapeGenerator<TInputType, TOutputType>::
SetWeightList(const WeightListType &list)
{
	m_WeightList = list;
}


// ------------------------------------------------------------------------
template<typename TInputType, typename TOutputType>
void
ShapeGenerator<TInputType, TOutputType>::
BeforeThreadedGenerateData()
{
	unsigned int numInputs = this->GetNumberOfInputs();
	if(numInputs != m_WeightList.size())
	{
		itkExceptionMacro(<< "Number of inputs doesn't match the number of weights");
	}
}



// ------------------------------------------------------------------------
template<typename TInputType, typename TOutputType>
void
ShapeGenerator<TInputType, TOutputType>::
ThreadedGenerateData(const OutputRegionType &region, itk::ThreadIdType id)
{
	// create the input iterators
	unsigned int numInputs = m_WeightList.size();
	typedef itk::ImageRegionConstIterator<InputType> InputIteratorType;
	std::vector<InputIteratorType> iterators;
	for(unsigned int i = 0; i < numInputs; i++)
	{
		InputIteratorType it(this->GetInput(i), region);
		iterators.push_back(it);
	}
	
	typedef itk::ImageRegionIterator<OutputType> OutputIteratorType;
	OutputIteratorType outIt(this->GetOutput(), region);
	while(!outIt.IsAtEnd())
	{
		OutputValueType out = itk::NumericTraits<OutputValueType>::ZeroValue();

		for(unsigned int i = 0; i < numInputs; i++)
		{
			InputValueType val = iterators[i].Get();
			out += static_cast<OutputValueType>(
					m_WeightList[i] * static_cast<double>(val));
			++(iterators[i]);
		}

		outIt.Set(out);
		
		++outIt;
	}
}


} /* manifold */ 

#endif
