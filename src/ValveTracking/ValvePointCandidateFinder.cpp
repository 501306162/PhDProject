#include "ValvePointCandidateFinder.h"
#include <itkRegionalMaximaImageFilter.h>

#include <itkImageRegionConstIterator.h>

namespace vt
{
bool sort_list(const ValvePointCandidateFinder::PointPairType &a, const ValvePointCandidateFinder::PointPairType &b)
{
	return (a.first > b.first);
}


// ------------------------------------------------------------------------
void ValvePointCandidateFinder::Compute()
{
	// find the maxima 
	typedef itk::RegionalMaximaImageFilter<RealType, LabelType> MaximaFinderType;
	MaximaFinderType::Pointer maximaFinder = MaximaFinderType::New();
	maximaFinder->SetInput(m_Input);
	maximaFinder->SetFullyConnected(true);
	maximaFinder->SetForegroundValue(255);
	maximaFinder->SetBackgroundValue(0);
	maximaFinder->Update();

	LabelType::Pointer maxima = maximaFinder->GetOutput();

	// iterate through the maxima and sort them by value
	itk::ImageRegionConstIterator<RealType> it(m_Input, m_Input->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<LabelType> lit(maxima, maxima->GetLargestPossibleRegion());

	while(!it.IsAtEnd())
	{
		if(lit.Get() == 255)
		{
			RealType::IndexType index = it.GetIndex();
			PointType p;
			m_Input->TransformIndexToPhysicalPoint(index, p);
			PointPairType pair;
			pair.first = it.Get();
			pair.second = p;

			m_Output.push_back(pair);
		}

		++it; ++lit;
	}


	// sort the list
	std::sort(m_Output.begin(), m_Output.end(), sort_list);
}

// ------------------------------------------------------------------------
void ValvePointCandidateFinder::GetOutput(PointListType &list, const unsigned int num)
{
	for(unsigned int i = 0; i < num; i++)
	{
		list.push_back(m_Output[i]);
	}
}

}
