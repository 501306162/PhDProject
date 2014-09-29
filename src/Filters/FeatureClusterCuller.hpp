#ifndef FEATURE_CLUSTER_CULLER_HPP
#define FEATURE_CLUSTER_CULLER_HPP

#include "FeatureClusterCuller.h"

namespace filter
{
// ------------------------------------------------------------------------
template<int TDimension>
FeatureClusterCuller<TDimension>::
FeatureClusterCuller()
{
	this->SetNumberOfRequiredInputs(0);
}

// ------------------------------------------------------------------------
template<int TDimension>
void
FeatureClusterCuller<TDimension>::
SetInput(const ClusterList &input)
{
	m_Input = input;
}

// ------------------------------------------------------------------------
template<int TDimension>
typename FeatureClusterCuller<TDimension>::ClusterList
FeatureClusterCuller<TDimension>::
GetOutput() const
{
	return m_Output;
}


// ------------------------------------------------------------------------
template<int TDimension>
void
FeatureClusterCuller<TDimension>::
Update()
{
	// first remove duplicates
	ClusterList deduped;
	DeDuplicate(m_Input, deduped);

	// remove ones that are members of a larger set
	for(unsigned int i = 0; i < deduped.size(); i++)
	{
		ClusterType &c1 = deduped[i];

		bool remove = false;
		for(unsigned int j = i; j < deduped.size(); j++)
		{
			if(i==j) continue;
			
			ClusterType &c2 = deduped[j];
			if(c2.size() > c1.size() && c2.contains(c1))
			{
				remove = true;
				break;
			}
		}

		if(!remove)
			m_Output.push_back(c1);

	}
}

// ------------------------------------------------------------------------
template<int TDimension>
void
FeatureClusterCuller<TDimension>::
DeDuplicate(ClusterList &input, ClusterList &output)
{
	for(unsigned int i = 0; i < input.size(); i++)
	{
		if(input[i].size() == 1) continue;

		bool duplicate = false;
		for(unsigned int j = i; j < input.size(); j++)
		{
			if(i == j) continue;
			duplicate = input[i].equals(input[j]);
			if(duplicate)
			{
				duplicate = true;
				break;
			}
		}

		if(!duplicate)
		{
			output.push_back(input[i]);
		}
	}

	
}


} /* filter */ 

#endif
