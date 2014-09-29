#ifndef FEATURE_CLUSTER_LEARNER_HPP
#define FEATURE_CLUSTER_LEARNER_HPP

#include "FeatureClusterLearner.h"

namespace filter
{
// ------------------------------------------------------------------------
template<int TDimension>
FeatureClusterLearner<TDimension>::
FeatureClusterLearner()
{
	this->SetNumberOfRequiredInputs(0);
	m_Gamma = 1.0;
	m_E = 1.0;
	m_EIncrement = 100;

}

// ------------------------------------------------------------------------
template<int TDimension>
void
FeatureClusterLearner<TDimension>::
SetInput(const ClusterType &input)
{
	m_Input = input;
	this->Modified();
}

// ------------------------------------------------------------------------
template<int TDimension>
void
FeatureClusterLearner<TDimension>::
SetFeatures(const FeatureList &features)
{
	m_Features = features;
}

// ------------------------------------------------------------------------
template<int TDimension>
void
FeatureClusterLearner<TDimension>::
SetDistanceMatrix(const MatrixType &mat)
{
	m_Distances = mat;
}

// ------------------------------------------------------------------------
template<int TDimension>
void
FeatureClusterLearner<TDimension>::
Update()
{
	// get the id and distances for this cluster
	m_ClusterId = m_Input.clusterId;
	unsigned int totalFeatureNumber = m_Features.size();

	// create the ids of the features
	m_NonMemberIndices.clear();
	for(unsigned int i = 0; i < totalFeatureNumber; i++)
	{
		m_NonMemberIndices.push_back(i);
	}

	// get the ids of the cluster members
	m_MemberIndices.clear();
	for(unsigned int i = 0; i < m_Input.clusterMembers.size(); i++)
	{
		int memberId = m_Input.clusterMembers[i].featureId;
		m_MemberIndices.push_back(memberId);		
		m_NonMemberIndices.erase(
				m_NonMemberIndices.begin()+(memberId-i));
	}

	m_E = 0.0;
	m_EIncrement = 100.0;
	m_E = Optimise(m_E, m_EIncrement);

	
	// sort out the output
	m_Output.clusterId = m_Input.clusterId;
	m_Output.keyFeature = m_Input.keyFeature;
	m_Output.clusterMembers.clear();
	
	for(unsigned int i = 0; i < m_MemberIndices.size(); i++)
	{
		if(m_Distances(m_ClusterId,m_MemberIndices[i]) <= m_E)
			m_Output.clusterMembers.push_back(m_Features[m_MemberIndices[i]]);
	}

}

// ------------------------------------------------------------------------
template<int TDimension>
double
FeatureClusterLearner<TDimension>::
Optimise(double e, double eInc)
{
	//std::cout << e << " " << eInc << std::endl;
	double gamma = ComputeGamma(e);
	
	// check for still in progress
	if(gamma >= m_Gamma)
	{
		e+=eInc;
		return Optimise(e, eInc);
	}


	// can Inc be reduced?
	if(eInc > 0.5)
	{
		e-=eInc;
		eInc /= 2;
		return Optimise(e, eInc);
	}


	//must have ended
	return e-eInc;		
}



// ------------------------------------------------------------------------
template<int TDimension>
typename FeatureClusterLearner<TDimension>::ClusterType
FeatureClusterLearner<TDimension>::
GetOutput() const
{
	return m_Output;
}


// ------------------------------------------------------------------------
template<int TDimension>
double
FeatureClusterLearner<TDimension>::
ComputeGamma(double e)
{
	double memberCount = 0.0;
	double nonMemberCount = 0.0;
	m_OutputMemberIds.clear();
	for(unsigned int i = 0; i < m_MemberIndices.size(); i++)
	{
		int ind = m_MemberIndices[i];
		if(m_Distances(m_ClusterId, ind) <= e)
		{
			m_OutputMemberIds.push_back(ind);
			memberCount+=1.0;
		}
	}


	for(unsigned int i = 0; i < m_NonMemberIndices.size(); i++)
	{
		int ind = m_NonMemberIndices[i];
		if(m_Distances(m_ClusterId, ind) <= e)
			nonMemberCount+=1.0;
	}

	return memberCount / nonMemberCount;
}
} /* filter */ 

#endif
