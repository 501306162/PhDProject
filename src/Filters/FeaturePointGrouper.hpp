#ifndef FEATURE_POINT_GROUPER_HPP
#define FEATURE_POINT_GROUPER_HPP

#include "FeaturePointGrouper.h"
#include "CommonDefinitions.h"

namespace filter
{
// ------------------------------------------------------------------------
template<int TDimension>
FeaturePointGrouper<TDimension>::
FeaturePointGrouper()
{
	this->SetNumberOfRequiredInputs(0);
	m_LocationThreshold = 1.0;
	m_AngleThreshold = 0.8;
	m_ScaleThreshold = 1.5;
	m_AdjustLocationByScale = true;
	m_GroupSizeThreshold = 0;

}

// ------------------------------------------------------------------------
template<int TDimension>
void 
FeaturePointGrouper<TDimension>::
SetInput(const FeatureListType &input)
{
	m_Input = input;
	this->Modified();
}

// ------------------------------------------------------------------------
template<int TDimension>
void 
FeaturePointGrouper<TDimension>::
Update()
{
	for(unsigned int i = 0; i < m_Input.size(); i++)
	{
		// create a new feature group
		FeatureGroup group;
		group.clusterId = i;
		group.keyFeature = m_Input[i];

		for(unsigned int j = 0; j < m_Input.size(); j++)
		{
			if(AreSimilarEnough(m_Input[i], m_Input[j]))
				group.clusterMembers.push_back(m_Input[j]);
		}

		if(group.size() > m_GroupSizeThreshold)
			m_Output.push_back(group);
		
	}
}

// ------------------------------------------------------------------------
template<int TDimension>
typename FeaturePointGrouper<TDimension>::FeatureGroupList
FeaturePointGrouper<TDimension>::
GetOutput() const
{
	return m_Output;
}


// ------------------------------------------------------------------------
template<int TDimension>
bool 
FeaturePointGrouper<TDimension>::
AreSimilarEnough(const FeatureType &f1, const FeatureType &f2)
{

	// compute the scale difference
	double scaleDiff = std::fabs(std::log(f2.keyPoint.scale / f1.keyPoint.scale));
	if(scaleDiff > m_ScaleThreshold)
		return false;


	// compute the angular difference
	double angleDiff = std::atan2(sin(f1.keyPoint.angle-f2.keyPoint.angle), cos(f1.keyPoint.angle-f2.keyPoint.angle));
	angleDiff = utils::RadiansToDegrees(utils::UnwrapAngle(angleDiff));
	if(angleDiff > 180)
		angleDiff = 360 - angleDiff;
	angleDiff = utils::DegreesToRadians(angleDiff);
	if(angleDiff > m_AngleThreshold)
		return false;


	// compute the location difference
	double locDiff = 0.0;
	for(unsigned int i = 0; i < TDimension; i++)
	{
		double diff = f1.keyPoint.location[i] - f2.keyPoint.location[i];
		locDiff += (diff*diff);				
	}
	locDiff = std::sqrt(locDiff);

	if(m_AdjustLocationByScale)
		locDiff /= f1.keyPoint.scale;

	if(locDiff > m_LocationThreshold)
		return false;


	return true;


}
} /* filter */ 

#endif
