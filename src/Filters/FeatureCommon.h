#ifndef FEATURE_COMMON_H
#define FEATURE_COMMON_H

#include <itkImage.h>

namespace filter
{

/**
 * Specification for a keypoint
 */
template<int TDimension>
struct OrientatedKeyPoint
{
	typedef itk::Point<double, TDimension> KeyPointLocation;
	KeyPointLocation location;
	double angle;
	double degrees;
	double scale;

	typedef std::vector<OrientatedKeyPoint> List;
};


/**
 * Specification of a histogram of gradients feature
 */
template<int TDimension>
struct HoGFeature
{
	typedef std::vector<std::vector<double> > HoGFeatureHistogramType;
	typedef OrientatedKeyPoint<TDimension> KeyPointType;
	KeyPointType keyPoint;
	int featureId;
	HoGFeatureHistogramType histogram;
	typedef std::vector<HoGFeature> List;


	bool operator< (const HoGFeature &rhs) const
	{
		return this->featureId < rhs.featureId;		
	}

	
	

};


template<int TDimension>
struct HoGFeatureCluster
{
	HoGFeatureCluster() : sorted(false) {}

	typedef HoGFeature<TDimension> FeatureType;
	typename FeatureType::List clusterMembers;
	FeatureType keyFeature;
	unsigned int clusterId;
	typedef std::vector<HoGFeatureCluster> List;
	bool sorted;

	void sort()
	{
		if(sorted) return;
		std::sort(clusterMembers.begin(), clusterMembers.end());		
		sorted=true;
	}

	bool contains(const HoGFeatureCluster &rhs) const
	{
		for(unsigned int i = 0; i < size(); i++)
		{
			if(rhs.keyFeature.featureId == clusterMembers[i].featureId)
				return true;			
		}

		return false;
	}

	bool equals(HoGFeatureCluster &rhs)
	{
		if(!rhs.sorted) rhs.sort();
		if(!sorted) sort();	

		if(rhs.size() != size()) return false;

		for(unsigned int i = 0; i < size(); i++)
		{
			if(clusterMembers[i].featureId != rhs.clusterMembers[i].featureId) 
				return false;
		}

		return true;
	}

	unsigned int size() const
	{
		return clusterMembers.size();
	}

	bool operator< (const HoGFeatureCluster &rhs) const
	{
		return this->size() > rhs.size();		
	}

};

} /* filter */ 
#endif
