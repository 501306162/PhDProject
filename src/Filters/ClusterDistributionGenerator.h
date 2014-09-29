#ifndef CLUSTER_DISTRIBUTION_GENERATOR_H
#define CLUSTER_DISTRIBUTION_GENERATOR_H

#include <itkProcessObject.h>
#include <FeatureCommon.h>

namespace filter
{
template<int TDimension>
class ClusterDistributionGenerator : public itk::ProcessObject
{
public:
	typedef ClusterDistributionGenerator Self;
	typedef itk::ProcessObject Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkTypeMacro(ClusterDistributionGenerator, ProcessObject);
	itkNewMacro(Self);

	typedef HoGFeatureCluster<TDimension> ClusterType;
	typedef typename ClusterType::FeatureType FeatureType;


	void SetInput(const ClusterType &cluster);
	void Update();

protected:
	ClusterDistributionGenerator();
	virtual ~ClusterDistributionGenerator() {}
	


private:
		
	ClusterType m_Cluster;
	
	ClusterDistributionGenerator(const Self&);
	void operator=(const Self&);
};

} /* filter */ 

#include "ClusterDistributionGenerator.hpp"

#endif
