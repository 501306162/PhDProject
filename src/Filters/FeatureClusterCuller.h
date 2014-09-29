#ifndef FEATURE_CLUSTER_CULLER_H
#define FEATURE_CLUSTER_CULLER_H

#include "FeatureCommon.h"
#include <itkProcessObject.h>

namespace filter
{
template<int TDimension>
class FeatureClusterCuller : public itk::ProcessObject
{
public:
	typedef FeatureClusterCuller Self;
	typedef itk::ProcessObject Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkTypeMacro(FeatureClusterCuller, ProcessObject);
	itkNewMacro(Self);


	typedef HoGFeatureCluster<TDimension> ClusterType;
	typedef typename ClusterType::List ClusterList;


	void SetInput(const ClusterList &input);
	ClusterList GetOutput() const;

	void Update();

	void DeDuplicate(ClusterList &input, ClusterList &output);

protected:
	FeatureClusterCuller();
	virtual ~FeatureClusterCuller() {}


private:
	
	ClusterList m_Input;
	ClusterList m_Output;


	FeatureClusterCuller(const Self&);
	void operator=(const Self&);

};
} /* filter */ 

#include "FeatureClusterCuller.hpp"

#endif
