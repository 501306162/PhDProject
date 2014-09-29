#ifndef FEATURE_CLUSTER_LEARNER_H
#define FEATURE_CLUSTER_LEARNER_H

#include <itkProcessObject.h>
#include "FeatureCommon.h"
#include "MatrixCommon.h"

namespace filter
{
template<int TDimension>
class FeatureClusterLearner : public itk::ProcessObject
{
public:
	typedef FeatureClusterLearner Self;
	typedef itk::ProcessObject Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkTypeMacro(FeatureClusterLearner, ProcessObject);
	itkNewMacro(Self);

	typedef HoGFeatureCluster<TDimension> ClusterType;
	typedef HoGFeature<TDimension> FeatureType;
	typedef typename FeatureType::List FeatureList;
	typedef utils::DoubleMatrixType MatrixType;


	void SetInput(const ClusterType &input);
	void SetFeatures(const FeatureList &features);
	void SetDistanceMatrix(const MatrixType &mat);




	itkSetMacro(EIncrement, double);
	itkSetMacro(Gamma, double);
	itkGetMacro(E, double);

	ClusterType GetOutput() const;

	void Update();

protected:
	FeatureClusterLearner();
	virtual ~FeatureClusterLearner() {}
	double ComputeGamma(double e);

	double Optimise(double e, double eInc);
	

private:

	double m_Gamma;
	double m_E;
	double m_EIncrement;
	ClusterType m_Input;
	ClusterType m_Output;
	FeatureList m_Features;
	MatrixType m_Distances;

	unsigned int m_ClusterId;
	std::vector<int> m_MemberIndices;
	std::vector<int> m_NonMemberIndices;

	std::vector<int> m_OutputMemberIds;

	FeatureClusterLearner(const Self&);
	void operator=(const Self&);
};
} /* filter */ 


#include "FeatureClusterLearner.hpp"

#endif
