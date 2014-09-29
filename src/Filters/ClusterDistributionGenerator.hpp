#ifndef CLUSTER_DISTRIBUTION_GENERATOR_HPP
#define CLUSTER_DISTRIBUTION_GENERATOR_HPP

#include "ClusterDistributionGenerator.h"
#include <itkVariableLengthVector.h>
#include <itkListSample.h>
#include <itkCovarianceSampleFilter.h>


#include <VonMisesDistribution.h>

namespace filter
{
// ------------------------------------------------------------------------
template<int TDimension>
ClusterDistributionGenerator<TDimension>::
ClusterDistributionGenerator()
{	
	this->SetNumberOfRequiredInputs(0);
}

// ------------------------------------------------------------------------
template<int TDimension>
void
ClusterDistributionGenerator<TDimension>::
SetInput(const ClusterType &cluster)
{
	m_Cluster = cluster;
}

// ------------------------------------------------------------------------
template<int TDimension>
void
ClusterDistributionGenerator<TDimension>::
Update()
{
	// first compute the mean and vairation of the locations
	unsigned int numFeatures = m_Cluster.size();
	unsigned int numHistBins = m_Cluster.keyFeature.histogram.size() * 
		m_Cluster.keyFeature.histogram.front().size();

	// set up the list samples
	typedef itk::VariableLengthVector<double> MeasurementVectorType;
	typedef itk::Statistics::ListSample<MeasurementVectorType> SampleType;
	
	
	typename SampleType::Pointer locationSamples = SampleType::New();
	typename SampleType::Pointer scaleSamples = SampleType::New();
	typename SampleType::Pointer histogramSamples = SampleType::New();

	locationSamples->SetMeasurementVectorSize(TDimension);
	scaleSamples->SetMeasurementVectorSize(1);
	histogramSamples->SetMeasurementVectorSize(numHistBins);

	std::cout << m_Cluster.keyFeature.keyPoint.angle << std::endl;
	std::cout << "----" << std::endl;


	std::vector<double> angles;

	for(unsigned int i = 0; i < numFeatures; i++)
	{
		FeatureType feature = m_Cluster.clusterMembers[i];

		MeasurementVectorType loc;
		loc.SetSize(TDimension);
		for(unsigned int j = 0; j < TDimension; j++)
		{
			loc[j] = feature.keyPoint.location[j];
		}

		locationSamples->PushBack(loc);


		angles.push_back(feature.keyPoint.angle);



		MeasurementVectorType scale;
		scale.SetSize(1);
		scale[0] = std::fabs(std::log(feature.keyPoint.scale));
		scaleSamples->PushBack(scale);


		int count = 0;
		MeasurementVectorType hist;
		hist.SetSize(numHistBins);
		for(unsigned int j = 0; j < feature.histogram.size(); j++)
		{
			for(unsigned int k = 0; k < feature.histogram.front().size(); k++)
			{
				hist[count] = feature.histogram[j][k];				
				count++;
			}
		}

		histogramSamples->PushBack(hist);

	}
	
	std::cout << locationSamples << std::endl;

	utils::VonMisesDistribution angleDist;
	angleDist.SetData(angles);
	angleDist.ComputeMean();
	angleDist.ComputeK();
	std::cout << angleDist.Evaluate(2.6) << std::endl;

	


}



} /* filter */ 
#endif
