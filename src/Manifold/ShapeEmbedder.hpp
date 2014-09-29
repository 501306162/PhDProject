#ifndef SHAPE_EMBEDDER_HPP
#define SHAPE_EMBEDDER_HPP

#include "ShapeEmbedder.h"
#include <itkSignedMaurerDistanceMapImageFilter.h>

namespace manifold
{
	
// ------------------------------------------------------------------------
bool my_compare(const std::pair<int, double> &a, const std::pair<int, double> &b)
{
	return a.second < b.second;
}


// ------------------------------------------------------------------------
template<int TDimension>
void
ShapeEmbedder<TDimension>::Update()
{
	// do some checks
	if(!m_Manifold)
	{
		itkExceptionMacro(<< "Manifold Has Not been set");
	}

	if(!m_DistanceMeasure)
	{
		itkExceptionMacro(<< "Distance Measure has not been set");
	}

	if(m_DistanceMaps.size() == 0)
	{
		itkExceptionMacro(<< "Distance Maps Have not been set");
	}


	if(m_Knn == 0)
	{
		itkExceptionMacro(<< "Knn has not been set");
	}

	// get the input and compute the distance map
	typename RealImageType::Pointer distanceMap = RealImageType::New();
	ComputeDistanceMap(itkDynamicCastInDebugMode<InputType*>(this->GetInput(0)), distanceMap);


	// compute the distance matrix and then get the embedding
	MatrixType distances;
	ComputeDsistanceMatrix(distanceMap, distances);

	MatrixType Y;
	MatrixType X;
	m_Manifold->GetEmbedding(X);
	m_Manifold->Transform(distances, Y);


	// get the knn
	IndexDistanceList knn;
	ComputeKnn(X, Y, knn);
	
	
	ProcessOutput(knn);
}



// ------------------------------------------------------------------------
template<int TDimension>
void 
ShapeEmbedder<TDimension>::
ProcessOutput(IndexDistanceList &input)
{

	double distanceSum = 0.0;
	
	for(unsigned int i = 0; i < input.size(); i++)
	{
		distanceSum += 1.0 / input[i].second;
		m_Indices.push_back(input[i].first);
	}

	for(unsigned int i = 0; i < input.size(); i++)
	{
		m_Weights.push_back( (1.0/input[i].second) / distanceSum );
	}
}


// ------------------------------------------------------------------------
template<int TDimension>
void 
ShapeEmbedder<TDimension>::
ComputeKnn(MatrixType &X, MatrixType &Y, IndexDistanceList &knn)
{
	// compute the distances to all the values
	IndexDistanceList sortList;

	for(unsigned int i = 0; i < X.rows(); i++)
	{
		double dist = (X.row(i)-Y.row(0)).norm();
		IndexDistancePair p(i, dist);
		sortList.push_back(p);
	}

	// sort the list
	std::sort(sortList.begin(), sortList.end(), my_compare);

	for(unsigned int i = 0; i < m_Knn; i++)
	{
		knn.push_back(sortList[i]);
	}
}

// ------------------------------------------------------------------------
template<int TDimension>
void 
ShapeEmbedder<TDimension>::
SetDistanceMaps(const DistanceMapList &maps)
{
	m_DistanceMaps = maps;
}


// ------------------------------------------------------------------------
template<int TDimension>
void 
ShapeEmbedder<TDimension>::
ComputeDsistanceMatrix(const RealImageType * input,
			MatrixType &mat)
{
	// compute the sub distance matrix
	mat = MatrixType::Zero(1, m_DistanceMaps.size());
	for(unsigned int i = 0; i < m_DistanceMaps.size(); i++)
	{
		m_DistanceMeasure->SetInput1(input);
		m_DistanceMeasure->SetInput2(m_DistanceMaps[i]);
		mat(0,i) = m_DistanceMeasure->GetValue();
	}


}



// ------------------------------------------------------------------------
template<int TDimension>
void 
ShapeEmbedder<TDimension>::
ComputeDistanceMap(const InputType * input, RealImagePointer &output)
{
	// compute the distance map
	typedef itk::SignedMaurerDistanceMapImageFilter<InputType, RealImageType> DistanceFilterType;
	typename DistanceFilterType::Pointer distanceFilter = DistanceFilterType::New();
	distanceFilter->SetInput(input);
	distanceFilter->SetUseImageSpacing(true);
	distanceFilter->SetSquaredDistance(false);
	distanceFilter->SetBackgroundValue(0);
	distanceFilter->Update();

	output = distanceFilter->GetOutput();

}

// ------------------------------------------------------------------------
template<int TDimension>
void 
ShapeEmbedder<TDimension>::
SetInput(const InputType * input)
{
	this->SetNthInput(0, const_cast<InputType *>(input));
}

// ------------------------------------------------------------------------
template<int TDimension>
ShapeEmbedder<TDimension>::
ShapeEmbedder()
{
	this->SetNumberOfRequiredInputs(1);
	m_Manifold = 0;
	m_DistanceMeasure = 0;
	m_Knn = 0;
}

} /* manifold */ 

#endif
