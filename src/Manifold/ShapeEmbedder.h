#ifndef SHAPE_EMBEDDER_H
#define SHAPE_EMBEDDER_H

#include <itkProcessObject.h>
#include <itkImage.h>


#include "Manifold.h"
#include "ImageToImageDistanceMeasure.h"

namespace manifold
{
template<int TDimension>
class ShapeEmbedder : public itk::ProcessObject
{
public:
	typedef ShapeEmbedder Self;
	typedef itk::ProcessObject Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkTypeMacro(ShapeEmbedder, ProcessObject);
	itkNewMacro(Self);

	typedef itk::Image<unsigned char, TDimension> InputType;

	typedef itk::Image<double, TDimension> RealImageType;
	typedef typename RealImageType::Pointer RealImagePointer;
	typedef std::vector<RealImagePointer> DistanceMapList;


	/** typedefs for the manifold and distance measure */
	typedef ImageToImageDistanceMeasure<RealImageType> DistanceMeasureType;
	typedef typename DistanceMeasureType::Pointer DistanceMeasurePointer;

	typedef std::pair<int, double> IndexDistancePair;
	typedef std::vector<IndexDistancePair> IndexDistanceList;
	typedef std::vector<double> WeightListType;
	typedef std::vector<int> IndexListType;


	typedef Manifold::Pointer ManifoldPointer;
	typedef Manifold::MatrixType MatrixType;

	void SetInput(const InputType * input);

	virtual void Update();

	itkSetObjectMacro(DistanceMeasure, DistanceMeasureType);
	itkSetObjectMacro(Manifold, Manifold);
	itkSetMacro(Knn, unsigned int);
	
	void SetDistanceMaps(const DistanceMapList &maps);
	IndexListType GetIndices() const { return m_Indices; }
	WeightListType GetWeights() const { return m_Weights; }

protected:
	ShapeEmbedder();
	virtual ~ShapeEmbedder() {}


	void ComputeDsistanceMatrix(const RealImageType * input,
			MatrixType &mat);

	// function to compute the distance map
	void ComputeDistanceMap(const InputType * input, 
			RealImagePointer &output);

	void ComputeKnn(MatrixType &X, MatrixType &Y, IndexDistanceList &knn);

	void ProcessOutput(IndexDistanceList &input);
	


private:

	unsigned int m_Knn;
	DistanceMeasurePointer m_DistanceMeasure;
	ManifoldPointer m_Manifold;
	DistanceMapList m_DistanceMaps;

	ShapeEmbedder(const Self&);
	void operator=(const Self&);


	WeightListType m_Weights;
	IndexListType m_Indices;

};
} /* manifold */ 


#include "ShapeEmbedder.hpp"

#endif
