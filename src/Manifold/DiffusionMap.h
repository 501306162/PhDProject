#ifndef DIFFUSION_MAP_H
#define DIFFUSION_MAP_H


#include "Manifold.h"

namespace manifold
{
class DiffusionMap : public Manifold
{
public:
	typedef DiffusionMap 					Self;
	typedef Manifold 						Superclass;
	typedef itk::SmartPointer<Self> 		Pointer;
	typedef itk::SmartPointer<const Self> 	ConstPointer;

	itkNewMacro(Self);
	itkTypeMacro(DiffusionMap, Manifold);

	void SetDistancePrecomputed(bool val);
	virtual void Fit(const MatrixType &X);
	virtual void Transform(const MatrixType &X, MatrixType &Y);
	virtual void SetKnn(int knn);
	void SetEpsilon(double eps);
	void SetComponents(int components);

	void GetEmbedding(MatrixType &E);
	

protected:
	DiffusionMap();

	void ComputeDistanceMatrix(const MatrixType &X, MatrixType &D);
	void MakeSymmetric(const MatrixType &W, MatrixType &W2);
	void ApplyKernelFunction(const MatrixType &M, MatrixType &W);
	void ApplyNormalisation(const MatrixType &W, MatrixType &P);
	void Decompose(const MatrixType &P, MatrixType &vecs, MatrixType &vals);
	void ComputeEmbedding(const MatrixType &vecs, const MatrixType &vals, MatrixType &E);

private:
	bool m_DistancePrecomputed;
	double m_Epsilon;
	double m_Alpha;
	MatrixType m_DistanceMatrix;
	int m_Knn;
	int m_Components;
	MatrixType m_Vecs;
	MatrixType m_Vals;
	MatrixType m_Embedding;
	MatrixType m_X;
	MatrixType m_W;
	double m_MeanDistance;

};


} /* manifold */ 

#endif
