#include "DiffusionMap.h"

#include <itkVector.h>
#include <itkVariableLengthVector.h>
#include <itkListSample.h>
#include <itkKdTreeGenerator.h>
#include <itkEuclideanDistanceMetric.h>


namespace manifold
{

// ------------------------------------------------------------------------
DiffusionMap::DiffusionMap()
{
	m_DistancePrecomputed = false;
	m_Knn = -1;
	m_Epsilon = 1.0;
	m_Components = 3;
}

// ------------------------------------------------------------------------
void DiffusionMap::Fit(const MatrixType &X)
{
	if(!m_DistancePrecomputed)
	{
		ComputeDistanceMatrix(X, m_DistanceMatrix);
	}
	else
	{
		m_DistanceMatrix = X;
	}

	m_X = X;

	m_MeanDistance = m_DistanceMatrix.mean();

	// make the matrix symmetric
	MatrixType symmetric;
	MakeSymmetric(m_DistanceMatrix, symmetric);


	// apply the gaussian kernel
	MatrixType W;
	ApplyKernelFunction(symmetric, W);
	
	
	m_W = W;
	
	MatrixType P;
	ApplyNormalisation(W, P);


	MatrixType Y;
	Decompose(P, m_Vecs, m_Vals);

	ComputeEmbedding(m_Vecs, m_Vals, m_Embedding);
	
}


// ------------------------------------------------------------------------
void DiffusionMap::GetEmbedding(MatrixType &E)
{
	E = m_Embedding;
}

// ------------------------------------------------------------------------
void DiffusionMap::ComputeEmbedding(const MatrixType &vecs, const MatrixType &vals, MatrixType &E)
{
	E = MatrixType::Zero(vecs.rows(), m_Components);
	for(int i = 0; i < m_Components; i++)
	{
		E.col(i) = vecs.col(i+1) * vals(1+i,0);
	}
}

// ------------------------------------------------------------------------
void DiffusionMap::Decompose(const MatrixType &P, MatrixType &vecs, MatrixType &vals)
{
	Eigen::EigenSolver<MatrixType> solver(P);
	if(solver.info() != Eigen::Success)
	{
		std::cout << "Oh Dear" << std::endl;
	}

	vecs = solver.eigenvectors().real();
	vals = solver.eigenvalues().real();

	
}


// ------------------------------------------------------------------------
void DiffusionMap::ApplyNormalisation(const MatrixType &W, MatrixType &P)
{

	MatrixType Q = W.rowwise().sum();
	MatrixType W2 = MatrixType::Zero(W.rows(),W.cols());
	for(int i = 0; i < W.rows(); i++)
	{
		for(int j = 0; j < W.rows(); j++)
		{
			double dev = pow(Q(i,0),-1.0) * pow(Q(j,0),-1.0);
			W2(i,j) = W(i,j) / dev;
		}
	}

	MatrixType D = W2.rowwise().sum();
	P = MatrixType::Zero(W.rows(), W.cols());
	
	for(int i = 0; i < W.rows(); i++)
	{
		for(int j = 0; j < W.rows(); j++)
		{
			P(i,j) = W2(i,j) / D(i,0);
		}
	}

	/*

	// first normalisation
	MatrixType D = W.rowwise().sum();

	for(int i = 0; i < D.rows(); i++)
	{
		if(D(i,0) != 0.0)
		{
			D2(i,i) = pow(D(i,0),-1.0);
		}
	}

	MatrixType K = (D2*W)*D2;
	
	// second normalisation
	D = K.rowwise().sum();
	for(int i = 0; i < D.rows(); i++)
	{
		if(D(i,0) != 0.0)
		{
			D2(i,i) = pow(D(i,0),-1.0);
		}
	}

	P = K*D2;
	*/
}

// ------------------------------------------------------------------------
void DiffusionMap::ApplyKernelFunction(const MatrixType &M, MatrixType &W)
{
	W = MatrixType::Zero(M.rows(), M.cols());
	for(unsigned int i = 0; i < M.rows(); i++)
	{
		for(unsigned int j = 0; j < M.cols(); j++)
		{
			double w = (M(i,j) * M(i,j)) / (2*m_MeanDistance * m_MeanDistance);			
			W(i,j) = exp(-w);
		}
	}
}

// ------------------------------------------------------------------------
void DiffusionMap::SetEpsilon(double eps)
{
	m_Epsilon = eps;
}

// ------------------------------------------------------------------------
void DiffusionMap::SetComponents(int components)
{
	m_Components = components;
}

// ------------------------------------------------------------------------
void DiffusionMap::MakeSymmetric(const MatrixType &W, MatrixType &WS)
{
	WS = W+W.transpose();
}


// ------------------------------------------------------------------------
void DiffusionMap::Transform(const MatrixType &X, MatrixType &Y)
{
	std::cout << X << std::endl;
	int numNewSamples = X.rows();
	
	if(X.cols() != m_X.cols())
	{
		std::cout << "New Data has different col number" << std::endl;
	}	

	int numSamples = m_X.rows();
	MatrixType W = MatrixType::Zero(numNewSamples, numSamples);

	// create a new P
	for(int i = 0; i < numNewSamples; i++)
	{
		MatrixType xi = X.row(i);
		for(int j = 0; j < numSamples; j++)
		{
			double dist = 0.0;
			MatrixType xj = m_X.row(j);
			if(m_DistancePrecomputed)
			{
				dist = X(i,j);								
			}
			else
			{
				dist = (xi-xj).norm();
			}
			double w = (dist*dist) / (2*m_MeanDistance * m_MeanDistance);
			W(i,j) = exp(-w);
		}
	}
	
	
	MatrixType Q = W.rowwise().sum();
	MatrixType Q2 = m_W.rowwise().sum();
	MatrixType W2 = MatrixType::Zero(W.rows(), W.cols());

	for(int i = 0; i < numNewSamples; i++)
	{
		for(int j = 0; j < numSamples; j++)
		{
			double dev = pow(Q(i,0), -1.0) * pow(Q2(j,0), -1.0);
			W2(i,j) = W(i,j) / dev;
		}
	}

	MatrixType D = W2.rowwise().sum();
	MatrixType P = MatrixType::Zero(W.rows(), W.cols());
	for(int i = 0; i < numNewSamples; i++)
	{
		for(int j = 0; j < numSamples; j++)
		{
			P(i,j) = W2(i,j) / D(i,0);
		}
	}

	Y = MatrixType::Zero(numNewSamples, m_Components);
	for( int m = 0; m < m_Components; m++ )
	{
		for( int i = 0; i < numNewSamples; i++ )
		{
			for( int j = 0; j < numSamples; j++ )
			{
				Y(i,m) += P(i,j)*m_Vecs(j,m+1);
			}		
		}
	}	
}


// ------------------------------------------------------------------------
void DiffusionMap::SetDistancePrecomputed(bool val)
{
	m_DistancePrecomputed = val;
}

void DiffusionMap::SetKnn(int knn)
{
	m_Knn = knn;
}


void DiffusionMap::ComputeDistanceMatrix(const MatrixType &X, MatrixType &D)
{
	// build the measurement vector
	unsigned int cols = X.cols();
	unsigned int rows = X.rows();
	typedef itk::VariableLengthVector<double> VectorType;
	typedef itk::Statistics::ListSample<VectorType> SampleType;

	// create the sample
	SampleType::Pointer sample = SampleType::New();
	sample->SetMeasurementVectorSize(cols);


	for(unsigned int i = 0; i < rows; i++)
	{
		VectorType mv;
		mv.SetSize(cols);
		for(unsigned int j = 0; j < cols; j++)
		{
			mv[j] = X(i,j);
		}

		sample->PushBack(mv);

	}

	// create the Kd tree
	typedef itk::Statistics::KdTreeGenerator<SampleType> TreeGenerator;
	TreeGenerator::Pointer generator = TreeGenerator::New();
	generator->SetSample(sample);
	generator->SetBucketSize(16);
	generator->Update();

	typedef TreeGenerator::KdTreeType TreeType;
	typedef TreeType::NearestNeighbors NeighboursType;
	typedef TreeType::KdTreeNodeType NodeType;
	TreeType::Pointer tree = generator->GetOutput();




	// find the neighbours
	unsigned int numNeighbours = rows;
	if(m_Knn > 0)
		numNeighbours = m_Knn;

	TreeType::InstanceIdentifierVectorType neighbours;


	// initialise the distance matrix and fill it up
	D = MatrixType::Zero(rows, rows);

	for(unsigned int i = 0; i < rows; i++)
	{
		VectorType mv = sample->GetMeasurementVector(i);
		tree->Search(mv, numNeighbours, neighbours);
		for(unsigned int j = 0; j < neighbours.size(); j++)
		{

			// compute the distances
			typedef itk::Statistics::EuclideanDistanceMetric<VectorType> DistanceMeasure;
			DistanceMeasure::Pointer d2 = DistanceMeasure::New();
			double dist = d2->Evaluate(mv, sample->GetMeasurementVector(neighbours[j]));

			// set the distance value
			D(i, neighbours[j]) = dist;
		}
	}


}

}
