#include "InitialParticleGenerator.h"

#include <itkNormalVariateGenerator.h>
#include <time.h>

namespace vt
{
// ------------------------------------------------------------------------
void InitialParticleGenerator::Generate()
{

	unsigned int minPoints = 10000;
	std::vector<std::string> names;
	std::vector<MatrixType> matrixes;

	std::map<std::string, ImagePointPair>::iterator mapIt = m_Data.begin();
	while(mapIt != m_Data.end())
	{
		MatrixType a, b;
		GetAlignedPoints3(mapIt->second.second.first, mapIt->second.first, a);
		GetAlignedPoints3(mapIt->second.second.second, mapIt->second.first, b);

		if(a.cols() < minPoints)
			minPoints = a.cols();

		names.push_back(mapIt->first);
		matrixes.push_back(a);
		matrixes.push_back(b);

		++mapIt;
	}



	const unsigned int numRows = 3*matrixes.size();
	const unsigned int numCols = minPoints;
	std::cout << numRows << " " << numCols << std::endl;

	// create the data matrix 
	MatrixType X(numRows, numCols);
	for(unsigned int i = 0; i < matrixes.size(); i++)
	{
		for(unsigned int j = 0; j < minPoints; j++)
		{
			for(unsigned int k = 0; k < 3; k++)
			{
				X(i*3+k,j) = matrixes[i](k,j);		
			}
		}
	}


	// compute the mean and covariance 
	MatrixType mean = X.rowwise().mean();

	// subtract the mean from the points and get covariance
	MatrixType centered = X.colwise() - X.rowwise().mean();
	MatrixType C = centered*centered.transpose() / double(X.cols() - 1);

	// get the eigen vectors 
	Eigen::SelfAdjointEigenSolver<MatrixType> solver(C);
	MatrixType rot = solver.eigenvectors();
	MatrixType scl = solver.eigenvalues();

	for(unsigned int i = 0; i < scl.rows(); i++)
	{
		scl(i,0) = std::sqrt(scl(i,0));
		if(std::isnan(scl(i,0))) scl(i,0) = 0.0;
	}

	typedef itk::Statistics::NormalVariateGenerator RNGType;
	RNGType::Pointer generator = RNGType::New();
	generator->Initialize((int) time(NULL));

	for(unsigned int i = 0; i < m_NumParticles; i++)
	{
		MatrixType sampleVec(numRows,1);
		for(unsigned int j = 0; j < numRows; j++)
		{
			 sampleVec(j,0) = generator->GetVariate() * 1.0 * scl(j,0);
		}		

		sampleVec = rot*sampleVec + mean;


		for(unsigned int n = 0; n < names.size(); n++)
		{
			Particle particle;
			for(unsigned int j = 0; j < 3; j++)
			{
				particle.p1[j] = sampleVec(n*6+j,0);
				particle.p2[j] = sampleVec(n*6+j+3,0);
			}			

			particle.length = particle.p1.EuclideanDistanceTo(particle.p2);
			particle.direction = particle.p2-particle.p1;

			m_MapOut[names[n]].push_back(particle);


		}
	}

	/*
	MatrixType a1, a2;
	GetAlignedPoints(m_Points1, a1);
	GetAlignedPoints(m_Points2, a2);

	PointList p1List, p2List;
	GetSamples(a1, p1List);
	GetSamples(a2, p2List);
	//GetSamples2(a1,a2,p1List,p2List);



	// combine the points to create the particles 
	for(unsigned int i = 0; i < m_NumParticles; i++)
	{
		Particle particle;
		particle.p1 = p1List[i];
		particle.p2 = p2List[i];
		particle.length = particle.p1.EuclideanDistanceTo(particle.p2);
		particle.direction = particle.p2-particle.p1;

		m_Output.push_back(particle);
	}
	*/

	/*
	MatrixType a1, a2;
	GetAlignedPoints(m_Points1, a1);
	GetAlignedPoints(m_Points2, a2);

	PointList p1List, p2List;
	GetSamples2(a1,a2,p1List,p2List);



	// combine the points to create the particles 
	for(unsigned int i = 0; i < m_NumParticles; i++)
	{
		Particle particle;
		particle.p1 = p1List[i];
		particle.p2 = p2List[i];
		particle.length = particle.p1.EuclideanDistanceTo(particle.p2);
		particle.direction = particle.p2-particle.p1;

		m_Output.push_back(particle);
	}
	*/
}

// ------------------------------------------------------------------------
void InitialParticleGenerator::GetSamples2(const MatrixType &aligned1, const MatrixType &aligned2, 
		PointList &samples1, PointList &samples2)
{
	MatrixType X(6, aligned1.cols());

	for(unsigned int i = 0; i < aligned1.cols(); i++)
	{
		for(unsigned int j = 0; j < 3; j++)
		{
			X(j,i) = aligned1(j,i);
			X(3+j,i) = aligned2(j,i);
		}
	}

	// compute the mean and covariance 
	MatrixType mean = X.rowwise().mean();
	std::cout << mean << std::endl;

	// subtract the mean from the points and get covariance
	MatrixType centered = X.colwise() - X.rowwise().mean();
	MatrixType C = centered*centered.transpose() / double(X.cols() - 1);

	// get the eigen vectors 
	Eigen::SelfAdjointEigenSolver<MatrixType> solver(C);
	MatrixType rot = solver.eigenvectors();
	MatrixType scl = solver.eigenvalues();
	std::cout << rot << std::endl;

	for(unsigned int i = 0; i < 6; i++)
	{
		scl(i,0) = std::sqrt(scl(i,0));
		if(std::isnan(scl(i,0))) scl(i,0) = 0.0;
	}


	

	typedef itk::Statistics::NormalVariateGenerator RNGType;
	RNGType::Pointer generator = RNGType::New();
	generator->Initialize((int) time(NULL));

	for(unsigned int i = 0; i < m_NumParticles; i++)
	{
		MatrixType sampleVec(6,1);
		for(unsigned int j = 0; j < 6; j++)
		{
			 sampleVec(j,0) = generator->GetVariate() * 1.0 * scl(j,0);
		}		

		sampleVec = rot*sampleVec + mean;

		PointType outP1, outP2;
		for(unsigned int j = 0; j < 3; j++)
		{
			outP1[j] = sampleVec(j,0);			
			outP2[j] = sampleVec(3+j,0);			
		}
		samples1.push_back(outP1);
		samples2.push_back(outP2);
	}
}


// ------------------------------------------------------------------------
void InitialParticleGenerator::GetSamples(const MatrixType &aligned, PointList &samples)
{
	// compute the mean and covariance 
	MatrixType mean = aligned.rowwise().mean();

	// subtract the mean from the points and get covariance
	MatrixType centered = aligned.colwise() - aligned.rowwise().mean();
	MatrixType C = centered*centered.transpose() / double(aligned.cols() - 1);

	// get the eigen vectors 
	Eigen::SelfAdjointEigenSolver<MatrixType> solver(C);
	MatrixType rot = solver.eigenvectors();
	MatrixType scl = solver.eigenvalues();
	std::cout << rot << std::endl;

	for(unsigned int i = 0; i < 3; i++)
	{
		scl(i,0) = std::sqrt(scl(i,0));
	}
	scl(0,0) = 0.0;

	

	typedef itk::Statistics::NormalVariateGenerator RNGType;
	RNGType::Pointer generator = RNGType::New();
	generator->Initialize((int) time(NULL));

	for(unsigned int i = 0; i < m_NumParticles; i++)
	{
		MatrixType sampleVec(3,1);
		for(unsigned int j = 0; j < 3; j++)
		{
			 sampleVec(j,0) = generator->GetVariate() * 1.0 * scl(j,0);
		}		

		sampleVec = rot*sampleVec + mean;

		PointType outP;
		for(unsigned int j = 0; j < 3; j++)
		{
			outP[j] = sampleVec(j,0);			
		}
		samples.push_back(outP);
	}
}

// ------------------------------------------------------------------------
void InitialParticleGenerator::GetAlignedPoints3(const MatrixType &input, ImageType::Pointer &image, MatrixType &aligned)
{
	aligned = MatrixType(3, input.rows());

	for(unsigned int i = 0; i < input.rows(); i++)
	{
		PointType p;
		for(unsigned int j = 0; j < 3; j++)
		{
			p[j] = input(i,j);
		}

		p = ProjectPoint2(m_Transform->TransformPoint(p), image);

		for(unsigned int j = 0; j < 3; j++)
		{
			aligned(j,i) = p[j];
		}
	}

}



// ------------------------------------------------------------------------
void InitialParticleGenerator::GetAlignedPoints(const MatrixType &input, MatrixType &aligned)
{

	aligned = MatrixType(3, input.rows());

	for(unsigned int i = 0; i < input.rows(); i++)
	{
		PointType p;
		for(unsigned int j = 0; j < 3; j++)
		{
			p[j] = input(i,j);
		}

		p = ProjectPoint(m_Transform->TransformPoint(p));

		for(unsigned int j = 0; j < 3; j++)
		{
			aligned(j,i) = p[j];
		}
	}
}


// ------------------------------------------------------------------------
void InitialParticleGenerator::GetAlignedPoints2(const MatrixType &input, MatrixType &aligned)
{

	aligned = MatrixType(3, input.rows());

	for(unsigned int i = 0; i < input.rows(); i++)
	{
		PointType p;
		for(unsigned int j = 0; j < 3; j++)
		{
			p[j] = input(i,j);
		}

		p = m_Transform->TransformPoint(p);

		for(unsigned int j = 0; j < 3; j++)
		{
			aligned(j,i) = p[j];
		}
	}
}

// ------------------------------------------------------------------------
InitialParticleGenerator::PointType InitialParticleGenerator::ProjectPoint2(const PointType &p, 
		const ImageType::Pointer &image)
{
	VectorType v1, normal, origin;
	for(unsigned int i = 0; i < 3; i++)
	{
		v1[i] = p[i];
		normal[i] = image->GetDirection()(i,2);
		origin[i] = image->GetOrigin()[i];
	}


	VectorType dir = v1-origin;
	double dist = dir * normal;
   	VectorType pos = v1 - (dist * normal);	

	PointType out;
	for(unsigned int i = 0; i < 3; i++)
	{
		out[i] = pos[i];
	}
	return out;
}

// ------------------------------------------------------------------------
InitialParticleGenerator::PointType InitialParticleGenerator::ProjectPoint(const PointType &p)
{
	VectorType v1, normal, origin;
	for(unsigned int i = 0; i < 3; i++)
	{
		v1[i] = p[i];
		normal[i] = m_Image->GetDirection()(i,2);
		origin[i] = m_Image->GetOrigin()[i];
	}

	VectorType dir = v1-origin;
	double dist = dir * normal;
   	VectorType pos = v1 - (dist * normal);	

	PointType out;
	for(unsigned int i = 0; i < 3; i++)
	{
		out[i] = pos[i];
	}
	return out;
}

}
