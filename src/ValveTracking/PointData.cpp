#include "PointData.h"
#include <Directory.h>
#include <ValveIO.h>
#include <QString>
#include <itkNormalVariateGenerator.h>


namespace vt
{
// ------------------------------------------------------------------------
PointData::Pointer PointData::Load(const std::string &folder, const unsigned int exclude, bool flat)
{
	PointData::Pointer obj = PointData::New();
	if(!flat)
		obj->LoadData(folder, exclude);
	else
		obj->LoadDataFlat(folder, exclude);
	return obj;
}

// ------------------------------------------------------------------------
void PointData::LoadDataFlat(const std::string &folder, const unsigned int exclude)
{
	m_Exclude = "d" + QString::number(exclude).toStdString() + ".txt";

	// get the subdirectories 
	utils::Directory::FilenamesType subDirs = utils::Directory::GetDirectories(folder);
	for(unsigned int i = 0; i < subDirs.size(); i++)
	{
		const std::string folderName = subDirs[i];
		const std::string valveType = utils::Directory::GetFileName(folderName);

		if(valveType.find("MV") == std::string::npos && valveType.find("TP") == std::string::npos) continue;
		std::cout << "Loading: " << valveType << std::endl;
			
		ValvePoints points;
		LoadFolderFlat(folderName, points);

		m_Data[valveType] = points;

		PointDistributionList pdList;
		ComputeDistributions(points, pdList);

		m_Distributions[valveType] = pdList;


	}



	// compute the covs and stuff
}

// ------------------------------------------------------------------------
void PointData::ComputeDistributions(const ValvePoints &points, PointDistributionList &dists)
{	
	for(unsigned int timeStep = 0; timeStep < points.size(); timeStep++)
	{
		MatrixType p1 = points[timeStep].first;
		MatrixType p2 = points[timeStep].second;	

		Distribution dist1, dist2;
		GetDist(p1, dist1);
		GetDist(p2, dist2);

		PointDistribution pd;
		pd.first = dist1;
		pd.second = dist2;

		dists.push_back(pd);

	}


}

// ------------------------------------------------------------------------
void PointData::GetDist(const MatrixType &points, Distribution &dist)
{
	// reshape the points
	MatrixType X = points.transpose();

	// compute the mean and covariance 
	MatrixType mean = X.rowwise().mean();

	// subtract the mean from the points and get covariance
	MatrixType centered = X.colwise() - X.rowwise().mean();
	MatrixType C = centered*centered.transpose() / double(X.cols() - 1);

	// get the eigen vectors 
	Eigen::SelfAdjointEigenSolver<MatrixType> solver(C);
	MatrixType rot = solver.eigenvectors();
	MatrixType scl = solver.eigenvalues();

	for(unsigned int i = 0; i < 2; i++)
	{
		scl(i,0) = std::sqrt(scl(i,0));
	}

	dist.scl = scl;
	dist.rot = rot;
	dist.mean = mean;
}

// ------------------------------------------------------------------------
void PointData::GetUpdate(const std::string &type, const unsigned int timeStep,
		PointType &p1, PointType &p2)
{
	PointDistribution pd = m_Distributions[type][timeStep];

	GetSample(pd.first, p1);
	GetSample(pd.second, p2);
}

// ------------------------------------------------------------------------
void PointData::GetSample(const Distribution & dist, PointType &sample)
{

	MatrixType sampleVec(2,1);
	for(unsigned int j = 0; j < 2; j++)
	{
		sampleVec(j,0) = m_Generator->GetVariate() * 3.0 * dist.scl(j,0);
	}		

	sampleVec = dist.rot*sampleVec + dist.mean;

	//std::cout << sampleVec << std::endl;

	for(unsigned int j = 0; j < 2; j++)
	{
		sample[j] = sampleVec(j,0);			
	}

	sample[2] = 0;
}


// ------------------------------------------------------------------------
void PointData::LoadData(const std::string &folder, const unsigned int exclude)
{
	m_Exclude = "d" + QString::number(exclude).toStdString() + ".txt";

	// get the subdirectories 
	utils::Directory::FilenamesType subDirs = utils::Directory::GetDirectories(folder);
	for(unsigned int i = 0; i < subDirs.size(); i++)
	{
		const std::string folderName = subDirs[i];
		const std::string valveType = utils::Directory::GetFileName(folderName);

		if(valveType.find("MV") == std::string::npos) continue;
		std::cout << "Loading: " << valveType << std::endl;
			
		ValvePoints points;
		LoadFolder(folderName, points);

		m_Data[valveType] = points;
	}
}

// ------------------------------------------------------------------------
void PointData::GetMembers(const std::string type, const unsigned int timeStep,
		const PointType &p1, const PointType &p2, MembershipType::Pointer &p1Member, 
		MembershipType::Pointer &p2Member)
{
	ValvePoints vpoints = m_Data[type];
	PointPairType points = vpoints[timeStep];

	// compute the transform 
	TransformType::Pointer transform = TransformType::New();
	TransformToOrigin(p1, p2, transform);

	GetMemberFunction(points.first, transform, p1Member);
	GetMemberFunction(points.second, transform, p2Member);
	

}

// ------------------------------------------------------------------------
void PointData::GetMemberFunction(const MatrixType &points, const TransformType::Pointer &transform,
		MembershipType::Pointer &member)
{
	SampleType::Pointer sample = SampleType::New();
	for(unsigned int i = 0; i < points.rows(); i++)
	{
		PointType p;
		for(unsigned int j = 0; j < 3; j++)
		{
			p[j] = points(i,j);
		}		

		VectorType mv;
		p = transform->GetInverseTransform()->TransformPoint(p);
		for(unsigned int j = 0; j < 3; j++)
		{
			mv[j] = p[j];
		}

		sample->PushBack(mv);
	}

	CovFilterType::Pointer covFilter = CovFilterType::New();
	covFilter->SetInput(sample);
	covFilter->Update();

	member = MembershipType::New();
	member->SetMean(covFilter->GetMean());
	member->SetCovariance(covFilter->GetCovarianceMatrix());
}




// ------------------------------------------------------------------------
void PointData::LoadFolderFlat(const std::string &folderName, ValvePoints &points)
{
	// extract the files in this folder
	utils::Directory::FilenamesType filenames = utils::Directory::GetFiles(folderName, ".txt");
	

	std::vector<ValveSequence<3>::Pointer> sequences;
	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		const std::string filename = filenames[i];

		// check for excluding
		if(filename.find(m_Exclude) != std::string::npos)
			continue;
		


		ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
		reader->SetFileName(filename);
		ValveSequence<3>::Pointer sequence = reader->GetOutput();

		if(sequence->GetNumberOfLines() != 25) continue;

		sequences.push_back(sequence);
				
	}



	// loop through all the time steps
	for(unsigned int timeStep = 0; timeStep < 25; timeStep++)
	{
		const unsigned int nextTimeStep = (timeStep+1)%25;		
		
		PointPairType vpoints;
		vpoints.first = MatrixType::Zero(sequences.size(), 2);
		vpoints.second = MatrixType::Zero(sequences.size(), 2); 


		for(unsigned int seqNum = 0; seqNum <sequences.size(); seqNum++)
		{
			ValveLine<3>::Pointer v1 = sequences[seqNum]->GetValveLine(timeStep);
			PointType p1 = v1->GetP1();
			PointType p2 = v1->GetP2();
			
			typedef ValveLine<3>::ImageType ImageType;
			ImageType::Pointer image = v1->GetImage();


			// check the trace 
			double trace = 0.0;
			for(unsigned int i = 0; i < 3; i++)
			{
				trace += image->GetDirection()(i,i);
			}
				
			VectorType x = p2-p1;
			x.Normalize();
			VectorType z;
			for(unsigned int i = 0; i < 3; i++)
			{
				z[i] = image->GetDirection()(i,2);
			}

			VectorType y = itk::CrossProduct(x,z);
			if(trace > 0) y = -y;
			
			ValveLine<3>::Pointer v2 = sequences[seqNum]->GetValveLine(nextTimeStep);
			PointType pp1 = v2->GetP1();
			PointType pp2 = v2->GetP2();

			vpoints.first(seqNum, 0) = x * (pp1-p1);				
			vpoints.first(seqNum, 1) = y * (pp1-p1);				
			vpoints.second(seqNum, 0) = x * (pp2-p2);				
			vpoints.second(seqNum, 1) = y * (pp2-p2);				
		}

		points.push_back(vpoints);
	}
}




// ------------------------------------------------------------------------
void PointData::LoadFolder(const std::string &folderName, ValvePoints &points)
{
	// extract the files in this folder
	utils::Directory::FilenamesType filenames = utils::Directory::GetFiles(folderName, ".txt");
	

	std::vector<ValveSequence<3>::Pointer> sequences;
	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		const std::string filename = filenames[i];

		// check for excluding
		if(filename.find(m_Exclude) != std::string::npos)
			continue;

		ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
		reader->SetFileName(filename);
		ValveSequence<3>::Pointer sequence = reader->GetOutput();
		sequences.push_back(sequence);
				
	}


	// loop through all the time steps
	for(unsigned int timeStep = 0; timeStep < 25; timeStep++)
	{
		const unsigned int nextTimeStep = (timeStep+1)%25;		
		
		PointPairType vpoints;
		vpoints.first = MatrixType::Zero(sequences.size(), 3);
		vpoints.second = MatrixType::Zero(sequences.size(), 3); 


		for(unsigned int seqNum = 0; seqNum <sequences.size(); seqNum++)
		{
			ValveLine<3>::Pointer v1 = sequences[seqNum]->GetValveLine(timeStep);
			TransformType::Pointer transform = TransformType::New();
			TransformToOrigin(v1->GetP1(), v1->GetP2(), transform);

			ValveLine<3>::Pointer v2 = sequences[seqNum]->GetValveLine(nextTimeStep);
			PointType p1 = transform->TransformPoint(v2->GetP1());
			PointType p2 = transform->TransformPoint(v2->GetP2());

			for(unsigned int i = 0; i < 3; i++)
			{
				vpoints.first(seqNum, i) = p1[i];				
				vpoints.second(seqNum, i) = p2[i];				
			}
		}

		points.push_back(vpoints);
	}
}

// ------------------------------------------------------------------------
void PointData::TransformToOrigin(const PointType &p1, const PointType &p2, TransformType::Pointer &transform)
{
	VectorType translation;
	for(unsigned int i = 0; i < 3; i++)
	{
		translation[i] = -p1[i];
	}

	// compute the scaling 
	double scale =  1.0 / (p2-p1).GetNorm();

	// compute the rotation
	VectorType direction = p2-p1;
	direction.Normalize();
	VectorType xAxis;
	xAxis.Fill(0);
	xAxis[0] = 1.0;
	TransformType::OutputVectorType axis = itk::CrossProduct(xAxis, direction);
	double angle = acos(xAxis*direction);

	transform->Translate(translation);
	transform->Scale(scale);
	transform->Rotate3D(-axis, angle);
}


}
