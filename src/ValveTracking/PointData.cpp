#include "PointData.h"
#include <Directory.h>
#include <ValveIO.h>
#include <QString>


namespace vt
{
// ------------------------------------------------------------------------
PointData::Pointer PointData::Load(const std::string &folder, const unsigned int exclude)
{
	PointData::Pointer obj = PointData::New();
	obj->LoadData(folder, exclude);
	return obj;
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
