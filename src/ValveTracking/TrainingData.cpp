#include "TrainingData.h"

#include <MatrixReader.h>
#include <itkPoint.h>

namespace vt 
{
// ------------------------------------------------------------------------
ValveTrainingData::Pointer ValveTrainingData::Load(const std::string &filename)
{
	ValveTrainingData::Pointer data = ValveTrainingData::New();
	data->LoadData(filename);
	data->SetNumTimeSteps(25);
	return data;
}


// ------------------------------------------------------------------------
void ValveTrainingData::GetTrainingPoints(const unsigned int exclude, const unsigned int timeStep, MatrixType &points)
{
	// get the output size of the matrixes	
	unsigned int outputCount = m_Points.rows() - OwnershipCount(m_PointOwners, exclude);
	unsigned int timed = TimeStepCount();
	unsigned int size = (int) ((double) outputCount / (double) timed);


	unsigned int featureLength = m_Points.cols();
	points = MatrixType(size, featureLength);

	unsigned int rowCount = 0;
	for(unsigned int i = 0; i < m_PointOwners.rows(); i++)
	{
		if(m_PointOwners(i,0) != (int) exclude && m_PointTimes(i,0) == (int) timeStep)
		{
			points.row(rowCount) = m_Points.row(i);
			rowCount++;
		}
	}

}

// ------------------------------------------------------------------------
void ValveTrainingData::GetTrainingPoints(const unsigned int exclude, const unsigned int timeStep, 
		std::vector<MatrixType> &points)
{
	std::vector<int> ownerList = OwnerList();
	
	typedef itk::Point<double, 3> PointType;
	typedef std::vector<MatrixType> PointList;
	typedef std::map<int, PointList> PointMap;

	PointMap tmpPoints;
	for(unsigned int i = 0; i < m_Points.rows(); i++)
	{
		if(m_PointOwners(i,0) != (int) exclude && m_PointTimes(i,0) == (int) timeStep)
		{
			int owner = m_PointOwners(i,0);
			tmpPoints[owner].push_back(m_Points.row(i));
		}
	}

	PointMap::iterator mapIt = tmpPoints.begin();
	while(mapIt != tmpPoints.end())
	{
		MatrixType pmat = MatrixType(3, mapIt->second.size());
		for(unsigned int i = 0; i < mapIt->second.size(); i++)
		{
			pmat.col(i) = mapIt->second[i].row(0).transpose();			
		}
		
		points.push_back(pmat);

		++mapIt;
	}


}

// ------------------------------------------------------------------------
void ValveTrainingData::GetTestData(const unsigned int include, const unsigned int timeStep, MatrixType &planes)
{
	unsigned int owned = OwnershipCount(m_Owners, include);
	unsigned int timed = TimeStepCount();

	unsigned int size = (int) ((double) owned / (double) timed);

	planes = MatrixType(size, m_Planes.cols());

	unsigned int count = 0;
	for(unsigned int i = 0; i < m_Owners.rows(); i++)
	{
		if(m_Owners(i,0) == (int) include && m_TimeSteps(i,0) == (int) timeStep)
		{
			planes.row(count) = m_Planes.row(i);
			count++;
		}
	}
}



// ------------------------------------------------------------------------
void ValveTrainingData::GetTrainingData(const unsigned int exclude, const unsigned int timeStep, MatrixType &planes)
{
	// get the output size of the matrixes	
	unsigned int outputCount = m_Planes.rows() - OwnershipCount(m_Owners, exclude);
	unsigned int timed = TimeStepCount();
	unsigned int size = (int) ((double) outputCount / (double) timed);


	unsigned int featureLength = m_Planes.cols();

	planes = MatrixType(size, featureLength);

	unsigned int rowCount = 0;
	for(unsigned int i = 0; i < m_Owners.rows(); i++)
	{
		if(m_Owners(i,0) != (int) exclude && m_TimeSteps(i,0) == (int) timeStep)
		{
			planes.row(rowCount) = m_Planes.row(i);
			rowCount++;
		}
	}
}


// ------------------------------------------------------------------------
std::vector<int> ValveTrainingData::OwnerList()
{
	std::vector<int> ownerList;
	std::map<int,int> checkList;	
	for(unsigned int i = 0; i < m_Owners.rows(); i++)
	{
		int owner = m_Owners(i,0);		
		if(checkList.count(owner) == 0)
		{
			checkList[owner] = owner;
			ownerList.push_back(owner);
		}
	}
	return ownerList;
}

// ------------------------------------------------------------------------
void ValveTrainingData::LoadData(const std::string &filename)
{
	utils::MatrixReader::Pointer reader = utils::MatrixReader::New();
	reader->SetFilename(filename);
	reader->Read();

	utils::MatrixDataSet::Pointer dataSet = reader->GetOutput();
	m_Planes = dataSet->doubleData["planes"];
	m_TimeSteps = dataSet->intData["time_steps"];
	m_Owners = dataSet->intData["owners"];
	m_Points = dataSet->doubleData["points"];
	m_PointOwners = dataSet->intData["point_owners"];
	m_PointTimes = dataSet->intData["point_times"];
}


// ------------------------------------------------------------------------
unsigned int ValveTrainingData::OwnershipCount(const IntMatrixType &mat, unsigned int ownerId)
{
	unsigned int rows = mat.rows();
	unsigned int count = 0;
	for(unsigned int i = 0; i < rows; i++)
	{
		if(mat(i,0) == (int) ownerId)
			count++;
	}

	return count;
}

// ------------------------------------------------------------------------
unsigned int ValveTrainingData::TimeStepCount()
{
	return m_NumTimeSteps;
}


//////////////////////////////////////////////////////////////////////////////




// ------------------------------------------------------------------------
TrainingData::Pointer TrainingData::Load(const std::string &filename)
{
	TrainingData::Pointer data = TrainingData::New();
	data->LoadData(filename);
	return data;
}


// ------------------------------------------------------------------------
void TrainingData::GetTestData(const unsigned int include, MatrixType &X, IntMatrixType &y)
{
	unsigned int size = OwnershipCount(include);
	X = MatrixType(size, m_X.cols());
	y = IntMatrixType(size,1);
	

	unsigned int count = 0;
	for(unsigned int i = 0; i < m_Owners.rows(); i++)
	{
		if(m_Owners(i,0) == (int) include)
		{
			X.row(count) = m_X.row(i);
			y.row(count) = m_y.row(i);
			count++;
		}
	}
}


// ------------------------------------------------------------------------
std::vector<int> TrainingData::OwnerList()
{
	std::vector<int> ownerList;
	std::map<int,int> checkList;	
	for(unsigned int i = 0; i < m_Owners.rows(); i++)
	{
		int owner = m_Owners(i,0);		
		if(checkList.count(owner) == 0)
		{
			checkList[owner] = owner;
			ownerList.push_back(owner);
		}

	}
	return ownerList;
}

// ------------------------------------------------------------------------
void TrainingData::LoadData(const std::string &filename)
{
	utils::MatrixReader::Pointer reader = utils::MatrixReader::New();
	reader->SetFilename(filename);
	reader->Read();

	utils::MatrixDataSet::Pointer dataSet = reader->GetOutput();
	m_X = dataSet->doubleData["X"];
	m_y = dataSet->intData["labels"];
	m_Owners = dataSet->intData["owners"];
}

// ------------------------------------------------------------------------
void TrainingData::GetTrainingData(const unsigned int exclude, MatrixType &X, IntMatrixType &y)
{
	// get the output size of the matrixes	
	unsigned int outputCount = m_X.rows() - OwnershipCount(exclude);
	unsigned int featureLength = m_X.cols();
	X = MatrixType(outputCount, featureLength);
	y = IntMatrixType(outputCount, 1);

	unsigned int rowCount = 0;
	for(unsigned int i = 0; i < m_Owners.rows(); i++)
	{
		if(m_Owners(i,0) != (int) exclude)
		{
			X.row(rowCount) = m_X.row(i);
			y.row(rowCount) = m_y.row(i);
			rowCount++;
		}
	}
}


// ------------------------------------------------------------------------
void TrainingData::GetPositiveTrainingData(const unsigned int exclude, MatrixType &X)
{
	unsigned int count = 0;
	for(unsigned int i = 0; i < m_y.rows(); i++)
	{
		if(m_y(i,0) == 1 && m_Owners(i,0) != (int) exclude)
			count++;		
	}

	X = MatrixType(count, m_X.cols());

	count = 0;
	for(unsigned int i = 0; i < m_y.rows(); i++)
	{
		if(m_y(i,0) == 1 && m_Owners(i,0) != (int) exclude)
		{
			MatrixType r = m_X.row(i).normalized();
			X.row(count) = r;
			count++;
		}
	}
}



// ------------------------------------------------------------------------
unsigned int TrainingData::OwnershipCount(unsigned int ownerId)
{
	unsigned int rows = m_Owners.rows();
	unsigned int count = 0;
	for(unsigned int i = 0; i < rows; i++)
	{
		if(m_Owners(i,0) == (int) ownerId)
			count++;
	}

	return count;
}

}
