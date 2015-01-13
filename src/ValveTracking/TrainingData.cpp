#include "TrainingData.h"

#include <MatrixReader.h>

namespace vt 
{
// ------------------------------------------------------------------------
ValveTrainingData::Pointer ValveTrainingData::Load(const std::string &filename)
{
	ValveTrainingData::Pointer data = ValveTrainingData::New();
	data->LoadData(filename);
	return data;
}


// ------------------------------------------------------------------------
void ValveTrainingData::GetTestData(const unsigned int include, const unsigned int timeStep, MatrixType &planes)
{
	unsigned int owned = OwnershipCount(include);
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
	unsigned int outputCount = m_Planes.rows() - OwnershipCount(exclude);
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
}


// ------------------------------------------------------------------------
unsigned int ValveTrainingData::OwnershipCount(unsigned int ownerId)
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
