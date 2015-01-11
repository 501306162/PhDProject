#include "TrainingData.h"

#include <MatrixReader.h>

namespace vt 
{
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
