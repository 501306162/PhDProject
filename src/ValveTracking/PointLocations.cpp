#include "PointLocations.h"

#include <MatrixReader.h>

namespace vt
{
// ------------------------------------------------------------------------
PointLocations::Pointer PointLocations::Load(const std::string &filename)
{
	PointLocations::Pointer obj = PointLocations::New();
	obj->LoadData(filename);
	return obj;
}

// ------------------------------------------------------------------------
void PointLocations::LoadData(const std::string &filename)
{
	utils::MatrixReader::Pointer reader = utils::MatrixReader::New();
	reader->SetFilename(filename);
	reader->Read();

	m_Data = reader->GetOutput()->doubleData;
}



// ------------------------------------------------------------------------
void  PointLocations::GetPoints(const int exclude, const std::string type, 
		unsigned int timeStep, MatrixType &p1, MatrixType &p2)
{
	MatrixType mat = m_Data[type];
	unsigned int count = 0;
	for(unsigned int i = 0; i < mat.rows(); i++)
	{
		if((int) mat(i,0) != exclude && (int) mat(i,1) == (int) timeStep)
			count++;
	}


	p1 = MatrixType(count, 3);
	p2 = MatrixType(count, 3);

	count = 0;
	for(unsigned int i = 0; i < mat.rows(); i++)
	{
		if((int) mat(i,0) != exclude && (int) mat(i,1) == (int) timeStep)
		{
			for(unsigned int j = 0; j < 3; j++)
			{
				p1(count, j) = mat(i, 2+j);
				p2(count, j) = mat(i, 5+j);
			}
			count++;
		}
	}
}


}
