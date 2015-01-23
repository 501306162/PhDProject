#include "LengthData.h"
#include <Directory.h>
#include <QString>
#include <MatrixReader.h>

namespace vt
{
// ------------------------------------------------------------------------
LengthData::Pointer LengthData::Load(const std::string &folderName)
{
	LengthData::Pointer obj = LengthData::New();
	obj->LoadData(folderName);
	return obj;

}

// ------------------------------------------------------------------------
void LengthData::LoadData(const std::string &folderName)
{
	// get the sub folders 
	utils::Directory::FilenamesType filenames = utils::Directory::GetFiles(folderName, ".hdf5");
	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		const std::string filename = filenames[i];
		QString fname = QString::fromStdString(utils::Directory::GetFileName(filename));
		fname = fname.replace(".hdf5","");
		const std::string typeName = fname.toStdString();

		utils::MatrixReader::Pointer reader = utils::MatrixReader::New();
		reader->SetFilename(filename);
		reader->Read();

		utils::MatrixDataSet::Pointer dset = reader->GetOutput();

		MatrixPair pair;
		pair.first = dset->doubleData["lengths"];
		pair.second = dset->intData["owners"];
		
		m_Data[typeName] = pair;
	}

}

// ------------------------------------------------------------------------
void LengthData::Prepare(const unsigned int exclude)
{
	// get the time steps 
	const unsigned int numTimeSteps = m_Data.begin()->second.first.cols();
	MatrixMap::iterator mapIt = m_Data.begin();
	while(mapIt != m_Data.end())
	{
		const std::string type = mapIt->first;
		MembershipListType memberList;

		for(unsigned int i = 0; i < numTimeSteps; i++)
		{
			MatrixType lengths = mapIt->second.first;
			IntMatrixType owners = mapIt->second.second;

			unsigned int rows = lengths.rows();

			SampleType::Pointer sample = SampleType::New();

			for(unsigned int j = 0; j < rows; j++)
			{
				// check the exclude 
				if(owners(j,0) != (int) exclude)
				{
					MeasurementType mv;
					mv[0] = lengths(j,i);
					sample->PushBack(mv);
				}
			}

			typedef itk::Statistics::CovarianceSampleFilter<SampleType> CovFilterType;
			CovFilterType::Pointer covFilter = CovFilterType::New();
			covFilter->SetInput(sample);
			covFilter->Update();

			MembershipType::Pointer member = MembershipType::New();
			member->SetMean(covFilter->GetMean());
			member->SetCovariance(covFilter->GetCovarianceMatrix());

			memberList.push_back(member);
			
			m_Functions[type] = memberList;
		}
	
		++mapIt;
	}


}

// ------------------------------------------------------------------------
double LengthData::Predict(const std::string &type, const unsigned int time, const double length)
{
	MembershipType::Pointer member = m_Functions[type][time];
	MeasurementType mv;
	mv[0] = length;
	return member->Evaluate(mv);
}



}
