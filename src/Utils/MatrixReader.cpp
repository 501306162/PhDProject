#include "MatrixReader.h"

namespace utils
{
// ------------------------------------------------------------------------
MatrixReader::MatrixReader()
{
	m_File = NULL;
	m_Filename = "";
}

// ------------------------------------------------------------------------
void MatrixReader::Update()
{
	Read();
}

// ------------------------------------------------------------------------
void MatrixReader::SetFilename(const std::string &filename)
{
	m_Filename = filename;
}

// ------------------------------------------------------------------------
void MatrixReader::Read()
{
	if(m_Filename.empty())
	{
		itkExceptionMacro(<< "Filename Has Not Been Set");
	}


	// try and open the file
	m_File = new H5::H5File(m_Filename.c_str(), H5F_ACC_RDONLY);

	std::vector<std::string> dsetNames;
	GetDataSetNames(dsetNames);
	if(dsetNames.size() == 0)
	{
		itkExceptionMacro(<< "No valid datasets found");

	}

	// initialise the output
	m_Output = MatrixDataSet::New();

	// loop through the names of the datasets and load the matrixes
	for(unsigned int i = 0; i < dsetNames.size(); i++)
	{
		// get the type of the matrix
		std::string name = dsetNames[i];
		H5::DataSet ds = m_File->openDataSet(name.c_str());	

		switch(ds.getTypeClass())
		{
			case H5T_FLOAT :
				{
				DoubleMatrixType mat;
				ReadMatrix(name, mat);
				m_Output->AddData(name, mat);
				break;
				}
			case H5T_INTEGER :
				{
				IntMatrixType mat;
				ReadMatrix(name, mat);
				m_Output->AddData(name, mat);
				break;
				}
			default :
				break;

		}
	}


	m_File->close();
	delete m_File;

}

// ------------------------------------------------------------------------
MatrixDataSet * MatrixReader::GetOutput() const
{
	return m_Output;
}


// ------------------------------------------------------------------------
void MatrixReader::ReadMatrix(const std::string &name, DoubleMatrixType &mat)
{
	// get the dataset and dataspace
	H5::DataSet ds = m_File->openDataSet(name.c_str());
	H5::DataSpace dspace = ds.getSpace();
	hsize_t dims[2];
	dspace.getSimpleExtentDims(dims);

	std::cout << dims[0] << " " << dims[1] << std::endl;


	// initialise the matrix
	mat = DoubleMatrixType::Zero(dims[0],dims[1]);
	ds.read(mat.data(), H5::PredType::NATIVE_DOUBLE);
}


// ------------------------------------------------------------------------
void MatrixReader::ReadMatrix(const std::string &name, IntMatrixType &mat)
{
	// get the dataset and dataspace
	H5::DataSet ds = m_File->openDataSet(name.c_str());
	H5::DataSpace dspace = ds.getSpace();
	hsize_t dims[2];
	dspace.getSimpleExtentDims(dims);

	// initialise the matrix
	mat = IntMatrixType::Zero(dims[0],dims[1]);
	ds.read(mat.data(), H5::PredType::NATIVE_INT);
}



// ------------------------------------------------------------------------
void MatrixReader::GetDataSetNames(std::vector<std::string> &dsetNames)
{
	unsigned int numberOfObjects = m_File->getNumObjs();
	for(unsigned int i = 0; i < numberOfObjects; i++)
	{	
		if(m_File->getObjTypeByIdx(i) == H5G_DATASET)
			dsetNames.push_back(m_File->getObjnameByIdx(i));
	}
}

} /* utils */ 
