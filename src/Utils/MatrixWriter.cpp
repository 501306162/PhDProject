#include "MatrixWriter.h"

#include <cstdio>
#include <H5Cpp.h>

namespace utils
{
// ------------------------------------------------------------------------
MatrixWriter::MatrixWriter()
{
	m_Filename = "";
	m_Output = 0;
	m_File = NULL;
	m_AddMatlabHeader = true;
}

// ------------------------------------------------------------------------
void MatrixWriter::SetInput(const MatrixDataSet * data)
{
	m_Output = data;
}

// ------------------------------------------------------------------------
void MatrixWriter::SetFilename(const std::string &filename)
{
	m_Filename = filename;
}

// ------------------------------------------------------------------------
void MatrixWriter::SetAddMatlabHeader(bool val)
{
	m_AddMatlabHeader = val;
}


// ------------------------------------------------------------------------
void MatrixWriter::Update()
{
	Write();
}

// ------------------------------------------------------------------------
void MatrixWriter::Write()
{
	// do some initial checks
	if(m_Filename.empty())
	{
		itkExceptionMacro(<< "Filename Has Not Been Set");
	}

	if(!m_Output)
	{
		itkExceptionMacro(<< "Output Has Not Been Set");
	}


	// create the hdf5 file object
	m_File = new H5::H5File(m_Filename.c_str(), H5F_ACC_TRUNC);

	// write the data to the output
	WriteData(m_Output->doubleData);
	WriteData(m_Output->intData);
	

	// close and delete the hdf5 file
	m_File->close();
	delete m_File;


	if(m_AddMatlabHeader)
		AddMatlabHeader();
}

// ------------------------------------------------------------------------
void MatrixWriter::AddMatlabHeader()
{

	// Prepend the header.
	FILE* fd = fopen(m_Filename.c_str(), "r+b");
	if (!fd)
	{
		itkExceptionMacro(<< "Couldn't Open Filename FOr Header Appending");
	}

	char header[512];
	memset(header, 0, sizeof(header));
	sprintf(header, "MATLAB 7.3 format file written by MatFile class, by Tim Hutt"); // Do not make this longer than 116 chars (I think)
	header[124] = 0;
	header[125] = 2;
	header[126] = 'I';
	header[127] = 'M';

	// Get file length.
	fseek(fd, 0, SEEK_END);
	long long length = ftell(fd);
	fseek(fd, 0, SEEK_SET);

	// TODO: Do this properly without reading entire file into memory.
	if (length > 1024L*1024L*1024L*10L) // 10 GB.
	{
		itkExceptionMacro(<< "File To Big To Write Header");
		fclose(fd);
	}

	unsigned char* buffer = new unsigned char[length];
	if (fread(buffer, 1, length, fd) != length)
	{
		delete[] buffer;
		fclose(fd);
	}

	fseek(fd, 0, SEEK_SET);


	fwrite(header, 1, sizeof(header), fd);
	fwrite(buffer, 1, length, fd);
	fclose(fd);

	delete[] buffer;

}

// ------------------------------------------------------------------------
void MatrixWriter::WriteData(const IntDataGroup &data)
{
	IntDataGroup::const_iterator it = data.begin();
	while(it != data.end())
	{
		std::string name = it->first;
		IntMatrixType mat = it->second;		
		WriteMatrix(name, mat);

		++it;
	}	
}


// ------------------------------------------------------------------------
void MatrixWriter::WriteMatrix(const std::string &name, const IntMatrixType &mat)
{
	// create the dataspace and data set
	H5::DataSpace * ds = CreateDataSpace(mat);
	H5::DataSet * dset = new H5::DataSet(m_File->createDataSet(
				name.c_str(), H5::PredType::NATIVE_INT, *ds));

	// write the data
	dset->write(mat.data(),H5::PredType::NATIVE_INT);
	dset->close();
	delete ds;
	delete dset;
}


// ------------------------------------------------------------------------
void MatrixWriter::WriteData(const DoubleDataGroup &data)
{
	DoubleDataGroup::const_iterator it = data.begin();
	while(it != data.end())
	{
		std::string name = it->first;
		DoubleMatrixType mat = it->second;		
		WriteMatrix(name, mat);

		++it;
	}	
}

// ------------------------------------------------------------------------
void MatrixWriter::WriteMatrix(const std::string &name, const DoubleMatrixType &mat)
{
	// create the dataspace and data set
	H5::DataSpace * ds = CreateDataSpace(mat);
	H5::DataSet * dset = new H5::DataSet(m_File->createDataSet(
				name.c_str(), H5::PredType::NATIVE_DOUBLE, *ds));

	// write the data
	dset->write(mat.data(),H5::PredType::NATIVE_DOUBLE);
	delete ds;
	delete dset;
}

// ------------------------------------------------------------------------
template<typename Derived>
H5::DataSpace * MatrixWriter::CreateDataSpace(const Eigen::MatrixBase<Derived> &mat)
{
	int rank = 2;
	hsize_t	dims[rank];
	dims[0] = mat.rows();
	dims[1] = mat.cols();	

	return new H5::DataSpace(rank, dims);
}

} /* 	utils */ 


