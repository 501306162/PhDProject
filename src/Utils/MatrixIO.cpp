#include "MatrixIO.h"

#include <cstdio>

#include <H5Cpp.h>

namespace utils
{

// ------------------------------------------------------------------------
MatrixIO::MatrixIO()
{
	m_UseMatHeaders = true;
	m_Filename = "";
}

// ------------------------------------------------------------------------
void MatrixIO::SetFilename(const std::string &filename)
{
	m_Filename = filename;
}

// ------------------------------------------------------------------------
void MatrixIO::AddInput(const std::string &name, const DoubleMatrixType &input)
{
	m_DataSet.doubleData[name] = input;
	m_DataSet.count++;
}


// ------------------------------------------------------------------------
void MatrixIO::AddInput(const std::string &name, const IntMatrixType &input)
{
	m_DataSet.intData[name] = input;
	m_DataSet.count++;
}

// ------------------------------------------------------------------------
void MatrixIO::Read() throw(itk::ExceptionObject)
{
	// do some initial checks
	if(m_Filename.empty())
	{
		itkExceptionMacro(<< "Filename Has Not Been Set");
	}	


	H5::H5File * file = new H5::H5File(m_Filename, H5F_ACC_RDONLY);
	unsigned int numObjects =  file->getNumObjs();

	std::vector<std::string> datasetNames;

	for(unsigned int i = 0; i < numObjects; i++)
	{
		// check if the object is a dataset 
		if(file->getObjTypeByIdx(i) == H5G_DATASET)
		{
			// do a switch on the data type
			H5::DataSet dset = file->openDataSet(file->getObjnameByIdx(i));

			// get the rank and dimensions
			H5::DataSpace dspace = dset.getSpace();
			hsize_t dims[2];
			int rank = dspace.getSimpleExtentNdims();
			int ndims = dspace.getSimpleExtentDims(dims,NULL);

			std::cout << rank << " " << dims[0] << " " << dims[1] << std::endl;


			if(dset.getTypeClass() == H5T_FLOAT)
			{
				// create the output array
				double * data = new double[dims[0]*dims[1]];
				dset.read(data, H5::PredType::NATIVE_DOUBLE);

				for(unsigned int i = 0; i < dims[0]; i++)
				{
					for(unsigned int j = 0; j < dims[1]; j++)
					{
						std::cout << data[i*dims[1]+j] << " ";
					}
					std::cout << std::endl;
				}
											
			}
			else if(dset.getTypeClass() == H5T_INTEGER)
			{

			}
			
		}

				
	}
	
	file->close();

}


// ------------------------------------------------------------------------
void MatrixIO::Write() throw (itk::ExceptionObject)
{
	// do some initial checks
	if(m_Filename.empty())
	{
		itkExceptionMacro(<< "Filename Has Not Been Set");
	}	

	if(m_DataSet.count == 0)
	{
		itkExceptionMacro(<< "Output Is Either Not set or is uninitialised");
	}


	// open the file
	Open(m_Filename, true);
	if(m_File < 0)
	{
		itkExceptionMacro(<< "Error Opening the output File");
	}

	// start writing the double buffers
	std::map<std::string, DoubleMatrixType>::iterator dIt = m_DataSet.doubleData.begin();
	while(dIt != m_DataSet.doubleData.end())
	{
		DoubleMatrixType mat = dIt->second;
		hsize_t dims[2];
		dims[0] = mat.rows();
		dims[1] = mat.cols();

		AddDataSet(dIt->first, 2, dims, H5T_NATIVE_DOUBLE, (void*) mat.data());
		++dIt;
	}



	// start writing the int buffers
	std::map<std::string, IntMatrixType>::iterator iIt = m_DataSet.intData.begin();
	while(dIt != m_DataSet.doubleData.end())
	{
		IntMatrixType mat = iIt->second;
		hsize_t dims[2];
		dims[0] = mat.rows();
		dims[1] = mat.cols();

		AddDataSet(iIt->first, 2, dims, H5T_NATIVE_INT, (void*) mat.data());
		++iIt;
	}


	Close();
	
	// set the header if we need it
	if(m_UseMatHeaders)
	{
		AddMatlabHeader();
	}
}


// ------------------------------------------------------------------------
void MatrixIO::AddMatlabHeader()
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
void MatrixIO::AddDataSet(const std::string &name, int rank, 
		const hsize_t * dims, hid_t type, const void * buffer)
{
	H5LTmake_dataset(m_File, name.c_str(), rank, dims, type, buffer);
}



// ------------------------------------------------------------------------
void MatrixIO::Open(const std::string &fn, bool create)
{
	if(create)
	{
		m_File = H5Fcreate(fn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);	
	}
	else
	{
		m_File = H5Fopen(fn.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	}
}


// ------------------------------------------------------------------------
void MatrixIO::Close()
{
	if(m_File >= 0)
	{
		H5Fclose(m_File);
		m_File = -1;
	}
}

// ------------------------------------------------------------------------
MatrixIO::~MatrixIO()
{
	Close();
}




} /* utils */ 
