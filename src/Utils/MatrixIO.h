#ifndef MATRIX_IO_H
#define MATRIX_IO_H

#include <itkObject.h>
#include <itkObjectFactory.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <Eigen/Dense>

#include <map>

namespace utils
{

class MatrixIO : public itk::Object
{
public:
	typedef MatrixIO						Self;
	typedef itk::Object 					Superclass;
	typedef itk::SmartPointer<Self>	    	Pointer;
	typedef itk::SmartPointer<const Self> 	ConstPointer;
	itkNewMacro(Self);
	itkTypeMacro(MatrixIO, itk::Object);


	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DoubleMatrixType;
	typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> IntMatrixType;

	typedef struct data_set_
	{
		std::map<std::string, DoubleMatrixType> doubleData;
		std::map<std::string, IntMatrixType> intData;
		int count;
	} DataSet;



	void SetFilename(const std::string &filename);
	void AddInput(const std::string &name, const DoubleMatrixType & mat);
	void AddInput(const std::string &name, const IntMatrixType & mat);
	void UseMatHeaders(bool val);

	void Write() throw(itk::ExceptionObject);
	void Read() throw(itk::ExceptionObject);

	virtual ~MatrixIO();

protected:
	MatrixIO();

private:
		
	void AddMatlabHeader();

	void Open(const std::string &fn, bool create);
	void Close();

	void AddDataSet(const std::string &name, int rank, const hsize_t * dims, hid_t type, const void * buffer);


	bool m_UseMatHeaders;
	DataSet m_DataSet;
	std::string m_Filename;


	hid_t m_File;
};

} /* utils */ 

#endif
