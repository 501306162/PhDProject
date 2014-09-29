#ifndef MATRIX_READER_H
#define MATRIX_READER_H

#include "MatrixCommon.h"

#include <itkObject.h>
#include <itkObjectFactory.h>


#include <H5Cpp.h>

namespace utils
{
class MatrixReader : public itk::Object
{
public:
	typedef MatrixReader 						Self;
	typedef itk::Object 						Superclass;
	typedef itk::SmartPointer<Self> 			Pointer;
	typedef itk::SmartPointer<const Self>		ConstPointer;
	itkTypeMacro(MatrixReader, itk::Object);
	itkNewMacro(Self);

	
	void SetInput(const MatrixDataSet * data);
	void SetFilename(const std::string &filename);

	MatrixDataSet * GetOutput() const;

	void Read();
	void Update();
protected:
	MatrixReader();
private:

	void ReadData(const DoubleDataGroup &data);
	void ReadData(const IntDataGroup &data);
	
	void ReadMatrix(const std::string &name, DoubleMatrixType &mat);
	void ReadMatrix(const std::string &name, IntMatrixType &mat);

	void GetDataSetNames(std::vector<std::string> &dsetNames);


	MatrixDataSet::Pointer m_Output;
	std::string m_Filename;


	H5::H5File * m_File;

};


} /* utils */ 

#endif
