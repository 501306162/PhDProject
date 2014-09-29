#ifndef MATRIX_WRITER_H
#define MATRIX_WRITER_H

#include "MatrixCommon.h"

#include <itkObject.h>
#include <itkObjectFactory.h>


#include <H5Cpp.h>

namespace utils
{
class MatrixWriter : public itk::Object
{
public:
	typedef MatrixWriter 						Self;
	typedef itk::Object 						Superclass;
	typedef itk::SmartPointer<Self> 			Pointer;
	typedef itk::SmartPointer<const Self>		ConstPointer;
	itkTypeMacro(MatrixWriter, itk::Object);
	itkNewMacro(Self);

	
	void SetInput(const MatrixDataSet * data);
	void SetFilename(const std::string &filename);
	void SetAddMatlabHeader(bool val);

	void Write();
	void Update();
protected:
	MatrixWriter();
private:

	void WriteData(const DoubleDataGroup &data);
	void WriteData(const IntDataGroup &data);
	
	void WriteMatrix(const std::string &name, const DoubleMatrixType &mat);
	void WriteMatrix(const std::string &name, const IntMatrixType &mat);


	template<typename Derived>
	H5::DataSpace * CreateDataSpace(const Eigen::MatrixBase<Derived> &mat);

	void AddMatlabHeader();

	const MatrixDataSet * m_Output;
	std::string m_Filename;
	bool m_AddMatlabHeader;


	H5::H5File * m_File;

};


} /* utils */ 

#endif
