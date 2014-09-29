#ifndef MATRIX_COMMON_H
#define MATRIX_COMMON_H

#include <Eigen/Dense>
#include <map>

#include <itkObject.h>
#include <itkObjectFactory.h>

namespace utils
{
/** typedefs for the types of matrixes */
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> DoubleMatrixType;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> IntMatrixType;

/** typedefs for the datasets */
typedef std::pair<std::string, DoubleMatrixType> DoubleDataSet;
typedef std::map<std::string, DoubleMatrixType> DoubleDataGroup;
typedef std::pair<std::string, IntMatrixType> IntDataSet;
typedef std::map<std::string, IntMatrixType> IntDataGroup;



class MatrixDataSet : public itk::Object
{
public :
	typedef MatrixDataSet 				Self;
	typedef itk::Object					Superclass;
	typedef itk::SmartPointer<Self>		Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkNewMacro(Self);
	itkTypeMacro(MatrixDataSet, itk::Object);


	MatrixDataSet()
	{
		this->numDataSets = 0;
	}
	
	void AddData(const std::string &name, const IntMatrixType &mat)
	{
		intData[name] = mat;
		numDataSets++;
	}

	void AddData(const std::string &name, const DoubleMatrixType &mat)
	{
		doubleData[name] = mat;
		numDataSets++;
	}

	DoubleDataGroup doubleData;
	IntDataGroup intData;
	int numDataSets;
};

} /* utils */ 

#endif
