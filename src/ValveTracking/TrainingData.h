#ifndef TRAINING_DATA_H
#define TRAINING_DATA_H

#include <itkObject.h>
#include <itkObjectFactory.h>
#include <MatrixCommon.h>

namespace vt 
{
class TrainingData : public itk::Object 
{
public:
	typedef TrainingData Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(TrainingData, Object);
	itkNewMacro(Self);

	typedef utils::DoubleMatrixType MatrixType;
	typedef utils::IntMatrixType IntMatrixType;

	static Pointer Load(const std::string &filename);

	void LoadData(const std::string &filename);

	void GetTrainingData(const unsigned int exclude, MatrixType &X, IntMatrixType &y);
	void GetTestData(const unsigned int include, MatrixType &X, IntMatrixType &y);
	std::vector<int> OwnerList();


protected:
	TrainingData() {}
	virtual ~TrainingData() {}


private:

	unsigned int OwnershipCount(unsigned int ownerId);

	MatrixType m_X;
	IntMatrixType m_y;
	IntMatrixType m_Owners;
	

	TrainingData(const Self&);
	void operator=(const Self&);
};

}

#endif
