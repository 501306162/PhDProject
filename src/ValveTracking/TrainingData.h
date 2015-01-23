#ifndef TRAINING_DATA_H
#define TRAINING_DATA_H

#include <itkObject.h>
#include <itkObjectFactory.h>
#include <MatrixCommon.h>

namespace vt 
{
class ValveTrainingData : public itk::Object 
{
public:
	typedef ValveTrainingData Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(ValveTrainingData, Object);
	itkNewMacro(Self);

	typedef utils::DoubleMatrixType MatrixType;
	typedef utils::IntMatrixType IntMatrixType;

	static Pointer Load(const std::string &filename);

	void LoadData(const std::string &filename);

	itkSetMacro(NumTimeSteps, unsigned int);

	void GetTrainingData(const unsigned int exclude, const unsigned int timeStep, MatrixType &planes);
	void GetTrainingPoints(const unsigned int exclude, const unsigned int timeStep, MatrixType &points);
	void GetTrainingPoints(const unsigned int exclude, const unsigned int timeStep, std::vector<MatrixType> &points);
	void GetTestData(const unsigned int include, const unsigned int timeStep, MatrixType &planes);
	std::vector<int> OwnerList();


protected:
	ValveTrainingData() {}
	virtual ~ValveTrainingData() {}


private:

	unsigned int OwnershipCount(const IntMatrixType &mat, unsigned int ownerId);
	unsigned int TimeStepCount();

	MatrixType m_Planes;
	IntMatrixType m_TimeSteps;
	unsigned int m_NumTimeSteps;
	IntMatrixType m_Owners;
	MatrixType m_Points;
	IntMatrixType m_PointOwners;
	IntMatrixType m_PointTimes;
	

	ValveTrainingData(const Self&);
	void operator=(const Self&);

};


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
	void GetPositiveTrainingData(const unsigned int exclude, MatrixType &X);
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
