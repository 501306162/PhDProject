#ifndef PATCH_TRAINING_DATA_H
#define PATCH_TRAINING_DATA_H

#include <itkObject.h>
#include <itkObjectFactory.h>
#include <MatrixCommon.h>

namespace vt
{
class PatchTrainingData : public itk::Object 
{
public:
	typedef PatchTrainingData Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(PatchTrainingData, Object);
	itkNewMacro(Self);

	typedef utils::DoubleMatrixType MatrixType;
	typedef utils::IntMatrixType IntMatrixType;

	typedef std::pair<MatrixType, IntMatrixType> MatrixOwnerSet;
	typedef std::vector<MatrixOwnerSet> MatrixPointSet;
   	typedef std::vector<MatrixPointSet> MatrixPointSetSequence;
	typedef std::map<std::string, MatrixPointSetSequence> DataMap;	

	static PatchTrainingData::Pointer Load(const std::string &configFilename);

	void LoadData(const std::string &configFilename);

	void GetAllData(const std::string &type, const unsigned int timeStep, 
		const unsigned int pointNum, MatrixType &X, IntMatrixType &y); 

	void GetTrainingData(const unsigned int exclude, 
		const std::string &type, const unsigned int timeStep, 
		const unsigned int pointNum, MatrixType &X, IntMatrixType &y);

	void GetPositiveTrainingData(const unsigned int exclude, 
		const std::string &type, const unsigned int timeStep, 
		const unsigned int pointNum, MatrixType &X);

	void GetNegativeTrainingData(const unsigned int exclude, 
		const std::string &type, const unsigned int timeStep, 
		const unsigned int pointNum, MatrixType &X);


protected:
	PatchTrainingData() {}
	virtual ~PatchTrainingData() {}

private:

	void extractData(const std::string &filename, MatrixPointSet &positives, MatrixPointSet &negatives);
	void extractMatrix(const MatrixOwnerSet &input, const unsigned int exclude,
		MatrixType &X);

	unsigned int ownerCount(const IntMatrixType &owners, const unsigned int search);

	PatchTrainingData(const Self&);
	void operator=(const Self&);

	DataMap m_PosData;
	DataMap m_NegData;
};



}

#endif
