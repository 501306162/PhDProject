#ifndef SVM_CLASSIFIER_H
#define SVM_CLASSIFIER_H

#include "TrainingData.h"
#include <svm.h>


namespace vt
{
class SVMClassifier : public itk::Object 
{
public:
	typedef SVMClassifier Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(SVMClassifier, Object);
	itkNewMacro(Self);

	typedef TrainingData::MatrixType MatrixType;
	typedef TrainingData::IntMatrixType IntMatrixType;

	typedef struct svm_model ModelType;
	typedef struct svm_problem ProblemType;
	typedef struct svm_node NodeType;
	typedef struct svm_parameter ParametersType;

	void Train(const MatrixType &X, const IntMatrixType &y);
	void PredictProbability(const MatrixType &X, IntMatrixType &classes, MatrixType &probs);


protected:
	SVMClassifier() : m_Problem(NULL) {}
	virtual ~SVMClassifier() {}

	

private:
	
	void BuildProblem(const MatrixType &X, const IntMatrixType &y, ProblemType * problem);
	void BuildNode(const MatrixType &row, NodeType * node);
	double Weight(const IntMatrixType &y, unsigned int label);
	unsigned int NonZeroNodes(const MatrixType &row);

	SVMClassifier(const Self&);
	void operator=(const Self&);

	ModelType * m_Model;
	ProblemType * m_Problem;
	

};


}

#endif
