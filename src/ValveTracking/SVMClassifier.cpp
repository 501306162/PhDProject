#include "SVMClassifier.h"

namespace vt 
{

void print_null(const char *s) {}

// ------------------------------------------------------------------------
void SVMClassifier::Train(const MatrixType &X, const IntMatrixType &y)
{
	svm_set_print_string_function(&print_null);


	if(m_Problem) delete m_Problem;

	//std::cout << "Building Problem" << std::endl;
	m_Problem = new ProblemType;
	BuildProblem(X, y, m_Problem);

	//std::cout << m_Problem->l << std::endl;

	// build the parameters 
	ParametersType * params = new ParametersType;
	params->svm_type = C_SVC;
	params->kernel_type = LINEAR;
	params->cache_size = 1000.00;
	params->C = 1;
	params->gamma = 1.0 / (double) X.cols();
	params->coef0 = 0;
	params->degree = 3;
	params->nu = 0.5;
	params->eps = 0.1;
	params->shrinking = 1;
	params->probability = 1;
	
	// compute the weights 
	params->nr_weight = 2;
	params->weight_label = new int[2];
	params->weight_label[0] = 0;
	params->weight_label[1] = 1;
	params->weight = new double[2];
	params->weight[0] = Weight(y,0);
	params->weight[1] = Weight(y,1);

	
	//std::cout << "Checking Parameters" << std::endl;
	//std::cout << svm_check_parameter(m_Problem, params) << std::endl;


	//std::cout << "Training..." << std::endl;
	m_Model = svm_train(m_Problem, params);

}

// ------------------------------------------------------------------------
void SVMClassifier::PredictProbability(const MatrixType &X, IntMatrixType &classes, MatrixType &probs)
{
	//std::cout << "Testing Values" << std::endl;
	// allocate the output 
	classes = IntMatrixType(X.rows(),1);
	probs = MatrixType(X.rows(),2);

	//std::cout << m_Model->l << std::endl;

	// create the nodes one by one
	for(unsigned int i = 0; i < X.rows(); i++)
	{
		unsigned int nonZeros = NonZeroNodes(X.row(i));
		NodeType * nodes = new NodeType[nonZeros+1];
		BuildNode(X.row(i), nodes);

		double *pr = new double[2];
		double cls = svm_predict_probability(m_Model,nodes,pr);

		classes(i,0) = (int) cls;
		probs(i,0) = pr[0];
		probs(i,1) = pr[1];
	}
}

// ------------------------------------------------------------------------
bool SVMClassifier::Save(const std::string &filename)
{
	if(svm_save_model(filename.c_str(), m_Model) == 0) return true;
	else return false;
}


// ------------------------------------------------------------------------
bool SVMClassifier::Load(const std::string &filename)
{
	m_Model = svm_load_model(filename.c_str());
	if(m_Model == NULL) return false;
	else return true;
}

// ------------------------------------------------------------------------
unsigned int SVMClassifier::NonZeroNodes(const MatrixType &row)
{
	// count the non zero values
	unsigned int nonZeros = 0;
	for(unsigned int j = 0; j < row.cols(); j++)
	{
		if(row(0,j) > 0.0)
			nonZeros++;
	}

	return nonZeros;
}

// ------------------------------------------------------------------------
double SVMClassifier::Weight(const IntMatrixType &y, unsigned int label)
{
	unsigned int count = 0;
	for(unsigned int i = 0; i < y.rows(); i++)
	{
		if(y(i,0) == (const int) label)
			count++;
	}
	return ((double) count) / ((double) y.rows());
}


// ------------------------------------------------------------------------
void SVMClassifier::BuildProblem(const MatrixType &X, const IntMatrixType &y, ProblemType * problem)
{
	problem->l = X.rows();
	problem->y = new double[problem->l];
	problem->x = new NodeType*[problem->l];

	for(int i = 0; i < problem->l; i++)
	{
		problem->y[i] = (double) y(i,0);		

		unsigned int nonZeros = NonZeroNodes(X.row(i));
		NodeType * nodes = new NodeType[nonZeros+1];
		BuildNode(X.row(i), nodes);
		problem->x[i] = nodes;
		
	}

}

// ------------------------------------------------------------------------
void SVMClassifier::BuildNode(const MatrixType &row, NodeType * nodes)
{
	// create the list of nodes
	unsigned int count = 0;
	for(unsigned int i = 0; i < row.cols(); i++)
	{
		if(row(0,i) > 0.0)
		{
			NodeType *n = new NodeType;
			n->index = i;
			n->value = row(0,i);
			nodes[count] = *n;
			count++;
		}
	}

	NodeType * endNode = new NodeType;
	endNode->index = -1;
	nodes[count] = *endNode;
}

}
