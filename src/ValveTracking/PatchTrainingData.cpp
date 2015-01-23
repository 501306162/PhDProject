#include "PatchTrainingData.h"
#include "ParameterHelpers.h"
#include <ConfigParser.h>
#include <MatrixReader.h>


namespace vt
{
// ------------------------------------------------------------------------
PatchTrainingData::Pointer PatchTrainingData::Load(const std::string &configFilename)
{
	PatchTrainingData::Pointer obj = PatchTrainingData::New();
	obj->LoadData(configFilename);
	return obj;
}

// ------------------------------------------------------------------------
void PatchTrainingData::LoadData(const std::string &configFilename)
{
	PatchExtractorParameters params(configFilename);

	// loop through the sub types and load the data
	for(unsigned int i = 0; i < params.valveSubTypes.size(); i++)
	{
		const std::string type = params.valveSubTypes[i];
		MatrixPointSetSequence posSequence, negSequence;

		std::cout << "Loading: " << type << std::endl;
		for(unsigned int t = 0; t < params.numberOfTimeSteps; t++)
		{
			const std::string fileName = params.getOutputFilename(t, type);
			MatrixPointSet pos, neg;
			extractData(fileName, pos, neg);
			std::cout << neg.front().first.rows() << std::endl;	

			posSequence.push_back(pos);
			negSequence.push_back(neg);
		}		

		m_PosData[type] = posSequence;
		m_NegData[type] = negSequence;
	}
}


// ------------------------------------------------------------------------
void PatchTrainingData::GetAllData(const std::string &type, const unsigned int timeStep, 
		const unsigned int pointNum, MatrixType &X, IntMatrixType &y)
{
	GetTrainingData(1000, type, timeStep, pointNum, X, y);
}

// ------------------------------------------------------------------------
void PatchTrainingData::extractData(const std::string &filename, MatrixPointSet &positives, MatrixPointSet &negatives)
{
	// load the data from the disk
	utils::MatrixReader::Pointer reader = utils::MatrixReader::New();
	reader->SetFilename(filename);
	reader->Read();
	utils::MatrixDataSet::Pointer dset = reader->GetOutput();

	

	MatrixOwnerSet p1Pos, p2Pos, p1Neg, p2Neg;
	p1Pos.first = dset->doubleData["p1-pos"];
	p1Pos.second = dset->intData["p1-pos-o"];

	p2Pos.first = dset->doubleData["p2-pos"];
	p2Pos.second = dset->intData["p2-pos-o"];

	positives.push_back(p1Pos);
	positives.push_back(p2Pos);

	p1Neg.first = dset->doubleData["p1-neg"];
	p1Neg.second = dset->intData["p1-neg-o"];

	p2Neg.first = dset->doubleData["p2-neg"];
	p2Neg.second = dset->intData["p2-neg-o"];

	negatives.push_back(p1Neg);
	negatives.push_back(p2Neg);
}

// ------------------------------------------------------------------------
void PatchTrainingData::GetTrainingData(const unsigned int exclude, 
		const std::string &type, const unsigned int timeStep, 
		const unsigned int pointNum, MatrixType &X, IntMatrixType &y)
{
	MatrixType posX, negX;
	GetPositiveTrainingData(exclude, type, timeStep, pointNum, posX);
	GetNegativeTrainingData(exclude, type, timeStep, pointNum, negX);

	const unsigned int totalOutput = posX.rows() + negX.rows();
	const unsigned int cols = posX.cols();

	X = MatrixType(totalOutput, cols);
	y = IntMatrixType(totalOutput, 1);

	X.block(0,0,negX.rows(), cols) = negX;
	X.block(negX.rows(),0,posX.rows(), cols) = posX;
	y.block(0,0,negX.rows(),1) = IntMatrixType::Zero(negX.rows(),1);
	y.block(negX.rows(),0,posX.rows(),1) = IntMatrixType::Ones(posX.rows(),1);
}


// ------------------------------------------------------------------------
void PatchTrainingData::GetNegativeTrainingData(const unsigned int exclude, 
		const std::string &type, const unsigned int timeStep, 
		const unsigned int pointNum, MatrixType &X)
{
	MatrixOwnerSet mat = m_NegData[type][timeStep][pointNum-1];
	extractMatrix(mat, exclude, X);
}


// ------------------------------------------------------------------------
void PatchTrainingData::GetPositiveTrainingData(const unsigned int exclude, 
		const std::string &type, const unsigned int timeStep, 
		const unsigned int pointNum, MatrixType &X)
{
	MatrixOwnerSet mat = m_PosData[type][timeStep][pointNum-1];
	extractMatrix(mat, exclude, X);
}


// ------------------------------------------------------------------------
unsigned int PatchTrainingData::ownerCount(const IntMatrixType &owners, const unsigned int search)
{
	for(unsigned int i = 0; i < owners.rows(); i++)
	{
		if(owners(i,0) == (int) search)
			return owners(i,2)-owners(i,1);		
	}

	return 0;
}

// ------------------------------------------------------------------------
void PatchTrainingData::extractMatrix(const MatrixOwnerSet &input, const unsigned int exclude,
		MatrixType &X)
{
	IntMatrixType owners = input.second;
	MatrixType mat = input.first;
	unsigned int rows = mat.rows() - ownerCount(owners, exclude);
	unsigned int cols = mat.cols();

	X = MatrixType(rows, cols);

	unsigned int outputCount = 0;
	for(unsigned int i = 0; i < owners.rows(); i++)
	{
		// check for the exclusion
		if(owners(i,0) == (int) exclude) continue;

		unsigned int startIdx = owners(i,1);		
		unsigned int endIdx = owners(i,2);		

		for(unsigned int j = startIdx; j < endIdx; j++)
		{
			X.row(outputCount) = mat.row(j);			
			outputCount++;
		}
	}
}

}
