#include <iostream>
#include <MatrixCommon.h>
#include <MatrixReader.h>
#include <ValveIO.h>
#include <TrainingData.h>
#include <SVMClassifier.h>


using namespace vt;
int main(int argc, char **argv)
{
	// load the data matrix 
	std::string dataMatrixFilename = argv[1];

	std::cout << "Loading Data" << std::endl;
	TrainingData::Pointer data = TrainingData::Load(dataMatrixFilename);
	std::vector<int> owners = data->OwnerList();
	std::sort(owners.begin(), owners.end());

	for(unsigned int i = 0; i < owners.size()-1; i++)
	{
		int testData = owners[i];
		
		typedef TrainingData::MatrixType MatrixType;
		typedef TrainingData::IntMatrixType IntMatrixType;
		MatrixType X, X_test;
		IntMatrixType y, y_test;
		data->GetTrainingData(testData,X,y);
		data->GetTestData(testData, X_test, y_test);


		//std::cout << "Training SVM" << std::endl;
		SVMClassifier::Pointer classifier = SVMClassifier::New();
		classifier->Train(X,y);

		MatrixType probs;
		IntMatrixType classes;

		unsigned int lastIndex = X_test.rows()-1;
		MatrixType xNew = X_test.block(lastIndex,0,1,X_test.cols());
		//std::cout << "Testing" << std::endl;
		classifier->PredictProbability(xNew, classes, probs);

		std::cout << "Index: " <<  testData <<  " Class: " << classes(0,0) <<  " Probs: " << probs(0,0) << " / " << probs(0,1) << std::endl;

	
	}



	return 0;
}
