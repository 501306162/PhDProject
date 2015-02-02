#include <iostream>
#include <PatchTrainingData.h>
#include <Directory.h>
#include <QString>
#include <ParameterHelpers.h>
#include <SVMClassifier.h>

using namespace vt;

typedef PatchExtractorParameters PatchParams;

typedef std::vector<SVMClassifier::Pointer> ClassifierList;
typedef std::map<std::string, ClassifierList> ClassifierMap;
void saveClassifiers(const PatchParams &params, const unsigned int exclude, const ClassifierMap &classifiers);
void loadClassifiers(const PatchParams & params, const unsigned int exclude, ClassifierMap &classifiers);
void loadOrTrainClassifiers(const PatchParams & params, const unsigned int exclude,
		const PatchTrainingData::Pointer &trainingData, ClassifierMap &classifiers);

typedef utils::DoubleMatrixType MatrixType;
typedef utils::IntMatrixType IntMatrixType;


int main(int argc, char ** argv)
{
	const std::string inputDirectory = argv[1];
	utils::Directory::FilenamesType dirs = utils::Directory::GetDirectories(inputDirectory);
	PatchExtractorParameters params(argv[2], true);


	// load the patch data 
	PatchTrainingData::Pointer patches = PatchTrainingData::Load(argv[2]);



	//for(params.timeStep = 0; params.timeStep < params.numberOfTimeSteps; params.timeStep++)
	//{
		//ClassifierMap classifiers;
		//loadOrTrainClassifiers(params, 1000, patches, classifiers);
	//}
	

	for(unsigned int i = 0; i < dirs.size(); i++)
	{

		const std::string patientDir = dirs[i];
		const unsigned int id = QString::fromStdString(utils::Directory::GetFileName(patientDir)).replace("d", "").toInt();
		std::cout << "Processing: " << id << std::endl;

		// loop through the time steps
		for(params.timeStep = 0; params.timeStep < params.numberOfTimeSteps; params.timeStep++)
		{
			ClassifierMap classifiers;
			loadOrTrainClassifiers(params, id, patches, classifiers);
		}
		

	}



	return 0;


}




// ------------------------------------------------------------------------
void traingClassifiers(const unsigned int exclude, 
		const PatchParams &params, 
		const PatchTrainingData::Pointer &trainingData,
	   	ClassifierMap &classifiers)
{
	for(unsigned int i = 0; i < params.valveSubTypes.size(); i++)
	{
		ClassifierList classifierList;

		
		std::string type = params.valveSubTypes[i];		
		for(unsigned int j = 0; j < 2; j++)
		{
			MatrixType X;
			IntMatrixType y;
			std::cout << "Training: " << type << " - " << j+1 << std::endl;
			trainingData->GetTrainingData(exclude, type, params.timeStep, j+1, X, y);


			SVMClassifier::Pointer classifier = SVMClassifier::New();
			classifier->Train(X,y);
			classifierList.push_back(classifier);
		}

		classifiers[type] = classifierList;
	}
}

// ------------------------------------------------------------------------
void loadOrTrainClassifiers(const PatchParams & params, const unsigned int exclude,
		const PatchTrainingData::Pointer &trainingData, ClassifierMap &classifiers)
{
	std::stringstream ss;
	ss << exclude << "-" << params.valveSubTypes[0] << "-" << params.timeStep << "-0.cls";
	std::string testFilename = utils::Directory::GetPath(params.outputDirectory, ss.str());
	std::cout << ss.str() << std::endl;
	if(utils::Directory::FileExists(testFilename))
	{
		loadClassifiers(params, exclude, classifiers);
	}
	else
	{
		traingClassifiers(exclude, params, trainingData, classifiers);
		saveClassifiers(params, exclude, classifiers);
	}
}

// ------------------------------------------------------------------------
void saveClassifiers(const PatchParams &params, const unsigned int exclude, const ClassifierMap &classifiers)
{
	ClassifierMap::const_iterator mapIt = classifiers.begin();
	while(mapIt != classifiers.end())
	{
		for(unsigned int i = 0; i < 2; i++)
		{
			std::stringstream ss;
			ss <<  exclude << "-" << mapIt->first << "-" << params.timeStep << "-" << i <<  ".cls";
			std::string filename = utils::Directory::GetPath(params.outputDirectory, ss.str());
			mapIt->second[i]->Save(filename);
		}
		++mapIt;
	}
}

// ------------------------------------------------------------------------
void loadClassifiers(const PatchParams &params, const unsigned int exclude, ClassifierMap &classifiers)
{
	for(unsigned int i = 0; i < params.valveSubTypes.size(); i++)
	{
		ClassifierList classifierList;
		std::string type = params.valveSubTypes[i];		
		for(unsigned int j = 0; j < 2; j++)
		{
			std::stringstream ss;
			ss << exclude << "-" << type << "-" << params.timeStep << "-" << j << ".cls";
			std::string filename = utils::Directory::GetPath(params.outputDirectory, ss.str());
			SVMClassifier::Pointer classifier = SVMClassifier::New();
			classifier->Load(filename);
			classifierList.push_back(classifier);
		}
		classifiers[type] = classifierList;
	}
}
