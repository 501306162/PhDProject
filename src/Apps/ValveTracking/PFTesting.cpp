#include <iostream>
#include <time.h>

#include "PFFunctions.h"

/*
#include <ResultViewer.h>
#include <chrono>
#include <random> 
#include <PointLocations.h>
#include <TestData.h>
#include <CommonDefinitions.h>
#include <PatchExtractor2.h>
#include <CMRFileExtractor.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPolyDataWriter.h>
#include <vtkVertexGlyphFilter.h>
#include <itkNormalVariateGenerator.h>
#include <InitialParticleGenerator.h>
#include <LengthData.h>
#include <ParameterHelpers.h>
#include "ExtractTrainingFeatures.h"
#include "FittingFunctions.h"
#include <itkGaussianMembershipFunction.h>
#include <itkCovarianceSampleFilter.h>
#include <itkSample.h>

using namespace vt;




typedef InitialParticleGenerator::Particle Particle;
typedef InitialParticleGenerator::ParticleList ParticleList;
typedef std::pair<double, Particle> PProbs;
typedef std::vector<PProbs> ParticleCostList;


bool prob_sort(const PProbs &a, const PProbs &b)
{
	return (a.first > b.first);
}

typedef itk::Point<double, 3> PointType;

typedef itk::Vector<double, 3> MVType;
typedef itk::Statistics::ListSample<MVType> SampleType;
typedef itk::Statistics::GaussianMembershipFunction<MVType> MembershipType;


typedef itk::Image<unsigned short, 3> ImageType;
PointType projectPoint(const PointType &p, const ImageType::Pointer &image);
typedef std::map<std::string, ParticleList> ParticleMap;
typedef ImageType::IndexType IndexType;

//struct ByIndex : public std::binary_function<IndexType, IndexType, bool>
//{
	//bool operator()(const IndexType &l, const IndexType &r) const
	//{
		
	//}

//}



typedef TestData::TransformType TransformType;
typedef PointLocations::MatrixType MatrixType;
void resample(const ParticleCostList &input, ParticleCostList &output);
void createInitialSamples(const MatrixType &points1, const ImageType::Pointer &image, 
		const TransformType::Pointer &transform, 
		unsigned int numSamples, std::vector<PointType> &samples);

void updateParticles(const ParticleCostList &input, const ImageType::Pointer &image, 
		PatchParams &params, PointData::Pointer &pointData, ParticleList &output);

void getOtherParticleDist(const ParticleMap &map, const std::string &type, MembershipType::Pointer &membership);

*/

int main(int argc, char ** argv)
{
	const std::string inputDataFolder = argv[1];
	const std::string parameterFilename = "/home/om0000/ValveTracking/TrainingData/conf6/extract_features.json";
	const std::string pointLocationsFilename = "/home/om0000/ValveTracking/TrainingData/PointLocations.hdf5";
	const std::string pointModelFolder = "/home/om0000/ValveTracking/SortedLines";
	const unsigned int numberOfParticles = 1000;


	// load the nessesary input data
	InputData inputData;
	loadInputData(inputDataFolder, inputData);
	PointLocations::Pointer pointLocations = PointLocations::Load(pointLocationsFilename);
	PatchParams params(parameterFilename, true);
	PointData::Pointer pointModel = PointData::Load(pointModelFolder, inputData.id, true);


	/*

	params.timeStep = 0;


	const unsigned int numParticles = 10000;

	// load the data 
	TestData::Pointer testData = TestData::Initialise(argv[1]);
	testData->Prepare(params.timeStep);

	// load the point data

	ParticleMap particleMap;
	for(unsigned int i = 0; i < types.size(); i++)
	{
		// get the images
		const std::string type = types[i];
		TestData::ImageGroup startingImages;
		testData->GetImageGroup(type, startingImages);

		// extract the points
		MatrixType points1, points2;
		pointLocations->GetPoints(testData->GetId(), "MV-" + type, params.timeStep, points1, points2);

		// generate the particles
		InitialParticleGenerator::Pointer generator = InitialParticleGenerator::New();
		generator->SetImage(startingImages.image);
		generator->SetNumParticles(numParticles);
		generator->SetTransform(testData->GetTransform());
		generator->SetPoints(points1, points2);
		generator->Generate();

		ParticleList particles = generator->GetOutput();
		particleMap[type] = particles;
	}

		
	for(unsigned int runs = 0; runs < 3; runs++)
	{
		

	for(params.timeStep = 0; params.timeStep < params.numberOfTimeSteps; params.timeStep++)
	{
		// load the classifiers
		PatchTrainingData::Pointer patchData = PatchTrainingData::New();
		ClassifierMap classifiers;
		loadOrTrainClassifiers(params, testData->GetId(), patchData, classifiers);
	

		testData->Prepare(params.timeStep);

		std::vector<std::vector<PointType> > dispPoints;

		// iterate through the map 
		ParticleMap::iterator mapIt = particleMap.begin();
		while(mapIt != particleMap.end())
		{
			const std::string type = mapIt->first;
			ParticleList &particles = mapIt->second;
			ClassifierList cls = classifiers["MV-"+type];

			MembershipType::Pointer otherMembers = MembershipType::New();
			getOtherParticleDist(particleMap, type, otherMembers);
			

			TestData::ImageGroup currentImages;
			testData->GetImageGroup(type, currentImages);

			typedef std::pair<std::pair<int, int>, double> MatrixIndexPair;
			typedef std::map<std::pair<int, int>, double> MatrixIndexMap;
			MatrixIndexMap ind1, ind2;

			// compute the particle cost
			ParticleCostList particleCosts;
			for(unsigned int i = 0; i < particles.size(); i++)
			{
				Particle &p = particles[i];
				MatrixType p1Feature, p2Feature;
				ValveLine<3>::Pointer valve = ValveLine<3>::New();
				valve->SetP1(p.p1);
				valve->SetP2(p.p2);
				valve->SetImage(currentImages.image);
				valve->UpdateIndexs();


				ImageType::IndexType pind1, pind2;
				currentImages.image->TransformPhysicalPointToIndex(p.p1, pind1);
				currentImages.image->TransformPhysicalPointToIndex(p.p2, pind2);

				std::pair<int, int> key1(pind1[0], pind1[1]);
				std::pair<int, int> key2(pind2[0], pind2[1]);

				double p1Prob, p2Prob;
				if(ind1.count(key1))
				{
					p1Prob = ind1[key1];
				}
				else
				{
					extractLBPFeature(params, valve, p.p1, p1Feature);
					IntMatrixType classes;
					MatrixType probs;
					cls[0]->PredictProbability(p1Feature, classes, probs);
					p1Prob = probs(0,1);
					ind1[key1] = p1Prob;
				}


				if(ind2.count(key2))
				{
					p2Prob = ind2[key2];
				}
				else
				{
					extractLBPFeature(params, valve, p.p2, p2Feature);
					IntMatrixType classes;
					MatrixType probs;
					cls[1]->PredictProbability(p2Feature, classes, probs);
					p2Prob = probs(0,1);
					ind2[key2] = p2Prob;
				}

				MVType center;
				for(unsigned int j = 0; j < 3; j++)
				{
					center[j] = p.p1[j] + (0.5 * (p.p2[j] - p.p1[j]));
				}

				double centerProb = otherMembers->Evaluate(center);

				double prob = p1Prob * p2Prob; // * centerProb;

				PProbs ou;
				ou.first = prob;
				ou.second = p;
				particleCosts.push_back(ou);
			}


			// resample the particles
			ParticleCostList resampled;
			resample(particleCosts, resampled);

			
			std::vector<PointType> subPoints;
			for(unsigned int pn = 0; pn < resampled.size(); pn++)
			{
				subPoints.push_back(resampled[pn].second.p1);
				subPoints.push_back(resampled[pn].second.p2);
			}

			dispPoints.push_back(subPoints);



			// set the particles for the next turn
			ParticleList updatedParticles;
			updateParticles(resampled, currentImages.image, params, pointModel, updatedParticles);

			particleMap[type] = updatedParticles;


			++mapIt;
		}



		


		std::vector<ImageType::Pointer> dispImages;
		testData->GetImages(dispImages);

		ResultViewer::Pointer viewer = ResultViewer::New();
		viewer->SetImages(dispImages);
		viewer->AddPoints(dispPoints);
		viewer->SetUp();

		std::stringstream ss1;
		ss1 << runs << "-" << params.timeStep << ".png";
		viewer->Save(ss1.str());

		



	}
	}


	
	*/
	return 0;
}

/*
// ------------------------------------------------------------------------
void updateParticles(const ParticleCostList &input, const ImageType::Pointer &image, 
		PatchParams &params, PointData::Pointer &pointData, ParticleList &output)
{
	for(unsigned int i = 0; i < input.size(); i++)
	{
		const Particle &particle = input[i].second;

		PointType up1, up2;
		pointData->GetUpdate("MV-4C", params.timeStep, up1, up2);

		PointData::VectorType x = particle.p2 - particle.p1;
		x.Normalize();

		PointData::VectorType z;
		for(unsigned int j = 0; j < 3; j++)
		{
			z[j] = image->GetDirection()(j,2);
		}
		

		double trace = 0.0;
		for(unsigned int j = 0; j < 3; j++)
		{
			trace += image->GetDirection()(j,j);
		}

		PointData::VectorType y = itk::CrossProduct(x,z);
		if(trace > 0) y = -y;

		//std::cout << up1 << " " << up2 << std::endl;
		Particle outP;
		for(unsigned int j = 0; j < 3; j++)
		{
			outP.p1[j] = particle.p1[j] + (up1[0] * x[j]);
			outP.p1[j] += (up1[1] * y[j]);
			outP.p2[j] = particle.p2[j] + (up2[0] * x[j]);
			outP.p2[j] += (up2[1] * y[j]);
		}

		outP.length = outP.p1.EuclideanDistanceTo(outP.p2);
		outP.direction = outP.p2-outP.p1;

		output.push_back(outP);
	}
}



// ------------------------------------------------------------------------
void resample(const ParticleCostList &input, ParticleCostList &output)
{
	unsigned int N = input.size();
	double W = 0.0;
	double wMax = 0.0;
	for(unsigned int i = 0; i < N; i++)
	{
		if(input[i].first > wMax)
			wMax = input[i].first;
		W += input[i].first;
	}

	// get a starting index 
    std::random_device rd;
    std::mt19937 gen(rd());
	typedef std::uniform_int_distribution<unsigned int> D;
	D d(0, N-1);

	typedef std::uniform_real_distribution<double> DR;
	DR dr(0.0, wMax);


	unsigned int index = d(gen);
	double beta = 0.0;
	for(unsigned int i = 0; i < N; i++)
	{
		beta += dr(gen);		
		while(input[index].first < beta)
		{
			beta -= input[index].first;
			index = (index+1) % N;
		}
		output.push_back(input[index]);
	}
}

// ------------------------------------------------------------------------
void getOtherParticleDist(const ParticleMap &map, const std::string &type, MembershipType::Pointer &membership)
{
	ParticleMap::const_iterator mapIt = map.begin();
	SampleType::Pointer sample = SampleType::New();
	while(mapIt != map.end())
	{
		if(mapIt->first != type)
		{
			for(unsigned int pnum = 0; pnum < mapIt->second.size(); pnum++)
			{
				const Particle &p = mapIt->second[pnum];
				MVType center;
				for(unsigned int i = 0; i < 3; i++)
				{
					center[i] = p.p1[i] + (0.5 * (p.p2[i]-p.p1[i]));
				}

				sample->PushBack(center);
				
			}
		}
		++mapIt;
	}


	typedef itk::Statistics::CovarianceSampleFilter<SampleType> CovFilterType;
	CovFilterType::Pointer covFilter = CovFilterType::New();
	covFilter->SetInput(sample);
	covFilter->Update();

	membership = MembershipType::New();
	membership->SetMean(covFilter->GetMean());
	membership->SetCovariance(covFilter->GetCovarianceMatrix());
}
*/

