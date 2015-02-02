#include "PFFunctions.h"
#include <Directory.h>
#include <ValveIO.h>


int main(int argc, char ** argv)
{
	const std::string inputDataFolder = argv[1];
	const std::string parameterFilename = "/home/om0000/ValveTracking/TrainingData/conf6/extract_features.json";
	const std::string pointLocationsFilename = "/home/om0000/ValveTracking/TrainingData/PointLocations.hdf5";
	const std::string pointModelFolder = "/home/om0000/ValveTracking/SortedLines";
	const unsigned int numberOfParticles = 10000;
	const unsigned int numberOfIterations =3;



	




	// load the nessesary input data
	InputData inputData;
	loadInputData(inputDataFolder, inputData);
	PointLocations::Pointer pointLocations = PointLocations::Load(pointLocationsFilename);
	PatchParams params(parameterFilename, true);
	PointData::Pointer pointModel = PointData::Load(pointModelFolder, inputData.id, true);


	// get the ground truth folders
	GTMap groundTruth;
	loadGroundTruth(inputData, groundTruth);






	// create the initial particles 
	InitialParticleGenerator::ParticleMap particleMap;
	createInitialParticles(inputData, pointLocations, numberOfParticles, particleMap);


	for(unsigned int iteration = 0; iteration < numberOfIterations; iteration++)
	{
		for(params.timeStep = 0; params.timeStep < params.numberOfTimeSteps; params.timeStep++)
		{
			std::cout << "Iteration: " << iteration << " time step: " << params.timeStep << std::endl;
			// load the classifier for this time step
			PatchTrainingData::Pointer patchData = PatchTrainingData::New();
			ClassifierMap classifiers;
			loadOrTrainClassifiers(params, 1000, patchData, classifiers);	



			PCostList pCostList;
			pCostList.resize(numberOfParticles);
			std::vector<std::vector<std::pair<PointType, PointType> > > viewList;
			std::vector<std::vector<double> > viewWeights;
			std::vector<std::pair<PointType, PointType> > viewGT;
			std::vector<std::pair<PointType, PointType> > viewRes;

			double c1sum = 0, c2sum = 0, c3sum = 0;
			

			// move through the particle map
			InitialParticleGenerator::ParticleMap::iterator mapIt = particleMap.begin();
			while(mapIt != particleMap.end())
			{
				const std::string type = mapIt->first;
				ImageType::Pointer image = inputData.valveTypeMap[type][params.timeStep];
				ParticleList &particles = mapIt->second;
				ClassifierList cls = classifiers[type];


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
					valve->SetImage(image);
					valve->UpdateIndexs();


					itk::ContinuousIndex<double, 3> pind1, pind2;
					image->TransformPhysicalPointToContinuousIndex(p.p1, pind1);
					image->TransformPhysicalPointToContinuousIndex(p.p2, pind2);

					std::pair<int, int> key1((int) (pind1[0] * 2.0) , (int) (2.0 *pind1[1]) );
					std::pair<int, int> key2((int) (pind2[0] * 2.0) , (int) (2.0 *pind2[1]) );

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

					double prob = p1Prob * p2Prob; // * centerProb;

					PProbs ou;
					ou.first = prob;
					ou.second = p;
					particleCosts.push_back(ou);
				}

			

				ParticleCostList resampled;
				resample(particleCosts, resampled);
				

				// get the costs
				PointType wp1, wp2;
				computeWeigthedMeanPoints(resampled, wp1, wp2);

				std::vector<double> costs;
				compareValve(groundTruth[type]->GetValveLine(params.timeStep), wp1, wp2, costs);
	
				std::cout << type << ": " << costs[0] << " " << costs[1] << " " << costs[2] << std::endl;
				c1sum+=costs[0];
				c2sum+=costs[1];
				c3sum+=costs[2];


				std::pair<PointType, PointType> gt,res;
				gt.first = groundTruth[type]->GetValveLine(params.timeStep)->GetP1();
				gt.second = groundTruth[type]->GetValveLine(params.timeStep)->GetP2();

				res.first = wp1;
				res.second = wp2;
				viewGT.push_back(gt);
				viewRes.push_back(res);



				std::vector<std::pair<PointType, PointType> > vlist;
				std::vector<double> vweights;
				
				for(unsigned int r = 0; r < resampled.size(); r++)
				{
					std::pair<PointType, PointType> vp;
					vp.first = resampled[r].second.p1;
					vp.second = resampled[r].second.p2;
					vweights.push_back(resampled[r].first);

					vlist.push_back(vp);
				}
				viewWeights.push_back(vweights);

				viewList.push_back(vlist);


				ParticleList updatedParticles;
				updateParticles(resampled, type, image, params, pointModel, updatedParticles);

				particleMap[type] = updatedParticles;

				++mapIt;
			}

			std::cout << "MSE: " << c3sum / static_cast<double>(inputData.valveTypes.size()) << std::endl;

			/*
			std::cout << pCostList.size() << std::endl;

			PCostList resampled;
			resample2(pCostList, resampled);


			std::vector<std::vector<std::pair<PointType, PointType> > > viewList;
			viewList.resize(inputData.valveTypeMap.size());
			for(unsigned int i = 0; i < resampled.size(); i++)
			{
				PCostMap::const_iterator mapIt = resampled[i].begin();
				int count = 0;
				while(mapIt != resampled[i].end())
				{
					std::pair<PointType, PointType> pt;
					pt.first = mapIt->second.second.p1;
					pt.second = mapIt->second.second.p2;
					viewList[count].push_back(pt);
					++mapIt; ++count;
				}
			}


			updateParticles2(resampled, inputData, params, pointModel, particleMap); 


			*/



			InitialParticleGenerator::ParticleMap::iterator mapIt2 = particleMap.begin();
			std::vector<ImageType::Pointer> images;
			std::vector<std::vector<std::pair<PointType, PointType> > > points;
			while(mapIt2 != particleMap.end())
			{
				images.push_back(inputData.valveTypeMap[mapIt2->first][params.timeStep]);
				++mapIt2;
			}

			std::stringstream ss;
			ss << "img" << iteration*params.numberOfTimeSteps+params.timeStep << ".png";

			ResultViewer::Pointer viewer = ResultViewer::New();
			viewer->SetImages(images);
			viewer->AddPoints(viewList);
			viewer->SetWeights(viewWeights);
			viewer->SetGT(viewGT);
			viewer->SetRes(viewRes);
			viewer->SetUp();
			viewer->Save(ss.str());
			//
			// c


		}		
	}





	return 0;
}
