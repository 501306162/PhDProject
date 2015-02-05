#include "PFFunctions.h"
#include <random>

#include <ValveOriginFinder.h>
#include <Directory.h>

#include <QString>
#include <ValveIO.h>

#include <vtkSmartPointer.h>
#include <vtkVertexGlyphFilter.h>

// ------------------------------------------------------------------------
void createInitialParticles(const InputData &data, PointLocations::Pointer &pointLocations, 
		const unsigned int numParticles, ParticleMap &particles)
{
	InitialParticleGenerator::Pointer particleGenerator = InitialParticleGenerator::New();
	particleGenerator->SetNumParticles(numParticles);
	particleGenerator->SetTransform(data.transform);


	// iterate through the image map 
	ValveImageMap::const_iterator mapIt = data.valveTypeMap.begin();
	while(mapIt != data.valveTypeMap.end())
	{
		// get the point locations
		InitialParticleGenerator::ImagePointPair points;
		points.first = mapIt->second[0];
		pointLocations->GetPoints(data.id, mapIt->first, 0, points.second.first, points.second.second);
		particleGenerator->AddInput(mapIt->first, points);
		

		++mapIt;

	}

	particleGenerator->Generate();
	particles = particleGenerator->GetMapOutput();
}

// ------------------------------------------------------------------------
void computeParticleCosts(const PatchParams &params, const ImageType::Pointer &image,
		const ClassifierList &classifiers, const ParticleList &particles, ParticleCostList &costs)
{
	typedef std::pair<std::pair<int, int>, double> MatrixIndexPair;
	typedef std::map<std::pair<int, int>, double> MatrixIndexMap;
	MatrixIndexMap ind1, ind2;

	// compute the particle cost
	ParticleCostList particleCosts;
	for(unsigned int i = 0; i < particles.size(); i++)
	{
		const Particle &p = particles.at(i);
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
			classifiers[0]->PredictProbability(p1Feature, classes, probs);
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
			classifiers[1]->PredictProbability(p2Feature, classes, probs);
			p2Prob = probs(0,1);
			ind2[key2] = p2Prob;
		}

		double prob = p1Prob * p2Prob; // * centerProb;

		PProbs ou;
		ou.first = prob;
		ou.second = p;
		costs.push_back(ou);
	}

}




// ------------------------------------------------------------------------
double prodCost(const PCostMap &map)
{
	double prod = 1.0;
	PCostMap::const_iterator mapIt = map.begin();
	while(mapIt != map.end())
	{
		prod *= mapIt->second.first;
		++mapIt;
	}

	return prod;
}

// ------------------------------------------------------------------------
void resample2(const PCostList &input, PCostList &output)
{
	std::vector<double> wcosts;
	unsigned int N = input.size();
	double W = 0.0;
	double wMax = 0.0;
	for(unsigned int i = 0; i < N; i++)
	{
		double w = prodCost(input[i]);
		wcosts.push_back(w);
		if(w > wMax)
			wMax = w;
		W += w;
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
		while(wcosts[index] < beta)
		{
			beta -= wcosts[index];
			index = (index+1) % N;
		}
		output.push_back(input[index]);
		//std::cout << output.size() << std::endl;
	}

}

// ------------------------------------------------------------------------
void loadGroundTruth(const InputData &input, GTMap &groundTruth)
{
	for(unsigned int i = 0; i < input.valveTypes.size(); i++)
	{
		std::string base = "/home/om0000/ValveTracking/SortedLines";		
		std::stringstream end;
		end << "d"  << input.id << ".txt";

		std::string fileName = utils::Directory::GetPath(
				utils::Directory::GetPath(base, input.valveTypes[i]), end.str());

		ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
		reader->SetFileName(fileName);
		groundTruth[input.valveTypes[i]] = reader->GetOutput();
	}
}

// ------------------------------------------------------------------------
void computeWeigthedMeanPoints(const ParticleCostList &input, PointType &p1, PointType &p2)
{
	double W = 0.0;
	for(unsigned int i = 0; i < input.size(); i++)
	{
		W += input[i].first;	
	}

	p1.Fill(0); p2.Fill(0);
	for(unsigned int i = 0; i < input.size(); i++)
	{
		for(unsigned int j = 0; j < 3; j++)
		{
			double w = input[i].first / W;
			p1[j] += input[i].second.p1[j] * w;
			p2[j] += input[i].second.p2[j] * w;
		}		
	}
}

// ------------------------------------------------------------------------
void compareValve(const ValveLine<3>::Pointer &input, const PointType &p1, const PointType &p2, std::vector<double> &dist)
{
	double c1 = input->GetP1().EuclideanDistanceTo(p1);
	double c2 = input->GetP2().EuclideanDistanceTo(p2);
	double c3 = c1 + c2;
	dist.push_back(c1);
	dist.push_back(c2);
	dist.push_back(c3 / 2.0);
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
void updateParticles(const ParticleCostList &input, const std::string type, const ImageType::Pointer &image, 
		PatchParams &params, PointData::Pointer &pointData, ParticleList &output)
{
	for(unsigned int i = 0; i < input.size(); i++)
	{
		const Particle &particle = input[i].second;

		PointType up1, up2;
		pointData->GetUpdate(type, params.timeStep, up1, up2);

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
void updateParticles2(const PCostList &input,  InputData &inputData, 
		PatchParams &params, PointData::Pointer &pointData, ParticleMap &output)
{
	for(unsigned int i = 0; i < input.size(); i++)
	{
		PCostMap::const_iterator mapIt = input[i].begin();
		while(mapIt != input[i].end())
		{
			const std::string type = mapIt->first;
			const Particle &particle = mapIt->second.second;
			const ImageType::Pointer image = inputData.valveTypeMap[type][params.timeStep];

			PointType up1, up2;
			pointData->GetUpdate(type, params.timeStep, up1, up2);

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

			output[type][i] = outP;

			++mapIt;
		}
	}
}





// ------------------------------------------------------------------------
void loadInputData(const std::string &folderName, InputData &inputData)
{
	std::cout << "--> Loading input data" << std::endl;
	inputData.extractor = CMRFileExtractor::New();
	inputData.extractor->SetFolderName(folderName);
	inputData.extractor->Extract();

	// check for the types of images availablefolders
	if(inputData.extractor->Has2CImage())
	{
		inputData.imageTypes.push_back("2C");
		inputData.valveTypes.push_back("MV-2C");
		for(unsigned int i = 0; i < 25; i++)
			inputData.valveTypeMap["MV-2C"].push_back(inputData.extractor->Get2CImage(i));
	}
	if(inputData.extractor->Has3CImage())
	{
		inputData.imageTypes.push_back("3C");
		inputData.valveTypes.push_back("MV-3C");
		for(unsigned int i = 0; i < 25; i++)
			inputData.valveTypeMap["MV-3C"].push_back(inputData.extractor->Get3CImage(i));
	}
	if(inputData.extractor->Has4CImage())
	{
		inputData.imageTypes.push_back("4C");
		inputData.valveTypes.push_back("MV-4C");
		inputData.valveTypes.push_back("TP-4C");
		for(unsigned int i = 0; i < 25; i++)
			inputData.valveTypeMap["MV-4C"].push_back(inputData.extractor->Get4CImage(i));
		for(unsigned int i = 0; i < 25; i++)
			inputData.valveTypeMap["TP-4C"].push_back(inputData.extractor->Get4CImage(i));
	}
	if(inputData.extractor->HasR2CImage())
	{
		inputData.imageTypes.push_back("R2C");
		inputData.valveTypes.push_back("TP-R2C");
		for(unsigned int i = 0; i < 25; i++)
			inputData.valveTypeMap["TP-R2C"].push_back(inputData.extractor->GetR2CImage(i));
	}

	


	// compute the transform 
	ValveOriginFinder::Pointer originFinder = ValveOriginFinder::New();
	originFinder->Set2CImage(inputData.extractor->Get2CImage(0));
	originFinder->Set3CImage(inputData.extractor->Get3CImage(0));
	originFinder->SetImageStack(inputData.extractor->GetStackImage(0));
	originFinder->Compute();

	// create the transform 
	inputData.transform = TransformType::New();
	inputData.transform->SetTranslation(originFinder->GetTranslation());
	inputData.transform->SetMatrix(originFinder->GetRotation());
	
	// get the id 
	inputData.id = QString::fromStdString(utils::Directory::GetFileName(folderName))
		.replace("d","").replace(".txt","").toInt();
}

// ------------------------------------------------------------------------
void SaveParticles(const InitialParticleGenerator::ParticleList &particles, const std::string &filename)
{


}


