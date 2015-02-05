#ifndef PF_FUNCTIONS_H
#define PF_FUNCTIONS_H

#include <iostream>
#include <vector>
#include <CMRFileExtractor.h>
#include <itkSimilarity3DTransform.h>
#include <PointLocations.h>
#include <InitialParticleGenerator.h>
#include <CommonDefinitions.h>
#include <ResultViewer.h>
#include <ValveLine.h>

#include "FittingFunctions.h"
#include "ExtractTrainingFeatures.h"

using namespace vt;


typedef itk::Similarity3DTransform<double> TransformType;
typedef InitialParticleGenerator::ParticleMap ParticleMap;
typedef std::map<std::string, std::vector<ImageType::Pointer> > ValveImageMap;

typedef std::map<std::string, ValveSequence<3>::Pointer> GTMap;

typedef InitialParticleGenerator::Particle Particle;
typedef InitialParticleGenerator::ParticleList ParticleList;
typedef std::pair<double, Particle> PProbs;
typedef std::vector<PProbs> ParticleCostList;

typedef std::map<std::string, PProbs> PCostMap;
typedef std::vector<PCostMap> PCostList;

// structure to hold the input data
typedef struct _in_data_
{
	std::vector<std::string> imageTypes;
	ValveImageMap valveTypeMap;
	std::vector<std::string> valveTypes;
	CMRFileExtractor::Pointer extractor;
	TransformType::Pointer transform;
	unsigned int id;


} InputData;

void loadGroundTruth(const InputData &input, GTMap &groundTruth);
void loadInputData(const std::string &folderName, InputData &inputData);
void createInitialParticles(const InputData &data, PointLocations::Pointer &pointLocations, 
		const unsigned int numParticles, ParticleMap &particles);

void SaveParticles(const InitialParticleGenerator::ParticleList &particles, const std::string &filename);
void resample(const ParticleCostList &input, ParticleCostList &output);
void resample2(const PCostList &input, PCostList &output);
void updateParticles(const ParticleCostList &input, const std::string type, const ImageType::Pointer &image, 
		PatchParams &params, PointData::Pointer &pointData, ParticleList &output);

void updateParticles2(const PCostList &input,   InputData &inputData, 
		PatchParams &params, PointData::Pointer &pointData, ParticleMap &output);
double prodCost(const PCostMap &map);
void computeWeigthedMeanPoints(const ParticleCostList &input, PointType &p1, PointType &p2);
void compareValve(const ValveLine<3>::Pointer &input, const PointType &p1, const PointType &p2, std::vector<double> &dist);

void computeParticleCosts(const PatchParams &params, const ImageType::Pointer &image,
		const ClassifierList &classifiers, const ParticleList &particles, ParticleCostList &costs); 

#endif
