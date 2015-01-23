#ifndef PARAMETER_HELPERS_H
#define PARAMETER_HELPERS_H

#include "ConfigParser.h"
#include <iostream>
#include <vector>

namespace vt 
{
class PatchExtractorParameters
{
public:
	unsigned int timeStep;
	std::string inputDirectory;
	std::string outputDirectory;
	std::string valveType;
	std::vector<std::string> valveSubDirectories;
	std::vector<std::string> valveSubTypes;
	std::string featureType;
	unsigned int lbpGridSize;
	unsigned int lbpNeighbors;
	double lbpRadius;
	unsigned int segmentationPatchSize;
	unsigned int negativeFeatureRegionSize;
	unsigned int featurePatchSize;
	double featurePatchDistance;
	double negativeIgnoreDistance;
	std::string outputFilenameTemplate;
	std::string name;
	unsigned int numberOfTimeSteps;
	unsigned int featureSize;

	PatchExtractorParameters() {}
	PatchExtractorParameters(const std::string &configFilename);
	std::string getOutputFilename(const unsigned int time, const std::string &type);
	void print();
};



}

#endif
