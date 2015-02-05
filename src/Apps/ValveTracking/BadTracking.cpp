#include <iostream>
#include <ValveIO.h>
#include <CMRFileExtractor.h>
#include <PatchExtractor2.h>

#include "PFFunctions.h"


using namespace vt;
typedef ValveSequence<3> Sequencetype;
typedef Sequencetype::ValveLineType ValveType;
typedef ValveType::ImageType ImageType;
typedef ValveType::PointType PointType;
typedef std::pair<PointType, PointType> PointPair;
typedef std::map<std::string, PointPair> PositionsMap;
typedef std::pair<ImageType::Pointer, ImageType::Pointer> PatchPair;
typedef std::map<std::string, PatchPair> PatchMap;


void extractPatch(const PointType &p1, const PointType &p2, const ImageType::Pointer &image, 
		const PointType &location, ImageType::Pointer &patch);


int main(int argc, char ** argv)
{
	// load the ground truth and the image data 
	const std::string testDataFolder = argv[1];
	InputData inputData;
	loadInputData(testDataFolder, inputData);
	

	// load hte ground truths for each valve type
	GTMap groundTruths;
	loadGroundTruth(inputData, groundTruths);



	PatchMap patchMap;
	PositionsMap positions;
	GTMap::iterator mapIt = groundTruths.begin();
	while(mapIt != groundTruths.end())
	{
		PointType p1 = mapIt->second->GetValveLine(0)->GetP1();
		PointType p2 = mapIt->second->GetValveLine(0)->GetP2();
		PointPair pp(p1,p2);
		positions[mapIt->first] = pp;

		PatchPair patches;
		extractPatch(p1, p2, mapIt->second->GetValveLine(0)->GetImage(), p1, patches.first);
		extractPatch(p1, p2, mapIt->second->GetValveLine(0)->GetImage(), p2, patches.second);

		patchMap[mapIt->first] = patches;

		++mapIt;
	}


	// start iterating through the images
	for(unsigned int timeStep = 1; timeStep < 25; timeStep++)
	{
		// iterate through the maps
		for(unsigned int typeId = 0; typeId < inputData.valveTypes.size(); typeId++)
		{
			const std::string type = inputData.valveTypes[typeId];
			PointPair pointPar = positions[type];
			PatchPair patches = patchMap[type];

			// compute the 

				



			
		}
			
	}


	return 0;
}


// ------------------------------------------------------------------------
void extractPatch(const PointType &p1, const PointType &p2, const ImageType::Pointer &image, 
		const PointType &location, ImageType::Pointer &patch)
{
	typedef PatchExtractor2<ImageType> ExtractorType;
	ExtractorType::Pointer extractor = ExtractorType::New();

	ExtractorType::VectorType line = p2-p1;
	line.Normalize();

	ExtractorType::DistanceType distance;
	distance.Fill(30);
	distance[2] = 0;

	ExtractorType::SizeType size;
	size.Fill(20);
	size[2] = 1;

	extractor->SetInput(image);
	extractor->SetSize(size);
	extractor->SetLine(line);
	extractor->SetDistance(distance);
	extractor->SetCenter(location);
	extractor->Update();

	patch = extractor->GetOutput();

}


