#include "ExtractTrainingFeatures.h"

#include <ConfigParser.h>
#include <Directory.h>
#include <QString>
#include <QFileInfo>
#include <FlipChecker.h>
#include <ValveNormaliser.h>
#include <PatchExtractor.h>
#include <SimpleMRFSegmenter.h>
#include <LBPFeatureExtractor.h>
#include <PatchExtractor2.h>

#include <itkBinaryContourImageFilter.h>

#include <ValveIO.h>



// ------------------------------------------------------------------------
void processDataSet(const Params &params, const unsigned int folderNumber, 
		const unsigned int time, utils::MatrixDataSet::Pointer &dataSet)
{
	// get the list of valve filenames in this subdirectory
	const std::string valveDirectory = params.valveSubDirectories[folderNumber]; 
	const std::string valveType = params.valveSubTypes[folderNumber];
	utils::Directory::FilenamesType inputFilenames = utils::Directory::GetFiles(valveDirectory, ".txt");
	std::sort(inputFilenames.begin(), inputFilenames.end(), fname_sort);


	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "Processing Valve Type: " << valveType << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;


	// loop through the set of valves in this directory
	unsigned int valvesToProcess = inputFilenames.size();
	unsigned int negativeP1FeatureCount = 0, positiveP1FeatureCount = 0;
	unsigned int negativeP2FeatureCount = 0, positiveP2FeatureCount = 0;


	// these will hold the information about which points belong to who
	OwnerLookupType positiveOwnerP1Lookup, negativeOwnerP1Lookup;
	OwnerLookupType positiveOwnerP2Lookup, negativeOwnerP2Lookup;

	MatrixListType neg1List, neg2List, pos1List, pos2List;
	std::vector<unsigned int> ownerIds;

	for(unsigned int instance = 0; instance < valvesToProcess; instance++)
	{
		const std::string instanceFilename = inputFilenames[instance];
		const unsigned int ownerId = getOwnerId(instanceFilename);
		ownerIds.push_back(ownerId);
		std::cout << "Extracting Features For: " << ownerId << std::endl;

		// extract the features for this valve at this time step
		ValveFeatures features;
		processInstance(params, ownerId, valveType, instanceFilename, time, features);


		OwnerIndicesType poisitiveP1OwnerIndices = updateOwnerLookupAndCounts(
				positiveP1FeatureCount, features.p1PositiveFeatures);
		positiveOwnerP1Lookup[ownerId] = poisitiveP1OwnerIndices;

		OwnerIndicesType poisitiveP2OwnerIndices = updateOwnerLookupAndCounts(
				positiveP2FeatureCount, features.p2PositiveFeatures);
		positiveOwnerP2Lookup[ownerId] = poisitiveP2OwnerIndices;

		OwnerIndicesType negativeP1OwnerIndices = updateOwnerLookupAndCounts(
				negativeP1FeatureCount, features.p1NegativeFeatures);
		negativeOwnerP1Lookup[ownerId] = negativeP1OwnerIndices;

		OwnerIndicesType negativeP2OwnerIndices = updateOwnerLookupAndCounts(
				negativeP2FeatureCount, features.p2NegativeFeatures);
		negativeOwnerP2Lookup[ownerId] = negativeP2OwnerIndices;


		neg1List.push_back(features.p1NegativeFeatures);
		neg2List.push_back(features.p2NegativeFeatures);
		pos1List.push_back(features.p1PositiveFeatures);
		pos2List.push_back(features.p2PositiveFeatures);

	}
	
	// build the output matrixes 
	MatrixType p1NegativeFeatures, p2NegativeFeatures, p1PositiveFeatures, p2PositiveFeatures;
	buildOutputMatrix(pos1List, positiveP1FeatureCount, p1PositiveFeatures);
	buildOutputMatrix(pos2List, positiveP2FeatureCount, p2PositiveFeatures);
	buildOutputMatrix(neg1List, negativeP1FeatureCount, p1NegativeFeatures);
	buildOutputMatrix(neg2List, negativeP2FeatureCount, p2NegativeFeatures);


	// build the owner matrixes
	IntMatrixType p1PosOwners(ownerIds.size(),3);
	IntMatrixType p2PosOwners(ownerIds.size(),3);
	IntMatrixType p1NegOwners(ownerIds.size(),3);
	IntMatrixType p2NegOwners(ownerIds.size(),3);
	for(unsigned int i = 0; i < ownerIds.size(); i++)
	{
		p1PosOwners(i,0) = ownerIds[i];
		p1PosOwners(i,1) = positiveOwnerP1Lookup[ownerIds[i]].first;	
		p1PosOwners(i,2) = positiveOwnerP1Lookup[ownerIds[i]].second;	

		p2PosOwners(i,0) = ownerIds[i];
		p2PosOwners(i,1) = positiveOwnerP2Lookup[ownerIds[i]].first;	
		p2PosOwners(i,2) = positiveOwnerP2Lookup[ownerIds[i]].second;	

		p1NegOwners(i,0) = ownerIds[i];
		p1NegOwners(i,1) = negativeOwnerP1Lookup[ownerIds[i]].first;	
		p1NegOwners(i,2) = negativeOwnerP1Lookup[ownerIds[i]].second;	

		p2NegOwners(i,0) = ownerIds[i];
		p2NegOwners(i,1) = negativeOwnerP2Lookup[ownerIds[i]].first;	
		p2NegOwners(i,2) = negativeOwnerP2Lookup[ownerIds[i]].second;	
	}



	dataSet->AddData("p1-neg", p1NegativeFeatures);
	dataSet->AddData("p2-neg", p2NegativeFeatures);
	dataSet->AddData("p1-pos", p1PositiveFeatures);
	dataSet->AddData("p2-pos", p2PositiveFeatures);

	dataSet->AddData("p1-neg-o", p1NegOwners);
	dataSet->AddData("p2-neg-o", p2NegOwners);
	dataSet->AddData("p1-pos-o", p1PosOwners);
	dataSet->AddData("p2-pos-o", p2PosOwners);
}

// ------------------------------------------------------------------------
void buildOutputMatrix(const MatrixListType &mats, const unsigned int count, MatrixType &output)
{
	output = MatrixType(count, mats.front().cols());
	unsigned int blockStart = 0;
	for(unsigned int i = 0; i < mats.size(); i++)
	{
		output.block(blockStart,0,mats[i].rows(), mats.front().cols()) = mats[i];
		blockStart += mats[i].rows();
	}
}



// ------------------------------------------------------------------------
OwnerIndicesType updateOwnerLookupAndCounts(unsigned int &count, const MatrixType &features)
{
	OwnerIndicesType indices;
	indices.first = count;
	count += features.rows();
	indices.second = count;
	return indices;
}



// ------------------------------------------------------------------------
void processInstance(const Params &params, const unsigned int id, 
		const std::string &type, const std::string &filename, 
		const unsigned int timeStep, ValveFeatures &features)
{

	// check if we are flipping the valve
	FlipChecker::Pointer checker = FlipChecker::New();
	bool flipImage = checker->FlipImage(type, id);
	bool flipPoints = checker->FlipPoints(type, id);
	std::cout << id << " " << type << " " << flipImage << " " << flipPoints << std::endl;

	// load and align the valve 
	ValveSequence<3>::Pointer sequence = loadValve(filename);
	ValveLine<3>::Pointer line = sequence->GetValveLine(timeStep);


	extractFeaturesForPoint(params, line, 1, features.p1PositiveFeatures, features.p1NegativeFeatures);
	extractFeaturesForPoint(params, line, 2, features.p2PositiveFeatures, features.p2NegativeFeatures);

}


void extractFeaturesForPoint(const Params &params, ValveLine<3>::Pointer &line,
		const unsigned int pointId,
		MatrixType &positiveFeatures, MatrixType &negativeFeatures)
{

	// compute the segmentation and extract the contours
	LabelType::Pointer segmentation = LabelType::New();
	computeSegmentation(line, segmentation, params.segmentationPatchSize, pointId);

	typedef itk::BinaryContourImageFilter<LabelType, LabelType> ContourFilterType;
	ContourFilterType::Pointer contourFilter = ContourFilterType::New();
	contourFilter->SetInput(segmentation);
	contourFilter->SetFullyConnected(false);
	contourFilter->SetBackgroundValue(0);
	contourFilter->SetForegroundValue(255);
	contourFilter->Update();
	LabelType::Pointer contours = contourFilter->GetOutput();


	// extract the indices that will be used for the negative samples
	typedef ImageType::RegionType RegionType;
	RegionType negativeRegion = getNegativeFeatureRegion(line, params.negativeFeatureRegionSize, pointId);
	IndexlistType negativeIndices;
	extractMaskIndices(contours, negativeRegion, line->GetIndex(pointId), params.negativeIgnoreDistance, negativeIndices);



	// extract the negative features
	MatrixType negFeatures;
	extractFeatures(params, line, segmentation, negativeIndices, negativeFeatures);




	// extract the positive features
	MatrixType positiveFeature;
	IndexlistType positiveIndices;
	IndexType positiveIndex;
	positiveIndex[0] = static_cast<int>(line->GetIndex(pointId)[0]);
	positiveIndex[1] = static_cast<int>(line->GetIndex(pointId)[1]);
	positiveIndex[2] = 0;

	positiveIndices.push_back(positiveIndex);
	
	MatrixType posFeatures;
	extractFeatures(params, line, segmentation, positiveIndices, positiveFeatures);


}


// ------------------------------------------------------------------------
void extractFeatures(const Params &params,
		const ValveLine<3>::Pointer &line, 
		const LabelType::Pointer &segmentation,
	   	const IndexlistType &indices, 
		MatrixType &features)
{
	unsigned int numFeatures = params.featureSize;
	features = MatrixType::Zero(indices.size(), numFeatures);

	for(unsigned int i = 0; i < indices.size(); i++)
	{
		IndexType ind = indices[i];
		ContIndexType contInd;
		PointType point;
		line->GetImage()->TransformIndexToPhysicalPoint(ind, point);


		MatrixType feature = MatrixType(1, params.featureSize);
		if(params.featureType == "LBP")
		{

			extractLBPFeature(params, line, point, feature);
		}

		features.row(i) = feature.row(0);
	}

}


// ------------------------------------------------------------------------
void extractLBPFeature(const Params & params, const ValveLine<3>::Pointer &line, 
		const PointType &loc, MatrixType &feature)
{
	typedef PatchExtractor2<ImageType> ExtractorType;
	ExtractorType::Pointer extractor = ExtractorType::New();

	ExtractorType::SizeType size;
	size.Fill(params.featurePatchSize);
	size[2] = 1;

	ExtractorType::DistanceType distance;
	distance.Fill(params.featurePatchDistance);
	distance[2] = 0.0;

	PointType p1 = line->GetP1();
	PointType p2 = line->GetP2();
	ExtractorType::VectorType vec = p2-p1;
	vec.Normalize();

	extractor->SetCenter(loc);
	extractor->SetSize(size);
	extractor->SetDistance(distance);
	extractor->SetLine(vec);
	extractor->SetInput(line->GetImage());
	extractor->Update();


	LBPFeatureExtractor::Pointer featureExtractor = LBPFeatureExtractor::New();
	featureExtractor->SetRadius(params.lbpRadius);
	featureExtractor->SetGridX(params.lbpGridSize);
	featureExtractor->SetGridY(params.lbpGridSize);
	featureExtractor->SetNumNeighbours(params.lbpNeighbors);
	featureExtractor->SetInput(extractor->GetOutput());
	featureExtractor->Extract(feature);
}


// ------------------------------------------------------------------------
void extractFeature(const ImageType::Pointer &mask, const ContIndexType &location, 
		unsigned int featureSize, MatrixType & feature)
{


	/*

	typedef itk::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(patch);
	writer->SetImageIO(itk::PNGImageIO::New());
	writer->SetFileName("patch.png");
	writer->Update();

	feature = MatrixType(1,featureSize*featureSize);
	itk::ImageRegionIterator<ImageType> it(patch, patch->GetLargestPossibleRegion());
	int count = 0;
	while(!it.IsAtEnd())
	{
		feature(0,count) = it.Get();	
		++it; ++count;
	}
	*/

	/*
	BinaryPatchFeatureExtractor::Pointer featureBuilder = BinaryPatchFeatureExtractor::New();
	featureBuilder->SetInput(patch);
	featureBuilder->Extract(feature);
	*/
}



// ------------------------------------------------------------------------
int getOwnerId(const std::string &filename)
{
	// get the id of the examination
	QFileInfo info(QString::fromStdString(filename));
	QString fname = info.fileName();
	fname = fname.replace(".txt","");
	unsigned int ownerId = fname.replace("d","").toInt();

	return ownerId;
}

// ------------------------------------------------------------------------
ValveSequence<3>::Pointer loadValve(const std::string &fname)
{
	// load the valve sequence
	ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
	reader->SetFileName(fname);
	ValveSequence<3>::Pointer sequence = reader->GetOutput();

	return sequence;
}



// ------------------------------------------------------------------------
void computeSegmentation(const ValveLine<3>::Pointer &input, LabelType::Pointer &segmentation, 
		unsigned int segmentationPatchSize, unsigned int pointNumber)
{
	// extract the mask that will be used for the segmentation
	ImagePatchExtractor::Pointer extractor = ImagePatchExtractor::New();
	ImagePatchExtractor::SizeType patchSize;
	patchSize.Fill(segmentationPatchSize);
	extractor->SetImage(input->GetImage());
	extractor->SetPatchSize(patchSize);

	if(pointNumber == 1)
		extractor->SetPatchCenter(input->GetInd1());
	else
		extractor->SetPatchCenter(input->GetInd2());

	LabelType::Pointer mask = extractor->ExtractMask();

	SimpleMRFSegmenter::Pointer segmenter = SimpleMRFSegmenter::New();
	segmenter->SetImage(input->GetImage());
	segmenter->SetMask(mask);
	segmenter->SetOutputValue(255);
	segmenter->SetSmoothnessCost(1.0);
	segmenter->Segment();

	segmentation = segmenter->GetOutput();
}

// ------------------------------------------------------------------------
void extractMaskIndices(const LabelType::Pointer &mask, RegionType &region, 
		const ContIndexType &positiveLocation, const double distance, IndexlistType &indices)
{
	typedef itk::ImageRegionConstIterator<LabelType> IteratorType;
	IteratorType maskIt(mask, region);

	while(!maskIt.IsAtEnd())
	{
		if(maskIt.Get() == 255)
		{
			IndexType testIndex = maskIt.GetIndex();

			// get the distance between the indices
			double sum = 0.0;
			for(unsigned int i = 0; i < 3; i++)
			{
				sum+= (positiveLocation[i]-testIndex[i])*(positiveLocation[i]-testIndex[i]);
			}
			sum = sqrt(sum);

			if(sum > distance)
				indices.push_back(maskIt.GetIndex());
		}
			
		++maskIt;
	}
}

// ------------------------------------------------------------------------
RegionType getNegativeFeatureRegion(const ValveLine<3>::Pointer &valve,
		const unsigned int size, const unsigned int pointToConsider)
{
	// get the negative feature extraction region
	ImagePatchExtractor::Pointer negativeRegionExtractor = ImagePatchExtractor::New();
	negativeRegionExtractor->SetImage(valve->GetImage());
	negativeRegionExtractor->SetPatchCenter(valve->GetIndex(pointToConsider));

	ImagePatchExtractor::SizeType negativePatchSize;
	negativePatchSize.Fill(size);
	negativeRegionExtractor->SetPatchSize(negativePatchSize);

	return negativeRegionExtractor->ExtractionRegion();
}





// ------------------------------------------------------------------------
bool fname_sort(const std::string &f1, const std::string &f2)
{
	QString s1 = QFileInfo(QString::fromStdString(f1)).fileName().replace("d","").replace(".txt","");
	QString s2 = QFileInfo(QString::fromStdString(f2)).fileName().replace("d","").replace(".txt","");

	return (s1.toInt() < s2.toInt());
}
