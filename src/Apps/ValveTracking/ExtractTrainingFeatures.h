#ifndef EXTRACT_TRAINING_FEATURES_H
#define EXTRACT_TRAINING_FEATURES_H

#include <iostream>
#include <vector>
#include <ValveLine.h>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <MatrixCommon.h>
#include <ParameterHelpers.h>

using namespace vt;


typedef PatchExtractorParameters Params;




typedef utils::DoubleMatrixType MatrixType;
typedef utils::IntMatrixType IntMatrixType;
typedef std::vector<MatrixType> MatrixListType;
typedef MatrixListType FeatureListType;


typedef struct valve_features_
{
	MatrixType p1PositiveFeatures;
	MatrixType p2PositiveFeatures;
	MatrixType p1NegativeFeatures;
	MatrixType p2NegativeFeatures;

} ValveFeatures;

typedef std::pair<int, int> OwnerIndicesType;
typedef std::map<int, OwnerIndicesType> OwnerLookupType;


typedef itk::Image<unsigned short, 3> ImageType;
typedef itk::Image<unsigned char, 3> LabelType;
typedef ImageType::PointType PointType;
typedef ImageType::RegionType RegionType;
typedef ImageType::IndexType IndexType;
typedef ValveLine<3>::ContIndexType ContIndexType;
typedef std::vector<IndexType> IndexlistType;


void extractFeatures(const Params &params,
		const ValveLine<3>::Pointer &line, 
		const LabelType::Pointer &segmentation,
	   	const IndexlistType &indices, 
		MatrixType &features);

OwnerIndicesType updateOwnerLookupAndCounts(unsigned int &count, const MatrixType &features);


void buildOutputMatrix(const MatrixListType &mats, const unsigned int count, MatrixType &output);

void extractFeature(const ImageType::Pointer &mask, const ContIndexType &location, 
		unsigned int featureSize, MatrixType & feature); 

void extractLBPFeature(const Params & params, const ValveLine<3>::Pointer &line, 
		const PointType &loc, MatrixType &feature);

void computeSegmentation(const ValveLine<3>::Pointer &input, LabelType::Pointer &segmentation, 
		unsigned int segmentationPatchSize, unsigned int pointNumber);

void extractMaskIndices(const LabelType::Pointer &mask, RegionType &region, 
		const ContIndexType &positiveLocation, const double distance, IndexlistType &indices);

RegionType getNegativeFeatureRegion(const ValveLine<3>::Pointer &valve,
		const unsigned int size, const unsigned int pointtoConsider);

void processDataSet(const Params &params, const unsigned int folderNumber, 
		const unsigned int time, utils::MatrixDataSet::Pointer &dataSet);

void processInstance(const Params &params, const unsigned int id, 
		const std::string &type, const std::string &filename, 
		const unsigned int timeStep, ValveFeatures &features);

void extractFeaturesForPoint(const Params &params, ValveLine<3>::Pointer &line,
		const unsigned int pointId,
		MatrixType &positiveFeatures, MatrixType &negativeFeatures);


int getOwnerId(const std::string &filename);

bool fname_sort(const std::string &f1, const std::string &f2);

ValveSequence<3>::Pointer loadValve(const std::string &fname);




#endif
