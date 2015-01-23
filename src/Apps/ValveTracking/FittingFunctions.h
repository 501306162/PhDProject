#ifndef FITTING_FUNCTIONS_H
#define FITTING_FUNCTIONS_H

#include <iostream>


#include <SVMClassifier.h>
#include <PatchTrainingData.h>
#include <ParameterHelpers.h>
#include <TestData.h>
#include <ValveLine.h>
#include <SimpleMRFSegmenter.h>
#include <PatchExtractor.h>
#include <LengthData.h>

#include <itkGaussianMembershipFunction.h>

#include <vtkBoundingBox.h>

#include <itkListSample.h>
#include <itkCovarianceSampleFilter.h>

#include "ExtractTrainingFeatures.h"

using namespace vt;


typedef ValveLine<3> ValveType;
typedef std::vector<SVMClassifier::Pointer> ClassifierList;
typedef std::map<std::string, ClassifierList> ClassifierMap;
typedef PatchExtractorParameters PatchParams;
typedef PatchTrainingData::MatrixType MatrixType;
typedef PatchTrainingData::IntMatrixType IntMatrixType;

typedef itk::Image<unsigned short, 3> ImageType;
typedef itk::Image<unsigned char, 3> LabelType;

typedef itk::Vector<double, 6> PlaneMeasurementType;
typedef itk::Statistics::ListSample<PlaneMeasurementType> PlaneSampleType;
typedef itk::Statistics::CovarianceSampleFilter<PlaneSampleType> PlaneCovarianceFilterType;
typedef itk::Statistics::GaussianMembershipFunction<PlaneMeasurementType> PlaneMembershipType;
typedef itk::Image<double, 3> RealImageType;

typedef struct opt_data 
{
	TestData::Pointer testData;
	ClassifierMap classifiers;
	PatchParams params;
	LengthData::Pointer lengths;
	PlaneCovarianceFilterType::Pointer planeMember;

	opt_data() {}

} OptData;

typedef Eigen::Vector3d VectorType;

void traingClassifiers(const unsigned int exclude, 
		const PatchParams &params, 
		const PatchTrainingData::Pointer &trainingData,
	   	ClassifierMap &classifiers);


void FlipImage(const LabelType::Pointer &input, LabelType::Pointer &output);
void FlipImage(const RealImageType::Pointer &input, RealImageType::Pointer &output);

void createLine(const TestData::LineType &input, ImageType::Pointer &image, ValveType::Pointer &output); 
void computeProbImage(const PatchParams &params, unsigned int pnum,
	   	unsigned int id, std::string type, ClassifierMap &classifiers, 
		const ValveType::Pointer &valve, const LabelType::Pointer &mask, RealImageType::Pointer &output);

void saveClassifiers(const PatchParams &params, const unsigned int exclude, const ClassifierMap &classifiers);
void loadClassifiers(const PatchParams & params, const unsigned int exclude, ClassifierMap &classifiers);
void loadOrTrainClassifiers(const PatchParams & params, const unsigned int exclude,
		const PatchTrainingData::Pointer &trainingData, ClassifierMap &classifiers);


void startingPlane(const std::vector<MatrixType> &pointsData, VectorType &point, VectorType &normal);

void computeBoundingBox(const MatrixType &points, const TestData::TransformType::Pointer &transform, 
		vtkBoundingBox &boundingBox);


void computePlaneCovariance(const MatrixType &data, PlaneCovarianceFilterType::Pointer &cov);

void computeSegmentation(const ValveLine<3>::Pointer &input, LabelType::Pointer &segmentation, 
		unsigned int segmentationPatchSize, unsigned int pointNumber);


void convert_x(const std::vector<double> &x, VectorType &point, VectorType &normal);
void convert_plane(const VectorType &point, const VectorType &normal, std::vector<double> &x);

double imageCost(PatchParams &params, TestData::Pointer &data, 
		VectorType &point, VectorType &normal, ClassifierMap &classifiers, 
		LengthData::Pointer &lengths, bool save=false);



int getBestLines(PatchParams &params, const std::string type, 
		TestData::Pointer &testData, TestData::LineTypeList &lines, 
		ClassifierMap &classifiers, LengthData::Pointer &lengths);

double f(const std::vector<double> &x, std::vector<double> &grad, void *f_data);
double f2(const std::vector<double> &x, std::vector<double> &grad, void *f_data);


void savePlane(const MatrixType &plane, std::string filename);
void savePlane(const VectorType &point, const VectorType &normal, std::string filename);

void createAlignedValve(const TestData::LineType &line, const ImageType::Pointer &image, 
		ValveLine<3>::Pointer &output, bool flipImage, bool flipPoints);


#endif
