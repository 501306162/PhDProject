#include <iostream>

#include <CandidateLineFinder.h>
#include <ResultViewer.h>
#include "ExtractTrainingFeatures.h"
#include "FittingFunctions.h"
#include <TestData.h>
#include <TrainingData.h>
#include <ValveLine.h>
#include <ValveNormaliser.h>
#include <FlipChecker.h>
#include <itkResampleImageFilter.h>
#include <ValveIO.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <CommonDefinitions.h>
#include <nlopt.hpp>
#include <itkThresholdImageFilter.h>
#include <ValvePointCandidateFinder.h>


#include <itkMaximumImageFilter.h>
#include <itkRegionalMaximaImageFilter.h>

#include <LengthData.h>

#include <itkBinaryContourImageFilter.h>

int main(int argc, char ** argv)
{
	const std::string patchParamsFilename = argv[1];
	const std::string testDataFilename = argv[2];
	const std::string valveDataFilename = "/home/om0000/ValveTracking/TrainingData/MV-Planes.hdf5";
	PatchParams params(patchParamsFilename);

	TestData::Pointer testData = TestData::Initialise(testDataFilename);
	LengthData::Pointer lengthData = LengthData::Load("/home/om0000/ValveTracking/TrainingData/Lengths");
	lengthData->Prepare(testData->GetId());

	std::cout << "Loading Plane Data" << std::endl;
	ValveTrainingData::Pointer valveData = ValveTrainingData::Load(valveDataFilename);


	std::cout << "Loading PatchData" << std::endl;
	PatchTrainingData::Pointer patchData = PatchTrainingData::New();

	VectorType startPoint, startNormal;


	for(params.timeStep = 0; params.timeStep < 1; params.timeStep++)
	{


		std::cout << "Working on: " << params.timeStep << std::endl;


		// load the input points and create the bounding box
		MatrixType pointsData;
		valveData->GetTrainingPoints(testData->GetId(), params.timeStep, pointsData);

		testData->Prepare(params.timeStep);

		vtkBoundingBox boundingBox;
		computeBoundingBox(pointsData, testData->GetTransform(), boundingBox);
		testData->SetBoundingBox(boundingBox);
		testData->ProcessImages();


		std::cout << "Loading / Training Classifiers" << std::endl;
		ClassifierMap classifiers;
		loadOrTrainClassifiers(params, 1000, patchData, classifiers);


		// compute the distribution for the plane and get the mean for initialisation
		MatrixType planeData;
		valveData->GetTrainingData(testData->GetId(), params.timeStep, planeData);
		PlaneCovarianceFilterType::Pointer planeDistribution = PlaneCovarianceFilterType::New();
		computePlaneCovariance(planeData, planeDistribution);


		TestData::VectorType point, normal;

		if(params.timeStep != 10000)
		{
			for(unsigned int i = 0; i < 3; i++)
			{
				point(i) = planeDistribution->GetMean()[i];
				normal(i) = planeDistribution->GetMean()[i+3];
			}
			//startingPlane(pointsData2, point, normal);
		}
		else
		{
			point = startPoint;
			normal = startNormal;
		}


		PointType itkp;
		itk::CovariantVector<double,3> itkn;
		for(unsigned int pn = 0; pn < 3; pn++)
		{
			itkp[pn] = point(pn);
			itkn[pn] = normal(pn);
		}

		itkp = testData->GetTransform()->TransformPoint(itkp);
		itkn = testData->GetTransform()->TransformCovariantVector(itkn);

		for(unsigned int pn = 0; pn < 3; pn++)
		{
			point(pn) = itkp[pn];
			normal(pn) = itkn[pn];
		}



		savePlane(point, normal, "plane.vtk");
		



		// get the best lines for each of the types
		TestData::LineGroup testLines;
		testData->GetLinesFromPlane(normal, point, testLines);

		TestData::LineGroup::iterator lineIt = testLines.begin();
		while(lineIt != testLines.end())
		{
			std::string lineType = lineIt->first;

			ValveType::Pointer bestLine = ValveType::New();
			TestData::ImageGroup images;
			testData->GetImageGroup(lineType, images);
			createLine(lineIt->second[0], images.image, bestLine);


			// get the candidate line
			typedef CandidateLineFinder::PointListType PointListType;
			std::vector<PointListType> candidatePoints;

			for(unsigned int pnum = 0; pnum < 2; pnum++)
			{

				RealImageType::Pointer prob = RealImageType::New();
				computeProbImage(params, pnum, testData->GetId(), lineType, 
						classifiers, bestLine, images.mask, prob);

				std::stringstream ss;
				ss << pnum << "-" << lineType;
				utils::RealVolumeIO::Write(ss.str()+"prob.nrrd", prob);
				utils::ImageVolumeIO::Write(ss.str()+"image.nrrd", images.image);
				utils::LabelVolumeIO::Write(ss.str()+"mask.nrrd", images.mask);




				ValvePointCandidateFinder::Pointer finder = ValvePointCandidateFinder::New();
				finder->SetInput(prob);
				finder->Compute();

				PointListType candidates;
				finder->GetOutput(candidates);
				candidatePoints.push_back(candidates);
			}

			CandidateLineFinder::Pointer lineFinder = CandidateLineFinder::New();
			lineFinder->SetType("MV-"+lineType);
			lineFinder->SetP1List(candidatePoints[0]);
			lineFinder->SetP2List(candidatePoints[1]);
			lineFinder->SetLengthData(lengthData);
			lineFinder->SetTimeStep(params.timeStep);
			lineFinder->SetOriginalLine(bestLine->GetP1(), bestLine->GetP2());
			lineFinder->Compute();

			CandidateLineFinder::OutputPointPair best = lineFinder->GetOutput(0);

			bestLine->SetP1(best.first);
			bestLine->SetP2(best.second);
			bestLine->UpdateIndexs();

			ValveWriter<3>::Pointer writer = ValveWriter<3>::New();
			writer->SetInput(bestLine);
			std::stringstream ss;
			ss << lineType << "-" << params.timeStep;
			writer->SetFileName(ss.str());
			writer->Write();



			++lineIt;

		}
	}






	return 0;
}
