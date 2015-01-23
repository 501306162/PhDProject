#include <iostream>

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
	//PatchTrainingData::Pointer patchData = PatchTrainingData::Load(patchParamsFilename);
	PatchTrainingData::Pointer patchData = PatchTrainingData::New();

	VectorType startPoint, startNormal;


	for(unsigned int runs = 0; runs < 1; runs++)
	{
		
		for(params.timeStep = 0; params.timeStep < params.numberOfTimeSteps; params.timeStep++)
		{


			std::cout << "Working on: " << params.timeStep << std::endl;


			// load the input points and create the bounding box
			MatrixType pointsData;
			valveData->GetTrainingPoints(testData->GetId(), params.timeStep, pointsData);
			
			std::vector<MatrixType> pointsData2;
			valveData->GetTrainingPoints(testData->GetId(), params.timeStep, pointsData2);


			vtkBoundingBox boundingBox;


			computeBoundingBox(pointsData, boundingBox);


			std::cout << "Loading / Training Classifiers" << std::endl;
			ClassifierMap classifiers;
			loadOrTrainClassifiers(params, 1000, patchData, classifiers);



			testData->Prepare(params.timeStep, boundingBox);
			testData->SaveImageGroup("2C", "MV-2C");


			// compute the distribution for the plane and get the mean for initialisation
			MatrixType planeData;
			valveData->GetTrainingData(testData->GetId(), params.timeStep, planeData);
			PlaneCovarianceFilterType::Pointer planeDistribution = PlaneCovarianceFilterType::New();
			computePlaneCovariance(planeData, planeDistribution);


			TestData::VectorType point, normal;

			if(params.timeStep == 0 && runs == 0)
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

			OptData * optData = new OptData;
			optData->testData = testData;
			optData->classifiers = classifiers;
			optData->params = params;
			optData->lengths = lengthData;
			optData->planeMember = planeDistribution;




			nlopt::opt opt(nlopt::LN_BOBYQA, 5);
			opt.set_min_objective(f, (void*) optData);
			opt.set_xtol_rel(1e-4);
			opt.set_maxeval(100);


			std::vector<double> initStep(5);
			double pointStep = atof(argv[3]);
			double normalStep = atof(argv[4]);
			for(unsigned int i = 0; i < 3; i++)
			{
				initStep[i] = pointStep;
			}

			initStep[3] = normalStep;
			initStep[4] = normalStep;

			opt.set_default_initial_step(initStep);

			std::vector<double> xStart;
			convert_plane(point, normal, xStart);

			double minf;
			nlopt::result result = opt.optimize(xStart,minf);
			std::cout << "Final: " << result << std::endl;

			MatrixType initialPlane(1,6), finalPlane(1,6);
			VectorType finalPoint, finalNormal;
			convert_x(xStart, finalPoint, finalNormal);


			imageCost(params, testData, finalPoint, finalNormal, classifiers, lengthData, true);


			for(unsigned int i = 0; i < 3; i++)
			{
				finalPlane(0,i) = finalPoint(i);
				finalPlane(0,i+3) = finalNormal(i);
			}


			std::cout << finalPlane << std::endl;
			for(unsigned int i = 0; i < 6; i++)
			{
				initialPlane(0,i) = planeDistribution->GetMean()[i];
			}

			startPoint = finalPoint;
			startNormal = finalNormal;


			MatrixType gt;
			valveData->GetTestData(testData->GetId(), params.timeStep, gt);
			VectorType gtPoint, gtNormal;
			for(unsigned int i = 0; i < 3; i++)
			{
				gtPoint(i) = gt(0,i);
				gtNormal(i) = gt(0,i+3);
			}



			std::vector<ImageType::Pointer> vimages;
			testData->GetImages(vimages);
			ResultViewer::Pointer viewer = ResultViewer::New();
			viewer->SetImages(vimages);
			viewer->SetBoundingBox(boundingBox);
			viewer->SetPlane(finalPoint, finalNormal);
			viewer->SetStartPlane(gtPoint, gtNormal);
			viewer->SetUp();

			std::stringstream ss; 
			ss << params.name << "/" << testData->GetId() << "-" << params.timeStep << "-image.png";

			viewer->Save(ss.str());
			//viewer->View();

			//testData->SaveImageGroup("2C", "2C");
			//testData->SaveImageGroup("3C", "3C");
			//testData->SaveImageGroup("4C", "4C");
			//savePlane(initialPlane, "initialPlane.vtk");
			//savePlane(finalPlane, "finalPlane.vtk");
			

		}

	}






	return 0;
}
