#include <iostream>
#include <itkRescaleIntensityImageFilter.h>
#include <Directory.h>
#include <itkImageMomentsCalculator.h>

#include <ComputeProbImage.h>
#include <PatchExtractor2.h>
#include <LineIntersectionFinder.h>
#include <itkNrrdImageIO.h>
#include <itkPNGImageIO.h>
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
#include <itkTimeProbe.h>
#include <itkAndImageFilter.h>

#include <PointData.h>

#include <itkMaximumImageFilter.h>
#include <itkRegionalMaximaImageFilter.h>
#include <ResultViewer.h>

#include <LengthData.h>

#include <itkBinaryContourImageFilter.h>

typedef itk::RGBPixel<unsigned char> RGBPixelType;
typedef CandidateLineFinder::OutputPointPair PointPairType;
typedef itk::Image<RGBPixelType, 3> RGBOutputType;

typedef std::map<std::string, PointPairType> PointPairMap;

void getRGB(const ImageType::Pointer &image, const LabelType::Pointer &mask,
		const RealImageType::Pointer &prob1, const RealImageType::Pointer &prob2,
		RGBOutputType::Pointer &output);

void computePlane(std::vector<PointType> &inputPoints, VectorType &normal, VectorType &point);
void computePlane(PointPairMap &inputPoints, VectorType &normal, VectorType &point);
void run(PatchParams &params, const std::string folder);
void computeTrackingPatch(const Params &params, const ImageType::Pointer &image, 
		const SVMClassifier::Pointer &cls,
		const unsigned int size, const PointPairType &line, PointType &center,
		OptData &data, RealImageType::Pointer &prob);


void getNextPoint(const Params &params, const ImageType::Pointer &image, 
		const SVMClassifier::Pointer &cls,
		const unsigned int size, const PointPairType &line, 
		OptData &data, const PointType &start, PointType &end);


void getBestPlane(const VectorType &startPoint, const VectorType &startNormal, 
	   OptData &data, VectorType &endPoint, VectorType &endNormal, PointPairMap &bestLines);	

void trackPoints(const PointPairMap &startPoints, const PatchParams &params, OptData &data, 
		PointPairMap &endPoints);

int main(int argc, char ** argv)
{
	const std::string patchParamsFilename = argv[1];
	PatchParams params(patchParamsFilename);

	run(params, argv[2]);
	return 0;
	
	utils::Directory::FilenamesType folders = utils::Directory::GetDirectories(argv[2]);



	for(unsigned int instance = 0; instance < folders.size(); instance++)
	{
		std::string folder = folders[instance];
		run(params, folder);
	}





	return 0;
}

// ------------------------------------------------------------------------
void run(PatchParams &params, const std::string folder)
{
	const std::string valveDataFilename = "/home/om0000/ValveTracking/TrainingData/MV-Planes.hdf5";
	std::string testDataFilename = folder;

	TestData::Pointer testData = TestData::Initialise(testDataFilename);
	LengthData::Pointer lengthData = LengthData::Load("/home/om0000/ValveTracking/TrainingData/Lengths");
	lengthData->Prepare(testData->GetId());
	PointData::Pointer pointModel = PointData::Load("/home/om0000/ValveTracking/SortedLines", testData->GetId());

	std::cout << "Loading Plane Data" << std::endl;
	ValveTrainingData::Pointer valveData = ValveTrainingData::Load(valveDataFilename);


	std::cout << "Loading PatchData" << std::endl;
	PatchTrainingData::Pointer patchData = PatchTrainingData::New();

	VectorType startPoint, startNormal;





	PointPairMap trackingPoints;
	for(params.timeStep = 0; params.timeStep < params.numberOfTimeSteps; params.timeStep++)
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


		// get the starting plane
		VectorType point, normal;
		meanPlane(planeDistribution, testData->GetTransform(), point, normal);
		savePlane(point, normal, "plane.vtk");


		OptData data;
		data.classifiers = classifiers;
		data.lengths = lengthData;
		data.testData = testData;
		data.planeMember = planeDistribution;
		data.params = params;
		data.pointData = pointModel;

		VectorType newPoint, newNormal;

		if(params.timeStep == 0)
		{
			getBestPlane(point, normal, data, newPoint, newNormal, trackingPoints);
		}
		else
		{
			PointPairMap endPoints;
			trackPoints(trackingPoints, params, data, endPoints);
			trackingPoints = endPoints;
		}

		VectorType np, nn;
		computePlane(trackingPoints, nn, np);


		ResultViewer::Pointer viewer = ResultViewer::New();
		viewer->SetPlane(np,nn);
		viewer->SetTransform(testData->GetTransform());
		viewer->SetStartPlane(point,normal);
		
		std::vector<ImageType::Pointer> viewImages;
		testData->GetImages(viewImages);
		viewer->SetImages(viewImages);
		viewer->SetBoundingBox(boundingBox);
		viewer->SetUp();
		std::stringstream ss;
		ss << testData->GetId() << "-" << params.timeStep << ".png";
		viewer->Save(ss.str());

			

	}

}


// ------------------------------------------------------------------------
void trackPoints(const PointPairMap &startPoints, const PatchParams &params, OptData &data, 
		PointPairMap &endPoints)
{
	// create the valve lines 
	PointPairMap::const_iterator mapIt = startPoints.begin();
	while(mapIt != startPoints.end())
	{
		const std::string lineType = mapIt->first;
		TestData::ImageGroup images;
		data.testData->GetImageGroup(lineType, images);
	

		PointType startP1 = mapIt->second.first;
		PointType startP2 = mapIt->second.second;

		PointData::MembershipType::Pointer pmem1, pmem2;
		data.pointData->GetMembers("MV-"+lineType, params.timeStep, startP1, startP2, pmem1, pmem2); 

		double diff = 1000;
		while(diff > 0.1)
		{

			


			RealImageType::Pointer prob1 = RealImageType::New();
			computeTrackingPatch(params, images.image, data.classifiers["MV-"+lineType][0],
					10, mapIt->second, startP1, data, prob1);

			RealImageType::Pointer prob2 = RealImageType::New();
			computeTrackingPatch(params, images.image, data.classifiers["MV-"+lineType][1],
					10, mapIt->second, startP2, data, prob2);


			itk::ImageRegionIterator<RealImageType> prob1It(prob1, prob1->GetLargestPossibleRegion());
			itk::ImageRegionIterator<RealImageType> prob2It(prob2, prob2->GetLargestPossibleRegion());
			while(!prob1It.IsAtEnd())
			{
				IndexType ind = prob1It.GetIndex();
				PointType point;
				prob1->TransformIndexToPhysicalPoint(ind, point);

				PointData::VectorType mv;
				for(unsigned int i = 0; i < 3; i++)
				{
					mv[i] = point[i];
				}

				double p1 = pmem1->Evaluate(mv);
				double p2 = pmem2->Evaluate(mv);

				prob1It.Set(prob1It.Get()*p1);
				prob2It.Set(prob2It.Get()*p2);

				++prob1It; ++prob2It;
			}



			typedef itk::ImageMomentsCalculator<RealImageType> MomentsType;
			MomentsType::Pointer moments = MomentsType::New();
			moments->SetImage(prob1);
			moments->Compute();

			PointType newP1 = moments->GetCenterOfGravity();
		

			MomentsType::Pointer moments2 = MomentsType::New();
			moments2->SetImage(prob2);
			moments2->Compute();

			PointType newP2 = moments2->GetCenterOfGravity();

		
			diff = (newP1-startP1).GetNorm() + (newP2-startP2).GetNorm();
			std::cout << "Diff: " << diff << std::endl;

			startP1 = newP1;
			startP2 = newP2;


			/*
			itk::ImageRegionIterator<RealImageType> probIt1(prob1, prob1->GetLargestPossibleRegion());
			itk::ImageRegionIterator<RealImageType> probIt2(prob2, prob2->GetLargestPossibleRegion());
			while(!probIt1.IsAtEnd())
			{
				IndexType index = probIt1.GetIndex();
				PointType point;
				prob1->TransformIndexToPhysicalPoint(index, point);

				PointData::VectorType mv;
				for(unsigned int i = 0; i < 3; i++)
				{
					mv[i] = point[i];
				}

				double pprob = pmem1->Evaluate(mv);
				double oprob = probIt.Get();
				double out = (pprob*3) * oprob;

				probIt.Set(out);

				++probIt;
			}

			utils::ImageVolumeIO::Write("image.nrrd", images.image);
			utils::RealVolumeIO::Write("prob2.nrrd", prob);
			*/


		}


		
		ValveLine<3>::Pointer v = ValveLine<3>::New();
		v->SetP1(startP1);
		v->SetP2(startP2);
		v->SetImage(images.image);

		std::stringstream ss;
		ss << lineType << "-" << params.timeStep;
		ValveWriter<3>::Pointer w = ValveWriter<3>::New();
		w->SetInput(v);
		w->SetFileName(ss.str());
		w->Write();

		PointPairType endP;
		endP.first = startP1;
		endP.second = startP2;

		endPoints[lineType] = endP;
		

		++mapIt;
	}


}

void getNextPoint(const Params &params, const ImageType::Pointer &image, 
		const SVMClassifier::Pointer &cls,
		const unsigned int size, const PointPairType &line, 
		OptData &data, const PointType &start, PointType &end)
{


}

// ------------------------------------------------------------------------
void computeTrackingPatch(const Params &params, const ImageType::Pointer &image, 
		const SVMClassifier::Pointer &cls,
		const unsigned int size, const PointPairType &line, PointType &center,
		OptData &data, RealImageType::Pointer &prob)
{
	ImagePatchExtractor::Pointer extractor = ImagePatchExtractor::New();
	extractor->SetMaskValue(255);
	extractor->SetImage(image);

	ImagePatchExtractor::SizeType psize;
	psize.Fill(size);
	psize[2] =1;
	extractor->SetPatchSize(psize);

	ContIndexType index;
	image->TransformPhysicalPointToContinuousIndex(center, index);
	extractor->SetPatchCenter(index);
	
	ComputeProbImage::Pointer computer = ComputeProbImage::New();
	computer->SetInput(image);
	computer->SetMask(extractor->ExtractMask());
	
	itk::Vector<double, 3> dir = (line.second - line.first);
	dir.Normalize();

	computer->SetLine(dir);
	computer->SetRadius(params.lbpRadius);
	computer->SetNeighbors(params.lbpNeighbors);
	computer->SetGridSize(params.lbpGridSize);
	computer->SetClassifier(cls);
	
	ComputeProbImage::SizeType patchSize;
	patchSize.Fill(params.featurePatchSize);
	patchSize[2] = 1;

	ComputeProbImage::DistanceType distance;
	distance.Fill(params.featurePatchDistance);
	distance[2] = 0.0;

	computer->SetPatchDistance(distance);
	computer->SetPatchSize(patchSize);
	computer->Update();

	prob = computer->GetOutput();

}

// ------------------------------------------------------------------------
void getBestPlane(const VectorType &startPoint, const VectorType &startNormal, 
	   OptData &data, VectorType &endPoint, VectorType &endNormal,
	   PointPairMap &bestLines)
{

	// get the best lines for each of the types
	TestData::LineGroup testLines;
	data.testData->GetLinesFromPlane(startNormal, startPoint, testLines);

	TestData::LineGroup::iterator lineIt = testLines.begin();

	std::vector<PointType> allPoints;	

	while(lineIt != testLines.end())
	{
		std::string lineType = lineIt->first;

		ValveType::Pointer bestLine = ValveType::New();
		TestData::ImageGroup images;
		data.testData->GetImageGroup(lineType, images);
		createLine(lineIt->second[0], images.image, bestLine);


		// get the candidate line
		typedef CandidateLineFinder::PointListType PointListType;
		std::vector<PointListType> candidatePoints;
		RealImageType::Pointer p1 = RealImageType::New();
		RealImageType::Pointer p2 = RealImageType::New();


		for(unsigned int pnum = 0; pnum < 2; pnum++)
		{

			RealImageType::Pointer prob = RealImageType::New();
			computeProbImage2(data.params, pnum, data.testData->GetId(), lineType, 
					data.classifiers, bestLine, images.mask, prob);

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
		lineFinder->SetLengthData(data.lengths);
		lineFinder->SetTimeStep(data.params.timeStep);
		lineFinder->SetOriginalLine(bestLine->GetP1(), bestLine->GetP2());
		lineFinder->Compute();

		bestLines[lineType] = lineFinder->GetOutput();

		std::vector<PointPairType> candididateLines;

		lineFinder->GetOutput(candididateLines, 1);

		for(unsigned int i = 0; i < candididateLines.size(); i++)
		{
			allPoints.push_back(candididateLines[i].first);
			allPoints.push_back(candididateLines[i].second);
		}

		/*
		   RGBOutputType::Pointer rgb = RGBOutputType::New();
		   getRGB(images.image, andFilter->GetOutput(), p1,p2, rgb);

		   typedef itk::ImageFileWriter<RGBOutputType> WriterType;
		   WriterType::Pointer rgbWriter = WriterType::New();
		   rgbWriter->SetImageIO(itk::NrrdImageIO::New());
		   rgbWriter->SetFileName(ss.str()+".nrrd");
		   rgbWriter->SetInput(rgb);
		   rgbWriter->Update();

*/




		/*


*/

		++lineIt;

	}

	computePlane(allPoints, endNormal, endPoint);

}

// ------------------------------------------------------------------------
void computePlane(PointPairMap &inputPoints, VectorType &normal, VectorType &point)
{
	std::vector<PointType> allPoints;
	PointPairMap::iterator mapIt = inputPoints.begin();
	while(mapIt != inputPoints.end())
	{
		allPoints.push_back(mapIt->second.first);
		allPoints.push_back(mapIt->second.second);
		mapIt++;
	}

	computePlane(allPoints, normal, point);


}


// ------------------------------------------------------------------------
void computePlane(std::vector<PointType> &inputPoints, VectorType &normal, VectorType &point)
{
	MatrixType P = MatrixType(3, inputPoints.size());
	for(unsigned int i = 0; i < inputPoints.size(); i++)
	{
		P(0,i) = inputPoints[i][0];		
		P(1,i) = inputPoints[i][1];		
		P(2,i) = inputPoints[i][2];		
	}

	MatrixType centroid = P.rowwise().mean();

	// subtract the mean from the points
	P = P.colwise() - P.rowwise().mean();
	MatrixType C = P*P.transpose();

	// get the eigen vectors
	Eigen::SelfAdjointEigenSolver<MatrixType> solver(C);

	normal = solver.eigenvectors().col(0);
	point = centroid;
}


// -----------------------------------------------------------------------
void getRGB(const ImageType::Pointer &image, const LabelType::Pointer &mask,
		const RealImageType::Pointer &prob1, const RealImageType::Pointer &prob2,
		RGBOutputType::Pointer &output)
{
	
	typedef itk::RescaleIntensityImageFilter<ImageType, LabelType> RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetInput(image);
	rescaler->SetOutputMaximum(255);
	rescaler->SetOutputMinimum(0);
	rescaler->Update();

	output->SetDirection(image->GetDirection());
	output->SetSpacing(image->GetSpacing());
	output->SetOrigin(image->GetOrigin());
	output->SetRegions(image->GetLargestPossibleRegion());
	output->Allocate();

	itk::ImageRegionConstIterator<LabelType> imIt(rescaler->GetOutput(), rescaler->GetOutput()->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<LabelType> maskIt(mask, mask->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<RealImageType> pIt(prob1, prob1->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<RealImageType> pIt2(prob2, prob2->GetLargestPossibleRegion());
	itk::ImageRegionIterator<RGBOutputType> outIt(output, output->GetLargestPossibleRegion());


	while(!outIt.IsAtEnd())
	{
		unsigned char im = imIt.Get();
		unsigned char ma = maskIt.Get();
		double p1 = pIt.Get();
		double p2 = pIt2.Get();

		RGBPixelType pix;
		for(unsigned int i = 0; i < 3; i++)
		{
			pix[i] = im;
		}

		if(ma == 255)
		{
			//pix[0] += 40;
		}

		pix[1] += (unsigned char) (p1*60);
		pix[0] += (unsigned char) (p2*60);

		//if(p1 > 0.1 || p2 > 0.1) pix[0] -= 40;


		outIt.Set(pix);


		++imIt; ++maskIt; ++pIt; ++pIt2; ++outIt;
	}


}




