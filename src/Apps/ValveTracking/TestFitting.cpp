#include <iostream>

#include <QFileInfo>
#include <QDir>
#include <FlipChecker.h>
#include <ValveIO.h>
#include <LBPFeatureExtractor.h>

#include <itkSimilarity3DTransform.h>
#include <itkResampleImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkBinaryContourImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <vtkLineSource.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <ValveNormaliser.h>
#include <PatchExtractor.h>
#include <BinaryPatchFeatureExtractor.h>
#include <itkGaussianMembershipFunction.h>
#include <itkDistanceToCentroidMembershipFunction.h>
#include <itkGaussianMembershipFunction.h>
#include <itkMatchCardinalityImageToImageMetric.h>


#include <CMRFileExtractor.h>
#include <ValvePlane.h>
#include <TrainingData.h>
#include <ValveOriginFinder.h>
#include <CommonDefinitions.h>
#include <SimpleMRFSegmenter.h>
#include <LineIntersectionFinder.h>

#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkPlaneSource.h>
#include <vtkSmartPointer.h>
#include <vtkBoundingBox.h>

#include <itkListSample.h>
#include <itkVector.h>
#include <itkCovarianceSampleFilter.h>

#include <nlopt.hpp>

#include <TestData.h>
#include <SVMClassifier.h>

#define DATA_PATH "/home/om0000/ValveTracking/"

using namespace vt;



typedef TestData::LabelType LabelType;
typedef TestData::ImageType ImageType;
typedef std::vector<TrainingData::Pointer> TrainingPointList;
typedef std::map<std::string, TrainingPointList> TrainingDataMap;
typedef utils::DoubleMatrixType MatrixType;
typedef itk::Vector<double, 6> PlaneMeasurementType;
typedef itk::Statistics::ListSample<PlaneMeasurementType> PlaneSampleType;
typedef itk::Vector<double, 256> PatchMeasurementType;
typedef itk::Statistics::ListSample<PatchMeasurementType> PatchSampleType;
typedef itk::Statistics::CovarianceSampleFilter<PlaneSampleType> PlaneCovarianceFilterType;
typedef itk::Statistics::CovarianceSampleFilter<PatchSampleType> PatchCovarianceFilterType;

typedef std::vector<PatchCovarianceFilterType::Pointer> PatchCovarianceFilterList;
typedef std::map<std::string, PatchCovarianceFilterList> PatchCovarianceMapType;

typedef std::vector<SVMClassifier::Pointer> ClassifierList;
typedef std::map<std::string, ClassifierList> ClassifierMap;

typedef LineIntersectionFinder::OutputLineType LineType;
typedef std::vector<LineType> LineTypeList;



typedef TestData::VectorType VectorType;
void convert_plane(const VectorType &point, const VectorType &plane, std::vector<double> &x);
void savePlane(const MatrixType &plane, std::string filename);
void savePlane(const VectorType &point, const VectorType &normal, std::string filename);
void saveLine(LineType &line, std::string filename);
void meanPlane(const ValveTrainingData::Pointer &data, const unsigned int exclude, MatrixType &plane);
void computeBoundingBox(const MatrixType &points, vtkBoundingBox &boundingBox);
void loadTrainingData(const unsigned int exclude, PatchCovarianceMapType &map, ValveTrainingData::Pointer &valveData);
void loadTrainingData2(const unsigned int exclude, ClassifierMap &map, ValveTrainingData::Pointer &valveData);
void createClassifier(const TrainingData::Pointer &training, const unsigned int exclude, SVMClassifier::Pointer &cls);
void convert_x(const std::vector<double> &x, VectorType &point, VectorType &normal);
void computePlaneCovariance(const MatrixType &data, PlaneCovarianceFilterType::Pointer &cov);
void computePatchCovariance(const MatrixType &data, PatchCovarianceFilterType::Pointer &cov);

void computeSegmentation(const ValveLine<3>::Pointer &input, LabelType::Pointer &segmentation, 
		unsigned int segmentationPatchSize, unsigned int pointNumber);
typedef ValveLine<3>::ContIndexType ContIndexType;
void extractFeature(const LabelType::Pointer &mask, const ContIndexType &location, 
		unsigned int featureSize, PatchMeasurementType & feature); 
void extractLBPFeature(const ImageType::Pointer &image,	const ContIndexType &location, 
		unsigned int featureSize, PatchMeasurementType & feature);

double imageCost(TestData::Pointer &data, VectorType &plane, VectorType &normal, PatchCovarianceMapType &patchCovariances);
double imageCost2(TestData::Pointer &data, VectorType &plane, VectorType &normal, ClassifierMap &classifiers);


typedef struct opt_data_
{
	TestData::Pointer data;
	VectorType normal;
	VectorType point;
	PatchCovarianceMapType patchMap;
	PlaneCovarianceFilterType::Pointer planeCov;
	ClassifierMap clsMap;
	VectorType originalNormal;
	VectorType originalPoint;
	int itCount;

} OptData;



double f(const std::vector<double> &x, std::vector<double> &grad, void *f_data);

int main(int argc, char ** argv)
{
	// load the input data
	const std::string inputDataFolder = argv[1];
	TestData::Pointer testData = TestData::Initialise(inputDataFolder);
	const unsigned int segmentationPatchSize = 40;
	const unsigned int featurePatchSize = 10;
	

	// load the training data
	PatchCovarianceMapType patchCovariances;
	ValveTrainingData::Pointer valveData;
	//loadTrainingData(testData->GetId(), patchCovariances, valveData);
	ClassifierMap clsMap;
	loadTrainingData2(testData->GetId(), clsMap, valveData);



	// compute the bounding box from the training points
	const unsigned int timeStep = 0;
	MatrixType trainingPoints;
	valveData->GetTrainingPoints(testData->GetId(), timeStep, trainingPoints);
	


	// find and set the bounding box and initialise the test data
	vtkBoundingBox boundingBox;
	computeBoundingBox(trainingPoints, boundingBox);
	testData->Prepare(timeStep, boundingBox);


	// compute the distribution for the plane and get the mean for initialisation
	MatrixType planeData;
	valveData->GetTrainingData(testData->GetId(), timeStep, planeData);
	PlaneCovarianceFilterType::Pointer planeDistribution = PlaneCovarianceFilterType::New();
	computePlaneCovariance(planeData, planeDistribution);

	TestData::VectorType point, normal;
	for(unsigned int i = 0; i < 3; i++)
	{
		point(i) = planeDistribution->GetMean()[i];
		normal(i) = planeDistribution->GetMean()[i+3];
	}


	OptData * optData = new OptData;
	optData->data = testData;
	//optData->patchMap = patchCovariancpatces;
	optData->clsMap = clsMap;
	optData->planeCov = planeDistribution;
	optData->originalNormal = normal;
	optData->originalPoint = point;
	optData->itCount = 0;
	//optData->checker = FlipChecker::New();


	/*
	int minIndex = 0;
	double minVal = std::numeric_limits<double>::max();
	for(int i = -50; i < 50; i++)
	{
		VectorType newPoint = point + (normal*i);
		double cost = imageCost2(testData, newPoint, normal, clsMap);
		if(cost < minVal)
		{
			minVal = cost;
			minIndex = i;
		}
	}

	VectorType finalPoint = point + (minIndex*normal);
	savePlane(point, normal, "initialPlane.vtk");
	savePlane(finalPoint, normal, "finalPlane.vtk");
	testData->SaveImageGroup("2C", "2C");
	testData->SaveImageGroup("3C", "3C");
	testData->SaveImageGroup("4C", "4C");

	*/
		
	//double dCost = imageCost(testData, point, normal, patchCovariances);
	//std::cout << dCost << std::endl;


	nlopt::opt opt(nlopt::LN_COBYLA, 5);
	opt.set_min_objective(f, (void*) optData);
	opt.set_xtol_rel(1e-4);
	opt.set_maxeval(100);

	
	std::vector<double> initStep(5);
	double pointStep = atof(argv[2]);
	double normalStep = atof(argv[3]);
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

	testData->SaveImageGroup("2C", "2C");
	testData->SaveImageGroup("3C", "3C");
	testData->SaveImageGroup("4C", "4C");
	savePlane(initialPlane, "initialPlane.vtk");
	savePlane(finalPlane, "finalPlane.vtk");


	return 0;
}

// ------------------------------------------------------------------------
void convert_x(const std::vector<double> &x, VectorType &point, VectorType &normal)
{
	for(unsigned int i = 0; i < 3; i++)
	{
		point(i) = x[i];
	}

	double theta = x[3];
	double rho = x[4];
	double r = 1.0;

	normal(0) = r * std::sin(theta) * std::cos(rho);
	normal(1) = r * std::sin(theta) * std::sin(rho);
	normal(2) = r * std::cos(theta);

}

// ------------------------------------------------------------------------
void convert_plane(const VectorType &point, const VectorType &normal, std::vector<double> &x)
{
	// convert the normal into spherical coordinates
	double theta = std::acos( normal(2) / normal.norm() );
	double rho = std::atan2(normal(1), normal(0));
	for(unsigned int i = 0; i < 3; i++)
	{
		x.push_back(point(i));
	}

	x.push_back(theta);
	x.push_back(rho);
}

// ------------------------------------------------------------------------
double f(const std::vector<double> &x, std::vector<double> &grad, void *f_data)
{
	OptData *data = (OptData*) f_data;
	
	convert_x(x, data->point, data->normal);

	PlaneMeasurementType mv;
	for(unsigned int i = 0; i < 3; i++)
	{
		mv[i] = data->point(i);
		mv[i+3] = data->normal(i);
	}

	for(unsigned int i = 0; i < 5; i++)
	{
		std::cout << x[i] << " ";
	}
	std::cout << "\n";


	std::stringstream ss;
	ss << data->itCount << ".vtk";
	data->itCount++;

	MatrixType op(1,6);
	for(unsigned int i = 0; i < 3; i++)
	{
		op(0,i) = data->point(i);
		op(0,i+3) = data->normal(i);
	}
	//savePlane(op, ss.str());

	grad.clear();

	typedef itk::Statistics::GaussianMembershipFunction<PlaneMeasurementType> MembershipType;
	MembershipType::Pointer member = MembershipType::New();
	member->SetMean(data->planeCov->GetMean());
	member->SetCovariance(data->planeCov->GetCovarianceMatrix());





	double val = imageCost2(data->data, data->point, data->normal, data->clsMap);
	double planeVal = -log(member->Evaluate(mv));
	double total = val;
	std::cout << "Index: " << data->itCount << std::endl;
	std::cout << "val: " << val << std::endl;
	std::cout << "Plane: " << planeVal << std::endl;
	std::cout << "Total: " << total << std::endl;
	std::cout << "---------------------------------" << std::endl;
	return total;

}


// ------------------------------------------------------------------------
double imageCost2(TestData::Pointer &data, VectorType &point, VectorType &normal, ClassifierMap &classifiers)
{
	typedef TestData::LineGroup LineGroup;
	LineGroup testLines;
	data->ExtractTestLines(normal, point, testLines);

	FlipChecker::Pointer checker = FlipChecker::New();

	LineGroup::iterator lineIt = testLines.begin();
	double totalCost = 0;
	while(lineIt != testLines.end())
	{
		double minValue = 1000;

		if(lineIt->second.size() == 0) {
			totalCost += minValue;
			lineIt++;
			continue;
		}

		for(unsigned int i = 0; i < lineIt->second.size(); i++)
		{
			// create a valve line from the data
			ValveLine<3>::Pointer vline = ValveLine<3>::New();
			vline->SetP1(lineIt->second[i].p1);
			vline->SetP2(lineIt->second[i].p2);
		

			TestData::ImageGroup images;
			data->GetImageGroup(lineIt->first, images);
			vline->SetImage(images.image);
			vline->UpdateIndexs();

			ValveNormaliser::Pointer normaliser = ValveNormaliser::New();
			normaliser->SetInput(vline);

			
			bool flipImage = checker->FlipImage("MV-" + lineIt->first, data->GetId());
			bool flipPoints = checker->FlipPoints("MV-" + lineIt->first, data->GetId());

			std::cout << "Flip: " << flipImage << " " << flipPoints << std::endl;
			normaliser->SetFlip(flipImage);
			normaliser->SetFlipPoints(flipPoints);
			normaliser->Normalise();
			ValveLine<3>::Pointer alignedValve = normaliser->GetOutput();

			ValveNormaliser::TransformType::Pointer transform = normaliser->GetTransform();
			typedef itk::ResampleImageFilter<LabelType, LabelType> ResamplerType;
			ResamplerType::Pointer resampler = ResamplerType::New();
			resampler->SetInput(images.seg);
			resampler->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<LabelType, double>::New());
			resampler->SetTransform(transform);
			resampler->SetOutputParametersFromImage(images.seg);
			resampler->Update();

			LabelType::Pointer segmentation = resampler->GetOutput();


			if(alignedValve->GetP1().EuclideanDistanceTo(alignedValve->GetP2()) < 0) continue;
			
			double lineSum = 0.0;

			ValveWriter<3>::Pointer writer = ValveWriter<3>::New();
			writer->SetInput(vline);
			writer->SetFileName(lineIt->first);
			writer->Write();

			savePlane(point, normal, "plane.vtk");

			for(unsigned int j = 0; j < 2; j++)
			{

				LabelType::Pointer segmentation1 = LabelType::New();
				//computeSegmentation(alignedValve, segmentation1, 40, j+1);
				



				PatchMeasurementType feature, featureTest;
				extractLBPFeature(alignedValve->GetImage(), alignedValve->GetIndex(j+1), 20, feature);
				//extractFeature(segmentation1, alignedValve->GetIndex(j+1), 10, featureTest);

				MatrixType test = MatrixType(1, feature.GetVectorDimension());
				for(unsigned int k = 0; k < feature.GetVectorDimension(); k++)
				{
					test(0,k) = feature[k];					
				}





				SVMClassifier::IntMatrixType classes;
				MatrixType probs;
				classifiers[lineIt->first][j]->PredictProbability(test, classes, probs);

				double prob = probs(0,1);
				double cost = -log(std::max(prob, 0.000000001));
				std::cout << j << " " << prob << " " << cost << std::endl;


				lineSum += cost;

			}

			std::cout << "Line Sum: " << lineSum << std::endl;

			if(lineSum < minValue)
				minValue = lineSum;
		}

		std::cout << "Min Value (" << lineIt->first  << "): " << minValue << std::endl;



		lineIt++;
		
		totalCost += minValue;

	}

	if(totalCost < 0.0) totalCost = 1000;

	return totalCost;

}




// ------------------------------------------------------------------------
void createClassifier(const TrainingData::Pointer &training, const unsigned int exclude, SVMClassifier::Pointer &cls)
{
	MatrixType X;
	SVMClassifier::IntMatrixType y;
	training->GetTrainingData(exclude, X,y);

	

	cls->Train(X,y);
}

void loadTrainingData2(const unsigned int exclude, ClassifierMap &map, ValveTrainingData::Pointer &valveData)
{
	TrainingData::Pointer ap1 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/LBP-MV-2C-P1.hdf5");
	SVMClassifier::Pointer ap1c = SVMClassifier::New();
	createClassifier(ap1, exclude, ap1c);
	
	TrainingData::Pointer ap2 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/LBP-MV-2C-P2.hdf5");
	SVMClassifier::Pointer ap2c = SVMClassifier::New();
	createClassifier(ap2, exclude, ap2c);

	ClassifierList l1;
	l1.push_back(ap1c);
	l1.push_back(ap2c);

	TrainingData::Pointer ap3 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/LBP-MV-3C-P1.hdf5");
	SVMClassifier::Pointer ap3c = SVMClassifier::New();
	createClassifier(ap3, exclude, ap3c);

	TrainingData::Pointer ap4 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/LBP-MV-3C-P2.hdf5");
	SVMClassifier::Pointer ap4c = SVMClassifier::New();
	createClassifier(ap4, exclude, ap4c);

	ClassifierList l2;
	l2.push_back(ap3c);
	l2.push_back(ap4c);

	TrainingData::Pointer ap5 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/LBP-MV-4C-P1.hdf5");
	SVMClassifier::Pointer ap5c = SVMClassifier::New();
	createClassifier(ap5, exclude, ap5c);

	TrainingData::Pointer ap6 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/LBP-MV-4C-P2.hdf5");
	SVMClassifier::Pointer ap6c = SVMClassifier::New();
	createClassifier(ap6, exclude, ap6c);

	ClassifierList l3;
	l3.push_back(ap5c);
	l3.push_back(ap6c);


	map["2C"] = l1;
	map["3C"] = l2;
	map["4C"] = l3;



	valveData = ValveTrainingData::Load("/home/om0000/ValveTracking/TrainingData/MV-Planes.hdf5");
	valveData->SetNumTimeSteps(25);

}

// ------------------------------------------------------------------------
void computeBoundingBox(const MatrixType &points, vtkBoundingBox &boundingBox)
{
	for(unsigned int i = 0; i < points.rows(); i++)
	{
		boundingBox.AddPoint(points(i,0), points(i,1), points(i,2));		
	}
	boundingBox.Inflate(20);
}



// ------------------------------------------------------------------------
void loadTrainingData(const unsigned int exclude, PatchCovarianceMapType &trainingData, ValveTrainingData::Pointer &valves)
{
	TrainingData::Pointer ap1 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/I-MV-2C-P1.hdf5");
	MatrixType ap1d;
	ap1->GetPositiveTrainingData(exclude, ap1d);
	PatchCovarianceFilterType::Pointer covap1 = PatchCovarianceFilterType::New();
	computePatchCovariance(ap1d, covap1);

	TrainingData::Pointer ap2 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/I-MV-2C-P2.hdf5");
	MatrixType ap2d;
	ap2->GetPositiveTrainingData(exclude, ap2d);
	PatchCovarianceFilterType::Pointer covap2 = PatchCovarianceFilterType::New();
	computePatchCovariance(ap2d, covap2);

	PatchCovarianceFilterList aList;
	aList.push_back(covap1);
	aList.push_back(covap2);
	trainingData["2C"] = aList;


	TrainingData::Pointer bp1 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/I-MV-3C-P1.hdf5");
	MatrixType bp1d;
	bp1->GetPositiveTrainingData(exclude, bp1d);
	PatchCovarianceFilterType::Pointer covbp1 = PatchCovarianceFilterType::New();
	computePatchCovariance(bp1d, covbp1);


	TrainingData::Pointer bp2 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/I-MV-3C-P2.hdf5");
	MatrixType bp2d;
	bp2->GetPositiveTrainingData(exclude, bp2d);
	PatchCovarianceFilterType::Pointer covbp2 = PatchCovarianceFilterType::New();
	computePatchCovariance(bp2d, covbp2);


	PatchCovarianceFilterList bList;
	bList.push_back(covbp1);
	bList.push_back(covbp2);
	trainingData["3C"] = bList;


	TrainingData::Pointer cp1 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/I-MV-4C-P1.hdf5");
	MatrixType cp1d;
	cp1->GetPositiveTrainingData(exclude, cp1d);
	PatchCovarianceFilterType::Pointer covcp1 = PatchCovarianceFilterType::New();
	computePatchCovariance(cp1d, covcp1);


	TrainingData::Pointer cp2 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/I-MV-4C-P2.hdf5");
	MatrixType cp2d;
	cp2->GetPositiveTrainingData(exclude, cp2d);
	PatchCovarianceFilterType::Pointer covcp2 = PatchCovarianceFilterType::New();
	computePatchCovariance(cp2d, covcp2);


	PatchCovarianceFilterList cList;
	cList.push_back(covcp1);
	cList.push_back(covcp2);
	trainingData["4C"] = cList;



	valves = ValveTrainingData::Load("/home/om0000/ValveTracking/TrainingData/MV-Planes.hdf5");
	valves->SetNumTimeSteps(25);


}

// ------------------------------------------------------------------------
void meanPlane(const ValveTrainingData::Pointer &data, const unsigned int exclude, MatrixType &plane)
{
	MatrixType allPlanes;
	data->GetTrainingData(exclude, 0, allPlanes);
	plane = allPlanes.colwise().mean();

}

// ------------------------------------------------------------------------
void savePlane(const VectorType &point, const VectorType &normal, std::string filename)
{
	MatrixType plane(1,6);
	for(unsigned int i = 0; i < 3; i++)
	{
		plane(0,i) = point(i);
		plane(0,i+3) = normal(i);
	}

	savePlane(plane, filename);

}


// ------------------------------------------------------------------------
void savePlane(const MatrixType &plane, std::string filename)
{
	vtkSmartPointer<vtkPlaneSource> source = vtkSmartPointer<vtkPlaneSource>::New();
	source->SetCenter(plane(0,0), plane(0,1), plane(0,2));
	source->SetNormal(plane(0,3), plane(0,4), plane(0,5));
	source->SetResolution(1,1);
	source->Update();

	vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
	poly = source->GetOutput();

	vtkSmartPointer<vtkPoints> p = poly->GetPoints();
	double cent[3];
	cent[0] = 0.0; cent[1] = 0.0; cent[2] = 0.0;

	for(unsigned int j = 0; j < p->GetNumberOfPoints(); j++)
	{
		cent[0] += p->GetPoint(j)[0] / (double) p->GetNumberOfPoints();
		cent[1] += p->GetPoint(j)[1] / (double) p->GetNumberOfPoints();
		cent[2] += p->GetPoint(j)[2] / (double) p->GetNumberOfPoints();
	}

	for(unsigned int j = 0; j < p->GetNumberOfPoints(); j++)
	{
		double pin[3];
		p->GetPoint(j,pin);		

		for(unsigned int k = 0; k < 3; k++)
		{
			pin[k] -= cent[k];
			pin[k] *= 100;
			pin[k] += cent[k];
		}

		p->SetPoint(j,pin);

	}


	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(poly);
	writer->Write();

}



// ------------------------------------------------------------------------
void saveLine(LineType &line, std::string filename)
{
	vtkSmartPointer<vtkLineSource> source = vtkSmartPointer<vtkLineSource>::New();
	source->SetPoint1(line.p1.GetDataPointer());
	source->SetPoint2(line.p2.GetDataPointer());
	source->Update();

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetInputData(source->GetOutput());
	writer->SetFileName(filename.c_str());
	writer->Write();

}

// ------------------------------------------------------------------------
void computePatchCovariance(const MatrixType &data, PatchCovarianceFilterType::Pointer &cov)
{
	PatchSampleType::Pointer sample = PatchSampleType::New();
	for(unsigned int i = 0; i < data.rows(); i++)
	{
		PatchMeasurementType mv;
		for(unsigned int j = 0; j < data.cols(); j++)
		{
			mv[j] = data(i,j);
		}
		
		sample->PushBack(mv);
	}

	cov = PatchCovarianceFilterType::New();
	cov->SetInput(sample);
	cov->Update();

}


// ------------------------------------------------------------------------
void computePlaneCovariance(const MatrixType &data, PlaneCovarianceFilterType::Pointer &cov)
{
	PlaneSampleType::Pointer sample = PlaneSampleType::New();
	for(unsigned int i = 0; i < data.rows(); i++)
	{
		PlaneMeasurementType mv;
		for(unsigned int j = 0; j < data.cols(); j++)
		{
			mv[j] = data(i,j);
		}
		
		sample->PushBack(mv);
	}

	cov = PlaneCovarianceFilterType::New();
	cov->SetInput(sample);
	cov->Update();
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
void extractImageFeature(const ImageType::Pointer &mask, const ContIndexType &location, 
		unsigned int featureSize, PatchMeasurementType & feature)
{
	ImagePatchExtractor::Pointer extractor = ImagePatchExtractor::New();
	extractor->SetImage(mask);
	extractor->SetPatchCenter(location);

	ImagePatchExtractor::SizeType size;
	size.Fill(featureSize);
	extractor->SetPatchSize(size);
	ImageType::Pointer patch = extractor->ExtractPatch();


	itk::ImageRegionIterator<ImageType> it(patch, patch->GetLargestPossibleRegion());
	int count = 0;
	while(!it.IsAtEnd())
	{
		feature[count] = it.Get();	
		++it; ++count;
	}

	feature.Normalize();
}

// ------------------------------------------------------------------------
void extractFeature(const LabelType::Pointer &mask, const ContIndexType &location, 
		unsigned int featureSize, PatchMeasurementType & feature)
{

	MaskPatchExtractor::Pointer extractor = MaskPatchExtractor::New();
	extractor->SetImage(mask);
	extractor->SetPatchCenter(location);

	ImagePatchExtractor::SizeType size;
	size.Fill(featureSize);
	extractor->SetPatchSize(size);
	LabelType::Pointer patch = extractor->ExtractPatch();


	itk::ImageRegionIterator<LabelType> it(patch, patch->GetLargestPossibleRegion());
	int count = 0;
	while(!it.IsAtEnd())
	{
		feature[count] = it.Get();	
		++it; ++count;
	}

	//feature.Normalize();
}


// ------------------------------------------------------------------------
void extractLBPFeature(const ImageType::Pointer &image,	const ContIndexType &location, 
		unsigned int featureSize, PatchMeasurementType & feature)
{
	ImagePatchExtractor::Pointer extractor = ImagePatchExtractor::New();
	extractor->SetImage(image);
	extractor->SetPatchCenter(location);

	ImagePatchExtractor::SizeType size;
	size.Fill(featureSize);
	extractor->SetPatchSize(size);
	ImageType::Pointer patch = extractor->ExtractPatch();


	LBPFeatureExtractor::Pointer featureExtractor = LBPFeatureExtractor::New();
	featureExtractor->SetInput(patch);
	MatrixType f;
	featureExtractor->Extract(f);

	for(unsigned int i = 0; i < f.cols(); i++)
	{
		feature[i] = f(0,i);
	}

}
