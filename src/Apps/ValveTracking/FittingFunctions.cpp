#include "FittingFunctions.h"
#include <FlipChecker.h>
#include <itkResampleImageFilter.h>
#include <ValveNormaliser.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkFlipImageFilter.h>
#include <itkPermuteAxesImageFilter.h>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataWriter.h>
#include <itkSmartPointer.h>
#include <vtkPlaneSource.h>
#include <ValveIO.h>
#include <CommonDefinitions.h>

#include <Directory.h>


// ------------------------------------------------------------------------
void createAlignedValve(const TestData::LineType &line, const ImageType::Pointer &image, 
		ValveLine<3>::Pointer &output, bool flipImage, bool flipPoints)
{
	if(!output) output = ValveLine<3>::New();

	ValveLine<3>::Pointer tmp = ValveLine<3>::New();
	tmp->SetP1(line.p1);
	tmp->SetP2(line.p2);
	tmp->SetImage(image);
	tmp->UpdateIndexs();


	ValveNormaliser::Pointer normaliser = ValveNormaliser::New();
	normaliser->SetInput(tmp);
	normaliser->SetFlip(flipImage);
	normaliser->SetFlipPoints(flipPoints);
	normaliser->Normalise();

	output = normaliser->GetOutput();

}

// ------------------------------------------------------------------------
void startingPlane(const std::vector<MatrixType> &pointsData, VectorType &point, VectorType &normal)
{
	// put into a big matrix 
	unsigned int count = 0;
	std::vector<MatrixType> flattened;
	for(unsigned int i = 0; i < pointsData.size(); i++)
	{
	
		std::cout << pointsData[i].cols() << std::endl;
		if(pointsData[i].cols() == 6)
		{
			MatrixType flat = MatrixType(1, 6*3);
			for(unsigned int j = 0; j < 3; j++)
			{
				for(unsigned int k = 0; k < 6; k++)
				{
					flat(0, k*3+j) = pointsData[i](k, j);
				}
				
			}
			flattened.push_back(flat);
			count++;
		}
		
	}



	MatrixType allPlanes = MatrixType(count, 6*3);
	for(unsigned int i = 0; i < count; i++)
	{
		allPlanes.row(i) = flattened[i].row(0);		
	}


	allPlanes = allPlanes.colwise() - allPlanes.rowwise().mean();
	
	// get the mean 
	MatrixType mean = allPlanes.colwise().mean();

	// reshape into P 
	MatrixType P = MatrixType(3,6);
	for(unsigned int i = 0; i < 3; i++)
	{
		for(unsigned int j = 0; j < 6; j++)
		{
			P(i,j) = mean(0,j*3+i);
		}
	}




	// subtract the mean from the points
	MatrixType centroid = P.rowwise().mean();
	P = P.colwise() - P.rowwise().mean();
	MatrixType C = P*P.transpose();

	// get the eigen vectors
	Eigen::SelfAdjointEigenSolver<MatrixType> solver(C);

	normal = solver.eigenvectors().col(0);
	point = centroid.transpose();
}

// ------------------------------------------------------------------------
int getBestLines(PatchParams &params, const std::string type,
		TestData::Pointer &testData, TestData::LineTypeList &lines, 
		ClassifierMap &classifiers, LengthData::Pointer &lengths)
{
	if(lines.size() == 0)
		return -1;
	if(lines.size() == 1)
		return 0;

	FlipChecker::Pointer checker = FlipChecker::New();
	bool flipImage = checker->FlipImage("MV-"+type, testData->GetId());
	bool flipPoints = checker->FlipPoints("MV-"+type, testData->GetId());

	if((type == "2C" || type == "3C") && flipImage)
	{
		flipPoints = !flipPoints;
	}	

	TestData::ImageGroup images;
	testData->GetImageGroup(type, images);

	double minValue = 1000;
	unsigned int minIndex = 0;
	for(unsigned int i = 0; i < lines.size(); i++)
	{
		// create a valve line from the data
		TestData::LineType line = lines[i];

		ValveLine<3>::Pointer alignedValve = ValveLine<3>::New();
		createAlignedValve(line, images.image, alignedValve, flipImage, flipPoints);

		double pointDist = alignedValve->GetP1().EuclideanDistanceTo(alignedValve->GetP2());
		double lengthProb = lengths->Predict("MV-" + type, params.timeStep, pointDist);
		double lengthCost = -log(std::max(lengthProb, 0.00000001));

		double lineSum = 0.0;
		for(unsigned int j = 0; j < 2; j++)
		{

			MatrixType feature;
			extractLBPFeature(params, alignedValve, alignedValve->GetIndex(j+1), feature);

			SVMClassifier::IntMatrixType classes;
			MatrixType probs;
			classifiers["MV-"+type][j]->PredictProbability(feature, classes, probs);

			double prob = probs(0,1);
			double cost = -log(std::max(prob, 0.000000001));

			lineSum += cost;

		}

		lineSum += lengthCost;

		if(lineSum < minValue)
		{
			minIndex = i;
			minValue = lineSum;
		}
	}


	return minIndex;
}

// ------------------------------------------------------------------------
void createLine(const TestData::LineType &input, ImageType::Pointer &image, ValveType::Pointer &output)
{
	output = ValveType::New();
	output->SetP1(input.p1);
	output->SetP2(input.p2);
	output->SetImage(image);
	output->UpdateIndexs();

}	

// ------------------------------------------------------------------------
double imageCost(PatchParams &params, TestData::Pointer &data, 
		VectorType &point, VectorType &normal, ClassifierMap &classifiers, 
		LengthData::Pointer &lengths, bool save)
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
		unsigned int minIndex = 0;

		if(lineIt->second.size() == 0) 
		{
			totalCost += minValue;
			lineIt++;
			continue;
		}

		const std::string type = lineIt->first;

		bool flipImage = checker->FlipImage("MV-" + lineIt->first, data->GetId());
		bool flipPoints = checker->FlipPoints("MV-" + lineIt->first, data->GetId());

		if((type == "2C" || type == "3C") && flipImage)
		{
			flipPoints = !flipPoints;
		}	

		TestData::ImageGroup images;
		data->GetImageGroup(type, images);

		for(unsigned int i = 0; i < lineIt->second.size(); i++)
		{
			// create a valve line from the data
			TestData::LineType line = lineIt->second[i];

			ValveLine<3>::Pointer alignedValve = ValveLine<3>::New();
			createAlignedValve(line, images.image, alignedValve, flipImage, flipPoints);

			double pointDist = alignedValve->GetP1().EuclideanDistanceTo(alignedValve->GetP2());
			double lengthProb = lengths->Predict("MV-" + type, params.timeStep, pointDist);
			double lengthCost = -log(std::max(lengthProb, 0.00000001));

			double lineSum = 0.0;
			for(unsigned int j = 0; j < 2; j++)
			{

				MatrixType feature;
				extractLBPFeature(params, alignedValve, alignedValve->GetPoint(j+1), feature);

				SVMClassifier::IntMatrixType classes;
				MatrixType probs;
				classifiers["MV-"+type][j]->PredictProbability(feature, classes, probs);

				double prob = probs(0,1);
				double cost = -log(std::max(prob, 0.000000001));

				lineSum += cost;

			}

			std::cout << "Costs: " << type << " " << lineSum << " " << lengthCost << std::endl;
			lineSum += lengthCost;

			if(lineSum < minValue)
			{
				minIndex = i;
				minValue = lineSum;
			}
		}

		// save the min index
		if(save)
		{
			TestData::LineType fline = lineIt->second[minIndex];
			ValveLine<3>::Pointer toSave = ValveLine<3>::New();
			createAlignedValve(fline, images.image, toSave, flipImage, flipPoints);
			ValveWriter<3>::Pointer writer = ValveWriter<3>::New();
			writer->SetFileName(type);
			writer->SetInput(toSave);
			writer->Write();
		}

		lineIt++;
		
		totalCost += minValue;

	}

	if(totalCost < 0.0) totalCost = 1000;

	return totalCost;

}



// ------------------------------------------------------------------------
void traingClassifiers(const unsigned int exclude, 
		const PatchParams &params, 
		const PatchTrainingData::Pointer &trainingData,
	   	ClassifierMap &classifiers)
{
	for(unsigned int i = 0; i < params.valveSubTypes.size(); i++)
	{
		ClassifierList classifierList;

		
		std::string type = params.valveSubTypes[i];		
		for(unsigned int j = 0; j < 2; j++)
		{
			MatrixType X;
			IntMatrixType y;
			std::cout << "Training: " << type << " - " << j+1 << std::endl;
			trainingData->GetTrainingData(exclude, type, params.timeStep, j+1, X, y);


			SVMClassifier::Pointer classifier = SVMClassifier::New();
			classifier->Train(X,y);
			classifierList.push_back(classifier);
		}

		classifiers[type] = classifierList;
	}
}

// ------------------------------------------------------------------------
void loadOrTrainClassifiers(const PatchParams & params, const unsigned int exclude,
		const PatchTrainingData::Pointer &trainingData, ClassifierMap &classifiers)
{
	std::stringstream ss;
	ss << exclude << "-" << params.valveSubTypes[0] << "-" << params.timeStep << "-0.cls";
	std::string testFilename = utils::Directory::GetPath(params.outputDirectory, ss.str());
	std::cout << ss.str() << std::endl;
	if(utils::Directory::FileExists(testFilename))
	{
		loadClassifiers(params, exclude, classifiers);
	}
	else
	{
		traingClassifiers(exclude, params, trainingData, classifiers);
		saveClassifiers(params, exclude, classifiers);
	}
}

// ------------------------------------------------------------------------
void saveClassifiers(const PatchParams &params, const unsigned int exclude, const ClassifierMap &classifiers)
{
	ClassifierMap::const_iterator mapIt = classifiers.begin();
	while(mapIt != classifiers.end())
	{
		for(unsigned int i = 0; i < 2; i++)
		{
			std::stringstream ss;
			ss <<  exclude << "-" << mapIt->first << "-" << params.timeStep << "-" << i <<  ".cls";
			std::string filename = utils::Directory::GetPath(params.outputDirectory, ss.str());
			mapIt->second[i]->Save(filename);
		}
		++mapIt;
	}
}

// ------------------------------------------------------------------------
void loadClassifiers(const PatchParams &params, const unsigned int exclude, ClassifierMap &classifiers)
{
	for(unsigned int i = 0; i < params.valveSubTypes.size(); i++)
	{
		ClassifierList classifierList;
		std::string type = params.valveSubTypes[i];		
		for(unsigned int j = 0; j < 2; j++)
		{
			std::stringstream ss;
			ss << exclude << "-" << type << "-" << params.timeStep << "-" << j << ".cls";
			std::string filename = utils::Directory::GetPath(params.outputDirectory, ss.str());
			SVMClassifier::Pointer classifier = SVMClassifier::New();
			classifier->Load(filename);
			classifierList.push_back(classifier);
		}
		classifiers[type] = classifierList;
	}
}




// ------------------------------------------------------------------------
void computeBoundingBox(const MatrixType &points, const TestData::TransformType::Pointer &transform, 
		vtkBoundingBox &boundingBox)
{

	for(unsigned int i = 0; i < points.rows(); i++)
	{

		PointType p;
		for(unsigned int j = 0; j < 3; j++)
		{
			p[j] = points(i,j);
		}


		boundingBox.AddPoint(p.GetDataPointer());		
	}
	boundingBox.Inflate(10);
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
double f2(const std::vector<double> &x, std::vector<double> &grad, void *f_data)
{
	OptData *data = (OptData*) f_data;

	// extract the points from x 
	for(unsigned int i = 0; i < 3; i++)
	{
		for(unsigned int j = 0; j < 3; j++)
		{
			
		}
	}


	VectorType point, normal;
	
	convert_x(x, point, normal);

	for(unsigned int i = 0; i < 5; i++)
	{
		std::cout << x[i] << " ";
	}
	std::cout << "\n";

	PlaneMeasurementType mv;
	for(unsigned int i = 0; i < 3; i++)
	{
		mv[i] = point(i);
		mv[i+3] = normal(i);
	}

	PlaneMembershipType::Pointer planeMember = PlaneMembershipType::New();
	planeMember->SetMean(data->planeMember->GetMean());
	planeMember->SetCovariance(data->planeMember->GetCovarianceMatrix());

	double planeProb = planeMember->Evaluate(mv);
	double planeCost = -log(std::max(planeProb, 0.00000001));


	grad.clear();


	double val = imageCost(data->params, data->testData, point, normal, data->classifiers, data->lengths, false);
	double total = val;
	std::cout << "Plane: " << planeCost << std::endl;
	std::cout << "val: " << val << std::endl;
	std::cout << "Total: " << total << std::endl;
	std::cout << "---------------------------------" << std::endl;
	return total;

}

// ------------------------------------------------------------------------
void computeProbImage(const PatchParams &params, unsigned int pnum,
	   	unsigned int id, std::string type, ClassifierMap &classifiers, 
		const ValveType::Pointer &valve, const LabelType::Pointer &mask, RealImageType::Pointer &output)
{
	output = RealImageType::New();
	output->SetDirection(mask->GetDirection());
	output->SetSpacing(mask->GetSpacing());
	output->SetOrigin(mask->GetOrigin());
	output->SetRegions(mask->GetLargestPossibleRegion());
	output->Allocate();
	output->FillBuffer(0);

	itk::ImageRegionIterator<LabelType> maskIt(mask, mask->GetLargestPossibleRegion());
	itk::ImageRegionIterator<RealImageType> probIt(output, output->GetLargestPossibleRegion());
	
	while(!maskIt.IsAtEnd())
	{
		if(maskIt.Get() == 255)
		{
			IndexType index = maskIt.GetIndex();
			PointType point;
			output->TransformIndexToPhysicalPoint(index, point);

			MatrixType feature;
			extractLBPFeature(params, valve, point, feature);

			MatrixType probs;
			IntMatrixType classes;
			classifiers["MV-"+type][pnum]->PredictProbability(feature, classes, probs);
			
			probIt.Set(probs(0,1));

		}

		
		++maskIt; ++probIt;
	}



}


// ------------------------------------------------------------------------
double f(const std::vector<double> &x, std::vector<double> &grad, void *f_data)
{
	OptData *data = (OptData*) f_data;

	VectorType point, normal;
	
	convert_x(x, point, normal);

	for(unsigned int i = 0; i < 5; i++)
	{
		std::cout << x[i] << " ";
	}
	std::cout << "\n";

	PlaneMeasurementType mv;
	for(unsigned int i = 0; i < 3; i++)
	{
		mv[i] = point(i);
		mv[i+3] = normal(i);
	}

	PlaneMembershipType::Pointer planeMember = PlaneMembershipType::New();
	planeMember->SetMean(data->planeMember->GetMean());
	planeMember->SetCovariance(data->planeMember->GetCovarianceMatrix());

	double planeProb = planeMember->Evaluate(mv);
	double planeCost = -log(std::max(planeProb, 0.00000001));


	grad.clear();


	double val = imageCost(data->params, data->testData, point, normal, data->classifiers, data->lengths, false);
	double total = val;
	std::cout << "Plane: " << planeCost << std::endl;
	std::cout << "val: " << val << std::endl;
	std::cout << "Total: " << total << std::endl;
	std::cout << "---------------------------------" << std::endl;
	return total;

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
void FlipImage(const LabelType::Pointer &input, LabelType::Pointer &output)
{
	if(!output) output = LabelType::New();

	ImageType::DirectionType outputDirection = input->GetDirection();

	typedef itk::PermuteAxesImageFilter<LabelType> FlipperType;
	FlipperType::PermuteOrderArrayType axes;
	axes[0] = 1;
	axes[1] = 0;
	axes[2] = 2;

	FlipperType::Pointer flipper = FlipperType::New();
	flipper->SetInput(input);
	flipper->SetOrder(axes);
	flipper->Update();

	output = flipper->GetOutput();
	output->SetDirection(outputDirection);

}


// ------------------------------------------------------------------------
void FlipImage(const RealImageType::Pointer &input, RealImageType::Pointer &output)
{
	if(!output) output = RealImageType::New();

	ImageType::DirectionType outputDirection = input->GetDirection();

	typedef itk::PermuteAxesImageFilter<RealImageType> FlipperType;
	FlipperType::PermuteOrderArrayType axes;
	axes[0] = 1;
	axes[1] = 0;
	axes[2] = 2;

	FlipperType::Pointer flipper = FlipperType::New();
	flipper->SetInput(input);
	flipper->SetOrder(axes);
	flipper->Update();

	output = flipper->GetOutput();
	output->SetDirection(outputDirection);

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


