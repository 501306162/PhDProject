#include <iostream>
#include <ValvePlane.h>
#include <Directory.h>
#include <QString>
#include <QDir>
#include <CMRFileExtractor.h>
#include <ValveOriginFinder.h>
#include <itkCenteredAffineTransform.h>
#include <itkResampleImageFilter.h>
#include <CommonDefinitions.h>
#include <itkSimilarity3DTransform.h>
#include <itkCovarianceSampleFilter.h>
#include <itkListSample.h>
#include <FlipChecker.h>
#include <MatrixCommon.h>
#include <MatrixWriter.h>

typedef utils::Directory::FilenamesType FilenamesType;


using namespace vt;
int main(int argc, char ** argv)
{
	const std::string inputDirectory = argv[1];
	const std::string dataDirectory = argv[2];


	ValvePlane::VectorType meanNormal;
	meanNormal.Fill(0);
	ValvePlane::PointType meanPoint;
	meanPoint.Fill(0);


	std::vector<int> ids;
	

	typedef itk::Vector<double, 6> MeasurementType;
	typedef itk::Statistics::ListSample<MeasurementType> SampleType;
	typedef itk::Statistics::CovarianceSampleFilter<SampleType> CovarianceFilterType;

	FilenamesType planeFilenames = utils::Directory::GetFiles(inputDirectory, ".txt");
	const unsigned int numberOfTimeSteps = 25;
	const unsigned int numberOfPlanes = planeFilenames.size();


	typedef utils::DoubleMatrixType MatrixType;
	typedef utils::IntMatrixType IntMatrixType;

	MatrixType planes = MatrixType::Zero(numberOfTimeSteps*numberOfPlanes, 6);
	IntMatrixType timeSteps = IntMatrixType::Zero(numberOfTimeSteps*numberOfPlanes, 1);
	IntMatrixType owners = IntMatrixType::Zero(numberOfTimeSteps*numberOfPlanes,1);
	
	typedef CMRFileExtractor::ImageType ImageType;
	typedef ImageType::PointType PointType;

	std::vector<PointType> allPoints;
	std::vector<unsigned int> pointOwners;
	std::vector<unsigned int> pointTimeSteps;

	// create the list of samples
	std::vector<SampleType::Pointer> samplesList;
	for(unsigned int i = 0; i < numberOfTimeSteps; i++) samplesList.push_back(SampleType::New());

	

	FlipChecker::Pointer checker = FlipChecker::New();
	for(unsigned int i = 0; i < planeFilenames.size(); i++)
	{

		// load the valve plane
		ValvePlaneSequence::Pointer sequence = ValvePlaneSequence::Load(planeFilenames[i]);
		bool needsFlip = checker->FlipImageFromFileName("MV-2C", planeFilenames[i]);

		// get the data directory
		std::string dataFolder = utils::Directory::GetPath(dataDirectory, sequence->GetName());

		// load up the data 
		CMRFileExtractor::Pointer extractor = CMRFileExtractor::New();
		extractor->SetFolderName(dataFolder);
		extractor->SetFlip(needsFlip);
		extractor->Extract();

		ValveOriginFinder::Pointer originFinder = ValveOriginFinder::New();
		originFinder->Set2CImage(extractor->Get2CImage(0));
		originFinder->Set3CImage(extractor->Get3CImage(0));
		originFinder->SetImageStack(extractor->GetStackImage(0));
		originFinder->Compute();

		typedef itk::Similarity3DTransform<double> TransformType;
		TransformType::Pointer transform = TransformType::New();
		transform->SetTranslation(originFinder->GetTranslation());
		transform->SetMatrix(originFinder->GetRotation());

		std::cout << sequence->GetName() << " " << sequence->GetNumberOfPlanes() << std::endl;
		
		for(unsigned int ts = 0; ts < numberOfTimeSteps; ts++)
		{
			ValvePlane::Pointer valve = sequence->GetValvePlane(ts);
			ValvePlane::VectorType normal = transform->GetInverseTransform()->TransformCovariantVector(valve->GetNormal());
			ValvePlane::PointType center = transform->GetInverseTransform()->TransformPoint(valve->GetCenter());


			ValvePlane::PointListType allPlanePoints = valve->GetAllPoints();
			for(unsigned int j = 0; j < allPlanePoints.size(); j++)
			{
				PointType newPoint = transform->GetInverseTransform()->TransformPoint(allPlanePoints[j]);
				allPoints.push_back(newPoint);
				pointOwners.push_back(QString::fromStdString(sequence->GetName()).replace("d","").toInt());
				pointTimeSteps.push_back(ts);
			}

			if(needsFlip)
			{
				normal = -normal;
			}

			for(unsigned int j = 0; j < 3; j++)
			{
				meanNormal[j] += (normal[j] / (double) planeFilenames.size());
				meanPoint[j] += (center[j] / (double) planeFilenames.size());
			}

			unsigned int matIndex = i*numberOfTimeSteps+ts;
			planes(matIndex,0) = center[0];
			planes(matIndex,1) = center[1];
			planes(matIndex,2) = center[2];
			planes(matIndex,3) = normal[0];
			planes(matIndex,4) = normal[1];
			planes(matIndex,5) = normal[2];

			timeSteps(matIndex,0) = ts;
			owners(matIndex,0) = QString::fromStdString(sequence->GetName()).replace("d","").toInt();
		}
	}

	// build the point matrixes
	MatrixType pointMatrix = MatrixType(allPoints.size(), 3);
	IntMatrixType pointOwnerMatrix = IntMatrixType(allPoints.size(), 1);
	IntMatrixType pointTimesMatrix = IntMatrixType(allPoints.size(), 1);
	for(unsigned int i = 0; i < allPoints.size(); i++)
	{
		pointMatrix(i,0) = allPoints[i][0];		
		pointMatrix(i,1) = allPoints[i][1];		
		pointMatrix(i,2) = allPoints[i][2];		
		pointOwnerMatrix(i,0) = pointOwners[i];
		pointTimesMatrix(i,0) = pointTimeSteps[i];
	}



	utils::MatrixDataSet::Pointer dataSet = utils::MatrixDataSet::New();
	dataSet->AddData("planes", planes);
	dataSet->AddData("time_steps", timeSteps);
	dataSet->AddData("owners", owners);
	dataSet->AddData("points", pointMatrix);
	dataSet->AddData("point_owners", pointOwnerMatrix);
	dataSet->AddData("point_times", pointTimesMatrix);

	utils::MatrixWriter::Pointer writer = utils::MatrixWriter::New();
	writer->SetInput(dataSet);
	writer->SetFilename(argv[3]);
	writer->Write();


	//for(unsigned int i = 0; i < numberOfTimeSteps; i++)
	//{
		//SampleType::Pointer samples = samplesList[i];

		//CovarianceFilterType::Pointer covFilter = CovarianceFilterType::New();
		//covFilter->SetInput(samples);
		//covFilter->Update();

		//MeasurementType mv = covFilter->GetMean();


		//vtkSmartPointer<vtkPlaneSource> source = vtkSmartPointer<vtkPlaneSource>::New();
		//source->SetCenter(mv[0], mv[1], mv[2]);
		//source->SetNormal(mv[3], mv[4], mv[5]);
		//source->SetResolution(1,1);
		//source->Update();

		//vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
		//poly = source->GetOutput();

		//vtkSmartPointer<vtkPoints> p = poly->GetPoints();
		//double cent[3];
		//cent[0] = 0.0; cent[1] = 0.0; cent[2] = 0.0;

		//for(unsigned int j = 0; j < p->GetNumberOfPoints(); j++)
		//{
			//cent[0] += p->GetPoint(j)[0] / (double) p->GetNumberOfPoints();
			//cent[1] += p->GetPoint(j)[1] / (double) p->GetNumberOfPoints();
			//cent[2] += p->GetPoint(j)[2] / (double) p->GetNumberOfPoints();
		//}

		//for(unsigned int j = 0; j < p->GetNumberOfPoints(); j++)
		//{
			//double pin[3];
			//p->GetPoint(j,pin);		

			//for(unsigned int k = 0; k < 3; k++)
			//{
				//pin[k] -= cent[k];
				//pin[k] *= 100;
				//pin[k] += cent[k];
			//}

			//p->SetPoint(j,pin);

		//}

		//std::stringstream ss;
		//ss << i << ".vtk";

		//vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
		//writer->SetFileName(ss.str().c_str());
		//writer->SetInputData(poly);
		//writer->Write();

	//}





	return 0;
}

/*
			typedef CMRFileExtractor::ImageType ImageType;
			ImageType::DirectionType dir1 = extractor->Get3CImage(ts)->GetDirection();
			ImageType::DirectionType dir2 = extractor->Get2CImage(ts)->GetDirection();
			ImageType::DirectionType newDir1, newDir2;
			for(unsigned int j = 0; j < 3; j++)
			{
				ValvePlane::VectorType v1;
				v1[0] = dir1(0,j);
				v1[1] = dir1(1,j);
				v1[2] = dir1(2,j);
				ValvePlane::VectorType v2 = transform->GetInverseTransform()->TransformCovariantVector(v1);
				newDir1(0,j) = v2[0];
				newDir1(1,j) = v2[1];
				newDir1(2,j) = v2[2];


				v1[0] = dir2(0,j);
				v1[1] = dir2(1,j);
				v1[2] = dir2(2,j);
				v2 = transform->GetInverseTransform()->TransformCovariantVector(v1);
				newDir2(0,j) = v2[0];
				newDir2(1,j) = v2[1];
				newDir2(2,j) = v2[2];

			}

			typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
			ResamplerType::Pointer resampler = ResamplerType::New();
			resampler->SetInput(extractor->Get2CImage(ts));
			resampler->SetTransform(transform);
			resampler->SetOutputParametersFromImage(extractor->Get2CImage(ts));
			resampler->SetOutputOrigin(transform->GetInverseTransform()->TransformPoint(extractor->Get2CImage(ts)->GetOrigin()));
			resampler->SetOutputDirection(newDir2);
			resampler->Update();

			utils::ImageVolumeIO::Write("image2C.nrrd", resampler->GetOutput());

			ResamplerType::Pointer resampler2 = ResamplerType::New();
			resampler2->SetInput(extractor->Get3CImage(ts));
			resampler2->SetTransform(transform);
			resampler2->SetOutputParametersFromImage(extractor->Get3CImage(ts));
			resampler2->SetOutputOrigin(transform->GetInverseTransform()->TransformPoint(extractor->Get3CImage(ts)->GetOrigin()));
			resampler2->SetOutputDirection(newDir1);
			resampler2->Update();

			utils::ImageVolumeIO::Write("image3C.nrrd", resampler2->GetOutput());
*/
