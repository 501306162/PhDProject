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
#include <vtkPlaneSource.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <itkSimilarity3DTransform.h>

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
	

	FilenamesType planeFilenames = utils::Directory::GetFiles(inputDirectory, ".txt");
	for(unsigned int i = 0; i < planeFilenames.size(); i++)
	{
		// load the valve plane
		ValvePlaneSequence::Pointer sequence = ValvePlaneSequence::Load(planeFilenames[i]);
		ValvePlane::Pointer valve = sequence->GetValvePlane(0);
		

		// get the data directory
		std::string dataFolder = utils::Directory::GetPath(dataDirectory, sequence->GetName());


		// load up the data 
		CMRFileExtractor::Pointer extractor = CMRFileExtractor::New();
		std::cout << dataFolder << std::endl;
		extractor->SetFolderName(dataFolder);
		extractor->SetDebug(true);
		//extractor->SetFlip(true);
		extractor->Extract();

		utils::ImageVolumeIO::Write("stack.nrrd", extractor->GetStackImage(0));
		utils::ImageVolumeIO::Write("2Co.nrrd", extractor->Get2CImage(0));
		utils::ImageVolumeIO::Write("3Co.nrrd", extractor->Get3CImage(0));


		ValveOriginFinder::Pointer originFinder = ValveOriginFinder::New();
		originFinder->Set2CImage(extractor->Get2CImage(0));
		originFinder->Set3CImage(extractor->Get3CImage(0));
		originFinder->SetImageStack(extractor->GetStackImage(0));
		originFinder->Compute();

		typedef itk::Similarity3DTransform<double> TransformType;
		TransformType::Pointer transform = TransformType::New();
		transform->SetTranslation(originFinder->GetTranslation());
		transform->SetMatrix(originFinder->GetRotation());


		
		std::cout << valve->GetNormal() << std::endl;
		std::cout << valve->GetCenter() << std::endl;
		//ValvePlane::VectorType normal = transform->GetInverseTransform()->TransformCovariantVector(valve->GetNormal());
		//ValvePlane::PointType center = transform->GetInverseTransform()->TransformPoint(valve->GetCenter());
		ValvePlane::VectorType normal = valve->GetNormal();
		ValvePlane::PointType center = valve->GetCenter();
		std::cout << normal << std::endl;
		std::cout << center << std::endl;

		for(unsigned int j = 0; j < 3; j++)
		{
			meanNormal[j] += (normal[j] / (double) planeFilenames.size());
			meanPoint[j] += (center[j] / (double) planeFilenames.size());
		}

		std::cout << sequence->GetName() << std::endl;

		typedef CMRFileExtractor::ImageType ImageType;
		utils::ImageVolumeIO::Write("orig.nrrd", extractor->Get2CImage(0));

		ImageType::DirectionType dir1 = extractor->Get3CImage(0)->GetDirection();
		ImageType::DirectionType dir2 = extractor->Get2CImage(0)->GetDirection();
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
		resampler->SetInput(extractor->Get2CImage(0));
		resampler->SetTransform(transform);
		resampler->SetOutputParametersFromImage(extractor->Get2CImage(0));
		resampler->SetOutputOrigin(transform->GetInverseTransform()->TransformPoint(extractor->Get2CImage(0)->GetOrigin()));
		resampler->SetOutputDirection(newDir2);
		resampler->Update();
		
		utils::ImageVolumeIO::Write("2C.nrrd", resampler->GetOutput());



		ResamplerType::Pointer resampler2 = ResamplerType::New();
		resampler2->SetInput(extractor->Get3CImage(0));
		resampler2->SetTransform(transform);
		resampler2->SetOutputParametersFromImage(extractor->Get3CImage(0));
		resampler2->SetOutputOrigin(transform->GetInverseTransform()->TransformPoint(extractor->Get3CImage(0)->GetOrigin()));
		resampler2->SetOutputDirection(newDir1);
		resampler2->Update();

		utils::ImageVolumeIO::Write("3C.nrrd", resampler2->GetOutput());
		
		vtkSmartPointer<vtkPlaneSource> source = vtkSmartPointer<vtkPlaneSource>::New();
		source->SetCenter(center[0], center[1], center[2]);
		source->SetNormal(normal[0], normal[1], normal[2]);
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

		std::stringstream ss;
		ss << i << ".vtk";

		vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
		writer->SetFileName(ss.str().c_str());
		writer->SetInputData(poly);
		writer->Write();

	}



	vtkSmartPointer<vtkPlaneSource> source = vtkSmartPointer<vtkPlaneSource>::New();
	source->SetCenter(meanPoint[0], meanPoint[1], meanPoint[2]);
	source->SetNormal(meanNormal[0], meanNormal[1], meanNormal[2]);
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
	writer->SetFileName("mean.vtk");
	writer->SetInputData(poly);
	writer->Write();



	return 0;
}


