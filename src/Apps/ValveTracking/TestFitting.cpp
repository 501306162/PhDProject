#include <iostream>

#include <QFileInfo>
#include <QDir>

#include <itkSimilarity3DTransform.h>
#include <itkResampleImageFilter.h>

#include <CMRFileExtractor.h>
#include <ValvePlane.h>
#include <TrainingData.h>
#include <ValveOriginFinder.h>
#include <CommonDefinitions.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkPlaneSource.h>
#include <vtkSmartPointer.h>

#define DATA_PATH "/home/om0000/ValveTracking/"

using namespace vt;



typedef std::vector<TrainingData::Pointer> TrainingPointList;
typedef std::map<std::string, TrainingPointList> TrainingDataMap;
typedef itk::Image<unsigned short, 3> ImageType;
typedef itk::Image<unsigned char, 3> LabelType;
typedef utils::DoubleMatrixType MatrixType;

typedef struct image_group_
{
	ImageType::Pointer image;
	LabelType::Pointer seg;
	LabelType::Pointer contours;
} ImageGroup;


typedef struct test_data_
{
	unsigned int id;
	ImageGroup c2Image;
	ImageGroup c3Image;
	ImageGroup c4Image;
	ImageType::Pointer stackImage;

} TestData;

typedef itk::Similarity3DTransform<double> TransformType;


void savePlane(const MatrixType &plane, std::string filename);
void loadTrainingData(TrainingDataMap &trainingData, ValveTrainingData::Pointer &valves);
void loadTestData(const std::string &folderName, TestData &testData);
void computeTransform(const TestData &testData, TransformType::Pointer &transform);
void applyTransformToTestData(const TestData &input, TransformType::Pointer &transform, TestData &output);
void applyTransform(const ImageType::Pointer &input, TransformType::Pointer &transform, ImageType::Pointer &output);
void meanPlane(const ValveTrainingData::Pointer &data, const unsigned int exclude, MatrixType &plane);

int main(int argc, char ** argv)
{
	// load the training data
	TrainingDataMap trainingData;
	ValveTrainingData::Pointer valveData;
	loadTrainingData(trainingData, valveData);


	// load the input data
	const std::string inputDataFolder = argv[1];
	TestData testData, alignedTestData;
	loadTestData(inputDataFolder, testData);
	
	// compute and apply the transform
	TransformType::Pointer transform = TransformType::New();
	computeTransform(testData, transform);
	applyTransformToTestData(testData, transform, alignedTestData);
	

	// compute the segmentation and contours


	MatrixType mean;
	meanPlane(valveData, testData.id, mean);
	savePlane(mean, "meanPlane.vtk");
	std::cout << mean << std::endl;








	return 0;
}

// ------------------------------------------------------------------------
void segmentTestData(TestData &data)
{

}


// ------------------------------------------------------------------------
void applyTransformToTestData(const TestData &input, TransformType::Pointer &transform, TestData &output)
{
	output.c2Image.image = ImageType::New();
	output.c3Image.image = ImageType::New();
	output.c4Image.image = ImageType::New();
	output.stackImage = ImageType::New();

	applyTransform(input.c2Image.image, transform, output.c2Image.image);
	applyTransform(input.c3Image.image, transform, output.c3Image.image);
	applyTransform(input.c4Image.image, transform, output.c4Image.image);
	applyTransform(input.stackImage, transform, output.stackImage);

	utils::ImageVolumeIO::Write("2C.nrrd", output.c2Image.image);
	utils::ImageVolumeIO::Write("3C.nrrd", output.c3Image.image);
	utils::ImageVolumeIO::Write("4C.nrrd", output.c4Image.image);
}

// ------------------------------------------------------------------------
void applyTransform(const ImageType::Pointer &input, TransformType::Pointer &transform, ImageType::Pointer &output)
{
	ImageType::DirectionType dir = input->GetDirection();
	ImageType::DirectionType newDir = transform->GetInverseMatrix() * dir;
	ImageType::PointType newOrigin = transform->GetInverseTransform()->TransformPoint(input->GetOrigin());


	typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
	ResamplerType::Pointer resampler = ResamplerType::New();
	resampler->SetInput(input);
	resampler->SetTransform(transform);
	resampler->SetOutputParametersFromImage(input);
	resampler->SetOutputOrigin(newOrigin);
	resampler->SetOutputDirection(newDir);
	resampler->Update();

	output = resampler->GetOutput();


}

// ------------------------------------------------------------------------
void computeTransform(const TestData &testData, TransformType::Pointer &transform)
{
	ValveOriginFinder::Pointer originFinder = ValveOriginFinder::New();
	originFinder->Set2CImage(testData.c2Image.image);
	originFinder->Set3CImage(testData.c3Image.image);
	originFinder->SetImageStack(testData.stackImage);
	originFinder->Compute();
	
	transform->SetTranslation(originFinder->GetTranslation());
	transform->SetMatrix(originFinder->GetRotation());
}



// ------------------------------------------------------------------------
void loadTestData(const std::string &folderName, TestData &testData)
{
	CMRFileExtractor::Pointer fileExtractor = CMRFileExtractor::New();
	fileExtractor->SetFolderName(folderName);
	fileExtractor->Extract();

	testData.c2Image.image = fileExtractor->Get2CImage(0);
	testData.c3Image.image = fileExtractor->Get3CImage(0);
	testData.c4Image.image = fileExtractor->Get4CImage(0);
	testData.stackImage = fileExtractor->GetStackImage(0);

	QFileInfo info(QString::fromStdString(folderName));
	testData.id = info.fileName().replace("d","").replace(".txt", "").toInt();
}



// ------------------------------------------------------------------------
void loadTrainingData(TrainingDataMap &trainingData, ValveTrainingData::Pointer &valves)
{
	TrainingData::Pointer ap1 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/MV-2C-P1.hdf5");
	TrainingData::Pointer ap2 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/MV-2C-P2.hdf5");
	TrainingData::Pointer bp1 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/MV-3C-P1.hdf5");
	TrainingData::Pointer bp2 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/MV-3C-P2.hdf5");
	TrainingData::Pointer cp1 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/MV-4C-P1.hdf5");
	TrainingData::Pointer cp2 = TrainingData::Load("/home/om0000/ValveTracking/TrainingData/MV-4C-P2.hdf5");

	trainingData["MV-2C"].push_back(ap1);
	trainingData["MV-2C"].push_back(ap2);
	trainingData["MV-3C"].push_back(bp1);
	trainingData["MV-3C"].push_back(bp2);
	trainingData["MV-4C"].push_back(cp1);
	trainingData["MV-4C"].push_back(cp2);

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

