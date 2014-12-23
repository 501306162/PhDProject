#include <iostream>

#include <CMRFileExtractor.h>
#include <Directory.h>
#include <ValveIO.h>
#include <CMRFileExtractor.h>
#include <ValveOriginFinder.h>
#include <itkSimilarity3DTransform.h>
#include <QString>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPolyDataWriter.h>
#include <vtkCellArray.h>
#include <vtkBoundingBox.h>
#include <vtkCubeSource.h>

#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkAppendPolyData.h>

#include <QDir>

using namespace vt;

typedef ValveOriginFinder::PointType PointType;
typedef ValveOriginFinder::RotationType RotationType;
typedef ValveOriginFinder::VectorType VectorType;


bool filename_sort(const std::string &f1, const std::string &f2);
std::string getDataFolder(const std::string &filename, std::string replace);
void computeTransform(const std::string &dataFilename, RotationType &rotation, VectorType &trans);


int main(int argc, char ** argv)
{
	// read the input files
	std::string inputFolder = argv[1];
	typedef utils::Directory Directory;
	Directory::FilenamesType filenames = Directory::GetFiles(inputFolder, "txt");



	// load the valves
	std::sort(filenames.begin(), filenames.end(), filename_sort);

	std::vector<PointType> p1List, p2List;
	std::vector<std::string> nameList;
	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		std::string dataFolder = getDataFolder(filenames[i], inputFolder);

		std::cout << "Folder: " << dataFolder  << std::endl;
		QDir dir(QString::fromStdString(dataFolder));
		if(!dir.exists()) continue;


		RotationType rotation;
		VectorType trans;
		computeTransform(dataFolder, rotation, trans);

		std::cout << "------------------------------------------" << std::endl;


		// load the valve sequence 
		ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
		reader->SetFileName(filenames[i]);
		ValveSequence<3>::Pointer seq = reader->GetOutput();

		

		PointType p1 = seq->GetValveLine(0)->GetP1();
		PointType p2 = seq->GetValveLine(0)->GetP2();

		typedef itk::Similarity3DTransform<double> TransformType;
		TransformType::Pointer transform = TransformType::New();
		transform->SetMatrix(rotation);
		transform->SetTranslation(trans);
		PointType pt1 = transform->GetInverseTransform()->TransformPoint(p1);
		PointType pt2 = transform->GetInverseTransform()->TransformPoint(p2);
		//PointType pt1 = transform->TransformPoint(p1);
		//PointType pt2 = transform->TransformPoint(p2);

		p1List.push_back(pt1);
		p2List.push_back(pt2);
		nameList.push_back(dataFolder);

	}
	
	PointType m1, m2;
	m1.Fill(0); m2.Fill(0);
	vtkBoundingBox bounding1, bounding2;
	for(unsigned int i = 0; i < p1List.size(); i++)
	{
		for(unsigned int j = 0; j < 3; j++)
		{
			m1[j] += p1List[i][j];
			m2[j] += p2List[i][j];

		}
		bounding1.AddPoint(p1List[i].GetDataPointer());		
		bounding1.AddPoint(p2List[i].GetDataPointer());		
		bounding2.AddPoint(p2List[i].GetDataPointer());		
	}

	bounding1.Inflate(5.0);

	for(unsigned int i = 0; i < 3; i++)
	{
		m1[i] /= p1List.size();
		m2[i] /= p1List.size();
	}

	for(unsigned int i = 0; i < p1List.size(); i++)
	{
		double dist = (p1List[i]-m1).GetNorm();
		double dist2 = (p2List[i]-m2).GetNorm();
		std::cout << i << " " << nameList[i] << " " << dist << " " << dist2 << std::endl;
	}


	vtkSmartPointer<vtkCubeSource> source = vtkSmartPointer<vtkCubeSource>::New();
	double bounds[6];
	bounding1.GetBounds(bounds);
	source->SetBounds(bounds);
	source->Update();


	vtkSmartPointer<vtkCubeSource> source2 = vtkSmartPointer<vtkCubeSource>::New();
	double bounds2[6];
	bounding2.GetBounds(bounds2);
	source2->SetBounds(bounds2);
	source2->Update();



	for(unsigned int i = 0; i < 1; i++)
	{
		std::string dataFolder = getDataFolder(filenames[i], inputFolder);
		RotationType rot;
		VectorType trans;
		computeTransform(argv[2], rot, trans);


		vtkSmartPointer<vtkTransform> transform = 
			vtkSmartPointer<vtkTransform>::New();

		//create the matrix
		vtkSmartPointer<vtkMatrix4x4> T = vtkSmartPointer<vtkMatrix4x4>::New();
		T->Zero();
		for(unsigned int i = 0; i < 3; i++)
		{
			for(unsigned int j = 0; j < 3; j++)
			{
				T->SetElement(i,j,rot(i,j));
			}
			T->SetElement(i,3,trans[i]);
		}
		T->SetElement(3,3,1);

		transform->SetMatrix(T);
		transform->Update();

		vtkSmartPointer<vtkTransformPolyDataFilter> filter = 
			vtkSmartPointer<vtkTransformPolyDataFilter>::New();
		filter->SetInputData(source->GetOutput());
		filter->SetTransform(transform);
		filter->Update();

		vtkSmartPointer<vtkTransformPolyDataFilter> filter2 = 
			vtkSmartPointer<vtkTransformPolyDataFilter>::New();
		filter2->SetInputData(source2->GetOutput());
		filter2->SetTransform(transform);
		filter2->Update();

		vtkSmartPointer<vtkAppendPolyData> appender = 
			vtkSmartPointer<vtkAppendPolyData>::New();
		appender->AddInputData(filter->GetOutput());
		appender->AddInputData(filter2->GetOutput());
		appender->Update();

		
		std::stringstream ss; 

		ss << i << ".vtk";
		vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
		writer->SetInputData(filter->GetOutput());
		writer->SetFileName(ss.str().c_str());
		writer->Update();


	}


	
	/*
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> vertices =
		vtkSmartPointer<vtkCellArray>::New();
	for(unsigned int i = 0; i < p1List.size(); i++)
	{
		vtkIdType pid[1];
		pid[0] = points->InsertNextPoint(p2List[i][0],p2List[i][1], p2List[i][2]);
		vertices->InsertNextCell(1,pid);

	}

	vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
	poly->SetPoints(points);
	poly->SetVerts(vertices);

	*/



	

	
	


	return 0;

}

// ------------------------------------------------------------------------
void computeTransform(const std::string &dataFilename, RotationType &rotation, VectorType &trans)
{
	CMRFileExtractor::Pointer fileExtractor = CMRFileExtractor::New();
	fileExtractor->SetFolderName(dataFilename);
	fileExtractor->Extract();



	// get the origin and the rotation
	ValveOriginFinder::Pointer originFinder = ValveOriginFinder::New();
	originFinder->Set2CImage(fileExtractor->Get2CImage(0));
	originFinder->Set3CImage(fileExtractor->Get3CImage(0));
	originFinder->SetImageStack(fileExtractor->GetStackImage(0));
	originFinder->Compute();



	PointType origin = originFinder->GetOrigin();
	rotation = originFinder->GetRotation();

	trans[0] = origin[0];
	trans[1] = origin[1];
	trans[2] = origin[2];
}



// ------------------------------------------------------------------------
std::string getDataFolder(const std::string &filename, std::string inputFolder)
{
	QString s1 = QString::fromStdString(filename);
	s1 = s1.replace(QString::fromStdString(inputFolder),"");
	s1 = s1.replace("/","");
	s1 = s1.replace(".txt","");

	QString par1 = "/Users/oliverferoze/ValveTracking/Testing/" + s1;
	return par1.toStdString();

}


// ------------------------------------------------------------------------
bool filename_sort(const std::string &f1, const std::string &f2)
{
	QString s1 = QString::fromStdString(f1);
	QString s2 = QString::fromStdString(f2);
	s1 = s1.replace("/Users/oliverferoze/ValveTracking/SortedLines/MV-2C/d","");
	s1 = s1.replace(".txt","");

	s2 = s2.replace("/Users/oliverferoze/ValveTracking/SortedLines/MV-2C/d","");
	s2 = s2.replace(".txt","");

	return (s1.toInt() < s2.toInt());
}

