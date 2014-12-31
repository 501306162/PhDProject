#include <iostream>

#include <CMRFileExtractor.h>
#include <Directory.h>
#include <ValveIO.h>
#include <ValveOriginFinder.h>
#include <itkSimilarity3DTransform.h>

#include <QString>
#include <QFileInfo>
#include <QVariant>
#include <QDir>
#include <QTextStream>

#include <qjson/serializer.h>

using namespace vt;

typedef ValveOriginFinder::PointType PointType;
typedef ValveOriginFinder::RotationType RotationType;
typedef ValveOriginFinder::VectorType VectorType;


bool filename_sort(const std::string &f1, const std::string &f2);
std::string getDataFolder(const std::string &filename, std::string replace, std::string dataFolder);
void computeTransform(const std::string &dataFilename, RotationType &rotation, VectorType &trans);
std::string getEntryName(const std::string & filePath);


int main(int argc, char ** argv)
{
	// read the input files
	std::string inputFolder = argv[1];
	std::string dFolder = argv[2];
	std::string outputFilename = argv[3];


	typedef utils::Directory Directory;
	Directory::FilenamesType filenames = Directory::GetFiles(inputFolder, "txt");

	// load the valves
	std::sort(filenames.begin(), filenames.end(), filename_sort);

	std::vector<PointType> p1List, p2List;
	std::vector<std::string> nameList;
	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		std::string dataFolder = getDataFolder(filenames[i], inputFolder, dFolder);

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

		p1List.push_back(pt1);
		p2List.push_back(pt2);
		nameList.push_back(dataFolder);

	}


	// create the output file
	QVariantMap outputList;
	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		std::string entryName = getEntryName(filenames[i]);
		QVariantList p1, p2;
		
		p1 << p1List[i][0];
		p1 << p1List[i][1];
		p1 << p1List[i][2];
		p2 << p2List[i][0];
		p2 << p2List[i][1];
		p2 << p2List[i][2];

		QVariantMap entry;
		entry["p1"] = p1;
		entry["p2"] = p2;
		
		outputList[QString::fromStdString(entryName)] = entry;
		
	}

	QFile outputFile(QString::fromStdString(outputFilename));
	outputFile.open(QIODevice::WriteOnly | QIODevice::Text);

	
	QJson::Serializer serialiser;
	bool ok;
	QString json = serialiser.serialize(outputList, &ok);
	
	QTextStream out(&outputFile);
	out << json;
	outputFile.close();

	return 0;

}

// ------------------------------------------------------------------------
std::string getEntryName(const std::string & filePath)
{
	QFileInfo f(QString::fromStdString(filePath));
	QString name = f.fileName();
	name = name.replace(".txt", "");
	return name.toStdString();
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
std::string getDataFolder(const std::string &filename, std::string inputFolder, std::string dataFolder)
{
	QString s1 = QString::fromStdString(filename);
	s1 = s1.replace(QString::fromStdString(inputFolder),"");
	s1 = s1.replace("/","");
	s1 = s1.replace(".txt","");

	QDir folder(QString::fromStdString(dataFolder));
	QString par1 = folder.absoluteFilePath(s1); 
	return par1.toStdString();

}


// ------------------------------------------------------------------------
bool filename_sort(const std::string &f1, const std::string &f2)
{
	QString s1 = QString::fromStdString(getEntryName(f1));
	QString s2 = QString::fromStdString(getEntryName(f2));

	s1 = s1.replace("d","");
	s2 = s2.replace("d","");

	return (s1.toInt() < s2.toInt());
}

/*
	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		std::string dataFolder = getDataFolder(filenames[i], inputFolder);
		std::cout << "Folder: " << dataFolder  << std::endl;
		QDir dir(QString::fromStdString(dataFolder));
		if(!dir.exists()) continue;


		RotationType rot;
		VectorType trans;
		std::cout << dataFolder << std::endl;
		computeTransform(dataFolder, rot, trans);


		vtkSmartPointer<vtkTransform> transform = 
			vtkSmartPointer<vtkTransform>::New();

		//create the matrix
		vtkSmartPointer<vtkMatrix4x4> T = vtkSmartPointer<vtkMatrix4x4>::New();
		T->Zero();
		for(unsigned int j = 0; j < 3; j++)
		{
			for(unsigned int k = 0; k < 3; k++)
			{
				T->SetElement(j,k,rot(j,k));
			}
			T->SetElement(j,3,trans[j]);
		}
		T->SetElement(3,3,1);

		transform->SetMatrix(T);
		transform->Update();
		
		// transform the points 
		vtkBoundingBox box;
		for(unsigned int j = 0; j < locations->GetNumberOfPoints(); j++)
		{
			double outputPoint[3];
			transform->TransformPoint(locations->GetPoint(j), outputPoint);
			box.AddPoint(outputPoint);
		}
		

		box.Inflate(0.5);

		vtkSmartPointer<vtkCubeSource> source = vtkSmartPointer<vtkCubeSource>::New();
		double bounds[6];
		box.GetBounds(bounds);
		source->SetBounds(bounds);
		source->Update();

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

		CMRFileExtractor::Pointer fileExtractor = CMRFileExtractor::New();
		fileExtractor->SetDebug(true);
		std::cout << "-----------------------------------------------" << std::endl;
		fileExtractor->SetFolderName(dataFolder);
		fileExtractor->Extract();

		typedef ValveSequence<3>::ImageType ImageType;
	

		vtkBoundingBox check;
		double b[6];
		filter->GetOutput()->GetBounds(b);
		check.SetBounds(b);

		
		// create the output image 
		ImageType::Pointer image = fileExtractor->Get3CImage(0);
		typedef itk::Image<double, 3> FloatType;
		FloatType::Pointer output = FloatType::New();
		output->SetOrigin(image->GetOrigin());
		output->SetDirection(image->GetDirection());
		output->SetSpacing(image->GetSpacing());
		output->SetRegions(image->GetLargestPossibleRegion());
		output->Allocate();
		output->FillBuffer(0.0);


		typedef itk::ImageRegionIterator<ImageType> IteratorType;
		typedef itk::ImageRegionIterator<FloatType> Iterator2Type;
		IteratorType inIt(image, image->GetLargestPossibleRegion());
		Iterator2Type outIt(output, output->GetLargestPossibleRegion());

		double h = atof(argv[3]);

		while(!inIt.IsAtEnd())
		{
			ImageType::IndexType index = inIt.GetIndex();
			ImageType::PointType point;
			image->TransformIndexToPhysicalPoint(index, point);

			// loop through the points
			vtkSmartPointer<vtkPoints> points = filter->GetOutput()->GetPoints();
			double kde = 0.0;
			for(unsigned int i = 0; i < points->GetNumberOfPoints(); i++)
			{
				double p[3];
				points->GetPoint(i,p);				

				typedef ImageType::PointType MVType;
				//typedef itk::Statistics::GaussianMembershipFunction<MVType> GaussianType;
				//GaussianType::Pointer gaussian = GaussianType::New();

				MVType mean;
				mean[0] = p[0];
				mean[1] = p[1];
				mean[2] = p[2];

				MVType t;
				t[0] = point[0];
				t[1] = point[1];
				t[2] = point[2];

				double dist = mean.EuclideanDistanceTo(t);

				itk::Statistics::GaussianDistribution::Pointer gaussian = itk::Statistics::GaussianDistribution::New();
				gaussian->SetMean(0.0);
				gaussian->SetVariance(1);
				kde += gaussian->EvaluatePDF(dist / h);

			}

			kde /= (((double) points->GetNumberOfPoints()) * h);


			
			outIt.Set(kde);
			if(check.ContainsPoint(point.GetDataPointer()))
			{
			//outIt.Set(inIt.Get());
			}
			else
			{
				//outIt.Set((unsigned short) (inIt.Get()/3));
			}

			++inIt; ++outIt;
		}
		
		std::stringstream ss2;
		ss2 << dir.dirName().toStdString() << ".nrrd";


		typedef itk::Image<unsigned char, 3> OutType;
		typedef itk::RescaleIntensityImageFilter<FloatType, OutType> RescalerType;
		RescalerType::Pointer rescaler = RescalerType::New();
		rescaler->SetInput(output);
		rescaler->SetOutputMaximum(255);
		rescaler->SetOutputMinimum(0);
		rescaler->Update();


		typedef itk::ImageFileWriter<FloatType> WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName(ss2.str());
		writer->SetImageIO(itk::NrrdImageIO::New());
		writer->SetInput(output);
		writer->Update();




		std::stringstream ss; 

		ss << dataFolder << ".vtk";
		vtkSmartPointer<vtkPolyDataWriter> writer2 = vtkSmartPointer<vtkPolyDataWriter>::New();
		writer2->SetInputData(filter->GetOutput());
		writer2->SetFileName(ss.str().c_str());
		writer2->Update();

		*/
