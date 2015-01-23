#include <iostream>

#include <Directory.h>
#include <ValveIO.h>
#include <QFileInfo>
#include <MatrixCommon.h>
#include <vtkPlaneSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <QVariant>
#include <QDir>

#include <QDebug>
#include <qjson/serializer.h>
#include <QFile>

using namespace vt;


typedef struct _v_instances
{
	std::string name;
	std::vector<std::string> filenames;
	std::string type;
	typedef std::vector<_v_instances> List;
} ValveInstance;
typedef std::map<std::string, ValveInstance::List> DataMap;
typedef std::map<std::string, std::vector<std::string> > FilenameMap;
typedef utils::DoubleMatrixType MatrixType;
typedef ValveSequence<3> SequenceType;
typedef SequenceType::ValveLineType LineType;
typedef LineType::PointType PointType;

void processInput(const std::string &inputDir, DataMap &output);
void getWantedDirectories(const std::string &inputDirs, 
		FilenameMap &dirs);
void getFiles(const std::vector<std::string> &dirs, FilenameMap &filenames);
void getStudyNames(const FilenameMap & filenames, std::vector<std::string> &studyNames);
void getDataWithEnoughPoints(const std::vector<std::string> &directories, const std::string type, std::vector<ValveInstance> &valves);
bool filenamesContainStudy(const std::string &study, const std::vector<std::string> &filenames, std::string &filename);
void computePlanes(ValveInstance &valve, QVariantList &output);
void computePlane(std::vector<PointType> &lines, MatrixType &normal, MatrixType &point);
QVariantList getPointData(const PointType &point);


int main(int argc, char **argv)
{
	const std::string inputDirectory = argv[1];
	const std::string outputDirectory = argv[2];
	DataMap data;
	processInput(inputDirectory, data);


	// process the valve points to create the planes
	DataMap::iterator mapIt = data.begin();
	while(mapIt != data.end())
	{

		// create the foler
		QDir outDir(QString::fromStdString(outputDirectory));
		QString subDir = QString::fromStdString(mapIt->first);
		QString newDir = outDir.absoluteFilePath(subDir); 
		outDir.mkpath(newDir);

		ValveInstance::List &valves = mapIt->second;
		for(unsigned int i = 0; i < valves.size(); i++)
		{
			
			MatrixType normals, points;
			QVariantList output;
			computePlanes(valves[i], output);			

			// create the filename 
			QString outputFilename = QDir(newDir).absoluteFilePath(QString::fromStdString(valves[i].name));
			QFile file(outputFilename);
			file.open(QIODevice::Text | QIODevice::WriteOnly);

			QJson::Serializer serialiser;
			serialiser.setIndentMode(QJson::IndentFull);
			bool ok;
			QString json = serialiser.serialize(output, &ok);
			QTextStream out(&file);

			out << json;
			file.close();

		}
		
		++mapIt;
	}

	


	return 0;
}

// ------------------------------------------------------------------------
void computePlanes(ValveInstance &valve, QVariantList &output)
{
	QVariantList valvesData;

	// load up the entries
	std::vector<SequenceType::Pointer> sequences;
	std::vector<PointType> points;
	for(unsigned int i = 0; i < valve.filenames.size(); i++)
	{
		ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
		reader->SetFileName(valve.filenames[i]);
		sequences.push_back(reader->GetOutput());
	}

	unsigned int timesteps = sequences.front()->GetNumberOfLines();
	for(unsigned int i = 0; i < timesteps; i++)
	{
		QVariantMap valveData;

		std::vector<PointType> inputPoints;
		QVariantList pointsData;

		for(unsigned int j = 0; j < sequences.size(); j++)
		{
			// assign the points
			QString filename = QString::fromStdString(valve.filenames[j]);
			filename = filename.replace("/" + QString::fromStdString(valve.name), "");
			QFileInfo info(filename);

			QVariantMap pointData;
			pointData["type"] = info.fileName();
			pointData["P1"] = getPointData(sequences[j]->GetValveLine(i)->GetP1());
			pointData["P2"] = getPointData(sequences[j]->GetValveLine(i)->GetP2());


			inputPoints.push_back(sequences[j]->GetValveLine(i)->GetP1());
			inputPoints.push_back(sequences[j]->GetValveLine(i)->GetP2());

			pointsData << pointData;
		}

		valveData["points"] = pointsData;

		MatrixType normal;
		MatrixType point;

		computePlane(inputPoints, normal, point);

		QVariantList normalData, centerData;
		for(unsigned int i = 0; i < 3; i++)
		{
			normalData << normal(i,0);
			centerData << point(i,0);
		}

		valveData["normal"] = normalData;
		valveData["center"] = centerData;
		
		valvesData << valveData;
	}

	output = valvesData;
}

// ------------------------------------------------------------------------
void computePlane(std::vector<PointType> &inputPoints, MatrixType &normal, MatrixType &point)
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

	normal = solver.eigenvectors().col(0).transpose();
	point = centroid.transpose();

	/*
	vtkSmartPointer<vtkPlaneSource> source = vtkSmartPointer<vtkPlaneSource>::New();
	source->SetCenter(centroid(0,0), centroid(1,0), centroid(2,0));
	source->SetNormal(normal(0,0), normal(1,0), normal(2,0));
	source->SetResolution(1,1);
	source->Update();

	vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
	poly = source->GetOutput();

	vtkSmartPointer<vtkPoints> p = poly->GetPoints();
	double cent[3];
	cent[0] = 0.0; cent[1] = 0.0; cent[2] = 0.0;

	for(unsigned int i = 0; i < p->GetNumberOfPoints(); i++)
	{
		cent[0] += p->GetPoint(i)[0] / (double) p->GetNumberOfPoints();
		cent[1] += p->GetPoint(i)[1] / (double) p->GetNumberOfPoints();
		cent[2] += p->GetPoint(i)[2] / (double) p->GetNumberOfPoints();
	}

	for(unsigned int i = 0; i < p->GetNumberOfPoints(); i++)
	{
		double pin[3];
		p->GetPoint(i,pin);		
		
		for(unsigned int j = 0; j < 3; j++)
		{
			pin[j] -= cent[j];
			pin[j] *= 100;
			pin[j] += cent[j];
		}

		p->SetPoint(i,pin);

	}
	*/
}



// ------------------------------------------------------------------------
void getDataWithEnoughPoints(const std::vector<std::string> &directories, const std::string type,
		std::vector<ValveInstance> &valves)
{
	FilenameMap filenames;
	getFiles(directories, filenames);

	std::vector<std::string> studyNames;
	getStudyNames(filenames, studyNames);


	// check for each study if there is enough points
	std::vector<ValveInstance> allValves;
	for(unsigned int i = 0; i < studyNames.size(); i++)
	{
		FilenameMap::iterator mapIt = filenames.begin();
		std::string studyName = studyNames[i];

		ValveInstance instance;
		instance.name = studyName;
		instance.type = type;

		while(mapIt != filenames.end())
		{
			std::string fname;
			if(filenamesContainStudy(studyName, mapIt->second, fname))
			{
				instance.filenames.push_back(fname);
			}
			++mapIt;
		}

		allValves.push_back(instance);
	}

	for(unsigned int i = 0; i < allValves.size(); i++)
	{
		if(allValves[i].filenames.size() >= 2)
			valves.push_back(allValves[i]);
	}
}






// ------------------------------------------------------------------------
void processInput(const std::string &inputDir, DataMap &output)
{
	FilenameMap directories;
	getWantedDirectories(inputDir, directories);

	// get the set of files for each of the directories 
	FilenameMap filenames;
	std::vector<ValveInstance> mvValves;
	getDataWithEnoughPoints(directories["MV"], "MV", mvValves);
	std::vector<ValveInstance> tpValves;
	getDataWithEnoughPoints(directories["TP"], "TP", tpValves);

	output["MV"] = mvValves;
	output["TP"] = tpValves;
}





// ------------------------------------------------------------------------
void getFiles(const std::vector<std::string> &dirs, FilenameMap &filenames)
{
	for(unsigned int i = 0; i < dirs.size(); i++)
	{
		utils::Directory::FilenamesType files = utils::Directory::GetFiles(dirs[i], ".txt");
		QFileInfo d(QString::fromStdString(dirs[i]));
		
		filenames[d.fileName().toStdString()] = files;
	}
}


// ------------------------------------------------------------------------
void getStudyNames(const FilenameMap & filenames, std::vector<std::string> &studyNames)
{
	std::map<std::string, int> checkMap;
	FilenameMap::const_iterator mapIt = filenames.begin();
	while(mapIt != filenames.end())
	{
		for(unsigned int i = 0; i < mapIt->second.size(); i++)
		{
			QFileInfo info(QString::fromStdString(mapIt->second[i]));
			std::string fpath = info.fileName().toStdString();
			if(checkMap.count(fpath) == 0)
			{
				checkMap[fpath] = 1;
				studyNames.push_back(fpath);
			}
		}
		++mapIt;
	}
}


// ------------------------------------------------------------------------
void getWantedDirectories(const std::string &inputDir, 
		 FilenameMap &dirs)
{
	utils::Directory::FilenamesType directories = utils::Directory::GetDirectories(inputDir);
	for(unsigned int i = 0; i < directories.size(); i++)
	{
		if(directories[i].find("MV") != std::string::npos)
			dirs["MV"].push_back(directories[i]);
		if(directories[i].find("TP") != std::string::npos)
			dirs["TP"].push_back(directories[i]);	
	}
}

// ------------------------------------------------------------------------
bool filenamesContainStudy(const std::string &study, const std::vector<std::string> &filenames, std::string &filename)
{
	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		if(filenames[i].find(study) != std::string::npos)
		{
			filename = filenames[i];
			return true;
		}
	}
	filename = "";
	return false;
}

// ------------------------------------------------------------------------
QVariantList getPointData(const PointType &point)
{
	QVariantList pointData;
	for(unsigned int i = 0; i < 3; i++)
	{
		pointData << point[i];	
	}

	return pointData;
}

