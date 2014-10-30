#include "io.h"


#include <qjson/serializer.h>
#include <qjson/parser.h>

#include <QFile>


// ------------------------------------------------------------------------
IO::IO(const std::string & filename, DataContainer * data)
{
	this->data = data;
	this->filename = filename;
}


// ------------------------------------------------------------------------
DataContainer * IO::Load(const std::string & filename)
{
	DataContainer * d = new DataContainer;

	IO io(filename, d);
	io.loadData();
	return d;
}


// ------------------------------------------------------------------------
bool IO::loadData()
{
	// try and open the file
	QFile file(QString::fromStdString(filename));
	if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;

	QTextStream in(&file);
	QJson::Parser parser;

	QVariantMap json = parser.parse(in.device()).toMap();

	std::string folderName = json["folder"].toString().toStdString();
	data->LoadData(folderName);
	loadImages(json);


	return true;
}


// ------------------------------------------------------------------------
void IO::loadImages(QVariantMap &input)
{
	QVariantList images = input["images"].toList();
	for(int i = 0; i < images.size(); i++)
	{
		QVariantMap imageData = images[i].toMap();

		QVariantList lines = imageData["lines"].toList();
		loadLines(lines, imageData["name"].toString().toStdString());
		 
	}

}

// ------------------------------------------------------------------------
void IO::loadLines(QVariantList &lines, const std::string &imageName)
{

	for (int i = 0; i < lines.size(); ++i)
	{
		QVariantMap line = lines[i].toMap();
		int slice = line["s"].toInt();
		int timeStep = line["t"].toInt();

		if(line.contains("MV"))
		{
			QVariantMap lmap = line["MV"].toMap();
			makeLine(lmap, imageName, Line::MV, timeStep, slice);			
		}
		if(line.contains("TP"))
		{
			QVariantMap lmap = line["TP"].toMap();
			makeLine(lmap, imageName, Line::TP, timeStep, slice);			
		}
		if(line.contains("AV"))
		{
			QVariantMap lmap = line["AV"].toMap();
			makeLine(lmap, imageName, Line::AV, timeStep, slice);			
		}

	}
}

// ------------------------------------------------------------------------
void IO::makeLine(QVariantMap &line, const std::string &imageName, Line::Type type, int t, int s)
{

	double * p1 = new double[3];
	double * p2 = new double[3];

	QVariantMap p1map = line["p1"].toMap();
	p1[0] = p1map["x"].toDouble();
	p1[1] = p1map["y"].toDouble();
	p1[2] = p1map["z"].toDouble();

	QVariantMap p2map = line["p2"].toMap();
	p2[0] = p2map["x"].toDouble();
	p2[1] = p2map["y"].toDouble();
	p2[2] = p2map["z"].toDouble();


	Line * lineOb = Line::NewLine(type, p1, p2);

	data->getInstance(imageName).images[t][s].lines[lineOb->getType()] = lineOb;
}

// ------------------------------------------------------------------------
bool IO::Save(const std::string &filename, DataContainer * data)
{
	IO io(filename, data);
	io.buildData();
	


	return io.isItOk();
}


// ------------------------------------------------------------------------
void IO::buildData()
{
	outputData.insert("folder", QString::fromStdString(data->getFolderName()));
	QVariantList imageList;
	for(unsigned int i = 0; i < data->numImages(); i++)
	{
		QVariantMap imageMap;
		buildImageData(i, imageMap);
		imageList << imageMap;		

	}

	outputData.insert("images", imageList);

	

	QJson::Serializer serialiser;

	serialiser.setIndentMode(QJson::IndentFull);
	QByteArray json = serialiser.serialize(outputData, &isOk);
	
	if(isOk)
	{
		QFile file(QString::fromStdString(filename));
		if(!file.open(QIODevice::WriteOnly | QIODevice::Text))
		{
			isOk = false;
			return;
		}
		QTextStream out(&file);
		out << QString::fromAscii(json);
		file.close();
		std::cout << "File written OK" << std::endl;
	}	
}

// ------------------------------------------------------------------------
void IO::buildImageData(unsigned int index, QVariantMap &dmap)
{
	DataInstance & instance = data->getInstance(index);
	dmap.insert("name", QString::fromStdString(instance.filename));
	QVariantList seqList;
	for(unsigned int i = 0; i < instance.images.size(); i++)
	{
		VolumeData &vol = instance.images[i];
		QVariantList volList;
		for(unsigned int j = 0; j < vol.size(); j++)
		{
			QVariantMap lineList;
			getLines(vol[j], lineList);

			if(lineList.size() > 0)
			{
				lineList.insert("t", i);
				lineList.insert("s", j);
				volList << lineList;
			}
		}

		if(volList.size() > 0)
			seqList << volList;
	}

	if(seqList.size() > 0)
		dmap.insert("lines", seqList);
}

// ------------------------------------------------------------------------
void IO::getLines(DataHolder &holder, QVariantMap &lines)
{
	Line::Map lineMap = holder.lines;
	Line::Map::iterator mapIt = lineMap.begin();
	while(mapIt != lineMap.end())
	{
		Line * line = mapIt->second;
		QVariantMap lmap;
		getLineData(line->getPoly()->GetPoints(), lmap);

		lines.insert(QString::fromStdString(line->getTypeString(line->getType())), lmap);
		++mapIt;
	}

}


// ------------------------------------------------------------------------
void IO::getLineData(vtkPoints * points, QVariantMap &pointData)
{
	QVariantMap p1;
	p1.insert("x", points->GetPoint(0)[0]);
	p1.insert("y", points->GetPoint(0)[1]);
	p1.insert("z", points->GetPoint(0)[2]);
	QVariantMap p2;
	p2.insert("x", points->GetPoint(1)[0]);
	p2.insert("y", points->GetPoint(1)[1]);
	p2.insert("z", points->GetPoint(1)[2]);

	pointData.insert("p1", p1);
	pointData.insert("p2", p2);
}











