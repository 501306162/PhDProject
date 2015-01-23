#include "ValvePlane.h"

#include <QFile>
#include <QFileInfo>
#include <QVariant>

#include <qjson/parser.h>


namespace vt 
{
// ------------------------------------------------------------------------
void ValvePlane::Initialise(const QVariantMap &planeData)
{
	// get the points
	QVariantList pointsData = planeData["points"].toList();
	for(int i = 0; i < pointsData.size(); i++)
	{
		QVariantMap pointData = pointsData[i].toMap();
		std::string type = pointData["type"].toString().toStdString();		

		PointType p1, p2;
		for(unsigned int j = 0; j < 3; j++)
		{
			p1[j] = pointData["P1"].toList()[j].toDouble();
			p2[j] = pointData["P2"].toList()[j].toDouble();
		}

		PointListType points;
		points.push_back(p1);
		points.push_back(p2);

		m_Points[type] = points;
	}

	for(unsigned int i = 0; i < 3; i++)
	{
		m_Center[i] = planeData["center"].toList()[i].toDouble();
		m_Normal[i] = planeData["normal"].toList()[i].toDouble();
	}
}


// ------------------------------------------------------------------------
ValvePlane::PointListType ValvePlane::GetAllPoints() const
{
	PointListType output;
	PointMapType::const_iterator mapIt = m_Points.begin();
	while(mapIt != m_Points.end())
	{
		for(unsigned int i = 0; i < mapIt->second.size(); i++)
		{
			output.push_back(mapIt->second[i]);			
		}
		++mapIt;
	}
	return output;
}

// ------------------------------------------------------------------------
ValvePlaneSequence::Pointer ValvePlaneSequence::Load(const std::string &filename)
{
	ValvePlaneSequence::Pointer plane = ValvePlaneSequence::New();
	plane->LoadFromFile(filename);
	return plane;

}

// ------------------------------------------------------------------------
void ValvePlaneSequence::LoadFromFile(const std::string &filename)
{
	QFile file(QString::fromStdString(filename));
	file.open(QIODevice::Text | QIODevice::ReadOnly);

	QFileInfo info(QString::fromStdString(filename));
	QString name = info.fileName();
	name = name.replace(".txt", "");

	m_Name = name.toStdString();

	bool ok;
	QJson::Parser parser;
	QVariantList planesData = parser.parse(&file, &ok).toList();
	file.close();

	unsigned int numPlanes = planesData.size();
	for(unsigned int i = 0; i < numPlanes; i++)
	{
		ValvePlane::Pointer plane = ValvePlane::New();
		plane->Initialise(planesData[i].toMap());		
		m_Planes.push_back(plane);
	}

	
}

}
