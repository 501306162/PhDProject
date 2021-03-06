#include "FlipChecker.h"
#include <QFile>
#include <QStringList>
#include <QString>
#include <QVariant>
#include <QFileInfo>
#include <qjson/parser.h>

namespace vt 
{
// ------------------------------------------------------------------------
FlipChecker::FlipChecker()
{
	// load the file
	QString flipFileName = "/home/om0000/Projects/PhDProjectNew/src/ValveTracking/ToFlip.json";
	QFile file(flipFileName);
	file.open(QIODevice::ReadOnly | QIODevice::Text);

	QJson::Parser parser;
	QVariantMap data = parser.parse(&file).toMap();
	BuildMap(data["coords"].toMap(), m_FlipMap);
	BuildMap(data["points"].toMap(), m_PointFlipMap);
	file.close();

}

// ------------------------------------------------------------------------
void FlipChecker::BuildMap(QVariantMap input, FlipMap &map)
{
	QStringList keys = input.keys();
	for (int i = 0; i < keys.size(); ++i)
	{
		std::string stdKey = keys[i].toStdString();
		QVariantList keyList = input[keys[i]].toList();
		FlipKeys flipKeys;

		for (int j = 0; j < keyList.size(); ++j)
		{
			map[stdKey].push_back(keyList[j].toInt());
		}
	}
}



// ------------------------------------------------------------------------
bool FlipChecker::FlipPoints(const std::string &type, const int number)
{
	if(m_PointFlipMap.count(type) == 0)
		return false;

	FlipKeys &list = m_PointFlipMap[type];
	for(unsigned int i = 0; i < list.size(); i++)
	{
		if(list[i] == number)
			return true;
	}

	return false;
}


// ------------------------------------------------------------------------
bool FlipChecker::FlipImage(const std::string &type, const int number)
{
	if(m_FlipMap.count(type) == 0)
		return false;

	FlipKeys &list = m_FlipMap[type];
	for(unsigned int i = 0; i < list.size(); i++)
	{
		if(list[i] == number)
		{	
			if(type == "MV-3C")
				return false;
			else 
				return true;
		}
	}

	if(type == "MV-3C")
		return true;
	else 
		return false;
}

// ------------------------------------------------------------------------
bool FlipChecker::FlipImageFromFileName(const std::string type, const std::string &filename)
{
	QFileInfo info(QString::fromStdString(filename));
	QString fname = info.fileName();
	fname = fname.replace("d","");
	fname = fname.replace(".txt","");

	return FlipImage(type, fname.toInt());
}

}
