#include "JsonReader.h"


#include <QFile>
#include <QVariant>
#include <qjson/parser.h>

namespace vt
{
// ------------------------------------------------------------------------
JsonReader::JsonReader()
{
	m_Filename = "";
}

// ------------------------------------------------------------------------
JsonReader::LineSet JsonReader::Read() const
{

	LineSet output;
	if(m_Filename.empty())
	{
		std::cout << "Filename is not set" << std::endl;
		return output;
	}


	QFile file(QString::fromStdString(m_Filename));
	if(!file.open(QIODevice::ReadOnly))
	{
		std::cout << "Couldn't open the file" << std::endl;
		exit(1);
	}


	QJson::Parser parser;
	bool ok;
	QVariant json = parser.parse(&file, &ok);
	file.close();

		
	if(!ok)
	{
		std::cout << "Couldn't parse the json" << std::endl;
		exit(1);
	}


	QVariantMap jsonMap = json.toMap();
	std::string folder = jsonMap["folder"].toString().toStdString();

	QVariantList images = jsonMap["images"].toList();
	for (int i = 0; i < images.size(); ++i)
	{
		// get the image details
		QVariantMap image = images[i].toMap();
		std::string imageType = image["image_type"].toString().toStdString();
		std::string imageName = image["name"].toString().toStdString();
		std::string imagePath = folder + "/" + imageName + ".nrrd";

		QVariantList lineList = image["lines"].toList();

		for (int i = 0; i < lineList.size(); ++i)
		{
			QVariantMap lineMap = lineList[i].toMap();
			int timeStep = lineMap["t"].toInt();

			QVariantMap::iterator mapIt = lineMap.begin();
			while(mapIt != lineMap.end())
			{
				if(mapIt.key() != "t" && mapIt.key() != "s")
				{
					// extract the line details
					std::string lineType = mapIt.key().toStdString();
					QVariantMap points = mapIt.value().toMap();
					
					Line newLine;
					newLine.t = timeStep;
					
					QVariantMap p1 = points["p1"].toMap();
					newLine.p1[0] = p1["x"].toDouble();
					newLine.p1[1] = p1["y"].toDouble();
					newLine.p1[2] = p1["z"].toDouble();


					QVariantMap p2 = points["p2"].toMap();
					newLine.p2[0] = p2["x"].toDouble();
					newLine.p2[1] = p2["y"].toDouble();
					newLine.p2[2] = p2["z"].toDouble();

					std::string key = lineType + ":" + imageType;

					output[key].lineType = lineType;
					output[key].imageType = imageType;
					output[key].imageFilename = imagePath;
					output[key].lines.push_back(newLine);

				}

				++mapIt;
			}

		}
	}


	return output;
}

}
