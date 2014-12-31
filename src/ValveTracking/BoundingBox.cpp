#include "BoundingBox.h"
#include <qjson/parser.h>
#include <QFile>
#include <QVariant>

namespace vt
{
// ------------------------------------------------------------------------
void BoundingBox::Load(const std::string &filename, const int exclude)
{
	//  open the file 
	QFile file(QString::fromStdString(filename));
	file.open(QIODevice::ReadOnly | QIODevice::Text);


	// parse the input file
	QJson::Parser parser;
	bool ok;
	QVariantMap pointMap = parser.parse(&file, &ok).toMap();
	file.close();

	std::cout << pointMap.size() << std::endl;


	// get the exclude key. If it is in the map, remove it
	QString excludeKey =  "d" + QString::number(exclude);
	if(pointMap.contains(excludeKey))
	{
		pointMap.remove(excludeKey);
	}

	// start reading the points from the list
	QVariantMap::iterator mapIt = pointMap.begin();
	PointContainerType::Pointer containerP1 = PointContainerType::New();
	PointContainerType::Pointer containerP2 = PointContainerType::New();
	while(mapIt != pointMap.end())
	{
		QVariantMap points = mapIt.value().toMap();
		QVariantList p1List = points["p1"].toList();
		QVariantList p2List = points["p2"].toList();

		PointType p1, p2;

		for(unsigned int i = 0; i < 3; i++)
		{
			p1[i] = p1List[i].toDouble();			
			p2[i] = p2List[i].toDouble();			
		}

		containerP1->push_back(p1);
		containerP2->push_back(p2);

		++mapIt;
	}

	SetBoundingBox(containerP1, m_BoundingBoxP1);
	SetBoundingBox(containerP2, m_BoundingBoxP2);
}

// ------------------------------------------------------------------------
void BoundingBox::SetBoundingBox(const PointContainerType::Pointer &points, BoundingBoxType::Pointer &box)
{
	// initialise if needed
	if(!box) box = BoundingBoxType::New();
	box->SetPoints(points);
	box->ComputeBoundingBox();
}

// ------------------------------------------------------------------------
void BoundingBox::TransformBoundingBox(const TransformType::Pointer &transform)
{
	const PointContainerType * inPoints1 = m_BoundingBoxP1->GetPoints();
	const PointContainerType * inPoints2 = m_BoundingBoxP2->GetPoints();
}


}
