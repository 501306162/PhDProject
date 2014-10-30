#ifndef IO_H
#define IO_H

#include <iostream>


#include <QVariant>
#include <QVariantList>
#include <QVariantMap>
#include <QDebug>
#include <QByteArray>

#include <vtkPoints.h>

#include "containers.h"

class IO
{
public:
	static DataContainer * Load(const std::string & filename);
	static bool Save(const std::string &filename, DataContainer * data);

	bool isItOk() { return isOk; }

private:

	IO(const std::string & filename, DataContainer * data);

	// writing methods 
	void buildData();
	void buildImageData(unsigned int index, QVariantMap &dmap);
	
	void getLines(DataHolder &holder, QVariantMap &lines);
	void getLineData(vtkPoints * points, QVariantMap &pointData);
	


	// loading methods 
	bool loadData();
	void loadImages(QVariantMap &input);
	void loadLines(QVariantList &lines, const std::string &imageName);
	void makeLine(QVariantMap &line, const std::string &imageName, Line::Type type, int t, int s);

	DataContainer * data;
	std::string filename;
	QVariantMap outputData;
	bool isOk;

};

#endif
