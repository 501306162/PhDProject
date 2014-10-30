#ifndef CONTAINERS_H
#define CONTAINERS_H

#include "common.h"
#include <vtkActor.h>

class DataContainer
{
public:


	DataContainer() {}
	virtual ~DataContainer() {}

	void LoadData(const std::string & folder);

	vtkImageSlice * getActor();
	vtkImageData * getVTKImage();
	std::vector<vtkActor*> getLines();
	Line::Map getLineData();


	void addLine(Line * line);
	void addLine(Line * line, int index, int time, int slice);
	void setViewedImage(unsigned int index);
	void setViewedTimeStep(unsigned int timestep);
	void setViewedSlice(unsigned int slice);
	void removeLine(Line::Type type);
	void propagateLine(Line::Type type);
	void lockLine(Line::Type type);

	std::string getFolderName() { return this->folderName; }

	unsigned int maxIndex();
	unsigned int maxTimeStep();
	unsigned int maxSlice();

	unsigned int numImages();
	unsigned int numTimeSteps(unsigned int index);
	unsigned int numSlices(unsigned int index, unsigned int timestep);

	std::string imageName(unsigned int index);

	DataInstance & getInstance(unsigned int index);
	DataInstance & getInstance(const std::string &filename);

	std::string filename(unsigned int index);


	bool linesAreLocked();

	DataHolder & getCurrentHolder();
private:


	std::string folderName;

	DataList data;


	unsigned int currentIndex;
	unsigned int currentSlice;
	unsigned int currentTimeStep;

};

#endif
