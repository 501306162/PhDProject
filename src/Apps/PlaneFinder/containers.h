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

	void addLine(Line * line);

	void setViewedImage(unsigned int index);
	void setViewedTimeStep(unsigned int timestep);
	void setViewedSlice(unsigned int slice);


	unsigned int maxIndex();
	unsigned int maxTimeStep();
	unsigned int maxSlice();

	unsigned int numImages();
	unsigned int numTimeSteps(unsigned int index);
	unsigned int numSlices(unsigned int index, unsigned int timestep);


	std::string filename(unsigned int index);

private:

	DataHolder & getCurrentHolder();


	DataList data;


	unsigned int currentIndex;
	unsigned int currentSlice;
	unsigned int currentTimeStep;

};

#endif
