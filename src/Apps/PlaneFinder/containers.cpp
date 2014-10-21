#include "containers.h"
#include "data_loader.h"

// ------------------------------------------------------------------------
void DataContainer::LoadData(const std::string & folder)
{
	DataLoader loader;
	loader.SetDataFolder(folder);
	loader.AddFilenameFilter("image*");
	loader.LoadData(data);

	currentIndex = 0;
	currentTimeStep = 0;
	currentSlice = 0;

}


// ------------------------------------------------------------------------
vtkImageSlice * DataContainer::getActor()
{
	return data[currentIndex].images[currentTimeStep][currentSlice].actor;
}

// ------------------------------------------------------------------------
vtkImageData * DataContainer::getVTKImage()
{
	return data[currentIndex].images[currentTimeStep][currentSlice].vtkImage;
}

// ------------------------------------------------------------------------
unsigned int DataContainer::maxIndex()
{
	return this->data.size() - 1;
}

// ------------------------------------------------------------------------
unsigned int DataContainer::maxTimeStep()
{
	return this->data[currentIndex].images.size() - 1;
}


// ------------------------------------------------------------------------
unsigned int DataContainer::maxSlice()
{
	return this->data[currentIndex].images[currentTimeStep].size() - 1;
}


// ------------------------------------------------------------------------
void DataContainer::setViewedImage(unsigned int index) 
{
	if(index > this->maxIndex())
		currentIndex = this->maxIndex();
	else
		currentIndex = index;
}

// ------------------------------------------------------------------------
unsigned int DataContainer::numImages()
{
	return data.size();
}

// ------------------------------------------------------------------------
unsigned int DataContainer::numTimeSteps(unsigned int index)
{
	if(index >= this->data.size())
		return 0;

	return data[index].images.size();
}


// ------------------------------------------------------------------------
std::string DataContainer::filename(unsigned int index)
{
	if(index >= this->data.size())
		return "";

	return data[index].filename;
}

// ------------------------------------------------------------------------
void DataContainer::addLine(Line *line)
{
	getCurrentHolder().lines.push_back(line);
}


// ------------------------------------------------------------------------
DataHolder & DataContainer::getCurrentHolder()
{
	return data[currentIndex].images[currentTimeStep][currentSlice];
}

// ------------------------------------------------------------------------
std::vector<vtkActor *> DataContainer::getLines()
{
	DataHolder dinst = getCurrentHolder();
	std::vector<vtkActor*> output;

	for(unsigned int i = 0; i < dinst.lines.size(); i++)
	{
		output.push_back(dinst.lines[i]->getActor());		
	}

	return output;
}

// ------------------------------------------------------------------------
unsigned int DataContainer::numSlices(unsigned int index, unsigned int timestep)
{
	if(index >= this->data.size())
		return 0;

	if(timestep >= this->data[index].images.size())
		return 0;

	return data[index].images[timestep].size();
}

// ------------------------------------------------------------------------
void DataContainer::setViewedTimeStep(unsigned int timestep) 
{
	if(timestep > this->maxTimeStep())
		currentTimeStep = this->maxTimeStep();
	else
		currentTimeStep = timestep;
}

// ------------------------------------------------------------------------
void DataContainer::setViewedSlice(unsigned int slice) 
{
	if(slice > this->maxSlice())
		currentSlice = this->maxSlice();
	else
		currentSlice = slice;
}
