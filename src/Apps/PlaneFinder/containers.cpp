#include "containers.h"
#include "data_loader.h"


#include <QFileInfo>


// ------------------------------------------------------------------------
DataContainer::DataContainer()
{
	saveName = "";
}


// ------------------------------------------------------------------------
bool DataContainer::hasSaveName()
{
	if(saveName.empty())
		return false;
	return true;
}

// ------------------------------------------------------------------------
DataInstance::ImageType DataContainer::imageType(unsigned int index)
{
	return data[index].imageType;
}

// ------------------------------------------------------------------------
std::string DataContainer::imageTypeString(unsigned int index)
{
	return data[index].getImageTypeString();
}

// ------------------------------------------------------------------------
void DataContainer::setImageType(unsigned int index, const std::string &type)
{
	data[index].imageType = DataInstance::getType(type);
}


// ------------------------------------------------------------------------
void DataContainer::LoadData(const std::string & folder)
{
	DataLoader loader;
	loader.SetDataFolder(folder);
	loader.AddFilenameFilter("*.nrrd");
	loader.LoadData(data);

	currentIndex = 0;
	currentTimeStep = 0;
	currentSlice = 0;


	QFileInfo info(QString::fromStdString(folder));
	this->folderName = folder;

}


// ------------------------------------------------------------------------
bool DataContainer::linesAreLocked()
{
	for(unsigned int i = 0; i < numImages(); i++)
	{
		for(unsigned int j = 0; j < numTimeSteps(i); j++)
		{
			for(unsigned int k = 0; k < numSlices(i,j); k++)
			{
				Line::Map &lines = data[i].images[j][k].lines;
				Line::Map::iterator mapIt = lines.begin();
				while(mapIt != lines.end())
				{
					if(!mapIt->second->isLocked())
						return false;
					++mapIt;
				}
			}
			
		}		
	}
	return true;
}

// ------------------------------------------------------------------------
std::string DataContainer::imageName(unsigned int index)
{
	return data[index].filename;
}


// ------------------------------------------------------------------------
DataInstance & DataContainer::getInstance(unsigned int index)
{
	return data[index];
}

// ------------------------------------------------------------------------
void DataContainer::lockLine(Line::Type type)
{
	DataHolder &holder = getCurrentHolder();
	Line * line = holder.lines[type];

	line->setLocked(!line->isLocked());
}

// ------------------------------------------------------------------------
void DataContainer::propagateLine(Line::Type lineType)
{
	DataHolder &holder = getCurrentHolder();
	Line * line = holder.lines[lineType];

	// loop through all timesteps of this sequence 
	DataInstance &sequence = this->data[currentIndex];
	for(unsigned int i = currentTimeStep+1; i < sequence.images.size(); i++)
	{
		// create and add a copy of the line
		Line * newLine = line->copy();	
		Line::Map &lmap = sequence.images[i][currentSlice].lines;
		
		if(lmap.count(newLine->getType()) > 0)
			delete lmap[newLine->getType()];

		sequence.images[i][currentSlice].lines[lineType] = newLine;
	}
	
}

// ------------------------------------------------------------------------
Line::Map DataContainer::getLineData()
{
	DataHolder holder = getCurrentHolder();
	return holder.lines;
}

// ------------------------------------------------------------------------
void DataContainer::removeLine(Line::Type lineType)
{
	DataHolder & holder = getCurrentHolder();
	Line * line = holder.lines[lineType];

	// loop through all timesteps of this sequence 
	DataInstance &sequence = this->data[currentIndex];
	for(unsigned int i = 0; i < sequence.images.size(); i++)
	{
		Line::Map &lmap = sequence.images[i][currentSlice].lines;
		
		if(lmap.count(line->getType()) > 0)
			lmap.erase(line->getType());
	}
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
	addLine(line, currentIndex, currentTimeStep, currentSlice);
}


// -----------------------------------------------------------------------
void DataContainer::addLine(Line * line, int index, int time, int slice)
{
	data[index].images[time][slice].lines[line->getType()] = line;
}


// ------------------------------------------------------------------------
DataHolder & DataContainer::getCurrentHolder()
{
	return data[currentIndex].images[currentTimeStep][currentSlice];
}


// ------------------------------------------------------------------------
DataInstance & DataContainer::getInstance(const std::string &filename)
{
	for(unsigned int i = 0; i < data.size(); i++)
	{
		if(filename.compare(data[i].filename) == 0)
			return data[i];		
	}

	return data.front();
}


// ------------------------------------------------------------------------
std::vector<vtkActor *> DataContainer::getLines()
{
	DataHolder dinst = getCurrentHolder();
	std::vector<vtkActor*> output;

	Line::Map::iterator it = dinst.lines.begin();
	while(it != dinst.lines.end())
	{
		output.push_back(it->second->getActor());		
		++it;
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
