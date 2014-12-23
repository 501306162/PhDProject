#ifndef IMAGE_EXPORTER_H
#define IMAGE_EXPORTER_H

#include "common.h"
#include "containers.h"


class ImageExporter
{
public:
	ImageExporter() {};
	void exportImages();
	
	void setData(DataContainer * data) { this->data = data; }
	void setFolder(const std::string &outputFolder) { this->folder = outputFolder; }

private:

	std::string folder;
	DataContainer * data;

};

#endif
