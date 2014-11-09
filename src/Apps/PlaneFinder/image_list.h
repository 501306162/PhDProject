#ifndef IMAGE_LIST_H
#define IMAGE_LIST_H

#include "common.h"
#include "containers.h"
#include <QTableWidget>

class ImageListDisplay 
{
public:

	std::string getImageTypeString();


	ImageListDisplay(DataContainer  * imageData);

	void setUp();

	QTableWidget * getImageList() { return imageList; }
	int selectedIndex();
	void setData(DataContainer * data);
	void showList();
	void setSelection(unsigned int row);


	

private:

	DataContainer * data;

	QTableWidget * imageList;

};

#endif
