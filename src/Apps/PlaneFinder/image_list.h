#ifndef IMAGE_LIST_H
#define IMAGE_LIST_H

#include "common.h"
#include "containers.h"
#include <QListWidget>

class ImageListDisplay 
{
public:
	ImageListDisplay(DataContainer  * imageData);

	void setUp();

	QListWidget * getImageList() { return imageList; }
	int selectedIndex();
	

private:
	QListWidget * imageList;

};

#endif
