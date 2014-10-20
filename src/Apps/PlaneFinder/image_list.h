#ifndef IMAGE_LIST_H
#define IMAGE_LIST_H

#include "common.h"
#include <QListWidget>

class ImageListDisplay 
{
public:
	ImageListDisplay(ImageDataList &imageData);

	void setUp();

	QListWidget * getImageList() { return imageList; }
	int selectedIndex();
	

private:
	ImageDataList images;
	QListWidget * imageList;

};

#endif
