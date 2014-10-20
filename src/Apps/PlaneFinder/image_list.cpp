#include "image_list.h"

#include <QListWidgetItem>

// ------------------------------------------------------------------------
ImageListDisplay::ImageListDisplay(ImageDataList &imageData)
{
	this->images = imageData;
}

// ------------------------------------------------------------------------
void ImageListDisplay::setUp()
{
	imageList = new QListWidget;
	for(unsigned int i = 0; i < this->images.size(); i++)
	{
		std::cout << this->images[i].filename << std::endl;
		new QListWidgetItem(
				QString::fromStdString(this->images[i].filename),
				this->imageList);
		
	}


	// set the selection procedure
	this->imageList->setSelectionMode(QAbstractItemView::SingleSelection);
	this->imageList->setItemSelected(this->imageList->item(0), true);
}


// ------------------------------------------------------------------------
int ImageListDisplay::selectedIndex()
{
	return this->imageList->row(this->imageList->selectedItems()[0]);	

}
