#include "image_list.h"

#include <QListWidgetItem>

// ------------------------------------------------------------------------
ImageListDisplay::ImageListDisplay(DataContainer *imageData)
{
	imageList = new QListWidget;
	showList(imageData);
}


// ------------------------------------------------------------------------
void ImageListDisplay::showList(DataContainer * imageData)
{
	if(imageData == NULL)
		return;

	for(unsigned int i = 0; i < imageData->numImages(); i++)
	{
		std::cout << imageData->filename(i) << std::endl;
		new QListWidgetItem(
				QString::fromStdString(imageData->filename(i)),
				this->imageList);
		
	}


	// set the selection procedure
	this->imageList->setSelectionMode(QAbstractItemView::SingleSelection);
	this->imageList->setItemSelected(this->imageList->item(0), true);
}

// ------------------------------------------------------------------------
void ImageListDisplay::setUp()
{
}


// ------------------------------------------------------------------------
int ImageListDisplay::selectedIndex()
{
	if(this->imageList->selectedItems().size() == 0)
		return -1;	
	return this->imageList->row(this->imageList->selectedItems()[0]);	

}
