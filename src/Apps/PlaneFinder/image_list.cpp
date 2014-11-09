#include "image_list.h"

#include <QtGui>

// ------------------------------------------------------------------------
ImageListDisplay::ImageListDisplay(DataContainer *imageData)
{
	imageList = new QTableWidget;

	QStringList headers;
	headers << "Name";
	headers << "Image Type";

	imageList->setColumnCount(2);
	imageList->setHorizontalHeaderLabels(headers);
	imageList->setShowGrid(false);
	imageList->setMinimumHeight(200);
	imageList->resizeColumnToContents(0);


	QHeaderView * view = imageList->horizontalHeader();
	view->setResizeMode(QHeaderView::ResizeToContents);

	QHeaderView * rowView = imageList->verticalHeader();
	rowView->setDefaultSectionSize(rowView->fontMetrics().height()+6);


	imageList->setSelectionBehavior(QAbstractItemView::SelectRows);
	imageList->setSelectionMode(QAbstractItemView::SingleSelection);


	setData(imageData);


}


// ------------------------------------------------------------------------
void ImageListDisplay::setSelection(unsigned int row)
{
	this->imageList->setRangeSelected(QTableWidgetSelectionRange(row,0,row,1), true);
}


// ------------------------------------------------------------------------
void ImageListDisplay::setData(DataContainer * data)
{
	this->data = data;
	showList();
	setSelection(0);
}


// ------------------------------------------------------------------------
void ImageListDisplay::showList()
{
	if(data == NULL)
		return;


	this->imageList->clear();

	this->imageList->setRowCount(data->numImages());

	for(unsigned int i = 0; i < data->numImages(); i++)
	{
		std::cout << data->filename(i) << std::endl;
		QTableWidgetItem * item1 = new QTableWidgetItem(
				QString::fromStdString(data->filename(i)));

		QFont font = item1->font();
		font.setPointSize(10);
		item1->setFont(font);

		QTableWidgetItem * item2 = new QTableWidgetItem(
				QString::fromStdString(data->imageTypeString(i)));
		item2->setFont(font);


		this->imageList->setItem(i,0, item1);
		this->imageList->setItem(i,1, item2);
		
		
	}


	// set the selection procedure
	//this->imageList->setSelectionMode(QAbstractItemView::SingleSelection);
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
