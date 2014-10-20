#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QVTKWidget.h>
#include <QtGui>

#include "common.h"
#include "image_list.h"
#include "image_viewer.h"

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow(ImageDataList &images);
	void setImageData(ImageDataList &images);
	virtual ~MainWindow();

protected:
	void closeEvent(QCloseEvent * event);
	void setTimeStep(unsigned int currentSelection);

protected slots:
	void imageSelectionChanged();

private:
	ImageDataList images;

	void createActions();
	void setUpSignals();
	QWidget * createLayout();
	QWidget * createImageList();
	QWidget * createButtonGroup();
	QWidget * createImageViewer();
	QWidget * createImageControls();


	ImageViewer * imageViewer;
	ImageListDisplay * imageList;


	// set of buttons 
	QPushButton * mvButton;
	QPushButton * tpButton;
	QPushButton * avButton;

	// sliders 
	QSlider * tSlider;
	QSlider * zSlider;


	unsigned int timeStep;


};



#endif
