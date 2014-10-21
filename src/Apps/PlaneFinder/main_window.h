#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QVTKWidget.h>
#include <QtGui>

#include "common.h"
#include "image_list.h"
#include "image_viewer.h"
#include "line.h"
#include "containers.h"

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow(DataContainer * data);
	virtual ~MainWindow();

protected:
	void closeEvent(QCloseEvent * event);
	void setTimeStep(unsigned int currentSelection);

protected slots:
	void imageSelectionChanged();
	void addLinePressed();


private:
	DataContainer  * data;
	vtkImageData * getCurrentVTKImage();

	QPushButton * getCheckedValveButton();

	void createActions();
	void setUpSignals();
	void updateSliders();
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

	QPushButton * addLineButton;
	QPushButton * removeLineButton;
	QPushButton * propagateLineButton;

	// sliders 
	QSlider * tSlider;
	QSlider * zSlider;

	QButtonGroup *valveGroup;





	unsigned int timeStep;


};



#endif
