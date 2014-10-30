#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QVTKWidget.h>
#include <QtGui>

#include "common.h"
#include "image_list.h"
#include "image_viewer.h"
#include "line.h"
#include "line_list.h"
#include "containers.h"
#include "key_controls.h"

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow(DataContainer * data);

	MainWindow();

	virtual ~MainWindow();



protected:
	void closeEvent(QCloseEvent * event);
	void setTimeStep(unsigned int currentSelection);

protected slots:
	void imageSelectionChanged();
	void sliderSelectionChanged();
	void addLinePressed();
	void removeLinePressed();
	void propagateLinePressed();
	void lineChanged();
	void tSliderRight();
	void tSliderLeft();
	void lockLine();
	void saveActionPressed();
	void newActionPressed();
	void loadActionPressed();



private:
	void setUpWindow();
	void initialiseFromData(DataContainer * data);


	DataContainer  * data;
	vtkImageData * getCurrentVTKImage();
	QPushButton * getCheckedValveButton();

	void updateAll(bool reset = false);


	bool okToSave();
	QString getSaveName();

	void createActions();
	void setUpSignals();
	void updateSliders();

	QWidget * createInfoPane();
	void setInfoPane();

	QWidget * createLayout();

	QWidget * createImageList();
	void setImageList();
	

	QWidget * createButtonGroup();
	QWidget * createImageViewer();
	QWidget * createLineList();
	QWidget * createImageControls();


	bool okToOverwrite();

	QString getNewFolder();
	QString getLoadFile();

	ImageViewer * imageViewer;
	ImageListDisplay * imageList;
	LineList * lineList;
	KeyControls * controls;

	QLabel * folderValue;


	QMenu * fileMenu;
	QAction * loadAction;
	QAction * newAction;
	QAction * saveAction;


	// set of buttons 
	QPushButton * mvButton;
	QPushButton * tpButton;
	QPushButton * avButton;

	QPushButton * addLineButton;
	QPushButton * removeLineButton;
	QPushButton * propagateLineButton;
	QPushButton * lockLineButton;


	// sliders 
	QSlider * tSlider;
	QSlider * zSlider;

	QButtonGroup *valveGroup;





	unsigned int timeStep;


};



#endif
