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

	void loadData(const std::string &filename);
	void saveData(const std::string &filename);
	void createNewSession(const std::string &folderName);


protected:
	void closeEvent(QCloseEvent * event);
	void setTimeStep(unsigned int currentSelection);

protected slots:
	void exportVideosPressed();
	void dicomParsed(QString &folder);
	void processDicomPressed();
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
	void saveAsActionPressed();
	void imageTypePressed();




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
	QWidget * createImageTypeControls();
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
	QAction * saveAsAction;
	QAction * processDicomAction;
	QAction * patchExtractorAction;
	QAction * exportVideos;


	// set of buttons 
	QPushButton * mvButton;
	QPushButton * tpButton;
	QPushButton * avButton;

	QComboBox * lineTypeCombo;

	QPushButton * addLineButton;
	QPushButton * removeLineButton;
	QPushButton * propagateLineButton;
	QPushButton * lockLineButton;

	QComboBox * imageTypeCombo;
	QPushButton * imageTypeButton;


	// sliders 
	QSlider * tSlider;
	QSlider * zSlider;

	QButtonGroup *valveGroup;





	unsigned int timeStep;


};



#endif
