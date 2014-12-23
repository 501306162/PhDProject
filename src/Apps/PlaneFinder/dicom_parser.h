#ifndef DICOM_PARSER_H
#define DICOM_PARSER_H

#include <QWidget>
#include <QTableWidget>
#include <QPushButton>
#include <QVTKWidget.h>
#include <vtkImageSlice.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <QTableWidgetItem>

#include "SeriesExtractor.h"

class DicomParser : public QWidget
{
	Q_OBJECT
public:
	DicomParser(QWidget * parent = 0);
	~DicomParser() {}

	void setFolderName(const QString &folderName);
	void setOutputFolder(const QString &folderName);
	void setUpDisplay();

	void getSelectedSeries(DicomSeriesList &series );
	void writeSeries(const DicomSeries & series, const QString &outputFolder);

signals:
	void saveDone(QString &outputFolder);

public slots:
	void selectItem(QTableWidgetItem * item);
	void okPressed();
	void cancelPressed();
	double posString(DicomImage &i);

private:

	void loadImages();

	void parseDicom();
	QWidget * createTable();
	QWidget * createViewer();

	QString folderName;
	QString outputFolder;
	QTableWidget * table;
	DicomSeriesList series;
	QVTKWidget * viewer;

	DicomSeriesList seriesAgain;

	vtkSmartPointer<vtkRenderer> renderer;

	std::vector<vtkSmartPointer<vtkImageSlice> > images;

	QPushButton * okButton;
	QPushButton * cancelButton;
};


#endif
