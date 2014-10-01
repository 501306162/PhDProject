#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QMainWindow>

#define MAINWINDOW_WIDTH 700
#define MAINWINDOW_HEIGHT 400


class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow();

	
protected:




private slots:
	void open();
	void save();
	void saveAs();


private:	

	void setUpMainWindow();

	QMenu *fileMenu;
	QAction *openAction;
	QAction *saveAction;
	QAction *saveAsAction;


};


#endif
