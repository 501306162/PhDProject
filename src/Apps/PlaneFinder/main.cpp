#include <iostream>

#include <QApplication>

#include "common.h"
#include "main_window.h"
#include "containers.h"

int main(int argc, char ** argv)
{
	std::string inputFolder = argv[1];


	DataContainer * container = new DataContainer;
	container->LoadData(inputFolder);


	QApplication app(argc, argv);

	MainWindow mainWindow(container);
	mainWindow.show();

	return app.exec();
}
