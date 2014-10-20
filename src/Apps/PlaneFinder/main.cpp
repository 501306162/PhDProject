#include <iostream>

#include <QApplication>

#include "common.h"
#include "data_loader.h"
#include "main_window.h"

int main(int argc, char ** argv)
{
	std::string inputFolder = argv[1];

	// load the input data
	DataLoader loader;
	loader.SetDataFolder(inputFolder);
	loader.AddFilenameFilter("image*");

	ImageDataList imageData;
	loader.LoadData(imageData);


	QApplication app(argc, argv);

	MainWindow mainWindow(imageData);
	mainWindow.show();

	return app.exec();
}
