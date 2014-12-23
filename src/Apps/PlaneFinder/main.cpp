#include <iostream>

#include <QApplication>

#include "common.h"
#include "main_window.h"
#include "containers.h"
#include "dicom_parser.h"

int main(int argc, char ** argv)
{
	QApplication app(argc, argv);
	/*
	DicomParser * parser = new DicomParser;


	parser->setFolderName(QString::fromStdString(argv[1]));
	parser->setUpDisplay();
	parser->show();
	*/

	MainWindow * window;

	if (argc == 3)
	{
		window = new MainWindow;
		window->loadData(argv[2]);
	}
	else if(argc > 1)
	{
		std::string inputFolder = argv[1];
		DataContainer * container = new DataContainer;
		container->LoadData(inputFolder);
		window = new MainWindow(container);
	}
	else
	{
		window =  new MainWindow;
	}

	window->show();

	return app.exec();


}
