#include <iostream>

#include <QApplication>

#include "common.h"
#include "main_window.h"
#include "containers.h"

int main(int argc, char ** argv)
{
	QApplication app(argc, argv);
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
