#include <QApplication>

#include "MainWindow.h"

int main(int argc, char ** argv)
{
	QApplication app(argc, argv);


	MainWindow window;
	window.show();

	return app.exec();

	/*

	std::string inputFolder = argv[1];
	std::cout << "Parsing DICOM files in folder: " << inputFolder << std::endl;


	SeriesExtractor extractor(inputFolder);

	std::vector<std::string> filters;
	filters.push_back("cine");
	filters.push_back("Cine");
	extractor.SetDescriptionFilters(filters);
	extractor.ExtractSeries();
	

	*/




	return 0;
}
