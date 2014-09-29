#include <iostream>

#include <gdcmReader.h>


#include "SeriesExtractor.h"


#include "common.h"

int main(int, char ** argv)
{
	std::string inputFolder = argv[1];
	std::cout << "Parsing DICOM files in folder: " << inputFolder << std::endl;


	SeriesExtractor extractor(inputFolder);
	extractor.ExtractSeries();
	



	return 0;
}
