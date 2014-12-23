#include <QApplication>

#include "MainWindow.h"
#include "SeriesExtractor.h"
#include "VolumeGrouper.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkGDCMImageIO.h>
#include <itkImageFileWriter.h>
#include <itkNrrdImageIO.h>
#include <itkJoinSeriesImageFilter.h>

#include <QString>

int main(int argc, char ** argv)
{
	/*
	QApplication app(argc, argv);
	


	MainWindow window;
	window.show();

	return app.exec();

	*/

	std::string inputFolder = argv[1];
	std::string outputFolder = argv[2];
	std::cout << "Parsing DICOM files in folder: " << inputFolder << std::endl;


	SeriesExtractor extractor(inputFolder);

	std::vector<std::string> filters;
	filters.push_back("cine");
	filters.push_back("Cine");
	extractor.SetDescriptionFilters(filters);
	extractor.ExtractSeries();


	DicomSeriesList series = extractor.GetOutput();
	std::cout << series.size() << std::endl;

	VolumeGrouper grouper;
	grouper.SetInput(series);
	grouper.Group();

	/*
	
	for(unsigned int i = 0; i < series.size(); i++)
	{
		typedef itk::Image<unsigned short, 3> SliceType;
		typedef itk::Image<unsigned short, 4> ImageType;

		std::cout << "Parsing: " << series[i].description << std::endl;


		typedef itk::JoinSeriesImageFilter<SliceType, ImageType> JoinerType;
		JoinerType::Pointer joiner = JoinerType::New();

		for(unsigned int j = 0; j < series[i].images.size(); j++)
		{
			typedef itk::ImageFileReader<SliceType> ReaderType;
			ReaderType::Pointer reader = ReaderType::New();
			reader->SetFileName(series[i].images[j].filename);
			reader->SetImageIO(itk::GDCMImageIO::New());			
			reader->Update();

			joiner->SetInput(j, reader->GetOutput());
		}

		joiner->Update();


		QString outputFilename = QString::fromStdString(outputFolder + "/image_" + series[i].description);
		outputFilename += "_" + QString::number(i) + ".nrrd";
		outputFilename = outputFilename.replace(" ", "_");
		
		typedef itk::ImageFileWriter<ImageType> WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetInput(joiner->GetOutput());
		writer->SetImageIO(itk::NrrdImageIO::New());
		writer->SetFileName(outputFilename.toStdString());
		writer->Write();


		
	}

	*/



	return 0;
}
