#include <iostream>
#include <Directory.h>
#include <JsonReader.h>
#include <JsonConverter.h>
#include <ValveIO.h>

#include <QFile>
#include <QVariant>
#include <qjson/parser.h>
#include <QDir>

#include <itkImage.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>



int main(int argc, char ** argv)
{
	// load the line data
	typedef utils::Directory::FilenamesType FilenamesType;
	FilenamesType filenames = utils::Directory::GetFiles(argv[1], ".json");

	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		std::string filename = filenames[i];		
		std::cout << filename << std::endl;

		QString inFilename = QString::fromStdString(filename);
		inFilename = inFilename.replace(QString::fromStdString(argv[1]),"");
		inFilename = inFilename.replace("/","");
		inFilename = inFilename.replace(".json","");
		std::cout << inFilename.toStdString() << std::endl;


		vt::JsonReader::Pointer reader = vt::JsonReader::New();
		reader->SetFilename(filename);
		typedef vt::JsonReader::LineSet LineSet;
		LineSet data = reader->Read();
		
		std::string partName = inFilename.toStdString();

		LineSet::iterator lineIt = data.begin();
		while(lineIt != data.end())
		{
			QString type = QString::fromStdString(lineIt->first);
			type = type.replace(":","-");

			QString output = QString::fromStdString(argv[2]);
			QDir outputDir(output);

			QString outputPath = outputDir.absoluteFilePath(type);
			std::cout << outputPath.toStdString() << std::endl;

			if(!outputDir.exists(outputPath))
				outputDir.mkdir(outputPath);
			

			QString fileName = QDir(outputPath).absoluteFilePath(QString::fromStdString(partName));



			vt::JsonConverter::Pointer converter = vt::JsonConverter::New();
			converter->SetInput(lineIt->second);
			converter->Convert();

			vt::ValveSequenceWriter<3>::Pointer writer = vt::ValveSequenceWriter<3>::New();
			writer->SetInput(converter->GetOutput());
			writer->SetFileName(fileName.toStdString());
			writer->Write();

			std::cout << lineIt->first << std::endl;


			++lineIt;
		}	






	}


	return 0;

}




