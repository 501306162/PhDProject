#include <iostream>
#include <ValveIO.h>
#include <Directory.h>
#include <qjson/parser.h>
#include <qjson/serializer.h>
#include <QDir>
#include <QFile>
#include <QVariant>
#include <QTextStream>

using namespace vt;

int main(int argc, char ** argv)
{
	// get the input filenames
	std::string inputFolder = argv[1];
	utils::Directory::FilenamesType filenames = utils::Directory::GetFiles(inputFolder, ".txt");

	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		// read the file 
		QFile * file = new QFile(QString::fromStdString(filenames[i]));
		file->open(QIODevice::ReadOnly | QIODevice::Text);


		QJson::Parser parser;
		bool ok;
		QVariantList list = parser.parse(file, &ok).toList();
		QVariantList outputList;
		file->close();
		delete file;

		for(int j = 0; j < list.size(); j++)
		{
			QString fname = list[j].toString();

			QDir input(fname);
			QDir output(QString::fromStdString(inputFolder));
			QString outputLine = output.absoluteFilePath(input.dirName());
			outputList.push_back(outputLine);
		}


		// save the output 
		QFile * outFile = new QFile(QString::fromStdString(filenames[i]));
		outFile->open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate);

		
		QJson::Serializer serialiser;
		bool ok2;
		QString json = serialiser.serialize(outputList, &ok2);
		QTextStream out(outFile);
		out << json;
		outFile->close();





	}





	return 0;
}
