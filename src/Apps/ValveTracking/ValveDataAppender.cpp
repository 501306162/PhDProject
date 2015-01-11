#include <iostream>

#include <Directory.h>
#include <ValveIO.h>
#include <ValveLine.h>

#include <QDir>
#include <QFileInfo>


using namespace vt;
int main(int argc, char ** argv)
{
	utils::Directory::FilenamesType filenames = utils::Directory::GetFiles(argv[1], ".txt");

	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
		reader->SetFileName(filenames[i]);
		ValveSequence<3>::Pointer sequence = reader->GetOutput();

		// create the output name
		QFileInfo info(QString::fromStdString(filenames[i]));
		QString part = info.fileName();
		part = part.replace(".txt","");

		QDir outputDir(QString::fromStdString(argv[2]));
		QString outputFile = outputDir.absoluteFilePath(part);

		ValveSequenceWriter<3>::Pointer writer = ValveSequenceWriter<3>::New();
		writer->SetInput(sequence);
		writer->SetFileName(outputFile.toStdString());
		writer->Write();
	}

	return 0;
}
