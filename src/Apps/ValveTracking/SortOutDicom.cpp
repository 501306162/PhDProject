#include <iostream>
#include <vector>

#include <QDir>

#include <SeriesExtractor.h>

void getSubFolders(const std::string &input, std::vector<std::string> &folders);


bool sort_filenames(const std::string &f1, const std::string &f2)
{
	QFileInfo inf1(QString::fromStdString(f1));
	QFileInfo inf2(QString::fromStdString(f2));


	int fn1 = inf1.fileName().replace("d","").toInt();
	int fn2 = inf2.fileName().replace("d","").toInt();

	return (fn1 < fn2);
	

}

int main(int argc, char ** argv)
{
	// get the input folders
	std::string inputFolder = argv[1];
	std::string outputFolder = argv[2];
	std::vector<std::string> folders;
	getSubFolders(inputFolder, folders);
	std::sort(folders.begin(), folders.end(), sort_filenames);

	int studies = 0;
	for(unsigned int i = 0; i < folders.size(); i++)
	{
		int studyCount = 0;
		SeriesExtractor extractor;
		extractor.SetDirectory(folders[i]);
		extractor.ExtractSeries();
		DicomSeriesList series = extractor.GetOutput();		
		std::map<std::string, int> studyMap;
		for(unsigned int j = 0; j < series.size(); j++)
		{
			std::string study = series[j].images.front().studyUID;			
			if(studyMap.count(study) == 0)
			{
				studies++;
				studyCount++;
			}
			studyMap[study] = 0;
		}

		std::cout << folders[i] << " " << studyCount << std::endl;

	}

	std::cout << studies << std::endl;
	

	return 0;
}

// ------------------------------------------------------------------------
void getSubFolders(const std::string &input, std::vector<std::string> &folders)
{
	QDir dir(QString::fromStdString(input));
	QStringList list = dir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);

	for(unsigned int i = 0; i < list.size(); i++)
	{
		std::string folderName = dir.absoluteFilePath(list[i]).toStdString();
		folders.push_back(folderName);
	}
}
