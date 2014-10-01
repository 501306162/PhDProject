#include <iostream>

#include <gdcmDirectory.h>
#include <gdcmTag.h>
#include <gdcmScanner.h>

#include <QFileInfo>
#include <QFile>

int main( int argc, char ** argv )
{
	std::string inputDirectory = argv[1];
	std::string outputDirectory = argv[2];


	// load up the files in the directory
	gdcm::Directory dir;
	dir.Load(argv[1], true);
	gdcm::Directory::FilenamesType filenames = dir.GetFilenames();
	

	// get the cine identifiers
	std::vector<std::string> filters;
	filters.push_back("Cine");
	filters.push_back("cine");


	// set up the scanning tag
	gdcm::Tag description(0x0008,0x103e);
	gdcm::Scanner scanner;
	scanner.AddTag(description);
	scanner.Scan(filenames);

	// loop through the images and check the values 
	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		const char * fname = filenames[i].c_str();
		if(!scanner.IsKey(fname)) continue;

		std::string desc = scanner.GetValue(fname, description);

		// check against the filter values 
		for(unsigned int j = 0; j < filters.size(); j++)
		{
			if(desc.find(filters[j])) break;

			QFileInfo file(fname);
			QString srcFile = QString::fromStdString(filenames[i]);
			QString trgFile = QString::fromStdString(outputDirectory);
			trgFile += "/" + file.baseName();

			if(QFile::exists(trgFile))
			{
				QFile::remove(trgFile);
			}
			QFile::copy(srcFile, trgFile);

		}
	} 
	
	return 0;

}



