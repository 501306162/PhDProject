#include <iostream>
#include <map>
#include <fstream>
#include <QString>

typedef std::map<int, std::map<int, int> > LookupMap;
void readLookUpFile(const std::string &filename, LookupMap &map);

int main(int argc, char ** argv)
{
	// load the registration and lookup file
	std::string regFile = argv[1];
	std::string lookupFilename = argv[2];
	
	
	// read the files
	LookupMap lookup;
	readLookUpFile(lookupFilename, lookup);


}



void readLookUpFile(const std::string &filename, LookupMap &map)
{
	std::ifstream input;
	input.open(filename.c_str());

	while(!input.eof())
	{
		std::string l;
		std::getline(input, l);

		if(l.empty()) continue;

		QString line = QString::fromStdString(l);

		std::cout << line.toStdString() << std::endl;

	}

}
