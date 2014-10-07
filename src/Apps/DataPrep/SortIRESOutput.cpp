#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <QString>
#include <QStringList>

typedef std::map<int, std::map<int, int> > LookupMap;
typedef std::vector<double> TranslationType;
typedef std::vector<double> RotationType;
typedef std::pair<TranslationType, RotationType> TransformationType;
typedef std::map<int, TransformationType> TransformationMap;


void readLookUpFile(const std::string &filename, LookupMap &map);
void sortRegistrationFile(const std::string &filename, LookupMap &lookup, TransformationMap &transforms);
void getImageIds(const QString &line, int &series, int &slice);
void getTranslation(const QString &line, TranslationType &trans);
void getRotation(const QString &line, RotationType &rot);



int main(int, char ** argv)
{
	// load the registration and lookup file
	std::string regFilename = argv[1];
	std::string lookupFilename = argv[2];
	std::string outputFilename = argv[3];

	// read the files
	LookupMap lookup;
	readLookUpFile(lookupFilename, lookup);
	
	// load the transformations
	TransformationMap transforms;
	sortRegistrationFile(regFilename, lookup, transforms);
	
	// save the output
	LookupMap::iterator mapIt = lookup.begin();


	std::ofstream output;
	output.open(outputFilename.c_str());

	

	while(mapIt != lookup.end())
	{

		std::map<int, int>::iterator mapIt2 = mapIt->second.begin();
		int seriesNumber = mapIt->first;

		while(mapIt2 != mapIt->second.end())
		{
			int sliceNumber = mapIt2->first;
			int seriesId = mapIt2->second;
			TransformationType trans = transforms[seriesId];

			std::stringstream ss;
			ss << seriesId << ":" << seriesNumber << ":" << sliceNumber << ":";
			for(unsigned int i = 0; i < 3; i++)
			{
				ss << trans.first[i] << ":";
			}

			for(unsigned int i = 0; i < 3; i++)
			{
				ss << trans.second[i] << ":";
			}

			std::string outputLine = ss.str();
			output << outputLine.substr(0,outputLine.size()-1) << std::endl;
			++mapIt2;
		}

		++mapIt;
	}


	output.close();
	return 0;
}

// ------------------------------------------------------------------------
void sortRegistrationFile(const std::string &filename, LookupMap &lookup, TransformationMap &transforms)
{
	// read the input file
	std::ifstream input;
	input.open(filename.c_str());
	if(!input.is_open())
	{
		std::cout << "couldn't open the registration file" << std::endl;
		exit(1);
	}

	// loop through the input 
	while(!input.eof())
	{
		std::string l;
		std::getline(input,l);

		if(l.empty()) continue;

		QString line = QString::fromStdString(l);
		line = line.replace("image (","");
		line = line.replace("): translation = ",""); 
		line = line.replace(" rotation = ","");
		line = line.replace("]","");

		QStringList parts = line.split("[");

		// find out the image ids 
		int series, slice, seriesNumber;
		TranslationType translation;
		RotationType rotation;
		getImageIds(parts[0], series, slice);
		seriesNumber = lookup[series][slice];
		getTranslation(parts[1], translation);
		getRotation(parts[2], rotation);


		transforms[seriesNumber] = TransformationType(translation, rotation);
	}

	input.close();


}

// ------------------------------------------------------------------------
void getTranslation(const QString &line, TranslationType &trans)
{
	QStringList parts = line.split(", ");
	trans.push_back(parts[0].toDouble());
	trans.push_back(parts[1].toDouble());
	trans.push_back(parts[2].toDouble());
}


// ------------------------------------------------------------------------
void getRotation(const QString &line, RotationType &rot)
{
	QStringList parts = line.split(", ");
	rot.push_back(parts[0].toDouble());
	rot.push_back(parts[1].toDouble());
	rot.push_back(parts[2].toDouble());
}





// ------------------------------------------------------------------------
void getImageIds(const QString &line, int &series, int &slice)
{
	QStringList parts = line.split(", ");
	series = parts[0].toInt();
	slice = parts[1].toInt();
}



// ------------------------------------------------------------------------
void readLookUpFile(const std::string &filename, LookupMap &map)
{
	// open the file
	std::ifstream input;
	input.open(filename.c_str());

	if(!input.is_open())
	{
		std::cout << "couldn't open the series lookup file" << std::endl;
		exit(1);
	}

	while(!input.eof())
	{
		std::string l;
		std::getline(input, l);

		if(l.empty()) continue;

		// extract the relevant information
		QString line = QString::fromStdString(l);
		line = line.replace("image (","");
		line = line.replace(") = ",",");

		QStringList tokens = line.split(",");
		if(tokens.size() != 3)
		{
			std::cout << "Number of numbers in the series file is not three, please check that it is valid" << std::endl;
			exit(1);
		}

		map[tokens[0].toInt()][tokens[1].toInt()]=tokens[2].toInt();
	}
	input.close();
}
