#include "functions.h"

#include <QString>
#include <QStringList>
#include <QXmlStreamReader>
#include <QFile>

#include <gdcmScanner.h>
#include <gdcmTag.h>

#include <itkImageFileReader.h>
#include <itkImage.h>
#include <itkGDCMImageIO.h>


// ------------------------------------------------------------------------
void readOptionsFile(const std::string &filename, SeriesTransform::Map &transforms)
{
	std::ifstream input;
	input.open(filename.c_str());
	if(!input.is_open()) 
	{
		std::cout << "Couldn't open options" << std::endl; 
		exit(1);
	}
	
	// read in the file and split out the values 
	while(!input.eof())
	{
		std::string l;
		std::getline(input, l);
		if(l.empty()) continue;

		QString line = QString::fromStdString(l);
		QStringList tokens = line.split(":");

		// build the transform object
		SeriesTransform transform;
		transform.dcmSeries = tokens[0].toInt();
		transform.series = tokens[1].toInt();
		transform.slice = tokens[2].toInt();

		transform.translation.push_back(tokens[3].toDouble());
		transform.translation.push_back(tokens[4].toDouble());
		transform.translation.push_back(tokens[5].toDouble());
		transform.rotation.push_back(tokens[6].toDouble());
		transform.rotation.push_back(tokens[7].toDouble());
		transform.rotation.push_back(tokens[8].toDouble());

		transforms[transform.dcmSeries] = transform;
	}
}


// ------------------------------------------------------------------------
void scanDicomFiles(const std::string &folder, gdcm::Directory::FilenamesType &filenames)
{
	gdcm::Directory dir;
	dir.Load(folder, true);
	filenames = dir.GetFilenames();
}


// ------------------------------------------------------------------------
void groupDicomFiles(const gdcm::Directory::FilenamesType &filenames,
	   	SeriesTransform::Map &transforms)
{
	// set up the dicom scanner
	gdcm::Scanner scanner;
	gdcm::Tag seriesNumberTag(0x0020,0x0011);
	gdcm::Tag instanceTag(0x0020,0x0013);
	gdcm::Tag numberOfImagesTag(0x0018,0x1090);
	scanner.AddTag(seriesNumberTag);
	scanner.AddTag(instanceTag);
	scanner.AddTag(numberOfImagesTag);
	scanner.Scan(filenames);

	gdcm::Scanner::MappingType mapping = scanner.GetMappings();

	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		const char * fname = filenames[i].c_str();

		// is the filename valid
		if(mapping.count(fname) > 0)
		{
			if(mapping[fname].count(seriesNumberTag) > 0 &&
					mapping[fname].count(instanceTag) > 0 &&
					mapping[fname].count(numberOfImagesTag) > 0)
			{

				// extract the dicom header information
				std::string seriesNumberStr, instanceNumberStr, numberOfImagesStr;
				seriesNumberStr = mapping[fname][seriesNumberTag];
				instanceNumberStr = mapping[fname][instanceTag];
				numberOfImagesStr = mapping[fname][numberOfImagesTag];

				unsigned int seriesNum = QString::fromStdString(seriesNumberStr).toInt();
				unsigned int instanceNum = QString::fromStdString(instanceNumberStr).toInt();
				unsigned int numberOfImages = QString::fromStdString(numberOfImagesStr).toInt();

				// resize the image list if needed
				if(transforms[seriesNum].imageFilenames.size() != numberOfImages)
				{
					transforms[seriesNum].imageFilenames.resize(numberOfImages);
				}

				// put the image file into the image list
				transforms[seriesNum].imageFilenames[instanceNum-1] = fname;
			}	
		}
	}
}

// ------------------------------------------------------------------------
void loadImageSeries(const gdcm::Directory::FilenamesType &filenames, 
		ImageList &images)
{
	typedef itk::GDCMImageIO IOType;
	typedef itk::ImageFileReader<ImageType> ReaderType;

	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetImageIO(IOType::New());
		reader->SetFileName(filenames[i]);
		
		try 
		{
			reader->Update();
		}
		catch(itk::ExceptionObject &e)
		{	
			std::cout << e << std::endl;
			exit(1);
		}

		images.push_back(reader->GetOutput());
	}
}

// ------------------------------------------------------------------------
void applyTransform(const ImageList &inputImages,
	   	const SeriesTransform &transform,
	   	ImageList &transformedImages)
{
	// loop through the images
	for(unsigned int i = 0; i < inputImages.size(); i++)
	{
		ImageType::Pointer inputImage = inputImages[i];


	}

}

// ------------------------------------------------------------------------
void readXMLValues(const std::string &input, OptionsData &options)
{

	// open the xml file
	QFile * file = new QFile(input.c_str());
	if(!file->open(QIODevice::ReadOnly | QIODevice::Text))
	{
		std::cout << "Couldn't open the xml file" << std::endl;
		exit(1);
	}

	// read the xml 
	QXmlStreamReader xml(file);

	while(!xml.atEnd())
	{
		QXmlStreamReader::TokenType token = xml.readNext();
		
		if(token == QXmlStreamReader::StartDocument) continue;
		if(token == QXmlStreamReader::StartElement)
		{
			if(xml.name() == "volume_size")
			{
				xml.readNext();

				// loop until the end element named volume size
				while(!(xml.tokenType() == QXmlStreamReader::EndElement &&
							xml.name() == "volume_size"))
				{
					// get the x y and z values
					if(xml.tokenType() == QXmlStreamReader::StartElement && xml.name() == "x")
					{
						xml.readNext();
						options.volumeSize[0] = xml.text().toString().toInt();
					}

					if(xml.tokenType() == QXmlStreamReader::StartElement && xml.name() == "y")
					{
						xml.readNext();
						options.volumeSize[1] = xml.text().toString().toInt();
					}

					if(xml.tokenType() == QXmlStreamReader::StartElement && xml.name() == "z")
					{
						xml.readNext();
						options.volumeSize[2] = xml.text().toString().toInt();
					}

					xml.readNext();
				}

			}

			if(xml.name() == "ROI_offset")
			{
				xml.readNext();

				// loop until the end element named volume size
				while(!(xml.tokenType() == QXmlStreamReader::EndElement &&
							xml.name() == "ROI_offset"))
				{
					// get the x y and z values
					if(xml.tokenType() == QXmlStreamReader::StartElement && xml.name() == "x")
					{
						xml.readNext();
						options.roiOffset[0] = xml.text().toString().toInt();
					}

					if(xml.tokenType() == QXmlStreamReader::StartElement && xml.name() == "y")
					{
						xml.readNext();
						options.roiOffset[1] = xml.text().toString().toInt();
					}

					if(xml.tokenType() == QXmlStreamReader::StartElement && xml.name() == "z")
					{
						xml.readNext();
						options.roiOffset[2] = xml.text().toString().toInt();
					}

					xml.readNext();
				}

			}
		}
	
	}
}

