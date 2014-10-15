#include <iostream>
#include <CommonDefinitions.h>

#include <itkRegionOfInterestImageFilter.h>
#include <itkImageRegionConstIterator.h>

#include <QString>
#include <QStringList>
#include <QDir>
#include <QFileInfo>

#include <map>
#include <vector>

#include <dcmtk/dcmdata/dcfilefo.h>
#include <dcmtk/dcmdata/dctk.h>
#include <dcmtk/config/osconfig.h>

typedef std::map<std::string, std::vector<std::string> > LookupMap;
typedef std::vector<std::string> FilenamesType;
void loadLookupFile(const std::string &filename, LookupMap &lookup);
void readFilenames(const std::string &inputDirectory, FilenamesType &imageFilenames);

int main(int, char ** argv)
{
	std::string inputFolder = argv[1];
	std::string outputFolder = argv[2];
	std::string lookupFilename = argv[3];


	// load the lookup file and read the image filenames
	LookupMap lookup;
	FilenamesType imageFilenames;
	readFilenames(inputFolder, imageFilenames);
	loadLookupFile(lookupFilename, lookup);

	// looop through the list of images
	for(unsigned int i = 0; i < imageFilenames.size(); i++)
	{
		std::string fname = imageFilenames[i];
		QFileInfo finfo(QString::fromStdString(fname));
		std::string basename = finfo.completeBaseName().toStdString() + ".nrrd";
		FilenamesType dicomFilenames = lookup[basename];

		// slice up the input image to recreate the output
		typedef utils::ImageVolume ImageType;
		ImageType::Pointer input = utils::ImageVolumeIO::Read(fname);

		const unsigned int slices = input->GetLargestPossibleRegion().GetSize()[2];

		for(unsigned int slice = 0; slice < slices; slice++)
		{
			typedef itk::RegionOfInterestImageFilter<ImageType, ImageType> ROIFilter;
			ROIFilter::Pointer roiFilter = ROIFilter::New();
			roiFilter->SetInput(input);

			ImageType::RegionType roi = input->GetLargestPossibleRegion();
			ImageType::IndexType roiIndex = roi.GetIndex();
			ImageType::SizeType roiSize = roi.GetSize();

			roiIndex[2] = slice;
			roiSize[2] = 1;
			roi.SetSize(roiSize);
			roi.SetIndex(roiIndex);

			roiFilter->SetRegionOfInterest(roi);
			roiFilter->Update();

			ImageType::Pointer imSlice = roiFilter->GetOutput();


			// load the dicom file
			std::string dicomFilename = dicomFilenames[slice];
			DcmFileFormat fileFormat;
			fileFormat.loadFile(dicomFilename.c_str());

			unsigned int numberOfPixels = imSlice->GetLargestPossibleRegion().GetNumberOfPixels();
			unsigned short * buffer = new unsigned short[numberOfPixels];
			fileFormat.getDataset()->putAndInsertUint16Array(DCM_PixelData, buffer, imSlice->GetLargestPossibleRegion().GetNumberOfPixels());

			// reset the pixel values
			itk::ImageRegionConstIterator<ImageType> it(imSlice, imSlice->GetLargestPossibleRegion());
			unsigned int count = 0;
			while(!it.IsAtEnd())
			{
				buffer[count] = it.Get();
				++it; ++count;
			}
			
			fileFormat.getDataset()->putAndInsertUint16Array(DCM_PixelData, buffer, imSlice->GetLargestPossibleRegion().GetNumberOfPixels());
			// create the output filename 
			std::stringstream ss;
			ss << outputFolder << "/dicom_" << i << "_" << slice << ".dcm";


			fileFormat.saveFile(ss.str().c_str());
		}
	}

	return 0;
}

// ------------------------------------------------------------------------
void readFilenames(const std::string &inputDirectory, FilenamesType &imageFilenames)
{
	QDir dir(QString::fromStdString(inputDirectory));
	QStringList filters;
	filters << "label*";
	QStringList filenames = dir.entryList(filters);

	for(int i = 0; i < filenames.size(); i++)
	{
		QString fname = dir.absoluteFilePath(filenames[i]);
		imageFilenames.push_back(fname.toStdString());
	}
}

// ------------------------------------------------------------------------
void loadLookupFile(const std::string &filename, LookupMap &lookup)
{
	std::ifstream file;
	file.open(filename.c_str());

	if(!file.is_open())
	{
		std::cout << "Couldn't open the lookup file" << std::endl;
		exit(1);
	}

	while(!file.eof())
	{
		std::string line;
		std::getline(file, line);
		if(line.empty()) continue;

		QString sline = QString::fromStdString(line);
		QStringList parts = sline.split(":");

		lookup[parts[0].toStdString()].push_back(parts[2].toStdString());
	}
}
