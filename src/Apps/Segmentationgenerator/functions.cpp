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
#include <itkPermuteAxesImageFilter.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_inverse.h>




// ------------------------------------------------------------------------
void computeTransform(const ImageType::Pointer &image, vnl_matrix<double> &transform)
{
	vnl_matrix<double> dirMat = image->GetDirection().GetVnlMatrix();
	transform.set_size(4,4);

	for(unsigned int i = 0; i < 3; i++)
	{
		for(unsigned int j = 0; j < 3; j++)
		{
			transform(i,j) = dirMat(i,j);
		}
	}


}

// ------------------------------------------------------------------------
void buildOutputVolumes(SeriesTransform::Map &transforms)
{
	// get the short axis series (we assume that these have more than one slice)
	ImageType::Pointer saVol = ImageType::New();
	int saSeriesNumber = getShortAxisSeriesNumber(transforms);

	// build the short axis volume
	ImageType::Pointer saVolume = ImageType::New();
	buildShortAxisVolume(transforms, saSeriesNumber, saVolume);
	saVolume->Print(std::cout);


}


// ------------------------------------------------------------------------
void buildShortAxisVolume(const SeriesTransform::Map &transforms, 
		const unsigned int seriesNumber, ImageType::Pointer &saVolume)
{
	// get the short axis transforms
	SeriesTransform::Map::const_iterator mapIt = transforms.begin();
	std::vector<SeriesTransform> saSlices;
	while(mapIt != transforms.end())
	{
		if(mapIt->second.series == seriesNumber)
		{
			unsigned int sliceNum = mapIt->second.slice;

			if(saSlices.size() < (sliceNum+1))
				saSlices.resize(sliceNum+1);

			saSlices[sliceNum] = mapIt->second;
		}
		
		++mapIt;
	}


	// get the dimensions of the output image
	ImageType::Pointer reference = saSlices[0].images[0];
	unsigned int x,y,z;
	x = reference->GetLargestPossibleRegion().GetSize()[0];
	y = reference->GetLargestPossibleRegion().GetSize()[1];
	z = saSlices.size();

	ImageType::RegionType region;
	ImageType::SizeType size;
	ImageType::IndexType index;
	size[0] = x;
	size[1] = y;
	size[2] = z;

	index.Fill(0);

	region.SetSize(size);
	region.SetIndex(index);
	
	// get the other parts
	ImageType::SpacingType spacing = reference->GetSpacing();
	spacing[2] = saSlices[0].sliceThickness;
	ImageType::DirectionType direction = reference->GetDirection();
	ImageType::PointType origin = reference->GetOrigin();


	saVolume->SetOrigin(origin);
	saVolume->SetDirection(direction);
	saVolume->SetSpacing(spacing);
	saVolume->SetRegions(region);
	saVolume->Allocate();
	saVolume->FillBuffer(0);
	
}

// ------------------------------------------------------------------------
int getShortAxisSeriesNumber(SeriesTransform::Map &transforms)
{
	SeriesTransform::Map::iterator mapIt = transforms.begin();
	std::map<int, int> seriesSliceMap;
	while(mapIt != transforms.end())
	{		
		if(seriesSliceMap.count(mapIt->second.series) == 0)
		{
			seriesSliceMap[mapIt->second.series] = 1;
		}
		else
		{
			seriesSliceMap[mapIt->second.series]++;
		}
		++mapIt;
	}

	std::map<int, int>::iterator sit = seriesSliceMap.begin();
	while(sit != seriesSliceMap.end())
	{
		if(sit->second > 1)
			return sit->first;

		++sit;
	}

	// dunno, return -1 as a guess
	return -1;
}


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
	gdcm::Tag seriesDescription(0x0008,0x103e);
	gdcm::Tag sliceThickness(0x0018,0x0050);
	scanner.AddTag(seriesNumberTag);
	scanner.AddTag(sliceThickness);
	scanner.AddTag(seriesDescription);
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
				std::string seriesNumberStr, instanceNumberStr, 
					numberOfImagesStr, descriptionStr, thicknessStr;
				seriesNumberStr = mapping[fname][seriesNumberTag];
				instanceNumberStr = mapping[fname][instanceTag];
				numberOfImagesStr = mapping[fname][numberOfImagesTag];
				descriptionStr = mapping[fname][seriesDescription];
				thicknessStr = mapping[fname][sliceThickness];



				unsigned int seriesNum = QString::fromStdString(seriesNumberStr).toInt();
				unsigned int instanceNum = QString::fromStdString(instanceNumberStr).toInt();
				unsigned int numberOfImages = QString::fromStdString(numberOfImagesStr).toInt();
				double thickness = QString::fromStdString(thicknessStr).toDouble();

				descriptionStr = QString::fromStdString(descriptionStr).replace(" ","").toStdString();


				if(transforms.count(seriesNum) == 0) continue;

				// resize the image list if needed
				if(transforms[seriesNum].imageFilenames.size() != numberOfImages)
				{
					transforms[seriesNum].imageFilenames.resize(numberOfImages);
				}

				// put the image file into the image list
				transforms[seriesNum].description = descriptionStr;
				transforms[seriesNum].sliceThickness = thickness;
				transforms[seriesNum].imageFilenames[instanceNum-1] = fname;
			}	
		}
	}
}

// ------------------------------------------------------------------------
void loadImageSeries(SeriesTransform &series,
		const std::vector<int> &instanceNumbers)
{
	// resize the image vector
	series.images.resize(series.imageFilenames.size());

	typedef itk::GDCMImageIO IOType;
	typedef itk::ImageFileReader<ImageType> ReaderType;

	for(unsigned int i = 0; i < instanceNumbers.size(); i++)
	{
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetImageIO(IOType::New());
		reader->SetFileName(series.imageFilenames[instanceNumbers[i]]);
		
		try 
		{
			reader->Update();
		}
		catch(itk::ExceptionObject &e)
		{	
			std::cout << e << std::endl;
			exit(1);
		}

		series.images[instanceNumbers[i]] = reader->GetOutput();
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

			if(xml.tokenType() == QXmlStreamReader::StartElement && xml.name() == "directory")
			{
				xml.readNext();
				options.dataDirectory = xml.text().toString().toStdString();
			}

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

// ------------------------------------------------------------------------
void getMinMaxBounds(const std::vector<BoundsType> &bounds, double * minMaxs)
{
	// initialise the min maxs 
	minMaxs[0] = 1000;
	minMaxs[1] = -1000;
	minMaxs[2] = 1000;
	minMaxs[3] = -1000;
	minMaxs[4] = 1000;
	minMaxs[5] = -1000;

	for(unsigned int i = 0; i < bounds.size(); i++)
	{
		for(unsigned int j = 0; j < 4; j++)
		{
			ImageType::PointType corner = bounds[i].corners[j];
			double x = corner[0];
			double y = corner[1];
			double z = corner[2];
			
			if(x < minMaxs[0]) minMaxs[0] = x;
			if(x > minMaxs[1]) minMaxs[1] = x;
			if(y < minMaxs[2]) minMaxs[2] = y;
			if(y > minMaxs[3]) minMaxs[3] = y;
			if(z < minMaxs[4]) minMaxs[4] = z;
			if(z > minMaxs[5]) minMaxs[5] = z;
		}	
	}
}

// ------------------------------------------------------------------------
void getImageBounds(const ImageType::Pointer &image, BoundsType &bounds)
{
	// get the image size
	ImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
	
	typedef ImageType::IndexType IndexType;
	IndexType ind1, ind2, ind3, ind4;
	ind1.Fill(0); ind2.Fill(0); ind3.Fill(0); ind4.Fill(0);

	ind2[0] = size[0]-1;
	ind3[1] = size[1]-1;
	ind4[0] = size[0]-1;
	ind4[1] = size[1]-1;

	image->TransformIndexToPhysicalPoint(ind1, bounds.corners[0]);
	image->TransformIndexToPhysicalPoint(ind2, bounds.corners[1]);
	image->TransformIndexToPhysicalPoint(ind3, bounds.corners[2]);
	image->TransformIndexToPhysicalPoint(ind4, bounds.corners[3]);

}

// ------------------------------------------------------------------------
void NormaliseImage(const ImageType::Pointer &input, 
		const ImageType::Pointer &reference,	
		ImageType::Pointer &output,
		SeriesTransform &trans)
{
	// flip the x and y axis 
	typedef itk::PermuteAxesImageFilter<ImageType> FlipperType;
	FlipperType::Pointer flipper = FlipperType::New();
	itk::FixedArray<unsigned int, 3> order;
	order[0] = 1;
	order[1] = 0;
	order[2] = 2;

	flipper->SetOrder(order);
	flipper->SetInput(input);
	flipper->Update();

	output = flipper->GetOutput();

	// get the reference offset
	vnl_vector<double> refOrigin(3);
	refOrigin[0] = reference->GetOrigin()[0];// + trans.translation[0];
	refOrigin[1] = reference->GetOrigin()[1];// + trans.translation[1];
	refOrigin[2] = reference->GetOrigin()[2];// + trans.translation[2];

	vnl_matrix<double> refRot = reference->GetDirection().GetVnlMatrix();
	vnl_matrix<double> refRotInv = vnl_inverse(refRot);
	vnl_vector<double> refOffset = refRotInv * refOrigin;

	// get the image origin
	vnl_vector<double> origin(3);
	origin[0] = output->GetOrigin()[0];
	origin[1] = output->GetOrigin()[1];
	origin[2] = output->GetOrigin()[2];
	
	// apply the rotation to the origin
	origin = refRotInv * origin;

	// apply the offset to the origin
	origin = origin - refOffset;

	// apply the scaling 
	origin /= reference->GetSpacing()[0];

	// put the values into the output image
	ImageType::PointType itkOrigin;
	itkOrigin[0] = origin[0];
	itkOrigin[1] = origin[1];
	itkOrigin[2] = origin[2];

	ImageType::SpacingType spacing;
	spacing.Fill(output->GetSpacing()[0] / reference->GetSpacing()[0]);

	// get the new direction
	ImageType::DirectionType dirOut = output->GetDirection();
	dirOut = reference->GetDirection().GetInverse() * dirOut.GetVnlMatrix();
	
	output->SetSpacing(spacing);
	output->SetDirection(dirOut);
	output->SetOrigin(itkOrigin);
}
