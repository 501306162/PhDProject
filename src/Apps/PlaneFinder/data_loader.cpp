#include "data_loader.h"

#include <QDir>
#include <QString>
#include <QStringList>

#include <itkNrrdImageIO.h>
#include <itkImageToVTKImageFilter.h>
#include <itkImageFileReader.h>
#include <itkFlipImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkRegionOfInterestImageFilter.h>

#include <vtkImageProperty.h>


// ------------------------------------------------------------------------
bool DataLoader::LoadData(DataList &imageData)
{
	FilenamesType filenames;
	GetFilenames(filenames);

	if(filenames.empty())
	{
		std::cout << "No valid image found, exiting" << std::endl;
		exit(1);
	}

	// loop through the filenames and load the images
	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		FilenameType fname = filenames[i];			
		DataInstance image;

		LoadImage(fname, image);
		image.filename = GetFileBasename(fname);

		imageData.push_back(image);
	}


	return true;

}


// ------------------------------------------------------------------------
void DataLoader::AddFilenameFilter(const std::string &filter)
{
	this->filters << QString::fromStdString(filter);
}


// ------------------------------------------------------------------------
std::string DataLoader::GetFileBasename(const std::string &filename)
{
	QFileInfo info(QString::fromStdString(filename));
	return info.baseName().toStdString();
}


// ------------------------------------------------------------------------
void DataLoader::LoadImage(const std::string &filename, DataInstance &image)
{
	typedef itk::Image<unsigned short, 4> ImageType;
	typedef itk::NrrdImageIO ImageIOType;
	typedef itk::ImageFileReader<ImageType> ReaderType;
	typedef itk::FlipImageFilter<ImageType> FlipperType;

	// load
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename);
	reader->SetImageIO(ImageIOType::New());

	// flip 
	FlipperType::Pointer flipper = FlipperType::New();
	FlipperType::FlipAxesArrayType axes;
	axes[0] = false;
	axes[1] = true;
	axes[2] = false;
	axes[3] = false;

	flipper->SetFlipAxes(axes);
	flipper->SetInput(reader->GetOutput());
	flipper->Update();

	ImageType::Pointer reference = flipper->GetOutput();
	ImageType::RegionType refRegion = reference->GetLargestPossibleRegion();
	const unsigned int timeSteps = refRegion.GetSize()[3];
	const unsigned int slices = refRegion.GetSize()[2];

	for(unsigned int i = 0; i < timeSteps; i++)
	{

		typedef itk::Image<unsigned short, 3> TimeStepType;
		typedef itk::ExtractImageFilter<ImageType, TimeStepType> ExtractorType;
		typedef itk::ImageToVTKImageFilter<TimeStepType> VTKFilterType;

		// create the extraction region 
		ImageType::RegionType exRegion;
		ImageType::SizeType exSize = refRegion.GetSize();
		exSize[3] = 0;
		ImageType::IndexType exIndex = refRegion.GetIndex();
		exIndex[3] = i;

		exRegion.SetSize(exSize);
		exRegion.SetIndex(exIndex);

		ExtractorType::Pointer extractor = ExtractorType::New();
		extractor->SetInput(reference);
		extractor->SetExtractionRegion(exRegion);
		extractor->SetDirectionCollapseToSubmatrix();
		extractor->Update();
		
		VolumeData imageVolume;



		for(unsigned int j = 0; j < slices; j++)
		{
			typedef itk::RegionOfInterestImageFilter<TimeStepType, TimeStepType> ROIFilter;
			ROIFilter::Pointer roiFilter = ROIFilter::New();
			roiFilter->SetInput(extractor->GetOutput());

			TimeStepType::RegionType roiRegion = extractor->GetOutput()->GetLargestPossibleRegion();
			TimeStepType::SizeType roiSize = roiRegion.GetSize();
			TimeStepType::IndexType roiIndex = roiRegion.GetIndex();
			roiSize[2] = 1;
			roiIndex[2] = j;
			roiRegion.SetSize(roiSize);
			roiRegion.SetIndex(roiIndex);

			roiFilter->SetRegionOfInterest(roiRegion);


			VTKFilterType::Pointer vtkFilter = VTKFilterType::New();
			vtkFilter->SetInput(roiFilter->GetOutput());
			vtkFilter->Update();

			// build the data
			DataHolder data;
			data.vtkImage = vtkSmartPointer<vtkImageData>::New();
			data.vtkImage->DeepCopy(vtkFilter->GetOutput());
			data.vtkImage->Print(std::cout);
			data.itkImage = roiFilter->GetOutput();
			data.timestep = i;
			data.slice = j;

			data.mapper = vtkSmartPointer<vtkImageSliceMapper>::New();
			data.mapper->SetInputData(data.vtkImage);

			data.actor = vtkSmartPointer<vtkImageSlice>::New();
			data.actor->SetMapper(data.mapper);

			data.actor->GetProperty()->SetColorLevel(177.48);
			data.actor->GetProperty()->SetColorWindow(512.04);


			imageVolume.push_back(data);
			
		}

		image.images.push_back(imageVolume);
	}
}




// ------------------------------------------------------------------------
void DataLoader::GetFilenames(FilenamesType &filenames)
{
	QDir dir(QString::fromStdString(folderName));
	//this->filters << "*.nrrd";
	QStringList files = dir.entryList(this->filters, QDir::Files);


	for(int i = 0; i < files.size(); i++)
	{
		QString qfullPath = dir.absoluteFilePath(files[i]);
		FilenameType fullPath = qfullPath.toStdString();
		filenames.push_back(fullPath);
	}
}
