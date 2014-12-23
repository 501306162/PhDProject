#include <iostream>

#include <CardiacDicomExtractor.h>
#include <ValveOriginFinder.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkNrrdImageIO.h>

#include <itkExtractImageFilter.h>


#include <CMRFileExtractor.h>

using namespace vt;

typedef itk::Image<unsigned short, 3> ImageType;

void LoadImage(const std::string &filename, ImageType::Pointer &image);

int main(int argc, char ** argv)
{

	CMRFileExtractor::Pointer fileExtractor = CMRFileExtractor::New();
	fileExtractor->SetFolderName(argv[1]);
	fileExtractor->Extract();

	/*
	std::string c2Filename = argv[1];
	std::string c3Filename = argv[2];
	std::string saStackFilename = argv[3];

	// load the images
	ImageType::Pointer c2Image = ImageType::New();
	ImageType::Pointer c3Image = ImageType::New();
	ImageType::Pointer stackImage = ImageType::New();
	LoadImage(c2Filename, c2Image);
	LoadImage(c3Filename, c3Image);
	LoadImage(saStackFilename, stackImage);

	*/

	ValveOriginFinder::Pointer originFinder = ValveOriginFinder::New();
	originFinder->Set2CImage(fileExtractor->Get2CImage(0));
	originFinder->Set3CImage(fileExtractor->Get3CImage(0));
	originFinder->SetImageStack(fileExtractor->GetStackImage(0));
	originFinder->Compute();






	return 0;
}

// ------------------------------------------------------------------------
void LoadImage(const std::string &filename, ImageType::Pointer &image)
{
	typedef itk::Image<unsigned short, 4> InImageType;
	typedef itk::ImageFileReader<InImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename);
	reader->SetImageIO(itk::NrrdImageIO::New());
	reader->Update();


	typedef itk::ExtractImageFilter<InImageType, ImageType> ExtractorType;
	ExtractorType::Pointer extractor = ExtractorType::New();
	extractor->SetInput(reader->GetOutput());

	InImageType::RegionType exRegion = reader->GetOutput()->GetLargestPossibleRegion();
	InImageType::SizeType exSize = exRegion.GetSize();
	InImageType::IndexType exIndex = exRegion.GetIndex();

	exSize[3] = 0;
	exIndex[3] = 0;
	exRegion.SetSize(exSize);
	exRegion.SetIndex(exIndex);

	extractor->SetExtractionRegion(exRegion);
	extractor->SetDirectionCollapseToSubmatrix();
	extractor->Update();

	image = extractor->GetOutput();



}

