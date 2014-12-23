#include <iostream>
#include <Directory.h>
#include <QString>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkNrrdImageIO.h>
#include <itkImageRegionIterator.h>
#include <itkImageFileWriter.h>


void filterFilenames(const std::vector<std::string> &filenames, std::string & filter, std::vector<std::string> &output);


typedef itk::Image<unsigned short, 4> ImageType;
bool imageSort(const ImageType::Pointer &im1, const ImageType::Pointer &im2)
{
	ImageType::PointType or1 = im1->GetOrigin();
	ImageType::PointType or2 = im2->GetOrigin();

	ImageType::DirectionType direction = im1->GetDirection();
	itk::Vector<double, 3> dir, v1, v2, normed;
   	for(unsigned int i = 0; i < 3; i++)
   	{
   		dir[i] = direction(i,2);
		v1[i] = or1[i];
		v2[i] = or2[i];
   	}

	std::cout << dir << std::endl;
	double mul1 = v1[0] / dir[0];
	double mul2 = v2[0] / dir[0];

	
	



	return (mul1<mul2);

}


int main(int argc, char ** argv)
{
	std::string folderName = argv[1];
	std::string filter = argv[2];

	typedef utils::Directory Directory;
	Directory::FilenamesType filenames = Directory::GetFiles(folderName, ".nrrd");

	Directory::FilenamesType goodFilenames;
	filterFilenames(filenames, filter, goodFilenames);

	// load the images
	std::vector<ImageType::Pointer> imageList;
	for(unsigned int i = 0; i < goodFilenames.size(); i++)
	{
		std::cout << goodFilenames[i] << std::endl;

		typedef itk::ImageFileReader<ImageType> ReaderType;

		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(goodFilenames[i]);
		reader->SetImageIO(itk::NrrdImageIO::New());
		reader->Update();

		imageList.push_back(reader->GetOutput());

	}

	std::sort(imageList.begin(), imageList.end(), imageSort);

	ImageType::Pointer outputImage = ImageType::New();
	ImageType::RegionType outputRegion;
	ImageType::SizeType outputSize = imageList.front()->GetLargestPossibleRegion().GetSize();
	outputSize[2] = imageList.size(); 
	ImageType::IndexType outputIndex;
	outputIndex.Fill(0);
	outputRegion.SetSize(outputSize);
	outputRegion.SetIndex(outputIndex);

	// compute the spacing 
	double meanSpacing = 0.0;
	for(unsigned int i = 1; i < imageList.size(); i++)
	{
		ImageType::Pointer im1 = imageList[i-1];
		ImageType::Pointer im2 = imageList[i];

		meanSpacing += (im1->GetOrigin()-im2->GetOrigin()).GetNorm();
	}

	meanSpacing /= (double) (imageList.size()-1);
	ImageType::SpacingType spacing = imageList.front()->GetSpacing();
	spacing[2] = meanSpacing;

	outputImage->SetRegions(outputRegion);
	outputImage->SetSpacing(spacing);
	outputImage->SetOrigin(imageList.front()->GetOrigin());
	outputImage->SetDirection(imageList.front()->GetDirection());
	outputImage->Allocate();

	itk::ImageRegionIterator<ImageType> outIt(outputImage, outputImage->GetLargestPossibleRegion());
	outIt.GoToBegin();

	while(!outIt.IsAtEnd())
	{
		ImageType::IndexType index = outIt.GetIndex();
		ImageType::IndexType testIndex = index;
		testIndex[2] = 0;

		outIt.Set(imageList[index[2]]->GetPixel(testIndex));


		++outIt;
	}

	std::string stackFilename = folderName + "/stack.nrrd";


	typedef itk::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(outputImage);
	writer->SetFileName(stackFilename);
	writer->SetImageIO(itk::NrrdImageIO::New());
	writer->Update();

	return 0;

}



// ------------------------------------------------------------------------
void filterFilenames(const std::vector<std::string> &filenames, 
		std::string & filter, std::vector<std::string> &output)
{
	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		QString test = QString::fromStdString(filenames[i]);
		QString query = QString::fromStdString(filter);
		if(test.contains(query))
			output.push_back(filenames[i]);
	}


}
