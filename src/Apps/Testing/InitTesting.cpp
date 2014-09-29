#include <iostream>

#include <itkImage.h>
#include <FilenamesReader.h>
#include <ImageIO.h>
#include <MatrixIO.h>
#include <MatrixReader.h>

#include <ImageToImageHeavisideDistanceMeasure.h>
#include <SimilarityMatrixBuilder.h>

#include <MatrixWriter.h>


int main(int, char *argv[])
{
	utils::FilenamesReader::FilenamesType fnames = utils::FilenamesReader::Read(argv[1], true);
	typedef itk::Image<unsigned char, 2> ImageType;
	typedef utils::ImageIO<ImageType> ImageIOType;

	typedef manifold::SimilarityMatrixBuilder<ImageType> MatrixBuilder;
	MatrixBuilder::Pointer builder = MatrixBuilder::New();


	for(unsigned int i = 0; i < fnames.size(); i++)
	{
		ImageType::Pointer image = ImageIOType::Read(fnames[i]);
		builder->AddImage(image);

	}


	
	typedef manifold::ImageToImageHeavisideDistanceMeasure<ImageType> DistanceMeasureType;
	DistanceMeasureType::Pointer measure = DistanceMeasureType::New();

	builder->SetDistanceMeasure(measure);
	builder->Update();

	utils::DoubleMatrixType mat = builder->GetOutput();

	utils::IntMatrixType dmat = utils::IntMatrixType::Zero(4, 2);
	
	int c = 0;
	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 2; j++)
		{
			dmat(i,j) = c;
			c++;
		}
	}

	std::cout << dmat << std::endl;


	utils::MatrixDataSet::Pointer dataSet = utils::MatrixDataSet::New();
	dataSet->AddData("data", mat);
	dataSet->AddData("int", dmat);

	typedef utils::MatrixWriter MatWriter;
	MatWriter::Pointer writer = MatWriter::New();
	writer->SetInput(dataSet);
	writer->SetFilename("test.mat");
	writer->Write();


	typedef utils::MatrixReader Reader;
	Reader::Pointer reader = Reader::New();
	reader->SetFilename("test.mat");
	reader->Read();


	utils::MatrixDataSet::Pointer output = reader->GetOutput();

	std::cout << output->intData["int"]<< std::endl;
	/*

	typedef utils::MatrixIO MatrixIOType;
	MatrixIOType::Pointer matrixIO = MatrixIOType::New();
	matrixIO->SetFilename("Test.mat");

	matrixIO->AddInput("data", mat);
	matrixIO->AddInput("data2", mat);
	matrixIO->Write();

	matrixIO->Read();

	*/

	return 0;
}
