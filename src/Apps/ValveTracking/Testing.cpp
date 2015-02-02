#include <iostream>
#include <MatrixCommon.h>
#include <MatrixReader.h>
#include <ValveIO.h>
#include <TrainingData.h>
#include <SVMClassifier.h>
#include <Directory.h>
#include <ParameterHelpers.h>
#include <itkPNGImageIO.h>
#include <itkImageFileWriter.h>
#include <QString>

#include <itkRescaleIntensityImageFilter.h>
#include <PatchExtractor2.h>

using namespace vt;
int main(int argc, char **argv)
{

	const std::string folderName = argv[1];
	utils::Directory::FilenamesType filenames = utils::Directory::GetFiles(folderName, ".txt");

	//PatchExtractorParameters params(argv[2]);


	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		const std::string filename = filenames[i];
		std::string fname = utils::Directory::GetFileName(filename);
		unsigned int id = QString::fromStdString(fname).replace("d","").replace(".txt","").toInt();


		ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
		reader->SetFileName(filenames[i]);
		ValveSequence<3>::Pointer sequence = reader->GetOutput();

		ValveLine<3>::Pointer line = sequence->GetValveLine(0);

		typedef ValveLine<3>::ImageType ImageType;

		typedef PatchExtractor2<ImageType> ExtractorType;
		ExtractorType::Pointer extractor = ExtractorType::New();

		ExtractorType::VectorType direction = line->GetP2()-line->GetP1();
		ExtractorType::SizeType size;
		size.Fill(300);
		size[2] = 1;

		ExtractorType::DistanceType distance;
		distance.Fill(200);
		distance[2] = 0;

		extractor->SetSize(size);
		extractor->SetDistance(distance);
		extractor->SetLine(direction);
		extractor->SetInput(line->GetImage());
		extractor->SetCenter(line->GetP1());

		ImageType::Pointer patch = extractor->GetOutput();

		typedef itk::Image<unsigned char, 3> OutputType;
		typedef itk::RescaleIntensityImageFilter<ImageType, OutputType> RescalerType;
		RescalerType::Pointer rescaler = RescalerType::New();
		rescaler->SetInput(patch);
		rescaler->SetOutputMaximum(255);
		rescaler->SetOutputMinimum(0);
		rescaler->Update();


		std::stringstream ss;
		ss << id << ".png";

		typedef itk::ImageFileWriter<OutputType> WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetInput(rescaler->GetOutput());
		writer->SetImageIO(itk::PNGImageIO::New());
		writer->SetFileName(ss.str());
		writer->Update();

	}


	return 0;
}
