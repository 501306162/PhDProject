#include <iostream>

#include <CommonDefinitions.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkNrrdImageIO.h>
#include <itkNiftiImageIO.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkJoinSeriesImageFilter.h>


#include <QDir>
#include <QStringList>
#include <QString>

bool file_sort(const QString &f1, const QString &f2)
{
	QString t1 = f1;
	QString t2 = f2;

	int idx1 = t1.remove("phiContour0_t=").remove("_final.nii").toInt();
	int idx2 = t2.remove("phiContour0_t=").remove("_final.nii").toInt();

	return (idx1 < idx2);

}	



int main(int, char ** argv)
{
	// load the input nifti image
	std::string inputFolderName = argv[1];
	std::string outputFilename = argv[2];


	// get the set of images from the input folder
	QStringList nameFilters;
	nameFilters << "phiContour*_t=*_final.nii";

	QDir inputFolder(QString::fromStdString(inputFolderName));
	QStringList inputFiles = inputFolder.entryList(nameFilters);
	unsigned int numTimeSteps = inputFiles.size();
	qSort(inputFiles.begin(), inputFiles.end(), file_sort);

	std::cout << "--> Found " << numTimeSteps << " segmentations" << std::endl;
	std::cout << "--> Loading images and thresholding " << std::endl;


	utils::LabelVolumeList labels;


	for(unsigned int i = 0; i < numTimeSteps; i++)
	{
		std::cout << "\t" << inputFiles[i].toStdString() << std::endl;

		typedef itk::Image<float,3> InputImageType;
		typedef itk::ImageFileReader<InputImageType> ReaderType;

		std::string filename = inputFolderName + "/" + inputFiles[i].toStdString(); 

		ReaderType::Pointer reader  = ReaderType::New();
		reader->SetFileName(filename);
		reader->SetImageIO(itk::NiftiImageIO::New());
		reader->Update();

		typedef itk::BinaryThresholdImageFilter<InputImageType, utils::LabelVolume> ThresholderType;
		ThresholderType::Pointer thresholder = ThresholderType::New();
		thresholder->SetInput(reader->GetOutput());
		thresholder->SetOutsideValue(0);
		thresholder->SetInsideValue(1);
		thresholder->SetLowerThreshold(0);
		thresholder->SetUpperThreshold(255);
		thresholder->Update();


		labels.push_back(thresholder->GetOutput());
	}


	
	// build the output image 
	typedef itk::Image<unsigned char, 4> OutputType;
	OutputType::Pointer output = OutputType::New();

	typedef itk::JoinSeriesImageFilter<utils::LabelVolume, OutputType> JoinerType;
	JoinerType::Pointer joiner = JoinerType::New();
	joiner->SetSpacing(1.0);
	joiner->SetOrigin(0);

	for(unsigned int i = 0; i < labels.size(); i++)
	{
		joiner->SetInput(i, labels[i]);
	}

	joiner->Update();
	output = joiner->GetOutput();



	// save the output
	typedef itk::ImageFileWriter<OutputType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(output);
	writer->SetFileName(outputFilename);
	writer->Update();

	/*

	utils::LabelVolume::Pointer ref = labels[0];

	OutputType::SpacingType outSpacing;
	OutputType::DirectionType outDirection;
	outDirection.SetIdentity();
	OutputType::PointType outOrigin;

	for(unsigned int i = 0; i < 3; i++)
	{
		outSpacing[i] = ref->GetSpacing()[i];
		outOrigin[i] = ref->GetOrigin()[i];
		

	
	}
	
	
	



	return 0;
	



	

	
   		
	utils::LabelVolumeIO::Write(outputFilename, thresholder->GetOutput());
	*/

	return 0;

}
