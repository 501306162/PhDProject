#include "ExtractTrainingFeatures.h"

#include <iostream>
#include <Directory.h>
#include <ValveIO.h>
#include <FlipChecker.h>
#include <ValveNormaliser.h>
#include <PatchExtractor.h>
#include <SimpleMRFSegmenter.h>
#include <BinaryPatchFeatureExtractor.h>
#include <MatrixWriter.h>
#include <CommonDefinitions.h>
#include <LBPFeatureExtractor.h>



#include <itkCastImageFilter.h>
#include <itkCenteredAffineTransform.h>
#include <itkResampleImageFilter.h>
#include <itkPNGImageIO.h>
#include <itkImageFileWriter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkPermuteAxesImageFilter.h>
#include <itkBinaryContourImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <ConfigParser.h>

#include <QFileInfo>
#include <QString>

using namespace vt;





void writePNG(const ImageType::Pointer &patch, const std::string filename);




int main(int argc, char **argv)
{
	// parse the input configuration
	Params params(argv[1], true);
	params.print();



	// the outer loop is for the time steps 
	for(unsigned int timeStep = 20; timeStep < params.numberOfTimeSteps; timeStep++)
	{
		// first we need to iterate through each of the input folders
		for(unsigned int folderNumber = 0; folderNumber < params.valveSubDirectories.size(); folderNumber++)
		{
			utils::MatrixDataSet::Pointer dataSet = utils::MatrixDataSet::New();
			const std::string type = params.valveSubTypes[folderNumber];
			processDataSet(params, folderNumber, timeStep, dataSet);


			// save the data
			std::string filename = params.getOutputFilename(timeStep, type);

			std::cout << "Saving To: " << filename << std::endl;
			utils::MatrixWriter::Pointer writer = utils::MatrixWriter::New();
			writer->SetInput(dataSet);
			writer->SetFilename(filename);
			writer->Write();
		}
	}
	
	return 0;


	/*
	utils::Directory::FilenamesType inputFilenames = utils::Directory::GetFiles(inputDirectory, ".txt");
	FlipChecker::Pointer flipChecker = FlipChecker::New();

	std::vector<MatrixType> negativeFeatures;
	std::vector<MatrixType> positiveFeatures;

	unsigned int negativeFeatureCount = 0;
	unsigned int positiveFeatureCount = 0;
	unsigned int samplesToProcess = inputFilenames.size();
	std::vector<unsigned int> featureOwnershipList;

	for(unsigned int i = 0; i < samplesToProcess; i++)
	{

				featureOwnershipList.push_back(ownerId);


		// load the valve sequence
		ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
		reader->SetFileName(inputFilenames[i]);
		ValveSequence<3>::Pointer sequence = reader->GetOutput();
		ValveLine<3>::Pointer line = sequence->GetValveLine(0);



		// normalise the line with a possible flip
		ValveNormaliser::Pointer normaliser = ValveNormaliser::New();
		normaliser->SetInput(line);
		bool toFlip = flipChecker->FlipImage(lookup, ownerId);
		bool toFlipPoints = flipChecker->FlipPoints(lookup, ownerId);
		if(lookup == "MV-3C") toFlip = !toFlip;

		normaliser->SetFlip(toFlip);
		normaliser->SetFlipPoints(toFlipPoints);
		normaliser->Normalise();
		ValveLine<3>::Pointer aligned = normaliser->GetOutput();

		std::cout << "Extracting: " << ownerId  << " "  << toFlip <<  " " << toFlipPoints << std::endl;

		ImagePatchExtractor::Pointer testExtractor = ImagePatchExtractor::New();
		testExtractor->SetImage(aligned->GetImage());
		testExtractor->SetPatchCenter(aligned->GetIndex(pointToConsider));

		ImagePatchExtractor::SizeType testPatchSize;
		testPatchSize.Fill(featureLength);
		testExtractor->SetPatchSize(testPatchSize);

		ImageType::Pointer testPatch = testExtractor->ExtractPatch();
		writePNG(testPatch, fname.toStdString()+".png");
		
	


		// compute the segmentation and extract the contours
		MaskType::Pointer seg = MaskType::New();
		computeSegmentation(aligned, seg, segmentationPatchSize, pointToConsider);

		typedef itk::BinaryContourImageFilter<MaskType, MaskType> ContourFilterType;
		ContourFilterType::Pointer contourFilter = ContourFilterType::New();
		contourFilter->SetInput(seg);
		contourFilter->SetFullyConnected(false);
		contourFilter->SetBackgroundValue(0);
		contourFilter->SetForegroundValue(255);
		contourFilter->Update();
		MaskType::Pointer contours = contourFilter->GetOutput();


		// extract the indices that will be used for the negative samples
		RegionType negativeRegion = getNegativeFeatureRegion(aligned, negativeFeatureRegionSize, pointToConsider);

		IndexlistType negativeIndices;
		extractMaskIndices(contours, negativeRegion, aligned->GetIndex(pointToConsider), negativeIgnoreDistance, negativeIndices);

		// extract the negative features
		MatrixType features;
		negativeFeatureCount += negativeIndices.size();
		extractFeatures(aligned->GetImage(), negativeIndices, featurePatchSize, features);
		negativeFeatures.push_back(features);


		// extract the positive features
		MatrixType positiveFeature;
		extractFeature(aligned->GetImage(), aligned->GetIndex(pointToConsider), featurePatchSize, positiveFeature);
		positiveFeatureCount++;
		positiveFeatures.push_back(positiveFeature);

	}
	

	unsigned int totalFeatureCount = positiveFeatureCount+negativeFeatureCount;
	MatrixType X = MatrixType::Zero(totalFeatureCount, featureLength);
	IntMatrixType labels = IntMatrixType(totalFeatureCount, 1);
	IntMatrixType ownership = IntMatrixType(totalFeatureCount,1);


	std::cout << "Adding negative features to matrix" << std::endl;

	// put all the negative features into the matrix
	unsigned int startRow = 0;
	for(unsigned int i = 0; i < negativeFeatures.size(); i++)
	{
		// assign the data
		MatrixType &n = negativeFeatures[i];
		unsigned int nSize = n.rows();
		X.block(startRow,0,nSize,featureLength) = n;


		// put the ownership values in
		IntMatrixType owner = IntMatrixType(nSize,1);
		owner.fill(featureOwnershipList[i]);
		ownership.block(startRow,0,nSize,1) = owner;

		// set the labels
		labels.block(startRow,0,nSize,1) = IntMatrixType::Zero(nSize,1);

		// update the start row 
		startRow+= nSize;
	}


	std::cout << "Adding positive features to matrix" << std::endl;

	// put the positive features in 
	for(unsigned int i = 0; i < positiveFeatures.size(); i++)
	{
		// assign the data
		MatrixType &n = positiveFeatures[i];
		unsigned int nSize = n.rows();
		X.block(startRow,0,nSize,featureLength) = n;


		// put the ownership values in
		IntMatrixType owner = IntMatrixType(nSize,1);
		owner.fill(featureOwnershipList[i]);
		ownership.block(startRow,0,nSize,1) = owner;

		// set the labels
		labels.block(startRow,0,nSize,1) = IntMatrixType::Ones(nSize,1);

		// update the start row 
		startRow+= nSize;
	}



	utils::MatrixDataSet::Pointer dataSet = utils::MatrixDataSet::New();
	dataSet->AddData("X", X);
	dataSet->AddData("owners", ownership);
	dataSet->AddData("labels", labels);

	utils::MatrixWriter::Pointer writer = utils::MatrixWriter::New();
	writer->SetInput(dataSet);
	writer->SetFilename(outputFileName);
	writer->Write();

	*/


	return 0;

}



// ------------------------------------------------------------------------
void writePNG(const ImageType::Pointer &patch, const std::string filename)
{
	typedef itk::Image<unsigned char, 3> PNGType;
	typedef itk::RescaleIntensityImageFilter<ImageType,PNGType> RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetInput(patch);
	rescaler->SetOutputMaximum(255);
	rescaler->SetOutputMinimum(0);
	

	typedef itk::ImageFileWriter<PNGType> TestWriterType;
	TestWriterType::Pointer testWriter = TestWriterType::New();
	testWriter->SetInput(rescaler->GetOutput());
	testWriter->SetImageIO(itk::PNGImageIO::New());
	testWriter->SetFileName(filename);
	//testWriter->Update();
	

}

