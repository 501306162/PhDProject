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



#include <itkCastImageFilter.h>
#include <itkCenteredAffineTransform.h>
#include <itkResampleImageFilter.h>
#include <itkPNGImageIO.h>
#include <itkImageFileWriter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkPermuteAxesImageFilter.h>
#include <itkBinaryContourImageFilter.h>
#include <itkImageRegionConstIterator.h>

#include <QFileInfo>
#include <QString>

using namespace vt;


typedef itk::Image<unsigned char, 3> MaskType;
typedef itk::Image<unsigned short, 3> ImageType;
typedef MaskType::IndexType IndexType;
typedef itk::ContinuousIndex<double, 3> ContIndexType;
typedef MaskType::RegionType RegionType;
typedef std::vector<IndexType> IndexlistType;
typedef Eigen::VectorXd FeatureType;
typedef utils::DoubleMatrixType MatrixType;
typedef utils::IntMatrixType IntMatrixType;


void computeSegmentation(const ValveLine<3>::Pointer &input, MaskType::Pointer &segmentation, 
		unsigned int segmentationPatchSize, unsigned int pointNumber);
void extractMaskIndices(const MaskType::Pointer &mask, RegionType &region, 
		const ContIndexType &positiveLocation, const double distance, IndexlistType &indices);
RegionType getNegativeFeatureRegion(const ValveLine<3>::Pointer &valve,
		const unsigned int size, const unsigned int pointtoConsider);
void extractFeatures(const MaskType::Pointer &seg, const IndexlistType &indices, 
		unsigned int featureSize, MatrixType &features);
void extractFeature(const MaskType::Pointer &mask, const ContIndexType &location, 
		unsigned int featureSize, MatrixType & feature); 
void writePNG(const ImageType::Pointer &patch, const std::string filename);


int main(int argc, char **argv)
{

	const std::string inputDirectory = argv[1];
	const std::string outputFileName = argv[2];
	const std::string lookup = argv[3];
	const unsigned int pointToConsider = atoi(argv[4]);

	const unsigned int segmentationPatchSize = 40;
	const unsigned int negativeFeatureRegionSize = 30;
	const unsigned int featurePatchSize = 10;
	const unsigned int featureLength = featurePatchSize*featurePatchSize;
	const double negativeIgnoreDistance = static_cast<double>(featurePatchSize) / 2.0;


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

		// get the id of the examination
		QFileInfo info(QString::fromStdString(inputFilenames[i]));
		QString fname = info.fileName();
		fname = fname.replace(".txt","");
		unsigned int ownerId = fname.replace("d","").toInt();
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
		testPatchSize.Fill(100);
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
		extractFeatures(seg, negativeIndices, featurePatchSize, features);
		negativeFeatures.push_back(features);


		// extract the positive features
		MatrixType positiveFeature;
		extractFeature(seg, aligned->GetIndex(pointToConsider), featurePatchSize, positiveFeature);
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



	return 0;

}

// ------------------------------------------------------------------------
void computeSegmentation(const ValveLine<3>::Pointer &input, MaskType::Pointer &segmentation, 
		unsigned int segmentationPatchSize, unsigned int pointNumber)
{

	// extract the mask that will be used for the segmentation
	ImagePatchExtractor::Pointer extractor = ImagePatchExtractor::New();
	ImagePatchExtractor::SizeType patchSize;
	patchSize.Fill(segmentationPatchSize);
	extractor->SetImage(input->GetImage());
	extractor->SetPatchSize(patchSize);

	if(pointNumber == 1)
		extractor->SetPatchCenter(input->GetInd1());
	else
		extractor->SetPatchCenter(input->GetInd2());

	MaskType::Pointer mask = extractor->ExtractMask();

	SimpleMRFSegmenter::Pointer segmenter = SimpleMRFSegmenter::New();
	segmenter->SetImage(input->GetImage());
	segmenter->SetMask(mask);
	segmenter->SetOutputValue(255);
	segmenter->SetSmoothnessCost(1.0);
	segmenter->Segment();

	segmentation = segmenter->GetOutput();
}

// ------------------------------------------------------------------------
void extractMaskIndices(const MaskType::Pointer &mask, RegionType &region, 
		const ContIndexType &positiveLocation, const double distance, IndexlistType &indices)
{
	typedef itk::ImageRegionConstIterator<MaskType> IteratorType;
	IteratorType maskIt(mask, region);

	while(!maskIt.IsAtEnd())
	{
		if(maskIt.Get() == 255)
		{
			IndexType testIndex = maskIt.GetIndex();

			// get the distance between the indices
			double sum = 0.0;
			for(unsigned int i = 0; i < 3; i++)
			{
				sum+= (positiveLocation[i]-testIndex[i])*(positiveLocation[i]-testIndex[i]);
			}
			sum = sqrt(sum);

			if(sum > distance)
				indices.push_back(maskIt.GetIndex());
		}
			
		++maskIt;
	}
}

// ------------------------------------------------------------------------
RegionType getNegativeFeatureRegion(const ValveLine<3>::Pointer &valve,
		const unsigned int size, const unsigned int pointToConsider)
{
	// get the negative feature extraction region
	ImagePatchExtractor::Pointer negativeRegionExtractor = ImagePatchExtractor::New();
	negativeRegionExtractor->SetImage(valve->GetImage());
	negativeRegionExtractor->SetPatchCenter(valve->GetIndex(pointToConsider));

	ImagePatchExtractor::SizeType negativePatchSize;
	negativePatchSize.Fill(size);
	negativeRegionExtractor->SetPatchSize(negativePatchSize);

	return negativeRegionExtractor->ExtractionRegion();
}

// ------------------------------------------------------------------------
void extractFeatures(const MaskType::Pointer &seg, const IndexlistType &indices, unsigned int featureSize, MatrixType &features)
{
	unsigned int numFeatures = featureSize*featureSize;
	features = MatrixType::Zero(indices.size(), numFeatures);

	for(unsigned int i = 0; i < indices.size(); i++)
	{
		IndexType ind = indices[i];
		ContIndexType contInd;
		for(unsigned int j = 0; j < 3; j++)
		{
			contInd[j] = ind[j];
		}

		MatrixType feature;
		extractFeature(seg, contInd, featureSize, feature);
		features.row(i) = feature.row(0);
	}

}


// ------------------------------------------------------------------------
void extractFeature(const MaskType::Pointer &mask, const ContIndexType &location, 
		unsigned int featureSize, MatrixType & feature)
{

	MaskPatchExtractor::Pointer extractor = MaskPatchExtractor::New();
	extractor->SetImage(mask);
	extractor->SetPatchCenter(location);

	MaskPatchExtractor::SizeType size;
	size.Fill(featureSize);
	extractor->SetPatchSize(size);
	MaskType::Pointer patch = extractor->ExtractPatch();


	typedef itk::ImageFileWriter<MaskType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(patch);
	writer->SetImageIO(itk::PNGImageIO::New());
	writer->SetFileName("patch.png");
	writer->Update();

	BinaryPatchFeatureExtractor::Pointer featureBuilder = BinaryPatchFeatureExtractor::New();
	featureBuilder->SetInput(patch);
	featureBuilder->Extract(feature);
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


