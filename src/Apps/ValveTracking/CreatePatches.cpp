#include <iostream>
#include <ParameterHelpers.h>
#include <Directory.h>
#include <ValveIO.h>
#include <PatchExtractor2.h>
#include <SimpleMRFSegmenter.h>
#include <PatchExtractor.h>
#include <itkBinaryContourImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkPNGImageIO.h>
#include <QString>
#include <itkImageRandomIteratorWithIndex.h>



using namespace vt;
typedef PatchExtractorParameters ParamsType;

int getInd(const std::string &filename);
int getInd2(const std::string &filename);
typedef ValveLine<3>::ImageType ImageType;
typedef itk::Image<unsigned char, 3> OutputType;
typedef itk::Image<unsigned char, 3> LabelType;
typedef ValveLine<3> ValveType;
typedef ValveType::PointType PointType;
typedef ImageType::RegionType RegionType;

void computeSegmentation(const ValveLine<3>::Pointer &input, LabelType::Pointer &segmentation, 
		unsigned int segmentationPatchSize, unsigned int pointNumber);
void saveFilenames(const std::vector<std::string> &filenames, const std::string &filename, bool pos=false);

OutputType::Pointer extractPatch(ParamsType &params, OutputType::Pointer &image, ValveLine<3>::Pointer &line, PointType &loc);
std::vector<PointType> getNegativeIndices(ParamsType &params, ValveLine<3>::Pointer &valve, const unsigned int pointId);
RegionType getNegativeFeatureRegion(const ValveLine<3>::Pointer &valve,
		const unsigned int size, const unsigned int pointtoConsider);
void extractMaskIndices(const LabelType::Pointer &mask, RegionType &region, 
		const PointType &positiveLocation, const double distance, std::vector<PointType> &points);


int main(int, char ** argv)
{
	// load the parameters 
	ParamsType params(argv[1]);
	const std::string directory = argv[2];
	const unsigned int timeStep = atoi(argv[3]);
	const unsigned int pointNumber = atoi(argv[4]);
	const std::string outputDirectory = argv[5];


	// create the output directory if it doesn exist
	utils::Directory::MkDir(outputDirectory);


	// get the files in the directory 
	utils::Directory::FilenamesType filenames = utils::Directory::GetFiles(directory, ".txt");
	std::vector<std::string> outputFilenames;
	std::vector<std::pair<int, std::vector<std::string> > > positiveFilenameMap;
	std::vector<std::pair<int, std::vector<std::string> > > negativeFilenameMap;
	
	for(unsigned int fnum = 0; fnum < filenames.size(); fnum++)
	{
		const std::string filename = filenames[fnum];
		int id = getInd(filename);

		// create the output directory 
		std::stringstream ssDir;
		ssDir << id;
		std::string outputSub = utils::Directory::GetPath(outputDirectory, ssDir.str());
		utils::Directory::MkDir(outputSub);


		// load the file
		ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
		reader->SetFileName(filename);		

		ValveSequence<3>::Pointer sequence = reader->GetOutput();
		ValveLine<3>::Pointer line = sequence->GetValveLine(timeStep);

		typedef itk::RescaleIntensityImageFilter<ImageType, OutputType> RescalerType;
		RescalerType::Pointer rescaler = RescalerType::New();
		rescaler->SetInput(line->GetImage());
		rescaler->SetOutputMaximum(255);
		rescaler->SetOutputMinimum(0);
		rescaler->Update();
		OutputType::Pointer rescaled = rescaler->GetOutput();



		// extract the positive patch
		PointType location = line->GetPoint(pointNumber);
		OutputType::Pointer patch = extractPatch(params, rescaled, line, location);


		// get the negative locations 
		std::vector<PointType> negativePoints = getNegativeIndices(params, line, pointNumber);
		std::vector<OutputType::Pointer> outputImages;
		std::vector<std::string> negFilenames, posFilenames;
		for(unsigned int ind = 0; ind < negativePoints.size(); ind++)
		{
			OutputType::Pointer negPatch = extractPatch(params, rescaled, line, negativePoints[ind]);
			outputImages.push_back(negPatch);
			std::stringstream ss;
			ss << "neg" << ind << ".png";
			std::string outputFilename = utils::Directory::GetPath(outputSub, ss.str());

			typedef itk::ImageFileWriter<OutputType> WriterType;
			WriterType::Pointer writer = WriterType::New();
			writer->SetInput(negPatch);
			writer->SetImageIO(itk::PNGImageIO::New());
			writer->SetFileName(outputFilename);
			negFilenames.push_back(outputFilename);
			writer->Update();

		}


		std::stringstream ss;
		ss << "pos" << 0 << ".png";
		std::string outputFilename = utils::Directory::GetPath(outputSub, ss.str());
		outputFilenames.push_back(outputFilename);


		typedef itk::ImageFileWriter<OutputType> WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetInput(patch);
		writer->SetImageIO(itk::PNGImageIO::New());
		writer->SetFileName(outputFilename);
		posFilenames.push_back(outputFilename);
		writer->Update();



		// save the filenames files 
		saveFilenames(posFilenames, utils::Directory::GetPath(outputSub, "pos.dat"));
		saveFilenames(negFilenames, utils::Directory::GetPath(outputSub, "neg.dat"));


		std::pair<int, std::vector<std::string> > currentPositiveFilenames(id, posFilenames);
		std::pair<int, std::vector<std::string> > currentNegativeFilenames(id, negFilenames);

		positiveFilenameMap.push_back(currentPositiveFilenames);
		negativeFilenameMap.push_back(currentNegativeFilenames);
	}


	// create the training dat files
	for(unsigned int fnum = 0; fnum < filenames.size(); fnum++)
	{
		const std::string filename = filenames[fnum];
		int id = getInd(filename);

		std::vector<std::string> posFilenames;
		std::vector<std::string> negFilenames;


		for(unsigned int i = 0; i < positiveFilenameMap.size(); i++)
		{
			if(id == positiveFilenameMap[i].first) continue;			

			for(unsigned int j = 0; j < positiveFilenameMap[i].second.size(); j++)
			{
				posFilenames.push_back(positiveFilenameMap[i].second[j]);
			}

			for(unsigned int j = 0; j < negativeFilenameMap[i].second.size(); j++)
			{
				negFilenames.push_back(negativeFilenameMap[i].second[j]);
			}
		}

		std::stringstream ss;
		ss << id;

		// save the filename lists
		saveFilenames(posFilenames, utils::Directory::GetPath(outputDirectory, "pos-" + ss.str() + ".dat"), true);
		saveFilenames(negFilenames, utils::Directory::GetPath(outputDirectory, "neg-" + ss.str() + ".dat"));
		
	}


	

}

// ------------------------------------------------------------------------
OutputType::Pointer extractPatch(ParamsType &params, OutputType::Pointer &image, ValveLine<3>::Pointer &line, PointType &loc)
{
	typedef PatchExtractor2<OutputType> ExtractorType;
	typedef ExtractorType::VectorType VectorType;
	typedef ExtractorType::DistanceType DistanceType;
	typedef ExtractorType::SizeType SizeType;

	SizeType size;
	size.Fill(params.featurePatchSize);
	size[2] = 1;

	DistanceType distance;
	distance.Fill(params.featurePatchDistance);
	distance[2] = 0.0;

	VectorType direction = line->GetP2()-line->GetP1();
	direction.Normalize();

	ExtractorType::Pointer extractor = ExtractorType::New();
	extractor->SetLine(direction);
	extractor->SetSize(size);
	extractor->SetDistance(distance);
	extractor->SetInput(image);
	extractor->SetCenter(loc);
	extractor->Update();

	return extractor->GetOutput();
}


// ------------------------------------------------------------------------
int getInd(const std::string &filename)
{
	return QString::fromStdString(utils::Directory::GetFileName(filename)).replace("d","").replace(".txt","").toInt();
}

int getInd2(const std::string &filename)
{

	return QString::fromStdString(utils::Directory::GetFileName(filename)).replace("p","").replace(".ng","").toInt();

}


// ------------------------------------------------------------------------
void computeSegmentation(const ValveLine<3>::Pointer &input, LabelType::Pointer &segmentation, 
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

	LabelType::Pointer mask = extractor->ExtractMask();

	SimpleMRFSegmenter::Pointer segmenter = SimpleMRFSegmenter::New();
	segmenter->SetImage(input->GetImage());
	segmenter->SetMask(mask);
	segmenter->SetOutputValue(255);
	segmenter->SetSmoothnessCost(1.0);
	segmenter->Segment();

	segmentation = segmenter->GetOutput();
}

// ------------------------------------------------------------------------
std::vector<PointType> getNegativeIndices(ParamsType &params, ValveLine<3>::Pointer &valve, const unsigned int pointId)
{
	// compute the segmentation and extract the contours
	LabelType::Pointer segmentation = LabelType::New();
	computeSegmentation(valve, segmentation, params.segmentationPatchSize, pointId);

	typedef itk::BinaryContourImageFilter<LabelType, LabelType> ContourFilterType;
	ContourFilterType::Pointer contourFilter = ContourFilterType::New();
	contourFilter->SetInput(segmentation);
	contourFilter->SetFullyConnected(false);
	contourFilter->SetBackgroundValue(0);
	contourFilter->SetForegroundValue(255);
	contourFilter->Update();
	LabelType::Pointer contours = contourFilter->GetOutput();


	// extract the indices that will be used for the negative samples
	typedef ImageType::RegionType RegionType;
	RegionType negativeRegion = getNegativeFeatureRegion(valve, params.negativeFeatureRegionSize, pointId);
	std::vector<PointType> negativeIndices;
	extractMaskIndices(contours, negativeRegion, valve->GetPoint(pointId), params.negativeIgnoreDistance, negativeIndices);

	return negativeIndices;
}


// ------------------------------------------------------------------------
void extractMaskIndices(const LabelType::Pointer &mask, RegionType &region, 
		const PointType &positiveLocation, const double distance, std::vector<PointType> &indices)
{
	typedef itk::ImageRandomIteratorWithIndex<LabelType> IteratorType;
	IteratorType maskIt(mask, region);

	unsigned int N = 400;
	maskIt.SetNumberOfSamples(N);
	while(!maskIt.IsAtEnd())
	{
		//if(maskIt.Get() == 255)
		//{
			ImageType::IndexType testIndex = maskIt.GetIndex();
			PointType point;
			mask->TransformIndexToPhysicalPoint(testIndex, point);

			// get the distance between the indices
			double dist = point.EuclideanDistanceTo(positiveLocation);

			if(dist > distance)
				indices.push_back(point);
		//}
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
void saveFilenames(const std::vector<std::string> &filenames, const std::string &filename, bool pos)
{
	std::ofstream file;
	file.open(filename.c_str());
	for(unsigned int i = 0; i < filenames.size(); i++)
	{
		if(pos)
		{
			file << filenames[i] << " " << 1 << " " << 0 << " " << 0 << " " << 20 << " " << 20 << "\n";
		}
		else
		{
			file << filenames[i] << "\n";
		}
	}
	file.close();
}



