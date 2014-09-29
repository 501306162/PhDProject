#include <iostream>

#include <ManifoldShapeModel.h>
#include <CommonDefinitions.h>
#include <FilenamesReader.h>
#include <DiffusionMap.h>
#include <ImageToImageHeavisideDistanceMeasure.h>
#include <SimilarityMatrixBuilder.h>
#include <ImageIO.h>
#include <MatrixReader.h>

#include <itkSimilarity3DTransform.h>
#include <MRFLabelAlignment.h>
#include <MRFGMMDataTerm.h>
#include <MRFSegmentation.h>
#include <MRFIsingSmoothnessTerm.h>
#include <MRFPriorGenerator.h>

#include <itkTimeProbe.h>


using namespace manifold;
using namespace utils;
int main(int argc, char *argv[])
{

	std::string labelFilenameFilename = argv[1];
	std::string distanceMapFilenamesFilename = argv[2];
	std::string imageFilenamesFilename = argv[3];
	std::string unalignedLabelFilenamesFilename = argv[4];
	std::string meanLabelFilename = argv[5];
	unsigned int imageToTest = atoi(argv[6]);	
	std::string gmmParamsFilename = argv[7];
	double priorSpreadValue = atof(argv[8]);
	double priorWeightValue = atof(argv[9]);
	double smoothnessWeightValue = atof(argv[10]);
	double entropyWeight = atof(argv[11]);



	// load up all the data
	FilenamesReader::FilenamesType labelFilenames = FilenamesReader::Read(labelFilenameFilename);
	FilenamesReader::FilenamesType distanceMapFilenames = FilenamesReader::Read(distanceMapFilenamesFilename);
	FilenamesReader::FilenamesType imageFilenames = FilenamesReader::Read(imageFilenamesFilename);
	FilenamesReader::FilenamesType unalignedLabelFilenames = FilenamesReader::Read(unalignedLabelFilenamesFilename);
	LabelVolume::Pointer           meanLabel = LabelVolumeIO::Read(meanLabelFilename);

	ImageVolume::Pointer testImage;
	LabelVolume::Pointer testLabel;

	std::vector<LabelVolume::Pointer> labels;
	std::vector<RealVolume::Pointer> maps;


	for(unsigned int i = 0; i < labelFilenames.size(); i++)
	{
		if(i == imageToTest)
		{
			testImage = ImageVolumeIO::Read(imageFilenames[i]);
			testLabel = LabelVolumeIO::Read(unalignedLabelFilenames[i]);
		}
		else
		{

			LabelVolume::Pointer label = LabelVolumeIO::Read(labelFilenames[i]);
			labels.push_back(label);
			RealVolume::Pointer map = RealVolumeIO::Read(distanceMapFilenames[i]);
			maps.push_back(map);

		}


	}


	// create the manifold
	typedef SimilarityMatrixBuilder<RealVolume> MatrixBuilder;
	MatrixBuilder::Pointer matBuilder = MatrixBuilder::New();
	matBuilder->SetImages(maps);


	std::cout << "Building Distance Matrix" << std::endl;
	typedef ImageToImageHeavisideDistanceMeasure<RealVolume> MeasureType;
	MeasureType::Pointer measure = MeasureType::New();
	measure->SetNumberOfThreads(8);
	measure->SetBoundryValue(0.0);
	matBuilder->SetDistanceMeasure(measure);
	matBuilder->SetIsSymmetrical(false);
	matBuilder->Update();
	



	

	std::cout << "Building Manifold" << std::endl;
	// build the manifold
	DiffusionMap::Pointer diffusionMap = DiffusionMap::New();
	diffusionMap->SetEpsilon(1);
	diffusionMap->SetComponents(3);
	diffusionMap->SetDistancePrecomputed(true);
	diffusionMap->Fit(matBuilder->GetOutput());




	typedef itk::Similarity3DTransform<double> TransformType;
	typedef ManifoldShapeModel<3, TransformType> ManifoldShapeModel;
	ManifoldShapeModel::Pointer shapeModel = ManifoldShapeModel::New();
	shapeModel->SetLabelData(labels);
	shapeModel->SetDistanceMaps(maps);
	shapeModel->SetMeanLabel(meanLabel);
	shapeModel->SetManifoldBuilder(diffusionMap);
	shapeModel->SetEntropyWeight(entropyWeight);
	shapeModel->SetDistanceMeasure(measure);




	// align the mean with the target
	typedef segmentation::MRFLabelAlignment<LabelVolume,LabelVolume, TransformType> AlignerType;
	AlignerType::Pointer aligner = AlignerType::New();
	aligner->SetFixedImage(testLabel);
	aligner->SetMovingImage(meanLabel);
	aligner->Update();
	
	AlignerType::ParametersType params = aligner->GetFinalParameters();
	TransformType::Pointer transform = aligner->GetTransform();


	// set up the data term
	MatrixReader::Pointer reader = MatrixReader::New();
	reader->SetFilename(gmmParamsFilename);
	reader->Update();
	MatrixDataSet::Pointer paramsDataSet = reader->GetOutput();




	typedef segmentation::MRFGMMDataTerm<ImageVolume::PixelType, double, 2> GMMDataTermType;
	GMMDataTermType::Pointer dataTerm = GMMDataTermType::New();
	
	typedef GMMDataTermType::InputType DataTermInputType;

	dataTerm->SetParameters(1, paramsDataSet->doubleData["PositiveParams"]);
	dataTerm->SetParameters(0, paramsDataSet->doubleData["NegativeParams"]);
	dataTerm->Initialise();


	// set up the smoothness term
	//
	//
	//
	typedef segmentation::MRFIsingSmoothnessTerm<ImageVolume::PixelType, double> SmoothnessTermType;
	SmoothnessTermType::Pointer smoothnessTerm = SmoothnessTermType::New();
	smoothnessTerm->SetSigma(10.0);


	LabelVolume::Pointer priorInitialiser = aligner->GetOutput();


	unsigned int count = 0;
	while(count < 3)
	{




		// create the prior image
		//
		//
		//
		typedef segmentation::MRFPriorGenerator<LabelVolume, RealVolume> PriorGeneratorType;
		PriorGeneratorType::Pointer priorGenerator = PriorGeneratorType::New();
		priorGenerator->SetInput(priorInitialiser);
		priorGenerator->SetSpread(priorSpreadValue);
		priorGenerator->Update();

		RealVolume::Pointer prior = priorGenerator->GetOutput();
		RealVolumeIO::Write("prior.nrrd", prior);

		typedef segmentation::MRFSegmentation<ImageVolume, RealVolume, LabelVolume, 2, double> SegmentationType;
		SegmentationType::Pointer seg = SegmentationType::New();
		seg->SetInput(testImage);
		seg->SetDataTerm(dataTerm);
		seg->SetPrior(prior);
		seg->SetSmoothnessTerm(smoothnessTerm);
		seg->SetNumberOfThreads(1);

		seg->SetPriorWeight(priorWeightValue);
		seg->SetSmoothnessWeight(smoothnessWeightValue);
		seg->Update();

		utils::LabelVolume::Pointer output = utils::LabelVolume::New();
		shapeModel->Optimise(seg->GetOutput(), testImage, transform, output);



		priorInitialiser = output;
		
		// now lets get to it
		utils::LabelVolumeIO::Write("seg.nrrd", seg->GetOutput());
		utils::LabelVolumeIO::Write("segOpt.nrrd", output);
		utils::ImageVolumeIO::Write("image.nrrd", testImage);


		count++;
	}













		



	





	return 0;
}
