#include <iostream>

#include <MRFGMMDataTerm.h>
#include <MatrixReader.h>
#include <ImageIO.h>

#include <CommonDefinitions.h>
#include <itkImageRegionIterator.h>

#include <MRFSegmentation.h>
#include <MRFIsingSmoothnessTerm.h>
#include <MRFImageLabeller.h>

#include <MRFPriorGenerator.h>
#include <MRFLabelAlignment.h>
#include <itkSimilarity3DTransform.h>

using namespace utils;
using namespace segmentation;
int main(int argc, char *argv[])
{
	
	// load an image
	ImageVolume::Pointer image = ImageVolumeIO::Read(argv[2]);
	LabelVolume::Pointer meanLabel = LabelVolumeIO::Read(argv[3]);
	LabelVolume::Pointer groundTruth = LabelVolumeIO::Read(argv[4]);



	// align the mean with the ground truth
	typedef itk::Similarity3DTransform<double> TransformType;
	typedef segmentation::MRFLabelAlignment<LabelVolume, LabelVolume, TransformType> LabelAlignerType;
	LabelAlignerType::Pointer aligner = LabelAlignerType::New();
	aligner->SetFixedImage(groundTruth);
	aligner->SetMovingImage(meanLabel);
	aligner->Update();
	



	
	// set up the data term
	MatrixReader::Pointer reader = MatrixReader::New();
	reader->SetFilename(argv[1]);
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



	// create the prior image
	//
	//
	//
	LabelVolume::Pointer label = LabelVolumeIO::Read(argv[3]);
	typedef segmentation::MRFPriorGenerator<LabelVolume, RealVolume> PriorGeneratorType;
	PriorGeneratorType::Pointer priorGenerator = PriorGeneratorType::New();
	priorGenerator->SetInput(aligner->GetOutput());
	priorGenerator->SetSpread(atof(argv[5]));
	priorGenerator->Update();

	RealVolume::Pointer prior = priorGenerator->GetOutput();
	RealVolumeIO::Write("prior.nrrd", prior);


	typedef segmentation::MRFSegmentation<ImageVolume, RealVolume, LabelVolume, 2, double> SegmentationType;
	SegmentationType::Pointer seg = SegmentationType::New();
	seg->SetInput(image);
	seg->SetDataTerm(dataTerm);
	seg->SetPrior(prior);
	seg->SetSmoothnessTerm(smoothnessTerm);
	seg->SetNumberOfThreads(1);

	seg->SetPriorWeight(atof(argv[6]));
	seg->SetSmoothnessWeight(atof(argv[7]));
	seg->Update();


	/*
	typedef segmentation::MRFImageLabeller<ImageVolume> Labeller;
	typedef Labeller::LabelSetType LabelSetType;
	LabelSetType labels;

	// create the label set
	itk::ImageRegionIterator<ImageVolume> it(image, image->GetLargestPossibleRegion());
	while(!it.IsAtEnd())
	{
		unsigned short val = it.Get();

		DataTermInputType inputParams(1);
		inputParams(0) = val;

		GMMDataTermType::OutputType output = dataTerm->Evaluate(inputParams);

		if(output[0] < output[1])
		{
			labels.push_back(1);			
		}
		else
		{
			labels.push_back(0);
		}

		++it;
	}
	

	Labeller::Pointer labeller = Labeller::New();
	labeller->UseImageParameters(image);
	labeller->SetLabelSet(labels);
	labeller->Update();

	*/

	ImageVolumeIO::Write("image.nrrd", image);
	LabelVolumeIO::Write("seg.nrrd", seg->GetOutput());




	return 0;
}

