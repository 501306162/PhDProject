#include <iostream>

#include <BoundingBox.h>
#include <CMRFileExtractor.h>
#include <ValveOriginFinder.h>

#include <itkSimilarity3DTransform.h>
#include <CommonDefinitions.h>


#include <SimpleMRFSegmenter.h>

using namespace vt;

int main(int argc, char ** argv)
{
	// load the image and the bounding box
	BoundingBox::Pointer boundingBox = BoundingBox::New();
	boundingBox->SetInfation(atof(argv[3]));
	boundingBox->Load(argv[1]);
	
	// load hte images and compute the reference coordinates
	CMRFileExtractor::Pointer extractor = CMRFileExtractor::New();
	extractor->SetFolderName(argv[2]);
	extractor->Extract();

	ValveOriginFinder::Pointer originFinder = ValveOriginFinder::New();
	originFinder->Set2CImage(extractor->Get2CImage(0));
	originFinder->Set3CImage(extractor->Get3CImage(0));
	originFinder->SetImageStack(extractor->GetStackImage(0));
	originFinder->Compute();


	// apply the transform to the bounding box
	typedef itk::Similarity3DTransform<double> TransformType;
	TransformType::Pointer transform = TransformType::New();
	transform->SetMatrix(originFinder->GetRotation());
	transform->SetTranslation(originFinder->GetTranslation());
	boundingBox->TransformBoundingBox(transform);


	BoundingBox::MaskType::Pointer mask = BoundingBox::MaskType::New();
	boundingBox->ComputeImageMask(extractor->Get2CImage(0), 1, mask);


	utils::LabelVolumeIO::Write("mask.nrrd", mask);
	utils::ImageVolumeIO::Write("image.nrrd", extractor->Get2CImage(0));


	SimpleMRFSegmenter::Pointer segmenter = SimpleMRFSegmenter::New();
	segmenter->SetImage(extractor->Get2CImage(0));
	segmenter->SetSmoothnessCost(atof(argv[4]));
	segmenter->SetMask(mask);
	segmenter->Segment();


	utils::LabelVolumeIO::Write("seg.nrrd", segmenter->GetOutput());


	return 0;
}
