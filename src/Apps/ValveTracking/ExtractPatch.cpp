#include <iostream>

#include <ValveIO.h>
#include <CommonDefinitions.h>
#include <PatchExtractor.h>
#include <CMRFileExtractor.h>
#include <itkResampleImageFilter.h>
#include <itkSimilarity3DTransform.h>
#include <itkCenteredAffineTransform.h>
#include <ValveOriginFinder.h>

using namespace vt;
int main(int argc, char ** argv)
{
	/*
	// load up a valve
	ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
	reader->SetFileName(argv[1]);
	ValveSequence<3>::Pointer seq = reader->GetOutput();
	ValveLine<3>::Pointer line = seq->GetValveLine(0);


	
	PatchExtractor::Pointer extractor = PatchExtractor::New();
	extractor->SetImage(line->GetImage());
	
	PatchExtractor::SizeType size;
	size[0]=50;
	size[1]=35;
	extractor->SetPatchSize(size);
	extractor->SetPatchCenter(line->GetInd1());

	utils::ImageVolumeIO::Write("test_patch.nrrd", extractor->ExtractPatch());
	utils::ImageVolumeIO::Write("test_image.nrrd", line->GetImage());
	utils::LabelVolumeIO::Write("test_mask.nrrd", extractor->ExtractMask());

	*/

	// load up a directory
	CMRFileExtractor::Pointer loader = CMRFileExtractor::New();
	loader->SetFolderName(argv[1]);
	loader->Extract();


	ValveOriginFinder::Pointer finder = ValveOriginFinder::New();
	finder->Set2CImage(loader->Get2CImage(0));
	finder->Set3CImage(loader->Get3CImage(0));
	finder->SetImageStack(loader->GetStackImage(0));
	finder->ComputeNormalise();

	typedef CMRFileExtractor::ImageType ImageType;
	typedef itk::Vector<double, 3> VectorType;
	ImageType::Pointer stack = loader->GetStackImage(0);
	ImageType::Pointer image = loader->Get2CImage(0);

	// tranlsation to the origin
	typedef itk::CenteredAffineTransform<double, 3> TransformType;
	TransformType::Pointer transform = TransformType::New();
	transform->SetCenter(finder->GetCenter());
	transform->Rotate3D(finder->GetAxis(), finder->GetAngle());

	typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
	ResamplerType::Pointer resampler = ResamplerType::New();
	resampler->SetInput(image);
	resampler->SetTransform(transform);
	resampler->SetOutputParametersFromImage(image);
	resampler->Update();


	utils::ImageVolumeIO::Write("rot.nrrd", resampler->GetOutput());
	utils::ImageVolumeIO::Write("image.nrrd", image);



	/*
	// load up a valve
	ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
	reader->SetFileName(argv[1]);
	ValveSequence<3>::Pointer seq = reader->GetOutput();
	ValveLine<3>::Pointer line = seq->GetValveLine(0);

	typedef ValveLine<3>::ImageType ImageType;
	typedef ImageType::PointType PointType;
	ImageType::Pointer image = line->GetImage();
	
	PointType p1 = line->GetP1();
	PointType p2 = line->GetP2();


	// tranlsation to the origin
	typedef itk::CenteredAffineTransform<double, 3> TransformType;
	TransformType::Pointer transform = TransformType::New();
	transform->SetCenter(p1);

	TransformType::OutputVectorType axis;
	for(unsigned int i = 0; i < 3; i++)
	{
		axis[i] = -image->GetDirection()(i,2);
	}

	itk::Vector<double, 3> vec1, vec2;
	for(unsigned int i = 0; i < 3; i++)
	{
		vec1[i] = p2[i]-p1[i];
		vec2[i] = image->GetDirection()(i,0);
	}

	vec1.Normalize();
	vec2.Normalize();

	double angle = acos(vec2*vec1);
	itk::Vector<double,3> axis2 = itk::CrossProduct(vec1,vec2);
	axis2.Normalize();

	std::cout << angle * (180 / M_PI) << std::endl;

	std::cout << axis << " " << axis2 << std::endl;


	transform->Rotate3D(axis, angle);

	typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
	ResamplerType::Pointer resampler = ResamplerType::New();
	resampler->SetInput(image);
	resampler->SetTransform(transform);
	resampler->SetOutputParametersFromImage(image);
	resampler->Update();


	utils::ImageVolumeIO::Write("rot.nrrd", resampler->GetOutput());





	vnl_matrix_fixed<double,3,3> originRotation = loader->GetStackImage(0)->GetDirection().GetVnlMatrix();
	for(unsigned int i = 0; i < 3; i++)
	{
		double z = originRotation(2, i);
		double y = originRotation(1, i);
		originRotation(2, i) = y;
		originRotation(1, i) = z;
	}

	std::cout << originRotation << std::endl;

	for(unsigned int i = 0; i < 3; i++)
	{
		translation[i] = loader->Get2CImage(0)->GetOrigin()[i];
	}

	transform->SetTranslation(translation);

	TransformType::MatrixType rot;
	rot =  loader->Get2CImage(0)->GetDirection() * loader->GetStackImage(0)->GetDirection().GetInverse();
	std::cout << rot << std::endl;
	transform->SetMatrix(rot);
	


	typedef CMRFileExtractor::ImageType ImageType;

	ImageType::Pointer input = loader->Get2CImage(0);
	ImageType::Pointer output = ImageType::New();
	output->SetSpacing(input->GetSpacing());
	output->SetRegions(input->GetLargestPossibleRegion());
	ImageType::DirectionType dir;
	dir.GetVnlMatrix() = originRotation;
	output->SetDirection(loader->GetStackImage(0)->GetDirection());
	output->Allocate();

	typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
	ResamplerType::Pointer resampler = ResamplerType::New();
	resampler->SetInput(loader->Get2CImage(0));
	resampler->SetTransform(transform);
	resampler->SetOutputParametersFromImage(output);
	resampler->Update();


	utils::ImageVolumeIO::Write("rot.nrrd", resampler->GetOutput());


	
	*/

	return 0;
}
