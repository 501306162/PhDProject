#include "TestData.h"

#include <QFileInfo>
#include <ValveOriginFinder.h>
#include <itkResampleImageFilter.h>
#include <SimpleMRFSegmenter.h>
#include <itkBinaryContourImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkNrrdImageIO.h>
#include <LineIntersectionFinder.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>


namespace vt
{
// ------------------------------------------------------------------------
TestData::Pointer TestData::Initialise(const std::string &directory)
{
	TestData::Pointer testData =TestData::New();
	testData->Load(directory);
	return testData;
}

// ------------------------------------------------------------------------
void TestData::Load(const std::string &directory)
{
	m_Extractor = CMRFileExtractor::New();
	m_Extractor->SetFolderName(directory);
	m_Extractor->Extract();

	QFileInfo info(QString::fromStdString(directory));
	m_Id = info.fileName().replace("d","").replace(".txt", "").toInt();
}


// ------------------------------------------------------------------------
void TestData::GetImages(std::vector<ImageType::Pointer> &images)
{
	TestGroup::iterator mapIt = m_AlignedGroup.begin();
	while(mapIt != m_AlignedGroup.end())
	{	
		images.push_back(mapIt->second.image);
		++mapIt;
	}
}

// ------------------------------------------------------------------------
void TestData::Prepare(unsigned int timestep)
{

	ImageGroup c2,c3,c4;
	c2.image = m_Extractor->Get2CImage(timestep);
	c3.image = m_Extractor->Get3CImage(timestep);
	c4.image = m_Extractor->Get4CImage(timestep);
	m_StackImage = m_Extractor->GetStackImage(timestep);

	m_InitialGroup["2C"] = c2;
	m_InitialGroup["3C"] = c3;
	m_InitialGroup["4C"] = c4;

	// compute the alignment of the initial group
	ComputeInitialTransform();
	//TransformInitialGroup();
	m_AlignedGroup = m_InitialGroup;
	
	//ProcessImages();

}

// ------------------------------------------------------------------------
void TestData::ProcessImages()
{
	TestGroup::iterator mapIt = m_AlignedGroup.begin();
	while(mapIt != m_AlignedGroup.end())
	{
		ComputeImageMask(mapIt->second);
		ComputeSegmentation(mapIt->second);
		ComputeContours(mapIt->second);
		mapIt++;
	}
}

// ------------------------------------------------------------------------
void TestData::ComputeImageMask(ImageGroup &images)
{
	images.mask = LabelType::New();
	images.mask->SetOrigin(images.image->GetOrigin());
	images.mask->SetDirection(images.image->GetDirection());
	images.mask->SetSpacing(images.image->GetSpacing());
	images.mask->SetRegions(images.image->GetLargestPossibleRegion());
	images.mask->Allocate();
	images.mask->FillBuffer(0);

	
	itk::ImageRegionIterator<LabelType> labelIt(images.mask, images.mask->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<ImageType> imIt(images.image, images.image->GetLargestPossibleRegion());



	while(!imIt.IsAtEnd())
	{
		ImageType::IndexType index = imIt.GetIndex();
		ImageType::PointType point;
		images.image->TransformIndexToPhysicalPoint(index, point);

		point = m_InitialTransform->GetInverseTransform()->TransformPoint(point);

		if(m_BoundingBox.ContainsPoint(point.GetDataPointer()))
			labelIt.Set(255);

		++imIt; ++labelIt;
	}

}

// ------------------------------------------------------------------------
void TestData::SaveImageGroup(std::string type, std::string name)
{
	ImageGroup group;
	GetImageGroup(type, group);
	SaveGroup(group, name);
}


// ------------------------------------------------------------------------
void TestData::GetImageGroup(const std::string type, ImageGroup &group)
{
	if(m_AlignedGroup.count(type) == 0)
		group = m_AlignedGroup["2C"];
	else
		group = m_AlignedGroup[type];
}



// ------------------------------------------------------------------------
void TestData::SaveGroup(const ImageGroup &group, std::string filename)
{
	if(group.image) SaveImage(group.image, filename+"-image.nrrd");
	if(group.seg) SaveLabel(group.seg, filename+"-seg.nrrd");
	if(group.contours) SaveLabel(group.contours, filename+"-contours.nrrd");
	if(group.mask) SaveLabel(group.mask, filename+"-mask.nrrd");
}


// ------------------------------------------------------------------------
void TestData::ComputeSegmentation(ImageGroup &images)
{
	SimpleMRFSegmenter::Pointer segmenter = SimpleMRFSegmenter::New();
	segmenter->SetImage(images.image);
	segmenter->SetMask(images.mask);
	segmenter->SetOutputValue(255);
	segmenter->SetSmoothnessCost(1.0);
	segmenter->SetMaskValue(255);
	segmenter->Segment();
	
	images.seg = segmenter->GetOutput();

}


// ------------------------------------------------------------------------
void TestData::ComputeSegmentation2(ImageGroup &images)
{
	typedef itk::OtsuThresholdImageFilter<ImageType, LabelType> ThresholderType;
	ThresholderType::Pointer segmenter = ThresholderType::New();
	segmenter->SetInput(images.image);
	segmenter->SetMaskImage(images.mask);
	segmenter->SetMaskValue(255);
	segmenter->Update();
	
	images.seg = segmenter->GetOutput();

}



// -----------------------------------------------------------------------
void TestData::ComputeContours(ImageGroup &image)
{
	typedef itk::BinaryContourImageFilter<LabelType, LabelType> ContourType;
	ContourType::Pointer contour = ContourType::New();
	contour->SetInput(image.seg);
	contour->SetFullyConnected(true);
	contour->SetBackgroundValue(0);
	contour->SetForegroundValue(255);
	contour->Update();

	image.contours = contour->GetOutput();
}



// ------------------------------------------------------------------------
void TestData::ComputeInitialTransform()
{
	ValveOriginFinder::Pointer originFinder = ValveOriginFinder::New();
	originFinder->Set2CImage(m_InitialGroup["2C"].image);
	originFinder->Set3CImage(m_InitialGroup["3C"].image);
	originFinder->SetImageStack(m_StackImage);
	originFinder->Compute();
	
	m_InitialTransform = TransformType::New();
	m_InitialTransform->SetTranslation(originFinder->GetTranslation());
	m_InitialTransform->SetMatrix(originFinder->GetRotation());
}

// ------------------------------------------------------------------------
void TestData::TransformInitialGroup()
{
	TestGroup::iterator mapIt = m_InitialGroup.begin();
	while(mapIt != m_InitialGroup.end())
	{
		ImageGroup newGroup;
		TransformImage(mapIt->second.image, m_InitialTransform, newGroup.image);
		m_AlignedGroup[mapIt->first] = newGroup;

		++mapIt;
	}
}

// ------------------------------------------------------------------------
void TestData::TransformImage(const ImageType::Pointer &input, 
		TransformType::Pointer &transform, ImageType::Pointer &output)
{
	ImageType::DirectionType dir = input->GetDirection();
	ImageType::DirectionType newDir = transform->GetInverseMatrix() * dir;
	ImageType::PointType newOrigin = transform->GetInverseTransform()->TransformPoint(input->GetOrigin());

	typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
	ResamplerType::Pointer resampler = ResamplerType::New();
	resampler->SetInput(input);
	resampler->SetTransform(transform);
	resampler->SetOutputParametersFromImage(input);
	resampler->SetOutputOrigin(newOrigin);
	resampler->SetOutputDirection(newDir);
	resampler->Update();

	output = resampler->GetOutput();
}

// ------------------------------------------------------------------------
void TestData::SaveLabel(const LabelType::Pointer &label, std::string filename)
{
	typedef itk::ImageFileWriter<LabelType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(label);
	writer->SetFileName(filename);
	writer->SetImageIO(itk::NrrdImageIO::New());
	writer->Update();
}

// ------------------------------------------------------------------------
void TestData::SaveImage(const ImageType::Pointer &label, std::string filename)
{
	typedef itk::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(label);
	writer->SetFileName(filename);
	writer->SetImageIO(itk::NrrdImageIO::New());
	writer->Update();
}

// ------------------------------------------------------------------------
void TestData::ExtractTestLines(const VectorType &normal, const VectorType &point, LineGroup &lines)
{
	TestGroup::iterator mapIt = m_AlignedGroup.begin();
	while(mapIt != m_AlignedGroup.end())
	{
		LineTypeList lineList;
		GetTestLines(mapIt->second, normal, point, lineList);

		lines[mapIt->first] = lineList;
		++mapIt;
	}

}

// ------------------------------------------------------------------------
void TestData::GetLinesFromPlane(const VectorType &normal, const VectorType &point, LineGroup &lines)
{
	TestGroup::iterator mapIt = m_AlignedGroup.begin();
	while(mapIt != m_AlignedGroup.end())
	{
		// find the intersection of the plane and the image
		LineIntersectionFinder::Pointer finder = LineIntersectionFinder::New();
		finder->SetImage(mapIt->second.image);
		finder->SetPlane(normal, point);
		finder->SetBoundingBox(m_BoundingBox);
		finder->Compute();

		LineTypeList lineList;
		lineList.push_back(finder->GetOutput());

		lines[mapIt->first] = lineList;
		++mapIt;
	}

}

// ------------------------------------------------------------------------
void TestData::GetTestLines(const ImageGroup &images, const VectorType &normal, const VectorType &point, LineTypeList &lines)
{
	// find the intersection of the plane and the image
	LineIntersectionFinder::Pointer finder = LineIntersectionFinder::New();
	finder->SetImage(images.image);
	finder->SetPlane(normal, point);
	finder->SetBoundingBox(m_BoundingBox);
	finder->Compute();

	LineIntersectionFinder::OutputLineType line = finder->GetOutput();

	// set up the interpolator
	typedef itk::NearestNeighborInterpolateImageFunction<LabelType> NNInterpolatorType;
	NNInterpolatorType::Pointer interpolator = NNInterpolatorType::New();
	interpolator->SetInputImage(images.seg);



	PointType start = line.p1;
	PointType end = line.p2;
	const unsigned int steps = 100;
	const double totalDist = (end-start).GetNorm();
	const double totalT = totalDist / line.direction.norm();
	const double incrementSize = totalT /  (double) steps;

	
	// get starting 
	bool up = false;
	if(interpolator->Evaluate(start) == 255)
		up = true;

	bool inLine = false;
	PointType segStart;

	for(unsigned int i = 0; i < steps; i++)
	{
		double t = i*incrementSize;	
		PointType newPoint = start;
		for(unsigned int j = 0; j < 3; j++)
		{
			newPoint[j] += line.direction(j) * t;	
		}

		if(!interpolator->IsInsideBuffer(newPoint)) continue;

		unsigned char value = interpolator->Evaluate(newPoint);
		if(up && value == 0)
		{
			up = false;
			segStart = newPoint;
			inLine = true;
		}

		if(!up && value == 255)
		{
			up = true;

			if(inLine)
			{
				LineType segment;
				segment.p1 = segStart;
				segment.p2 = newPoint;
				segment.direction = line.direction;

				lines.push_back(segment);
				inLine = false;
			}
		}
	}
}

}
