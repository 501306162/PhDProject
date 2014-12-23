#include "JsonConverter.h"

#include <itkImageFileReader.h>
#include <itkNrrdImageIO.h>
#include <itkImageToVTKImageFilter.h>
#include <itkFlipImageFilter.h>
#include <itkRegionOfInterestImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkImageFileWriter.h>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>


namespace vt
{
// ------------------------------------------------------------------------
void JsonConverter::Convert()
{

	m_Output = ValveSequence<3>::New();

	unsigned int timeSteps = m_Input.lines.size();
	
	m_VtkImage = vtkSmartPointer<vtkImageData>::New();
	m_Reference = ReferenceType::New();

	LoadImages(m_Input.imageFilename, m_Reference, m_VtkImage);

	for(unsigned int i = 0; i < timeSteps; i++)
	{
		ReferenceType::RegionType refRegion = m_Reference->GetLargestPossibleRegion();


		typedef itk::Image<unsigned short, 3> TimeStepType;
		typedef itk::ExtractImageFilter<ReferenceType, TimeStepType> ExtractorType;
		typedef itk::ImageToVTKImageFilter<TimeStepType> VTKFilterType;

		// create the extraction region 
		ReferenceType::RegionType exRegion;
		ReferenceType::SizeType exSize = refRegion.GetSize();
		exSize[3] = 0;
		ReferenceType::IndexType exIndex = refRegion.GetIndex();
		exIndex[3] = i;

		exRegion.SetSize(exSize);
		exRegion.SetIndex(exIndex);

		ExtractorType::Pointer extractor = ExtractorType::New();
		extractor->SetInput(m_Reference);
		extractor->SetExtractionRegion(exRegion);
		extractor->SetDirectionCollapseToIdentity();
		extractor->Update();


		typedef itk::RegionOfInterestImageFilter<ImageType, ImageType> ROIFilter;
		ROIFilter::Pointer roiFilter = ROIFilter::New();
		roiFilter->SetInput(extractor->GetOutput());

		TimeStepType::RegionType roiRegion = extractor->GetOutput()->GetLargestPossibleRegion();
		TimeStepType::SizeType roiSize = roiRegion.GetSize();
		TimeStepType::IndexType roiIndex = roiRegion.GetIndex();
		roiSize[2] = 1;
		roiIndex[2] = 0;
		roiRegion.SetSize(roiSize);
		roiRegion.SetIndex(roiIndex);

		roiFilter->SetRegionOfInterest(roiRegion);
		roiFilter->Update();

		// convert the lines
		ValveLine<3>::Pointer valve = ValveLine<3>::New();

		ImageType::Pointer image = roiFilter->GetOutput();
		ConvertLine(m_Input.lines[i], image, valve);
		FinishUp(i, valve);

		m_Output->AddValveLine(valve);

	}
	std::cout << m_Output->GetNumberOfLines() << std::endl;
}

// ------------------------------------------------------------------------
void JsonConverter::FinishUp(const unsigned int &timeStep, ValveLine<3>::Pointer &valve)
{
	ReferenceType::RegionType refRegion = m_Reference->GetLargestPossibleRegion();

	typedef itk::ExtractImageFilter<ReferenceType, ImageType> ExtractorType;

	// create the extraction region 
	ReferenceType::RegionType exRegion;
	ReferenceType::SizeType exSize = refRegion.GetSize();
	exSize[3] = 0;
	ReferenceType::IndexType exIndex = refRegion.GetIndex();
	exIndex[3] = timeStep;

	exRegion.SetSize(exSize);
	exRegion.SetIndex(exIndex);

	ExtractorType::Pointer extractor = ExtractorType::New();
	extractor->SetInput(m_Reference);
	extractor->SetExtractionRegion(exRegion);
	extractor->SetDirectionCollapseToSubmatrix();
	extractor->Update();

	ImageType::Pointer image = extractor->GetOutput();


	ImageType::PointType p1,p2;

	image->TransformContinuousIndexToPhysicalPoint(valve->GetInd1(), p1);
	image->TransformContinuousIndexToPhysicalPoint(valve->GetInd2(), p2);

	valve->SetP1(p1);
	valve->SetP2(p2);
	valve->SetImage(image);

}



// ------------------------------------------------------------------------
void JsonConverter::ConvertLine(const Line & line, ImageType::Pointer &image, ValveLine<3>::Pointer & output)
{
	double p1[3];
	double p2[3];

	ConvertPoint(image, line.p1, p1);
	ConvertPoint(image, line.p2, p2);

	PointType ip1, ip2;
	for(unsigned int i = 0; i < 3; i++)
	{
		ip1[i] = p1[i];
		ip2[i] = p2[i];
	}

	typedef itk::ContinuousIndex<double, 3> IndexType;
	IndexType index1, index2;
	image->TransformPhysicalPointToContinuousIndex(ip1, index1);
	image->TransformPhysicalPointToContinuousIndex(ip2, index2);

	output->SetInd1(index1);
	output->SetInd2(index2);

}



// ------------------------------------------------------------------------
void JsonConverter::ConvertPoint(ImageType::Pointer &image, const double * input, double * output)
{
	double * vtkOrigin = m_VtkImage->GetOrigin();

	// get end point 
	//
	ImageType::PointType point;
	ImageType::IndexType index;
	index.Fill(0);
	index[1] = image->GetLargestPossibleRegion().GetSize()[1]-1;
	image->TransformIndexToPhysicalPoint(index, point);


	double xOff = input[0] - vtkOrigin[0];
	output[0] = xOff + image->GetOrigin()[0];

	double yOff = input[1] - vtkOrigin[1];
	output[1] = point[1] - yOff;

	output[2] = point[2];

}


// ------------------------------------------------------------------------
void JsonConverter::LoadImages(const std::string &filename, 
		ReferenceType::Pointer &ref, vtkImageData* image)
{
	typedef itk::ImageFileReader<ReferenceType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename);
	reader->SetImageIO(itk::NrrdImageIO::New());
	reader->Update();
	
	ref = reader->GetOutput();


	// flip 
	typedef itk::FlipImageFilter<ReferenceType> FlipperType;

	FlipperType::Pointer flipper = FlipperType::New();
	FlipperType::FlipAxesArrayType axes;
	axes[0] = false;
	axes[1] = true;
	axes[2] = false;
	axes[3] = false;

	flipper->SetFlipAxes(axes);
	flipper->SetInput(ref);
	flipper->Update();


	ReferenceType::Pointer reference = flipper->GetOutput();
	ReferenceType::RegionType refRegion = reference->GetLargestPossibleRegion();
	const unsigned int slices = refRegion.GetSize()[2];


	typedef itk::Image<unsigned short, 3> TimeStepType;
	typedef itk::ExtractImageFilter<ReferenceType, TimeStepType> ExtractorType;
	typedef itk::ImageToVTKImageFilter<TimeStepType> VTKFilterType;

	// create the extraction region 
	ReferenceType::RegionType exRegion;
	ReferenceType::SizeType exSize = refRegion.GetSize();
	exSize[3] = 0;
	ReferenceType::IndexType exIndex = refRegion.GetIndex();
	exIndex[3] = 0;

	exRegion.SetSize(exSize);
	exRegion.SetIndex(exIndex);

	ExtractorType::Pointer extractor = ExtractorType::New();
	extractor->SetInput(reference);
	extractor->SetExtractionRegion(exRegion);
	extractor->SetDirectionCollapseToSubmatrix();
	extractor->Update();


	for(unsigned int j = 0; j < slices; j++)
	{
		typedef itk::RegionOfInterestImageFilter<TimeStepType, TimeStepType> ROIFilter;
		ROIFilter::Pointer roiFilter = ROIFilter::New();
		roiFilter->SetInput(extractor->GetOutput());

		TimeStepType::RegionType roiRegion = extractor->GetOutput()->GetLargestPossibleRegion();
		TimeStepType::SizeType roiSize = roiRegion.GetSize();
		TimeStepType::IndexType roiIndex = roiRegion.GetIndex();
		roiSize[2] = 1;
		roiIndex[2] = j;
		roiRegion.SetSize(roiSize);
		roiRegion.SetIndex(roiIndex);

		roiFilter->SetRegionOfInterest(roiRegion);


		VTKFilterType::Pointer vtkFilter = VTKFilterType::New();
		vtkFilter->SetInput(roiFilter->GetOutput());
		vtkFilter->Update();
		image->DeepCopy(vtkFilter->GetOutput());

		break;
	}
}




}



