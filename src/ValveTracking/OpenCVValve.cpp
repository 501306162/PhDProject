#include "OpenCVValve.h"

#include <itkCastImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>

#include <opencv2/highgui/highgui.hpp>

#include <itkImageFileWriter.h>
#include <itkPNGImageIO.h>
#include <itkImageRegionConstIterator.h>


namespace vt
{
// ------------------------------------------------------------------------
void OpenCVValve::InitialiseFromValve(const ValveLine<2>::Pointer &input)
{
	ConvertImage(input->GetImage(), m_Image);

	// set the points as the continuous indexs
	m_P1.x = input->GetInd1()[0];
	m_P1.y = input->GetInd1()[1];
	m_P2.x = input->GetInd2()[0];
	m_P2.y = input->GetInd2()[1];
}

// ------------------------------------------------------------------------
void OpenCVValve::ConvertImage(const ImageType::Pointer &input, MatPtr &mat)
{
	// cast the image to uchar 
	typedef itk::Image<unsigned char, 2> OutputImageType;
	typedef itk::RescaleIntensityImageFilter<ImageType, OutputImageType> CasterType;
	CasterType::Pointer caster = CasterType::New();
	caster->SetOutputMaximum(255);
	caster->SetOutputMinimum(0);
	caster->SetInput(input);
	caster->Update();


	OutputImageType::Pointer output = caster->GetOutput();


	typedef itk::ImageFileWriter<OutputImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetImageIO(itk::PNGImageIO::New());
	writer->SetInput(output);
	writer->SetFileName("test.png");
	writer->Update();



	ImageType::SizeType size = input->GetLargestPossibleRegion().GetSize();
	unsigned int rows = size[1];
	unsigned int cols = size[0];

	mat = new MatType(rows,cols, CV_8UC1);


	itk::ImageRegionConstIterator<OutputImageType> it(output, output->GetLargestPossibleRegion());
	it.GoToBegin();
	while(!it.IsAtEnd())
	{
		OutputImageType::IndexType index = it.GetIndex();
		unsigned char val = it.Get();
		mat->at<unsigned char>(cv::Point(index[0], index[1])) = val;
		++it;
	}
}

}
