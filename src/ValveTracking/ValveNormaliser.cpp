#include "ValveNormaliser.h"
#include <itkPermuteAxesImageFilter.h>
#include <itkCenteredAffineTransform.h>
#include <itkResampleImageFilter.h>


namespace vt
{
// ------------------------------------------------------------------------
ValveNormaliser::ValveNormaliser() 
{
	m_FlipPoints = false;
	m_Flip = false;
}

// ------------------------------------------------------------------------
void ValveNormaliser::Normalise()
{
	ValveType::Pointer workingValve = ValveType::New();
	if(m_Flip)
		FlipValve(m_Valve, workingValve);
	else
		workingValve = m_Valve;

	if(m_FlipPoints)
		FlipPoints(workingValve, workingValve);


	AlignValve(workingValve, m_Output);
}

// ------------------------------------------------------------------------
void ValveNormaliser::UnNormalise()
{
	ImageType::Pointer image = m_Valve->GetImage();

	typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
	ResamplerType::Pointer resampler = ResamplerType::New();
	resampler->SetInput(image);
	resampler->SetTransform(m_Transform->GetInverseTransform());
	resampler->SetOutputParametersFromImage(image);
	resampler->Update();

	m_Output = ValveType::New();
	m_Output->SetImage(resampler->GetOutput());
	m_Output->SetP1(m_Transform->TransformPoint(m_Valve->GetP1()));
	m_Output->SetP2(m_Transform->TransformPoint(m_Valve->GetP2()));
	m_Output->UpdateIndexs();

	if(m_FlipPoints)
		FlipPoints(m_Output, m_Output);

	if(m_Flip)
		FlipValve(m_Output, m_Output);

}


// ------------------------------------------------------------------------
void ValveNormaliser::FlipPoints(const ValveType::Pointer &input, ValveType::Pointer &output)
{
	ValveType::PointType tmp1, tmp2;
	tmp1 = input->GetP1();
	tmp2 = input->GetP2();


	output->SetImage(input->GetImage());
	output->SetP1(tmp2);
	output->SetP2(tmp1);
	output->UpdateIndexs();

}


// ------------------------------------------------------------------------
void ValveNormaliser::FlipValve(const ValveType::Pointer &input, ValveType::Pointer &output)
{
	if(!output) output = ValveType::New();

	ImageType::DirectionType outputDirection = input->GetImage()->GetDirection();

	typedef itk::PermuteAxesImageFilter<ImageType> FlipperType;
	FlipperType::PermuteOrderArrayType axes;
	axes[0] = 1;
	axes[1] = 0;
	axes[2] = 2;

	FlipperType::Pointer flipper = FlipperType::New();
	flipper->SetInput(input->GetImage());
	flipper->SetOrder(axes);
	flipper->Update();

	ImageType::Pointer outputImage = flipper->GetOutput();
	outputImage->SetDirection(outputDirection);

	output->SetImage(outputImage);
	
	// flip the points as well
	ContIndexType ind1 = input->GetInd1();
	ContIndexType ind2 = input->GetInd2();

	double tmp1, tmp2;
	tmp1 = ind1[0];
	tmp2 = ind2[0];

	ind1[0] = ind1[1];
	ind1[1] = tmp1;

	ind2[0] = ind2[1];
	ind2[1] = tmp2;

	output->SetInd1(ind1);
	output->SetInd2(ind2);
	output->UpdatePoints();
}

// ------------------------------------------------------------------------
void ValveNormaliser::AlignValve(const ValveType::Pointer &input, ValveType::Pointer &output)
{
	if(!output) output = ValveType::New();

	ImageType::Pointer image = input->GetImage();
	PointType p1 = input->GetP1();
	PointType p2 = input->GetP2();


	// tranlsation to the origin
	m_Transform = TransformType::New();
	m_Transform->SetCenter(p1);

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

	m_Transform->Rotate3D(axis, angle);

	typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
	ResamplerType::Pointer resampler = ResamplerType::New();
	resampler->SetInput(image);
	resampler->SetTransform(m_Transform);
	resampler->SetOutputParametersFromImage(image);
	resampler->Update();


	// create the output 
	if(!output) output = ValveLine<3>::New();
	output->SetImage(resampler->GetOutput());
	output->SetP1(m_Transform->TransformPoint(input->GetP1()));
	output->SetP2(m_Transform->GetInverseTransform()->TransformPoint(input->GetP2()));
	output->UpdateIndexs();
}

}
