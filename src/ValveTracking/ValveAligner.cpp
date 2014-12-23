#include "ValveAligner.h"
#include <itkLandmarkBasedTransformInitializer.h>
#include <itkResampleImageFilter.h>

namespace vt 
{
// ------------------------------------------------------------------------
void ValveAlignementTransformExtractor::Extract()
{
	typedef itk::LandmarkBasedTransformInitializer<TransformType, 
			ValveType::ImageType, ValveType::ImageType> InitialiserType;


	typedef InitialiserType::LandmarkPointContainer ContainerType;
	typedef InitialiserType::LandmarkPointType LandmarkType;

	ContainerType fixedLandmarks;
	ContainerType movingLandmarks;

	fixedLandmarks.push_back(m_FixedValve->GetP1());
	fixedLandmarks.push_back(m_FixedValve->GetP2());

	movingLandmarks.push_back(m_MovingValve->GetP1());
	movingLandmarks.push_back(m_MovingValve->GetP2());


	InitialiserType::Pointer initialiser = InitialiserType::New();
	initialiser->SetFixedLandmarks(fixedLandmarks);
	initialiser->SetMovingLandmarks(movingLandmarks);

	m_Transform = TransformType::New();
	m_Transform->SetIdentity();
	initialiser->SetTransform(m_Transform);
	initialiser->InitializeTransform();
}



// ------------------------------------------------------------------------
void ValveLineResampler::Resample()
{
	typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
	ResamplerType::Pointer resampler = ResamplerType::New();
	resampler->SetInput(m_InputValve->GetImage());
	resampler->SetUseReferenceImage(true);
	resampler->SetOutputParametersFromImage(m_ReferenceValve->GetImage());
	resampler->SetTransform(m_Transform);
	resampler->SetDefaultPixelValue(0);

	resampler->Update();


	// create the output valve
	m_Output = ValveType::New();
	m_Output->SetImage(resampler->GetOutput());
	



	// transform the points
	m_Output->SetP1(m_Transform->GetInverseTransform()->TransformPoint(m_InputValve->GetP1()));
	m_Output->SetP2(m_Transform->GetInverseTransform()->TransformPoint(m_InputValve->GetP2()));
	m_Output->UpdateIndexs();
}

// ------------------------------------------------------------------------
void ValveSequenceAligner::Align()
{
	// compute the transform 
	ValveAlignementTransformExtractor::Pointer transformComputer = ValveAlignementTransformExtractor::New();
	transformComputer->SetFixedValve(m_FixedSequence->GetValveLine(0));
	transformComputer->SetMovingValve(m_MovingSequence->GetValveLine(0));
	transformComputer->Extract();

	TransformType::Pointer transform = transformComputer->GetTransform();

	m_Output = SequenceType::New();
	
	for(unsigned int i = 0; i < m_FixedSequence->GetNumberOfLines(); i++)
	{
		ValveLineResampler::Pointer resampler = ValveLineResampler::New();
		resampler->SetTransform(transform);
		resampler->SetInputValve(m_MovingSequence->GetValveLine(i));
		resampler->SetReferenceValve(m_FixedSequence->GetValveLine(0));
		resampler->Resample();

		m_Output->AddValveLine(resampler->GetOutput());
	}
}


}
