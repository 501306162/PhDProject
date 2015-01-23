#include "ValveLine.h"

#include <vtkLineSource.h>
#include <itkPermuteAxesImageFilter.h>

namespace vt {
// ------------------------------------------------------------------------
template<int VDimensions>
void ValveLine<VDimensions>::UpdateIndexs()
{
	m_Image->TransformPhysicalPointToContinuousIndex(m_P1, m_Ind1);
	m_Image->TransformPhysicalPointToContinuousIndex(m_P2, m_Ind2);
}


// ------------------------------------------------------------------------
template<int VDimensions>
void ValveLine<VDimensions>::UpdatePoints()
{
	m_Image->TransformContinuousIndexToPhysicalPoint(m_Ind1, m_P1);
	m_Image->TransformContinuousIndexToPhysicalPoint(m_Ind2, m_P2);
}

// ------------------------------------------------------------------------
template<int VDimensions>
vtkSmartPointer<vtkPolyData> 
ValveLine<VDimensions>::GetPolyData() const
{

	vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();

	if(VDimensions == 3)
	{
		lineSource->SetPoint1(m_P1[0], m_P1[1], m_P1[2]);
		lineSource->SetPoint2(m_P2[0], m_P2[1], m_P2[2]);
	}
	else
	{
		lineSource->SetPoint1(m_P1[0], m_P1[1], 0.0);
		lineSource->SetPoint2(m_P2[0], m_P2[1], 0.0);
	}

	lineSource->Update();
	return lineSource->GetOutput();
}


// ------------------------------------------------------------------------
template<int VDimensions>
void ValveLine<VDimensions>::FlipPoints()
{
	PointType tmp1, tmp2;
	tmp1 = GetP1();
	tmp2 = GetP2();

	SetP1(tmp2);
	SetP2(tmp1);
	UpdateIndexs();

}

// ------------------------------------------------------------------------
template<int VDimensions>
void ValveLine<VDimensions>::FlipImage()
{
	typename ImageType::DirectionType outputDirection = m_Image->GetDirection();

	typedef itk::PermuteAxesImageFilter<ImageType> FlipperType;
	typename FlipperType::PermuteOrderArrayType axes;
	axes[0] = 1;
	axes[1] = 0;
	axes[2] = 2;

	typename FlipperType::Pointer flipper = FlipperType::New();
	flipper->SetInput(m_Image);
	flipper->SetOrder(axes);
	flipper->Update();

	typename ImageType::Pointer outputImage = flipper->GetOutput();
	outputImage->SetDirection(outputDirection);

	
	// flip the points as well
	ContIndexType ind1 = GetInd1();
	ContIndexType ind2 = GetInd2();

	double tmp1, tmp2;
	tmp1 = ind1[0];
	tmp2 = ind2[0];

	ind1[0] = ind1[1];
	ind1[1] = tmp1;

	ind2[0] = ind2[1];
	ind2[1] = tmp2;

	SetImage(outputImage);
	SetInd1(ind1);
	SetInd2(ind2);
	UpdatePoints();
}

}

template class vt::ValveLine<2>;
template class vt::ValveLine<3>;
template class vt::ValveSequence<2>;
template class vt::ValveSequence<3>;
