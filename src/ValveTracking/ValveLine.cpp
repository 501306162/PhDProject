#include "ValveLine.h"

#include <vtkLineSource.h>


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


}

template class vt::ValveLine<2>;
template class vt::ValveLine<3>;
template class vt::ValveSequence<2>;
template class vt::ValveSequence<3>;
