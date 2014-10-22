#include "point_picker.h"

#include <vtkRenderWindowInteractor.h>
#include <vtkPropPicker.h>
#include <vtkRenderer.h>
#include <vtkPointData.h>
#include <vtkProp.h>
#include <vtkCell.h>

// ------------------------------------------------------------------------
void PointPicker::OnLeftButtonDown()
{
	int * clickPos = this->GetInteractor()->GetEventPosition();

	// get the xyz at this location
	vtkSmartPointer<vtkPropPicker> picker = 
		vtkSmartPointer<vtkPropPicker>::New();
	picker->Pick(clickPos[0], clickPos[1], 0, this->currentRenderer);





	// trigger the parent event 
	vtkInteractorStyleImage::OnLeftButtonDown();

}

vtkStandardNewMacro(PointPicker);
