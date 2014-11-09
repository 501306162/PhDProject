#include "point_picker.h"

#include <vtkRenderWindowInteractor.h>
#include <vtkPropPicker.h>
#include <vtkRenderer.h>
#include <vtkPointData.h>
#include <vtkProp.h>
#include <vtkCell.h>
#include <vtkProperty.h>


// ------------------------------------------------------------------------
PointPicker::PointPicker()
{
	pickTolerance = 0.01;


	pointPicker = vtkSmartPointer<vtkPointPicker>::New();
	cellPicker = vtkSmartPointer<vtkCellPicker>::New();
	pointPicker->SetTolerance(pickTolerance);


	// set up the ghost point
	ghostPoints = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkPoints>::New();
	ghostPoints->InsertNextPoint(0,0,0);

	ghostPoly = vtkSmartPointer<vtkPolyData>::New();
	ghostPoly->SetPoints(ghostPoints);

	ghostFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
	ghostFilter->SetInputData(ghostPoly);
	ghostFilter->Update();

	ghostMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	ghostMapper->SetInputConnection(ghostFilter->GetOutputPort());

	ghostActor = vtkSmartPointer<vtkActor>::New();
	ghostActor->SetMapper(ghostMapper);
	//ghostActor->VisibilityOff();
	ghostActor->GetProperty()->SetPointSize(10);
	ghostActor->GetProperty()->SetColor(1,0,0);

	isMoving = false;

}



// ------------------------------------------------------------------------
void PointPicker::OnLeftButtonDown()
{
	int * clickPos = this->GetInteractor()->GetEventPosition();

	// get the xyz at this location
	pointPicker->Pick(clickPos[0], clickPos[1], 0, this->currentRenderer);
	double * pos = pointPicker->GetPickPosition();

	if(PointAtPickLocation())
	{
		this->isMoving = true;
		this->pickedPoint = pointPicker->GetPointId();
		std::cout << "Moving point: " << this->pickedPoint << std::endl;
		this->pickedData = dynamic_cast<vtkPolyData*>(pointPicker->GetDataSet());

		SetPointLocation();
	}

	

	// trigger the parent event 
	//vtkInteractorStyleImage::OnLeftButtonDown();
}

// ------------------------------------------------------------------------
void PointPicker::SetPointLocation()
{
	int * clickPos = this->GetInteractor()->GetEventPosition();
	std::cout << clickPos[0] << " " << clickPos[1] << " " << clickPos[2] << std::endl;

	// get the xyz at this location
	int success = cellPicker->Pick(clickPos[0], clickPos[1], 0, this->currentRenderer);
	std::cout << success << std::endl;

	double * pos = cellPicker->GetPickPosition();

	std::cout << "Settting location of point " << this->pickedPoint << ": " << 
		pos[0] << " " << pos[1] << " " << pos[2] << std::endl;

	pickedData->GetPoints()->SetPoint(this->pickedPoint, pos[0], pos[1], pos[2]);
	pickedData->Modified();
	window->Render();

}

// ------------------------------------------------------------------------
bool PointPicker::PointAtPickLocation()
{
	int pickedId = pointPicker->GetPointId();
	if(pickedId >= 0 && pickedId < 2)
	{
		return true;
	}

	return false;
}

// ------------------------------------------------------------------------
void PointPicker::OnMouseMove()
{
	// if not moving we can do default behavoiur
	if(!isMoving)
	{
		return;
	}

	SetPointLocation();

	

}


// ------------------------------------------------------------------------
void PointPicker::OnLeftButtonUp()
{
	isMoving = false;
	this->pointPicker->SetTolerance(pickTolerance);
}


// ------------------------------------------------------------------------

vtkStandardNewMacro(PointPicker);
