#ifndef POINT_PICKER_H
#define POINT_PICKER_H

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageSlice.h>
#include <vtkImageData.h>
#include <vtkPointPicker.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkRenderWindow.h>
#include <vtkCellPicker.h>

#include "common.h"

class PointPicker : public vtkInteractorStyleImage 
{
public:
	static PointPicker * New();
	PointPicker();
	vtkTypeMacro(PointPicker, vtkInteractorStyleImage);

	virtual void OnLeftButtonDown();
	virtual void OnMouseMove();
	virtual void OnLeftButtonUp();

	void SetRenderer(vtkRenderer * renderer) { this->currentRenderer = renderer; }
	void SetRenderWindow(vtkRenderWindow * window) { this->window = window; }
	void SetData(DataHolder * data) { this->data = data; }
	void SetPickedData(vtkPolyData * data) { pickedData = data; }
	

private:
	vtkRenderer * currentRenderer;
	DataHolder * data;

	bool PointAtPickLocation();
   	void SetPointLocation();	
	vtkSmartPointer<vtkPointPicker> pointPicker;
	vtkSmartPointer<vtkPoints> ghostPoints;

	vtkSmartPointer<vtkPolyDataMapper> ghostMapper;
	vtkSmartPointer<vtkActor> ghostActor;
	vtkSmartPointer<vtkPolyData> ghostPoly;
	vtkSmartPointer<vtkVertexGlyphFilter> ghostFilter;
	vtkSmartPointer<vtkCellPicker> cellPicker;

	vtkRenderWindow * window;

	bool isMoving;
	int pickedPoint;
	vtkPolyData * pickedData;

};


#endif
