#ifndef POINT_PICKER_H
#define POINT_PICKER_H

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageSlice.h>
#include <vtkImageData.h>

#include "common.h"

class PointPicker : public vtkInteractorStyleImage 
{
public:
	static PointPicker * New();
	vtkTypeMacro(PointPicker, vtkInteractorStyleImage);

	virtual void OnLeftButtonDown();

	void SetRenderer(vtkRenderer * renderer) { this->currentRenderer = renderer; }
	void SetData(DataHolder * data) { this->data = data; }

private:
	vtkRenderer * currentRenderer;
	DataHolder * data;
};


#endif
