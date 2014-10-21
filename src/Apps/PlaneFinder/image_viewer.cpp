#include "image_viewer.h"

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkImageProperty.h>


// ------------------------------------------------------------------------
ImageViewer::ImageViewer(DataContainer * data)
{
	this->data = data;

	this->widget = new QVTKWidget;

	style = vtkSmartPointer<PointPicker>::New();
	style->SetInteractionModeToImage2D();

	renderer = vtkSmartPointer<vtkRenderer>::New();
	widget->GetInteractor()->SetInteractorStyle(style);
	widget->GetRenderWindow()->AddRenderer(renderer);

}

// ------------------------------------------------------------------------
void ImageViewer::updateImage()
{
	renderer->RemoveAllViewProps();
	renderer->AddActor(data->getActor());

	std::vector<vtkActor*> lines = data->getLines();
	for(unsigned int i = 0; i < lines.size(); i++)
	{
		renderer->AddActor(lines[i]);
	}

	renderer->ResetCamera();
	style->SetRenderer(renderer);
	style->SetActor(data->getActor());
	style->SetImage(data->getVTKImage());
	widget->GetRenderWindow()->Render();
}






