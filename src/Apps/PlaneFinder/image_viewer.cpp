#include "image_viewer.h"

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkImageProperty.h>


// ------------------------------------------------------------------------
ImageViewer::ImageViewer(DataContainer * data)
{
	setData(data);

	this->widget = new QVTKWidget;

	style = vtkSmartPointer<PointPicker>::New();
	style->SetInteractionModeToImage2D();

	renderer = vtkSmartPointer<vtkRenderer>::New();
	widget->GetInteractor()->SetInteractorStyle(style);
	widget->GetRenderWindow()->AddRenderer(renderer);


}

// ------------------------------------------------------------------------
void ImageViewer::setData(DataContainer * data)
{
	this->data = data;
}


// ------------------------------------------------------------------------
void ImageViewer::updateImage(bool reset)
{
	renderer->RemoveAllViewProps();

	for(unsigned int i = 0; i < displayedLines.size(); i++)
	{
		std::cout << displayedLines.size() << std::endl;
		renderer->AddActor(displayedLines[i]);
	}

	renderer->AddActor(data->getActor());
	
	if(reset) resetCamera();
	style->SetRenderer(renderer);
	style->SetData(&data->getCurrentHolder());
	style->SetRenderWindow(widget->GetRenderWindow());
	widget->GetRenderWindow()->Render();
}

// ------------------------------------------------------------------------
void ImageViewer::resetCamera()
{
	renderer->ResetCamera();
}


// ------------------------------------------------------------------------
void ImageViewer::showAllLines()
{
	displayedLines = data->getLines();	
}

// ------------------------------------------------------------------------
void ImageViewer::setLineToDisplay(Line::Type type)
{
	displayedLines.clear();
	if(data->getLineData().count(type) > 0)
		displayedLines.push_back(data->getLineData()[type]->getActor());
}




