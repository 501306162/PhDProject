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
		displayedLines[i]->SetPosition(0,0,0.1);
		renderer->AddActor(displayedLines[i]);
	}

	renderer->AddActor(data->getActor());
	
	style->SetRenderer(renderer);
	style->SetData(&data->getCurrentHolder());
	style->SetRenderWindow(widget->GetRenderWindow());
	if(reset) resetCamera();
	widget->GetRenderWindow()->Render();
	
}

// ------------------------------------------------------------------------
void ImageViewer::resetCamera()
{
	std::cout << "Reseting camera" << std::endl;
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

	displayedLines = data->getLines();
	return;

	displayedLines.clear();
	if(data->getLineData().count(type) > 0)
	{
		displayedLines.push_back(data->getLineData()[type]->getActor());
	}
}




