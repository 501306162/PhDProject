#include "image_viewer.h"

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>


// ------------------------------------------------------------------------
ImageViewer::ImageViewer(ImageDataList &images)
{
	this->images = images;
	createViewerData();

	this->widget = new QVTKWidget;

	style = vtkSmartPointer<vtkInteractorStyleImage>::New();
	style->SetInteractionModeToImage2D();

	renderer = vtkSmartPointer<vtkRenderer>::New();
	setViewedImage(0);
	widget->GetInteractor()->SetInteractorStyle(style);
	widget->GetRenderWindow()->AddRenderer(renderer);
}

// ------------------------------------------------------------------------
void ImageViewer::setViewedImage(unsigned int index)
{
	renderer->RemoveAllViewProps();
	renderer->AddActor(this->viewerData[index].imageSlices[0]);
	renderer->ResetCamera();
	widget->GetRenderWindow()->Render();
		
}

// ------------------------------------------------------------------------
void ImageViewer::createViewerData()
{
	for(unsigned int i = 0; i < this->images.size(); i++)
	{
		ImageData &image = this->images[i];
		ViewerData vData;
		for(unsigned int j = 0; j < image.images.size(); j++)
		{
			vtkSmartPointer<vtkImageSliceMapper> mapper = 
				vtkSmartPointer<vtkImageSliceMapper>::New();
			mapper->SetInputData(image.images[j]);	
			vData.imageMappers.push_back(mapper);

			vtkSmartPointer<vtkImageSlice> slice = 
				vtkSmartPointer<vtkImageSlice>::New();
			slice->SetMapper(mapper);
			vData.imageSlices.push_back(slice);
		}

		this->viewerData.push_back(vData);
	}
}


