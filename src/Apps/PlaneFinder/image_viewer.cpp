#include "image_viewer.h"

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkImageProperty.h>


// ------------------------------------------------------------------------
ImageViewer::ImageViewer(ImageDataList &images)
{
	this->images = images;
	currentIndex = 0;
	currentSlice = 0;
	currentTimeStep = 0;

	createViewerData();


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
	renderer->AddActor(this->viewerData[currentIndex].actors[currentTimeStep][currentSlice]);
	renderer->ResetCamera();
	style->SetRenderer(renderer);
	style->SetActor(this->viewerData[currentIndex].actors[currentTimeStep][currentSlice]);
	style->SetImage(this->images[currentIndex].images[currentTimeStep][currentSlice]);
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
			ViewerData::ImageVolumeMappers mappers;
			ViewerData::ImageVolumeActors actors;

			for(unsigned int k = 0; k < image.images[j].size(); k++)
			{
				vtkSmartPointer<vtkImageSliceMapper> mapper = 
					vtkSmartPointer<vtkImageSliceMapper>::New();
				mapper->SetInputData(image.images[j][k]);	
				mappers.push_back(mapper);

				vtkSmartPointer<vtkImageSlice> slice = 
					vtkSmartPointer<vtkImageSlice>::New();
				slice->SetMapper(mapper);
				slice->GetProperty()->SetColorLevel(177.48);
				slice->GetProperty()->SetColorWindow(512.04);
				actors.push_back(slice);
			}

			vData.mappers.push_back(mappers);
			vData.actors.push_back(actors);
		}

		this->viewerData.push_back(vData);
	}
}


// ------------------------------------------------------------------------
unsigned int ImageViewer::maxIndex()
{
	return this->viewerData.size() - 1;
}

// ------------------------------------------------------------------------
unsigned int ImageViewer::maxTimeStep()
{
	return this->viewerData[currentIndex].actors.size() - 1;
}


// ------------------------------------------------------------------------
unsigned int ImageViewer::maxSlice()
{
	return this->viewerData[currentIndex].actors[currentTimeStep].size() - 1;
}


// ------------------------------------------------------------------------
void ImageViewer::setViewedImage(unsigned int index) 
{
	if(index > this->maxIndex())
		currentIndex = this->maxIndex();
	else
		currentIndex = index;
}


// ------------------------------------------------------------------------
void ImageViewer::setViewedTimeStep(unsigned int timestep) 
{
	if(timestep > this->maxTimeStep())
		currentTimeStep = this->maxTimeStep();
	else
		currentTimeStep = timestep;
}

// ------------------------------------------------------------------------
void ImageViewer::setViewedSlice(unsigned int slice) 
{
	if(slice > this->maxSlice())
		currentSlice = this->maxSlice();
	else
		currentSlice = slice;
}
