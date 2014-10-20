#ifndef IMAGE_VIEWER_H
#define IMAGE_VIEWER_H

#include "common.h"

#include <QVTKWidget.h>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkImageActor.h>
#include <vtkImageSlice.h>
#include <vtkImageSliceMapper.h>
#include <vtkRenderer.h>
#include <vtkInteractorStyleImage.h>

class ImageViewer 
{
public:

	typedef struct viewer_data_
	{
		std::vector<vtkSmartPointer< vtkImageSlice > > imageSlices;
		std::vector<vtkSmartPointer< vtkImageSliceMapper > > imageMappers;
	} ViewerData;

	ImageViewer(ImageDataList &imageData);


	void setViewedImage(unsigned int index);

	QVTKWidget * getWidget() { return widget; }

protected:
	void createViewerData();

private:
	ImageDataList images;
	std::vector<ViewerData> viewerData;
	QVTKWidget * widget;
	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkInteractorStyleImage> style;
};

#endif
