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

#include "point_picker.h"
#include "containers.h"

class ImageViewer 
{
public:

	typedef struct viewer_data_
	{
		typedef std::vector<vtkSmartPointer<vtkImageSlice> > ImageVolumeActors;
		typedef std::vector<ImageVolumeActors> ImageSequenceActors;
		ImageSequenceActors actors;

		typedef std::vector<vtkSmartPointer<vtkImageSliceMapper> > ImageVolumeMappers;
		typedef std::vector<ImageVolumeMappers> ImageSequenceMappers;
		ImageSequenceMappers mappers;

	} ViewerData;

	ImageViewer(DataContainer * data);

	void resetCamera();
	void setData(DataContainer * data);


	void updateImage(bool reset = false);
	void showAllLines();
	void setLineToDisplay(Line::Type type);

	QVTKWidget * getWidget() { return widget; }

protected:
	void createViewerData();

private:
	DataContainer * data;
	QVTKWidget * widget;
	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<PointPicker> style;
	std::vector<vtkActor*> displayedLines;

};

#endif
