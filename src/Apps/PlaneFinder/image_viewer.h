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

	ImageViewer(ImageDataList &imageData);


	void setViewedImage(unsigned int index);
	void setViewedTimeStep(unsigned int timestep);
	void setViewedSlice(unsigned int slice);
	void updateImage();
	unsigned int maxIndex();
	unsigned int maxTimeStep();
	unsigned int maxSlice();

	QVTKWidget * getWidget() { return widget; }

protected:
	void createViewerData();

private:
	ImageDataList images;
	std::vector<ViewerData> viewerData;
	QVTKWidget * widget;
	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<PointPicker> style;

	unsigned int currentIndex;
	unsigned int currentSlice;
	unsigned int currentTimeStep;
};

#endif
