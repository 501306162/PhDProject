#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <iostream>

#include <QMainWindow>
#include <QVTKWidget.h>
#include <QListWidget>
#include <QPushButton>

#include <itkImage.h>
#include <itkSimilarity2DTransform.h>

#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkImageData.h>
#include <vtkImageSliceMapper.h>
#include <vtkImageSlice.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor2D.h>


#include <ValveLine.h>
#include <ValveIO.h>


using namespace vt;

class MainWindow : public QMainWindow
{
	Q_OBJECT
public:
	
	typedef ValveSequence<3>::ImageType ImageType;
	typedef ImageType::PointType PointType;


	typedef itk::Image<unsigned short, 2> ImageSliceType;

	typedef struct display_line_
	{
		ValveLine<3>::Pointer valve;
		vtkSmartPointer<vtkPolyData> poly;
		vtkSmartPointer<vtkPolyDataMapper> polyMapper;
		vtkSmartPointer<vtkActor> polyActor;

		vtkSmartPointer<vtkImageData> image;
		vtkSmartPointer<vtkImageSliceMapper> imageMapper;
		vtkSmartPointer<vtkImageSlice> imageActor;
		vtkSmartPointer<vtkActor2D> labelActor;

	} DisplayLine;

	typedef std::vector<DisplayLine> LineSet;
	typedef std::vector<std::string> FilenamesType;
	typedef itk::Similarity2DTransform<double> TransformType;

	void showImages(const unsigned int index);

	MainWindow(const std::string &folder);
	MainWindow();
	~MainWindow() {}

	void init();

protected slots:
	void flipImagePressed();
	void flipLinePressed();
	void savePressed();

	void changeImage();

protected:
	void layout();
	void setUpViewer();
	void setUpList();
	void populateList();


	void loadData();
	void transformData();
	void transformLineSet(LineSet &lineSet);
	void getDataFiles(const std::string &folder, FilenamesType &files);
	//void getPoints(const std::string &filename, vtkImageData* image, Line &l, PointType &p1, PointType &p2);
	void loadVtkImage(const std::string &filename, int timeStep, vtkImageData * image);
	void convertPoint(vtkImageData * image, ImageType::Pointer &itkImage, double * in, double * out);
	void buildDisplayLine(const ValveSequence<3>::Pointer &valve, DisplayLine &line); 

	void computeTransform(ImageSliceType::PointType &p1x, ImageSliceType::PointType &p1y, 
		ImageSliceType::PointType &p2x, ImageSliceType::PointType &p2y,
		TransformType::Pointer &transform);


	void transformImage(ImageSliceType::Pointer &input, ImageSliceType::Pointer &ref, 
			TransformType::Pointer &transform, ImageSliceType::Pointer &output);
	
	void flattenImage(const ImageType::Pointer &input, ImageSliceType::Pointer &slice);

	void createVtkImage(const ImageSliceType::Pointer &input, vtkSmartPointer<vtkImageData> &output);
	void compute2DPoints(const ValveLine<3>::Pointer &valve, 
			ImageSliceType::PointType &op1, ImageSliceType::PointType &op2, 
			ImageSliceType::Pointer &image);


private:
	std::string dataFolder;
	LineSet data;

	QVTKWidget * viewer;
	vtkSmartPointer<vtkRenderer> renderer1;
	vtkSmartPointer<vtkRenderer> renderer2;

	QListWidget *nameList;

	std::string currentName;
	unsigned int currentIndex;
	
	QPushButton * flipLine;
	QPushButton * save;
	QPushButton * flipImage;
	
	std::vector<ValveSequence<3>::Pointer> inputValves;
	std::vector<DisplayLine> displayValves;


};


#endif
