
#include <QApplication>
#include <QVTKWidget.h>

#include "main_window.h"






// ------------------------------------------------------------------------
int main(int argc, char ** argv)
{

	QApplication app(argc, argv);

	MainWindow window(argv[1]);
	window.show();
	return app.exec();

	

	/*


	std::string dataFolder = argv[1];
	std::vector<std::string> filenames;
	getDataFiles(dataFolder, filenames);
	




	LineData::iterator lineIt = data.begin();
	while(lineIt != data.end())
	{
		std::cout << lineIt->firjjst << " " << lineIt->second.size() << std::endl;
		//pixelAtPoint(lineIt->second.front().lines.front().p1,lineIt->second.front().lines.front().p2, lineIt->second.front().imageFilename);
		checkLineSet(lineIt->second);
		++lineIt;
	}
	*/


	return 0;
}

/*

// ------------------------------------------------------------------------
void checkLineSet(LineSet &lines)
{
	typedef ImageType::PointType PointType;

	
	exit(1);
	
}

// ------------------------------------------------------------------------
void getPoints(const std::string &filename, Line &line, ImageType::PointType &p1, ImageType::PointType &p2)
{


}

// ------------------------------------------------------------------------
void view(vtkImageData * im1, vtkImageData * im2)
{
	
	vtkSmartPointer<vtkRenderer> ren1 = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkImageSliceMapper> map1 = vtkSmartPointer<vtkImageSliceMapper>::New();
	vtkSmartPointer<vtkImageSlice> actor1 = vtkSmartPointer<vtkImageSlice>::New();
	//ren1->SetViewport(0,0,0.5,1.0);
	map1->SetInputData(im1);
	actor1->SetMapper(map1);
	ren1->AddViewProp(actor1);


	vtkSmartPointer<vtkRenderer> ren2 = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkImageSliceMapper> map2 = vtkSmartPointer<vtkImageSliceMapper>::New();
	vtkSmartPointer<vtkImageSlice> actor2 = vtkSmartPointer<vtkImageSlice>::New();
	window->AddRenderer(ren2);
	ren2->SetViewport(0.5,0,1.0,1.0);
	map2->SetInputData(im2);
	actor2->SetMapper(map2);
	ren2->AddActor(actor2);
	ren2->ResetCamera();


	vtkSmartPointer<vtkRenderWindow> window = vtkSmartPointer<vtkRenderWindow>::New();
	window->AddRenderer(ren1);


	vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
	style->SetInteractionModeToImage2D();

	interactor->SetInteractorStyle(style);
	interactor->SetRenderWindow(window);
	interactor->Initialize();
	std::cout << "Hello" << std::endl;

	interactor->Start();


}


// ------------------------------------------------------------------------
unsigned short pixelAtPoint(double * p, double *p2, std::string &filename)
{
	// load the vtk image
	vtkSmartPointer<vtkImageData> vtkImage = vtkSmartPointer<vtkImageData>::New();
	loadVtkImage(filename, vtkImage);



	// load the itk image
	typedef itk::Image<unsigned short, 3> ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();

	itk::NrrdImageIO::Pointer imageIO = itk::NrrdImageIO::New();

	reader->SetImageIO(imageIO);
	reader->SetFileName(filename);
	reader->Update();
	ImageType::Pointer image = reader->GetOutput();


	double p1out[3];
	double p2out[3];

	convertPoint(vtkImage, image, p, p1out);
	convertPoint(vtkImage, image, p2, p2out);


	vtkSmartPointer<vtkLineSource> line = vtkSmartPointer<vtkLineSource>::New();
	line->SetPoint1(p1out);
	line->SetPoint2(p2out);
	line->Update();

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetInputData(line->GetOutput());
	writer->SetFileName("line.vtk");
	writer->Update();
	writer->Write();


	utils::ImageVolumeIO::Write("image.nrrd", image);

	std::cout << filename << std::endl;

	return 0;
	return 0;

}

// ------------------------------------------------------------------------
void convertPoint(vtkImageData * vtkImage, ImageType::Pointer &image, double * input, double * output)
{
}



// ------------------------------------------------------------------------
void getDataFiles(const std::string &folder, std::vector<std::string> &filenames)
{
}

// ------------------------------------------------------------------------
void loadVtkImage(const std::string &filename, vtkImageData * image)
{


}


*/
