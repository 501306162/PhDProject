#include "main_window.h"

#include <QtGui>
#include <QtCore>

#include <qjson/parser.h>

#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkImageProperty.h>
#include <vtkInteractorStyleImage.h>

#include <CommonDefinitions.h>

#include <itkSimilarity2DTransform.h>
#include <itkLandmarkBasedTransformInitializer.h>
#include <itkImageToVTKImageFilter.h>
#include <itkRegionOfInterestImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkFlipImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkNrrdImageIO.h>
#include <itkPermuteAxesImageFilter.h>

#include <Directory.h>

// ------------------------------------------------------------------------
MainWindow::MainWindow()
{
	this->init();
}


// ------------------------------------------------------------------------
MainWindow::MainWindow(const std::string &folder) : 
	dataFolder(folder)
{
	this->init();
}	


// ------------------------------------------------------------------------
void MainWindow::savePressed()
{
	QString dir = QFileDialog::getExistingDirectory(this, "Get Dir", "/Users/oliverferoze/Uni/Data");

	for(unsigned int i = 0; i < inputValves.size(); i++)
	{

		ValveSequence<3>::Pointer v = inputValves[i];
		ValveSequenceWriter<3>::Pointer writer = ValveSequenceWriter<3>::New();
		QDir dirPath(dir);
		QString fpath = dirPath.absoluteFilePath(nameList->item(i)->text().replace(".txt", ""));

		writer->SetFileName(fpath.toStdString());
		writer->SetInput(v);
		writer->Write();
	}

	

}

// ------------------------------------------------------------------------
void MainWindow::init()
{
	layout();
	loadData();

	this->nameList->setItemSelected(this->nameList->item(0), true);
	showImages(0);
}



// ------------------------------------------------------------------------
void MainWindow::changeImage()
{
	int row = this->nameList->row(this->nameList->selectedItems().first());
	showImages(row);

}


// -----------------------------------------------------------------------
void MainWindow::layout()
{
	QVBoxLayout * vlayout = new QVBoxLayout;
	QHBoxLayout * hlayout = new QHBoxLayout;
	QHBoxLayout * buttonLayout = new QHBoxLayout;
	setUpViewer();
	setUpList();
	

	hlayout->addWidget(this->nameList);
	hlayout->addWidget(this->viewer);
	QWidget * hwidget = new QWidget(this);
	hwidget->setLayout(hlayout);

	vlayout->addWidget(hwidget);


	this->flipImage = new QPushButton("Flip Image", this);
	this->flipLine = new QPushButton("Flip Line", this);
	this->save = new QPushButton("Save", this);
	buttonLayout->addWidget(this->flipImage);
	buttonLayout->addWidget(this->flipLine);
	buttonLayout->addWidget(this->save);


	connect(this->flipImage, SIGNAL(pressed()), this, SLOT(flipImagePressed()));
	connect(this->flipLine, SIGNAL(pressed()), this, SLOT(flipLinePressed()));
	connect(this->save, SIGNAL(pressed()), this, SLOT(savePressed()));
	connect(this->nameList, SIGNAL(itemSelectionChanged()), this, SLOT(changeImage()));

	QWidget * bwidget = new QWidget(this);
	bwidget->setLayout(buttonLayout);

	vlayout->addWidget(bwidget);
	
	QWidget * central = new QWidget(this);
	central->setLayout(vlayout);
	this->setCentralWidget(central);
}

// ------------------------------------------------------------------------
void MainWindow::setUpList()
{
	this->nameList = new QListWidget(this);
}


// ------------------------------------------------------------------------
void MainWindow::flipLinePressed()
{
	int row = nameList->row(nameList->selectedItems().front());
	ValveSequence<3>::Pointer v1 = inputValves[row];
	for(unsigned int i = 0; i < v1->GetNumberOfLines(); i++)
	{
		ValveLine<3>::PointType temp = v1->GetValveLine(i)->GetP1();		
		v1->GetValveLine(i)->SetP1(v1->GetValveLine(i)->GetP2());
		v1->GetValveLine(i)->SetP2(temp);
		v1->GetValveLine(i)->UpdateIndexs();
	}

	DisplayLine newLine;
	buildDisplayLine(v1, newLine);
	displayValves[row] = newLine;


	showImages(row);




}

// ------------------------------------------------------------------------
void MainWindow::flipImagePressed()
{
	int row = nameList->row(nameList->selectedItems().front());
	ValveSequence<3>::Pointer v1 = inputValves[row];
	for(unsigned int i = 0; i < v1->GetNumberOfLines(); i++)
	{
		
		// flip the points as well
		ValveLine<3>::ContIndexType ind1 = v1->GetValveLine(i)->GetInd1();
		ValveLine<3>::ContIndexType ind2 = v1->GetValveLine(i)->GetInd2();

		double tmp1, tmp2;
		tmp1 = ind1[0];
		tmp2 = ind2[0];

		ind1[0] = ind1[1];
		ind1[1] = tmp1;

		ind2[0] = ind2[1];
		ind2[1] = tmp2;

		v1->GetValveLine(i)->SetInd1(ind1);
		v1->GetValveLine(i)->SetInd2(ind2);
		v1->GetValveLine(i)->UpdatePoints();



		ImageType::Pointer im = v1->GetValveLine(i)->GetImage();		
		typedef itk::PermuteAxesImageFilter<ImageType> PermType;
		PermType::Pointer perm = PermType::New();
		perm->SetInput(im);


		itk::FixedArray<unsigned int,  3> order;
		order[0] = 1;
		order[1] = 0;
		order[2] = 2;

		perm->SetOrder(order);
		perm->Update();

		ImageType::Pointer newImage = perm->GetOutput();
		newImage->SetDirection(im->GetDirection());

		v1->GetValveLine(i)->SetImage(newImage);


		
	}

	DisplayLine newLine;
	buildDisplayLine(v1, newLine);
	displayValves[row] = newLine;

	showImages(row);
}



// ------------------------------------------------------------------------
void MainWindow::showImages(const unsigned int index)
{
	DisplayLine &line = displayValves[index];

	this->renderer1->RemoveAllViewProps();
	this->renderer2->RemoveAllViewProps();
	this->renderer1->AddActor(displayValves[0].imageActor);
	this->renderer1->AddActor(displayValves[0].polyActor);
	this->renderer2->AddActor(line.imageActor);
	this->renderer2->AddActor(line.polyActor);
	this->renderer1->ResetCamera();
	this->renderer2->SetActiveCamera(this->renderer1->GetActiveCamera());

	this->viewer->GetRenderWindow()->Render();

}


// ------------------------------------------------------------------------
void MainWindow::setUpViewer()
{
	this->viewer = new QVTKWidget(this);
	this->viewer->setMinimumWidth(600);
	this->viewer->setMinimumHeight(400);

	this->renderer1 = vtkSmartPointer<vtkRenderer>::New();
	this->renderer1->SetViewport(0.0,0.0,0.5,1.0);
	this->renderer2 = vtkSmartPointer<vtkRenderer>::New();
	this->renderer2->SetViewport(0.5,0.0,1.0,1.0);

	this->viewer->GetRenderWindow()->AddRenderer(this->renderer1);
	this->viewer->GetRenderWindow()->AddRenderer(this->renderer2);

	vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
	this->viewer->GetInteractor()->SetInteractorStyle(style);

	
}


// ------------------------------------------------------------------------
void MainWindow::computeTransform(ImageSliceType::PointType &p1x, ImageSliceType::PointType &p1y, 
		ImageSliceType::PointType &p2x, ImageSliceType::PointType &p2y,
		TransformType::Pointer &transform)
{
	typedef itk::LandmarkBasedTransformInitializer<TransformType, 
			ImageSliceType, ImageSliceType> InitialiserType;
	typedef InitialiserType::LandmarkPointContainer ContainerType;
	typedef InitialiserType::LandmarkPointType LandmarkType;

	ContainerType fixedLandmarks;
	ContainerType movingLandmarks;

	fixedLandmarks.push_back(p1x);
	fixedLandmarks.push_back(p1y);

	movingLandmarks.push_back(p2x);
	movingLandmarks.push_back(p2y);


	InitialiserType::Pointer initialiser = InitialiserType::New();
	initialiser->SetFixedLandmarks(fixedLandmarks);
	initialiser->SetMovingLandmarks(movingLandmarks);

	transform->SetIdentity();
	initialiser->SetTransform(transform);
	initialiser->InitializeTransform();
}


// ------------------------------------------------------------------------
void MainWindow::loadData()
{
	FilenamesType dataFiles;
	this->getDataFiles(this->dataFolder, dataFiles);

	for(unsigned int i = 0; i < dataFiles.size(); i++)
	{
		std::cout << dataFiles[i] << std::endl;

		ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
		reader->SetFileName(dataFiles[i]);

		QFileInfo file(QString::fromStdString(dataFiles[i]));
		QListWidgetItem * name = new QListWidgetItem(file.fileName());
		this->nameList->addItem(name);
			
		ValveSequence<3>::Pointer valve = reader->GetOutput();
		inputValves.push_back(valve);

		DisplayLine displayValve;
		buildDisplayLine(valve, displayValve);

		displayValves.push_back(displayValve);
	}

}

// ------------------------------------------------------------------------
void MainWindow::buildDisplayLine(const ValveSequence<3>::Pointer &valve, DisplayLine &line)
{
	ValveLine<3>::Pointer v1 = inputValves.front()->GetValveLine(0);
	ValveLine<3>::Pointer v2 = valve->GetValveLine(0);

	// create the 2d images
	ImageSliceType::Pointer im1 = ImageSliceType::New();
	ImageSliceType::Pointer im2 = ImageSliceType::New();
	flattenImage(v1->GetImage(), im1);
	flattenImage(v2->GetImage(), im2);


	// create the 2D points
	ImageSliceType::PointType p1x, p1y, p2x, p2y;
	compute2DPoints(v1, p1x, p1y, im1);
	compute2DPoints(v2, p2x, p2y, im2);


	std::cout << p1x << " " << p1y << " " << p2x << " " << p2y << std::endl;


	// compute the transform
	TransformType::Pointer transform = TransformType::New();
	computeTransform(p1x, p1y, p2x, p2y, transform);

	//transform->Print(std::cout);
	
	// apply the transform
	ImageSliceType::Pointer output = ImageSliceType::New();
	transformImage(im2, im1, transform, output);

	//utils::ImageIO<ImageSliceType>::Write("test.nrrd", output);
	//exit(1);

	line.image = vtkSmartPointer<vtkImageData>::New();
	createVtkImage(output, line.image);

	line.imageMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
	line.imageMapper->SetInputData(line.image);
	
	line.imageActor = vtkSmartPointer<vtkImageSlice>::New();
	line.imageActor->SetMapper(line.imageMapper);


	line.imageActor->GetProperty()->SetColorLevel(177.48);
	line.imageActor->GetProperty()->SetColorWindow(512.04);



	// create the line
	line.poly = valve->GetValveLine(0)->GetPolyData();
	line.polyMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	line.polyMapper->SetInputData(line.poly);
	line.polyActor = vtkSmartPointer<vtkActor>::New();
	line.polyActor->SetMapper(line.polyMapper);
	


}	


// ------------------------------------------------------------------------
void MainWindow::compute2DPoints(const ValveLine<3>::Pointer &valve, 
			ImageSliceType::PointType &op1, ImageSliceType::PointType &op2, 
			ImageSliceType::Pointer &image)
{
	typedef itk::ContinuousIndex<double, 2> IndType;
	IndType ind1, ind2;
	for(unsigned int i = 0; i < 2; i++)
	{
		ind1[i] = valve->GetInd1()[i];
		ind2[i] = valve->GetInd2()[i];
	}

	image->TransformContinuousIndexToPhysicalPoint(ind1, op1);
	image->TransformContinuousIndexToPhysicalPoint(ind2, op2);
}



// ------------------------------------------------------------------------
void MainWindow::createVtkImage(const ImageSliceType::Pointer &input, 
		vtkSmartPointer<vtkImageData> &output)
{
	typedef itk::FlipImageFilter<ImageSliceType> FlipperType;
	FlipperType::Pointer flipper = FlipperType::New();
	flipper->SetInput(input);

	FlipperType::FlipAxesArrayType axes;
	axes[0] = false;
	axes[1] = true;

	flipper->SetFlipAxes(axes);

	typedef itk::ImageToVTKImageFilter<ImageSliceType> ConverterType;
	ConverterType::Pointer converter = ConverterType::New();
	converter->SetInput(flipper->GetOutput());
	converter->Update();

	output->DeepCopy(converter->GetOutput());
}


// ------------------------------------------------------------------------
void MainWindow::flattenImage(const ImageType::Pointer &input, ImageSliceType::Pointer &slice)
{
	typedef itk::ExtractImageFilter<ImageType, ImageSliceType> ExtractorType;
	ExtractorType::Pointer extractor = ExtractorType::New();

	ImageType::RegionType exRegion;
	ImageType::SizeType exSize = input->GetLargestPossibleRegion().GetSize();
	ImageType::IndexType exIndex = input->GetLargestPossibleRegion().GetIndex();

	exSize[2] = 0;
	exIndex[2] = 0;
	exRegion.SetSize(exSize);
	exRegion.SetIndex(exIndex);

	extractor->SetExtractionRegion(exRegion);
	extractor->SetInput(input);
	extractor->SetDirectionCollapseToIdentity();

	extractor->Update();
	slice = extractor->GetOutput();



}


// ------------------------------------------------------------------------
void MainWindow::transformImage(ImageSliceType::Pointer &input, ImageSliceType::Pointer &ref, 
			TransformType::Pointer &transform, ImageSliceType::Pointer &output)
{
	typedef itk::ResampleImageFilter<ImageSliceType, ImageSliceType> ResamplerType;
	ResamplerType::Pointer resampler = ResamplerType::New();
	resampler->SetInput(input);
	resampler->SetTransform(transform);
	resampler->SetUseReferenceImage(true);
	resampler->SetOutputParametersFromImage(ref);

	resampler->Update();

	output = resampler->GetOutput();
}



// ------------------------------------------------------------------------
void MainWindow::getDataFiles(const std::string &folder, FilenamesType &files)
{
	utils::Directory::Pointer dir = utils::Directory::New();
	dir->SetDirectory(folder);
	dir->SetExtension(".txt");
	files = dir->GetOutput();
}


