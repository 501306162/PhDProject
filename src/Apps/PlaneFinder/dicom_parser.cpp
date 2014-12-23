#include "dicom_parser.h"


#include <QtGui>

#include <vtkImageData.h>
#include <vtkImageSliceMapper.h>
#include <vtkRenderer.h>
#include <vtkDICOMImageReader.h>
#include <vtkRenderWindow.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageProperty.h>


#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkGDCMImageIO.h>
#include <itkNrrdImageIO.h>
#include <itkJoinSeriesImageFilter.h>
#include <itkImageToVTKImageFilter.h>

// ------------------------------------------------------------------------
bool image_sort(const DicomImage &im1, const DicomImage &im2)
{
	return (im1.triggerTime < im2.triggerTime);
}


bool series_sort(const DicomSeries &s1, const DicomSeries &s2)
{

	return(std::atoi(s1.number.c_str()) < std::atoi(s2.number.c_str()));

	std::vector<double> pos1 = s1.images.front().imagePosition;
	std::vector<double> pos2 = s2.images.front().imagePosition;

	if(pos1.size() == 0)
		return false;

	std::vector<double> diff;
	double d1 = 0, d2 = 0;
	for(unsigned int i = 0; i < 3; i++)
	{
		d1 += (pos1[i]*pos1[i]);
		d2 += (pos2[i]*pos2[i]);
	}

	d1 = sqrt(d1);
	d2 = sqrt(d2);


	return (d1<d2);

}



// ------------------------------------------------------------------------
DicomParser::DicomParser(QWidget * parent) : QWidget(parent)
{
	outputFolder = "";
	folderName = "";
	
}

// ------------------------------------------------------------------------
void DicomParser::setFolderName(const QString &folderName)
{
	this->folderName = folderName;
}

// ------------------------------------------------------------------------
void DicomParser::setOutputFolder(const QString &folderName)
{
	this->outputFolder = folderName;
}

// ------------------------------------------------------------------------
void DicomParser::setUpDisplay()
{

	parseDicom();
	loadImages();


	QWidget * w1 = new QWidget(this);
	QHBoxLayout * l1 = new QHBoxLayout;
	l1->addWidget(createTable());
	l1->addWidget(createViewer());
	w1->setLayout(l1);
	

	
	okButton = new QPushButton("OK", this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(okPressed()));
	
	cancelButton = new QPushButton("Cancel", this);
	connect(cancelButton, SIGNAL(clicked()), this, SLOT(cancelPressed()));

	QHBoxLayout * l2 = new QHBoxLayout;
	l2->addWidget(cancelButton);
	l2->addWidget(okButton);

	QWidget * w2 = new QWidget(this);
	w2->setLayout(l2);

	QVBoxLayout * l3 = new QVBoxLayout;
	l3->addWidget(w1);
	l3->addWidget(w2);

	this->setLayout(l3);

}

// ------------------------------------------------------------------------
double DicomParser::posString(DicomImage &im)
{
	std::stringstream ss;
	double dist = 0;
	for(unsigned int i = 0; i < 3; i++)
	{
		dist += (im.imagePosition[i] * im.imagePosition[i]);
	}

	dist = sqrt(dist);

	return dist;
}


// ------------------------------------------------------------------------
void DicomParser::getSelectedSeries(DicomSeriesList & selectedSeries)
{
	for(unsigned int i = 0; i < seriesAgain.size(); i++)
	{
		QTableWidgetItem * item = table->item(i,4);
		if(item->checkState() == Qt::Checked)
			selectedSeries.push_back(seriesAgain[i]);		
	}
}

// ------------------------------------------------------------------------
void DicomParser::okPressed()
{
	if(outputFolder.isEmpty())
	{
		outputFolder = QFileDialog::getExistingDirectory(this, 
				"Choose output directory","/Users/oliverferoze/Uni/Data/ValveTracking/Testing");
	}

	if(outputFolder.isEmpty())
		return;


	DicomSeriesList outputSeries;
	getSelectedSeries(outputSeries);

	for(unsigned int i = 0; i < outputSeries.size(); i++)
	{
		writeSeries(outputSeries[i], outputFolder);
	}

	emit saveDone(outputFolder);
	this->close();

}

// ------------------------------------------------------------------------
void DicomParser::writeSeries(const DicomSeries &series, const QString &outputFolder)
{
	// create the series name
	QString name = QString::fromStdString(series.description);
	name.replace(" ","_");
	name += "_" + QString::number(series.images.front().seriesNumber);
	name += ".nrrd";
	
	QDir path(outputFolder);
	QString fullPath = path.absoluteFilePath(name);
	

	std::vector<DicomImage> images = series.images;
	std::sort(images.begin(), images.end());

	// write and build the output images 
	typedef itk::Image<unsigned short, 3> ImageType;
	typedef itk::Image<unsigned short, 4> OutputImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;
	typedef itk::JoinSeriesImageFilter<ImageType, OutputImageType> JoinerType;

	JoinerType::Pointer joiner = JoinerType::New();
	ImageType::Pointer orig;
	for(unsigned int i = 0; i < images.size(); i++)
	{
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(images[i].filename);
		std::cout << images[i].filename << std::endl;
		reader->SetImageIO(itk::GDCMImageIO::New());
		reader->Update();

		ImageType::Pointer im = reader->GetOutput();
		if(i == 0) orig = im;


		im->SetOrigin(orig->GetOrigin());
		im->SetDirection(orig->GetDirection());
		im->SetSpacing(orig->GetSpacing());

		joiner->SetInput(i, reader->GetOutput());
	}


	std::cout << joiner->GetOutput()->GetDirection() << std::endl;


	typedef itk::ImageFileWriter<OutputImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(joiner->GetOutput());
	writer->SetFileName(fullPath.toStdString());
	writer->SetImageIO(itk::NrrdImageIO::New());
	writer->Update();

}


// ------------------------------------------------------------------------
void DicomParser::cancelPressed()
{
	this->close();

}

// ------------------------------------------------------------------------
QWidget * DicomParser::createTable()
{
	// create the table and the headers
	table = new QTableWidget(this);

	QStringList headers;
	headers << "Series Name";
	headers << "Number";
	headers << "Number";
	headers << "Study";
	headers << "Include";

	table->setRowCount(this->series.size());
	table->setColumnCount(5);
	table->setHorizontalHeaderLabels(headers);
	table->horizontalHeader()->setResizeMode(QHeaderView::Stretch);
	table->setSelectionMode(QAbstractItemView::SingleSelection);
	table->setSelectionBehavior(QAbstractItemView::SelectRows);
	table->setShowGrid(false);



	double start = posString(seriesAgain.front().images.front());
	double last = 0.0;
	for(unsigned int i = 0; i < seriesAgain.size(); i++)
	{
		QString name = QString::fromStdString(seriesAgain[i].description);
		QTableWidgetItem * item = new QTableWidgetItem(name);		
		item->setFlags(item->flags()^Qt::ItemIsEditable);
		table->setItem(i,0,item);

		double pos = posString(seriesAgain[i].images[0]);
		double diff = pos - start;
		double comp = pos - last;
		last = pos;
		std::stringstream ssp;
		ssp << diff  << " - " << comp;
		QTableWidgetItem * number = new QTableWidgetItem(QString::number(seriesAgain[i].images[0].seriesNumber));
		QTableWidgetItem * num = new QTableWidgetItem(QString::fromStdString(ssp.str()));
		QTableWidgetItem * study = new QTableWidgetItem(QString::fromStdString(seriesAgain[i].images[0].studyUID));

		num->setFlags(item->flags());
		table->setItem(i,3, study);
		table->setItem(i,2,num);
		table->setItem(i,1,number);


		QTableWidgetItem * check = new QTableWidgetItem;
		check->setCheckState(Qt::Unchecked);
		check->setFlags(item->flags() ^ Qt::ItemIsUserCheckable);
		table->setItem(i,4, check);
	}

	table->resizeColumnToContents(3);

	return table;
}


// ------------------------------------------------------------------------
QWidget * DicomParser::createViewer()
{

	viewer = new QVTKWidget(this);
	viewer->setMinimumSize(150,150);


	renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(images.front());
	viewer->GetRenderWindow()->AddRenderer(renderer);

	vtkSmartPointer<vtkInteractorStyleImage> style = 
		vtkSmartPointer<vtkInteractorStyleImage>::New();
	style->SetInteractionModeToImage2D();

	viewer->GetInteractor()->SetInteractorStyle(style);
	viewer->GetRenderWindow()->Render();

	connect(table, SIGNAL(itemClicked(QTableWidgetItem*)), this, SLOT(selectItem(QTableWidgetItem*)));

	return viewer;
}


// ------------------------------------------------------------------------
void DicomParser::selectItem(QTableWidgetItem * item)
{
	// update the the image
	int index = item->row();
	
	renderer->RemoveAllViewProps();
	renderer->AddActor(images[index]);
	renderer->ResetCamera();

	viewer->GetRenderWindow()->Render();

	if(item->column() != 4) 
		return;

	// check the button 
	QTableWidgetItem * checkItem = table->item(index, 4);
	if(checkItem->checkState() == Qt::Checked)
		checkItem->setCheckState(Qt::Unchecked);
	else 
		checkItem->setCheckState(Qt::Checked);

}

// ------------------------------------------------------------------------
void DicomParser::loadImages()
{
	// group the series into sa and not sa 
	DicomSeriesList saList;
	DicomSeriesList otherList;

	for(unsigned int i = 0; i < series.size(); i++)
	{
		std::string description = series[i].description;
		QString desc = QString::fromStdString(description);
		
		if(desc.contains("Cine_SA"))
			saList.push_back(series[i]);
		else if(desc.contains("cine_SA"))
			saList.push_back(series[i]);
		else if(desc.contains("Cine_sa"))
			saList.push_back(series[i]);
		else if(desc.contains("Cine_ipat SA"))
			saList.push_back(series[i]);
		//else if(desc.contains("sa"))
			//saList.push_back(series[i]);
		else
			otherList.push_back(series[i]);
	}

	// sort the sa list based on position
	std::sort(saList.begin(), saList.end(), series_sort);
	std::sort(otherList.begin(), otherList.end(), series_sort);


	for(unsigned int i = 0; i < saList.size(); i++)
	{
		seriesAgain.push_back(saList[i]);
	}




	for(unsigned int i = 0; i < otherList.size(); i++)
	{
		seriesAgain.push_back(otherList[i]);
	}


	for(unsigned int i = 0; i < seriesAgain.size(); i++)
	{
		std::string filename = seriesAgain[i].images.front().filename;

		typedef itk::Image<unsigned short,3> ImageType;
		typedef itk::ImageFileReader<ImageType> ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetImageIO(itk::GDCMImageIO::New());
		reader->SetFileName(filename);

		typedef itk::ImageToVTKImageFilter<ImageType> VTKFilterType;
		VTKFilterType::Pointer vtkFilter = VTKFilterType::New();
		vtkFilter->SetInput(reader->GetOutput());
		vtkFilter->Update();


		vtkSmartPointer<vtkImageData> image = 
			vtkSmartPointer<vtkImageData>::New();
		image->DeepCopy(vtkFilter->GetOutput());
		image->Print(std::cout);
		

		vtkSmartPointer<vtkImageSliceMapper> mapper =
			vtkSmartPointer<vtkImageSliceMapper>::New();
		mapper->SetInputData(image);

		vtkSmartPointer<vtkImageSlice> actor = 
			vtkSmartPointer<vtkImageSlice>::New();
		actor->SetMapper(mapper);


		actor->GetProperty()->SetColorLevel(177.48);
		actor->GetProperty()->SetColorWindow(512.04);


		images.push_back(actor);


		
	}

}




// ------------------------------------------------------------------------
void DicomParser::parseDicom()
{
	SeriesExtractor extractor;

	std::string fname = folderName.toStdString();

	extractor.SetDirectory(fname);
	extractor.ExtractSeries();

	series = extractor.GetOutput();

}

