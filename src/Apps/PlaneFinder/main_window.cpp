#include "main_window.h"

#include <QtCore>

// ------------------------------------------------------------------------
MainWindow::MainWindow(ImageDataList &imageData)
{
	// set up the window
	setImageData(imageData);
	createActions();
	setCentralWidget(createLayout());
	setUpSignals();

}


// ------------------------------------------------------------------------
void MainWindow::setUpSignals()
{
	connect(this->imageList->getImageList(), SIGNAL(itemSelectionChanged()),
		this, SLOT(imageSelectionChanged()));

}

// ------------------------------------------------------------------------
void MainWindow::imageSelectionChanged()
{
	imageViewer->setViewedImage( imageList->selectedIndex());


}

// ------------------------------------------------------------------------
void MainWindow::setTimeStep(unsigned int currentSelection)
{

}



// ------------------------------------------------------------------------
QWidget * MainWindow::createImageControls()
{
	QGridLayout *layout = new QGridLayout;
	QLabel * tLabel = new QLabel("Time Step: ");
	QLabel * zLabel = new QLabel("Slice: ");
	tSlider = new QSlider;
	tSlider->setOrientation(Qt::Horizontal);
	zSlider = new QSlider;
	zSlider->setOrientation(Qt::Horizontal);
	
	layout->addWidget(zLabel, 0, 0);
	layout->addWidget(tLabel, 1, 0);
	layout->addWidget(zSlider, 0, 1);
	layout->addWidget(tSlider, 1, 1);

	QWidget * widget = new QWidget;
	widget->setLayout(layout);
	return widget;

}

// ------------------------------------------------------------------------
QWidget * MainWindow::createLayout()
{
	QWidget * widget = new QWidget();
	QWidget * leftWidget = new QWidget();
	QVBoxLayout * leftLayout = new QVBoxLayout();

	// on the left we have the image list and the buttons
	leftLayout->addWidget(createImageList());
	leftLayout->addWidget(createImageControls());
	leftLayout->addWidget(createButtonGroup());
	leftWidget->setLayout(leftLayout);
	leftWidget->setMaximumWidth(300);



	QHBoxLayout * mainLayout = new QHBoxLayout;
	mainLayout->addWidget(leftWidget);
	mainLayout->addWidget(createImageViewer());
	

	widget->setLayout(mainLayout);

	return widget;

}


// ------------------------------------------------------------------------
QWidget * MainWindow::createButtonGroup()
{
	QGroupBox * buttonBox = new QGroupBox;
	mvButton = new QPushButton("Mitral Valve");
	tpButton = new QPushButton("Tricuspid Valve");
	avButton = new QPushButton("Auortic Valve");

	QVBoxLayout * layout = new QVBoxLayout;
	layout->addWidget(mvButton);
	layout->addWidget(tpButton);
	layout->addWidget(avButton);
	
	buttonBox->setLayout(layout);
	return buttonBox;
}

// ------------------------------------------------------------------------
void MainWindow::closeEvent(QCloseEvent * event)
{

}

// ------------------------------------------------------------------------
void MainWindow::setImageData(ImageDataList &imageData)
{
	this->images = imageData;
}



// ------------------------------------------------------------------------
void MainWindow::createActions()
{
	
}

// ------------------------------------------------------------------------
QWidget * MainWindow::createImageList()
{
	this->imageList = new ImageListDisplay(this->images);
	this->imageList->setUp();
	return this->imageList->getImageList();
}

// ------------------------------------------------------------------------
QWidget * MainWindow::createImageViewer()
{
	imageViewer = new ImageViewer(this->images);
	imageViewer->getWidget()->setMinimumSize(QSize(500,500));
	return imageViewer->getWidget();
}


// ------------------------------------------------------------------------
MainWindow::~MainWindow()
{
	delete this->imageList;
	delete this->imageViewer;
}
