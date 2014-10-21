#include "main_window.h"

#include <QtCore>


// ------------------------------------------------------------------------
MainWindow::MainWindow(DataContainer *data) : data(data)
{
	// set up the window
	createActions();
	setCentralWidget(createLayout());
	setUpSignals();
}


// ------------------------------------------------------------------------
void MainWindow::setUpSignals()
{
	// image control slots
	connect(this->imageList->getImageList(), SIGNAL(itemSelectionChanged()),
		this, SLOT(imageSelectionChanged()));
	connect(this->zSlider, SIGNAL(valueChanged(int)),
		this, SLOT(imageSelectionChanged()));
	connect(this->tSlider, SIGNAL(valueChanged(int)),
		this, SLOT(imageSelectionChanged()));


	// button slots 
	connect(this->addLineButton, SIGNAL(pressed()), this, SLOT(addLinePressed()));

	//connect(this->mvButton, SIGNAL(pressed()), this, SLOT(toggleValveType()));
	//connect(this->tpButton, SIGNAL(pressed()), this, SLOT(toggleValveType()));
	//connect(this->avButton, SIGNAL(pressed()), this, SLOT(toggleValveType()));
}


// ------------------------------------------------------------------------
void MainWindow::addLinePressed()
{
	Line * line = Line::NewLine(data->getVTKImage(), 
			Line::getType(valveGroup->checkedId()));
	data->addLine(line);
	imageViewer->updateImage();
}



// ------------------------------------------------------------------------
QPushButton * MainWindow::getCheckedValveButton()
{
	if(mvButton->isChecked())
		return mvButton;
	if(avButton->isChecked())
		return avButton;
	if(tpButton->isChecked())
		return tpButton;
}

// ------------------------------------------------------------------------
void MainWindow::imageSelectionChanged()
{
	data->setViewedImage(imageList->selectedIndex());
	updateSliders();
	data->setViewedTimeStep(tSlider->value());
	data->setViewedSlice(zSlider->value());
	imageViewer->updateImage();

}

// ------------------------------------------------------------------------
void MainWindow::updateSliders()
{
	unsigned int maxSlice = data->maxSlice();
	unsigned int maxTimeStep = data->maxTimeStep();
	unsigned int currentSlice = zSlider->value();
	unsigned int currentTimeStep = tSlider->value();

	if(currentSlice > maxSlice)
		currentSlice = maxSlice;
	if(currentTimeStep > maxTimeStep)
		currentTimeStep = maxTimeStep;

	tSlider->setMaximum(maxTimeStep);
	tSlider->setValue(currentTimeStep);
	zSlider->setMaximum(maxSlice);
	zSlider->setValue(currentSlice);


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
	addLineButton = new QPushButton("Add Line");
	removeLineButton = new QPushButton("Remove Line");
	propagateLineButton = new QPushButton("Propagate Line");
	
	mvButton = new QPushButton("Mitral Valve");
	tpButton = new QPushButton("Tricuspid Valve");
	avButton = new QPushButton("Auortic Valve");
	mvButton->setCheckable(true);
	tpButton->setCheckable(true);
	avButton->setCheckable(true);

	mvButton->setChecked(true);

	valveGroup = new QButtonGroup;
	valveGroup->addButton(mvButton,0);
	valveGroup->addButton(tpButton,1);
	valveGroup->addButton(avButton,2);
	valveGroup->setExclusive(true);


	QVBoxLayout * layout = new QVBoxLayout;
	layout->addWidget(addLineButton);
	layout->addWidget(removeLineButton);
	layout->addWidget(propagateLineButton);

	QFrame * line = new QFrame;
	line->setFrameShape(QFrame::HLine);
	line->setFrameShadow(QFrame::Sunken);
	layout->addWidget(line);


	layout->addWidget(mvButton);
	layout->addWidget(tpButton);
	layout->addWidget(avButton);
	
	buttonBox->setLayout(layout);
	return buttonBox;
}

// ------------------------------------------------------------------------
void MainWindow::closeEvent(QCloseEvent * event)
{
	delete data;
}

// ------------------------------------------------------------------------
void MainWindow::createActions()
{
	
}

// ------------------------------------------------------------------------
QWidget * MainWindow::createImageList()
{
	this->imageList = new ImageListDisplay(this->data);
	this->imageList->setUp();
	return this->imageList->getImageList();
}

// ------------------------------------------------------------------------
QWidget * MainWindow::createImageViewer()
{
	imageViewer = new ImageViewer(this->data);
	imageViewer->getWidget()->setMinimumSize(QSize(500,500));
	imageSelectionChanged();
	return imageViewer->getWidget();
}


// ------------------------------------------------------------------------
MainWindow::~MainWindow()
{
	delete this->imageList;
	delete this->imageViewer;
}
