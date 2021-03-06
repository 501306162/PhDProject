#include "main_window.h"


// ------------------------------------------------------------------------
MainWindow::MainWindow(DataContainer *data)
{
	setUpWindow();
	initialiseFromData(data);
}

// ------------------------------------------------------------------------
MainWindow::MainWindow() : data(NULL)
{
	setUpWindow();
}



// ------------------------------------------------------------------------
void MainWindow::createActions()
{
	fileMenu = this->menuBar()->addMenu(tr("&File"));

	newAction = new QAction(tr("&New"), this);
	newAction->setShortcut(QKeySequence::New);
	newAction->setToolTip("Create a new session");
	fileMenu->addAction(newAction);


	loadAction = new QAction(tr("&Load"), this);
	loadAction->setShortcut(QKeySequence::Open);
	loadAction->setToolTip("Load a new set of lines");
	fileMenu->addAction(loadAction);


	saveAction = new QAction(tr("&Save"), this);
	saveAction->setShortcut(QKeySequence::Save);
	saveAction->setToolTip("Save the current session");
	fileMenu->addAction(saveAction);

	saveAsAction = new QAction(tr("&Save As"), this);
	saveAsAction->setShortcut(QKeySequence::SaveAs);
	saveAsAction->setToolTip("Save the current Session as");
	fileMenu->addAction(saveAsAction);
	fileMenu->addSeparator();
		
	processDicomAction = new QAction(tr("Process Dicom"), this);
	processDicomAction->setShortcut(QKeySequence::Print);
	processDicomAction->setToolTip("Process a dicom examination");
	fileMenu->addAction(processDicomAction);


	patchExtractorAction = new QAction(tr("Extract Patches"), this);
	patchExtractorAction->setToolTip("Extract the patches");
	fileMenu->addAction(patchExtractorAction);


	exportVideos = new QAction(tr("Extract Videos"), this);
	fileMenu->addAction(exportVideos);
	
}


// ------------------------------------------------------------------------
void MainWindow::setUpWindow()
{
	createActions();
	QWidget * centralWidget = createLayout();
	setCentralWidget(centralWidget);
	setUpSignals();
}


// ------------------------------------------------------------------------
void MainWindow::setUpSignals()
{
	// image control slots
	connect(this->imageList->getImageList(), SIGNAL(itemSelectionChanged()),
		this, SLOT(imageSelectionChanged()));
	connect(this->zSlider, SIGNAL(valueChanged(int)),
		this, SLOT(sliderSelectionChanged()));
	connect(this->tSlider, SIGNAL(valueChanged(int)),
		this, SLOT(sliderSelectionChanged()));

	// button slots 
	connect(this->addLineButton, SIGNAL(pressed()), this, SLOT(addLinePressed()));
	connect(this->removeLineButton, SIGNAL(pressed()), this, SLOT(removeLinePressed()));
	connect(this->propagateLineButton, SIGNAL(pressed()), this, SLOT(propagateLinePressed()));
	connect(this->lockLineButton, SIGNAL(pressed()), this, SLOT(lockLine()));
	connect(this->lineList->getWidget(), SIGNAL(itemSelectionChanged()), 
			this, SLOT(lineChanged()));
	connect(this->controls, SIGNAL(rightPressed()), this, SLOT(tSliderRight()));
	connect(this->controls, SIGNAL(leftPressed()), this, SLOT(tSliderLeft()));
	connect(this->controls, SIGNAL(lockPressed()), this, SLOT(lockLine()));
	connect(this->controls, SIGNAL(propagatePressed()), this, SLOT(propagateLinePressed()));
	connect(this->imageTypeButton, SIGNAL(pressed()), this, SLOT(imageTypePressed()));

	// connect the load and save 
	connect(this->saveAction, SIGNAL(triggered()), this, SLOT(saveActionPressed()));
	connect(this->newAction, SIGNAL(triggered()), this, SLOT(newActionPressed()));
	connect(this->loadAction, SIGNAL(triggered()), this, SLOT(loadActionPressed()));
	connect(this->saveAsAction, SIGNAL(triggered()), this, SLOT(saveAsActionPressed()));
	connect(this->processDicomAction, SIGNAL(triggered()), this, SLOT(processDicomPressed()));
	connect(this->exportVideos, SIGNAL(triggered()), this, SLOT(exportVideosPressed()));
}



// ------------------------------------------------------------------------
QWidget * MainWindow::createInfoPane()
{


	QGroupBox * box = new QGroupBox("Data Info", this);
	QLabel * folderName = new QLabel("Folder: ", this);
	folderValue = new QLabel("", this);

	QGridLayout * layout = new QGridLayout;
	layout->addWidget(folderName, 1, 0);
	layout->addWidget(folderValue, 1, 1);
	box->setLayout(layout);

	return box;
}



// ------------------------------------------------------------------------
QWidget * MainWindow::createImageControls()
{
	QGridLayout *layout = new QGridLayout;
	QLabel * tLabel = new QLabel("Time Step: ", this);
	QLabel * zLabel = new QLabel("Slice: ", this);
	tSlider = new QSlider(this);
	tSlider->setOrientation(Qt::Horizontal);
	zSlider = new QSlider(this);
	zSlider->setOrientation(Qt::Horizontal);
	
	layout->addWidget(zLabel, 0, 0);
	layout->addWidget(tLabel, 1, 0);
	layout->addWidget(zSlider, 0, 1);
	layout->addWidget(tSlider, 1, 1);

	QWidget * widget = new QWidget(this);
	widget->setLayout(layout);
	return widget;

}

// ------------------------------------------------------------------------
QWidget * MainWindow::createLineList()
{
	lineList = new LineList(data);
	return lineList->getWidget();
}

// ------------------------------------------------------------------------
QWidget * MainWindow::createLayout()
{
	QWidget * leftWidget = new QWidget(this);
	QVBoxLayout * leftLayout = new QVBoxLayout();

	// on the left we have the image list and the buttons
	QWidget * viewWidget = createImageViewer();

	leftLayout->addWidget(createImageList());
	leftLayout->addWidget(createImageTypeControls());
	leftLayout->addWidget(createLineList());
	leftLayout->addWidget(createImageControls());
	leftLayout->addWidget(createButtonGroup());
	leftWidget->setLayout(leftLayout);
	leftWidget->setMaximumWidth(400);


	QWidget * rightWidget = new QWidget(this);
	QVBoxLayout * rightLayout = new QVBoxLayout;
	rightLayout->addWidget(viewWidget);
	rightLayout->addWidget(createInfoPane());
	rightWidget->setLayout(rightLayout);


	QHBoxLayout * mainLayout = new QHBoxLayout;
	mainLayout->addWidget(leftWidget);
	mainLayout->addWidget(rightWidget);


	QWidget * widget = new QWidget(this);
	widget->setLayout(mainLayout);

	return widget;

}

// ------------------------------------------------------------------------
QWidget * MainWindow::createImageTypeControls()
{
	QWidget * widget = new QWidget(this);
	QHBoxLayout * layout = new QHBoxLayout;

	imageTypeCombo = new QComboBox(this);
	imageTypeCombo->addItems(DataInstance::getTypes());
	layout->addWidget(imageTypeCombo);

	imageTypeButton = new QPushButton("Apply", this);
	layout->addWidget(imageTypeButton);

	widget->setLayout(layout);

	return widget;


}


// ------------------------------------------------------------------------
QWidget * MainWindow::createButtonGroup()
{
	QGroupBox * buttonBox = new QGroupBox(this);
	addLineButton = new QPushButton("Add Line", this);
	propagateLineButton = new QPushButton("Propagate Line", this);
	lockLineButton = new QPushButton("Lock Line", this);
	removeLineButton = new QPushButton("Remove Line", this);
	
	lineTypeCombo = new QComboBox(this);
	lineTypeCombo->addItem("Mitral Valve");
	lineTypeCombo->addItem("Tricuspid Valve");
	lineTypeCombo->addItem("Aortic Valve");
	lineTypeCombo->addItem("Pulmanary Valve");


	mvButton = new QPushButton("Mitral Valve", this);
	tpButton = new QPushButton("Tricuspid Valve", this);
	avButton = new QPushButton("Auortic Valve", this);
	mvButton->setCheckable(true);
	tpButton->setCheckable(true);
	avButton->setCheckable(true);

	mvButton->setChecked(true);

	/*
	valveGroup = new QButtonGroup(this);
	
	valveGroup->addButton(mvButton,0);
	valveGroup->addButton(tpButton,1);
	valveGroup->addButton(avButton,2);
	valveGroup->setExclusive(true);
	*/


	QVBoxLayout * layout = new QVBoxLayout;
	layout->addWidget(addLineButton);
	layout->addWidget(lineTypeCombo);
	layout->addWidget(propagateLineButton);
	layout->addWidget(lockLineButton);
	layout->addWidget(removeLineButton);

	//QFrame * line = new QFrame(this);
	//line->setFrameShape(QFrame::HLine);
	//line->setFrameShadow(QFrame::Sunken);
	//layout->addWidget(line);


	//layout->addWidget(mvButton);
	//layout->addWidget(tpButton);
	//layout->addWidget(avButton);
	
	buttonBox->setLayout(layout);
	return buttonBox;
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
	controls = new KeyControls;
	imageViewer->getWidget()->installEventFilter(controls);
	return imageViewer->getWidget();
}






