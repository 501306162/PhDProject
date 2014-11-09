#include "main_window.h"

#include <QtCore>
#include "io.h"

// ------------------------------------------------------------------------
void MainWindow::saveActionPressed()
{
	if(okToSave())
	{
		if(!data->hasSaveName())
			return saveAsActionPressed();

		std::string filename = data->getSaveName();
		saveData(filename);
	}
}

// ------------------------------------------------------------------------
void MainWindow::saveData(const std::string &filename)
{
	bool ok = IO::Save(filename, data);
	if(!ok)
	{
		std::cout << "File couldn't be saved" << std::endl;
	}
}

// ------------------------------------------------------------------------
void MainWindow::saveAsActionPressed()
{
	if(okToSave())
	{
		QString filename = getSaveName();
		if(filename.isEmpty())
			return;

		saveData(filename.toStdString());
	}
}

// ------------------------------------------------------------------------
void MainWindow::loadData(const std::string &filename)
{
	DataContainer * data = IO::Load(filename);
	if(data != NULL)
		initialiseFromData(data);
}


// ------------------------------------------------------------------------
QString MainWindow::getSaveName()
{
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
			"/Users/oliverferoze/Uni/Data", tr("Json Files (*.json *.txt)"));

	return fileName;

}


// ------------------------------------------------------------------------
void MainWindow::imageTypePressed()
{
	QString type = imageTypeCombo->currentText();
	unsigned int selectedIndex = imageList->selectedIndex();
	data->setImageType(selectedIndex, type.toStdString());
	imageList->showList();
	imageList->setSelection(selectedIndex);

}

// ------------------------------------------------------------------------
bool MainWindow::okToSave()
{
	
	return true;
	//if(!data->linesAreLocked())
	//{
		//int ret = QMessageBox::question(this, tr("Hi"),
				//tr("There are unlocked lines in this session, These will be discarded"),
				//QMessageBox::Ok | QMessageBox::Cancel);

		//if(ret == QMessageBox::Ok)
			//return true;
		//else
			//return false;
	//}
	//else
	//{
		//return true;
	//}
}


// ------------------------------------------------------------------------
void MainWindow::tSliderRight()
{
	int val = tSlider->value();
	if(val == tSlider->maximum()) 
		tSlider->setValue(tSlider->minimum());
	else
		tSlider->setValue(val+1);
}

// ------------------------------------------------------------------------
void MainWindow::tSliderLeft()
{
	int val = tSlider->value();
	if(val == tSlider->minimum()) 
		tSlider->setValue(tSlider->maximum());
	else
		tSlider->setValue(val-1);
}


// ------------------------------------------------------------------------
void MainWindow::lineChanged()
{
	if(!lineList->lineSelected())
	{
		imageViewer->showAllLines();
	}
	else
	{
		Line::Type type = lineList->getSelectedLineType();
		lineList->setCurrentIndex();
		imageViewer->setLineToDisplay(type);
	}

	imageViewer->updateImage();

}

// ------------------------------------------------------------------------
void MainWindow::propagateLinePressed()
{
	if(!lineList->lineSelected())
		return; 

	data->propagateLine(lineList->getSelectedLineType());
}

// ------------------------------------------------------------------------
void MainWindow::removeLinePressed()
{
	if(!lineList->lineSelected())
		return; 

	// get the line 
	data->removeLine(lineList->getSelectedLineType());
	lineList->updateLines();
	imageViewer->updateImage();
}


// ------------------------------------------------------------------------
void MainWindow::addLinePressed()
{
	Line * line = Line::NewLine(data->getVTKImage(), 
			Line::getTypeEnum(lineTypeCombo->currentIndex()));
	data->addLine(line);
	data->propagateLine(line->getType());
	imageViewer->showAllLines();
	imageViewer->updateImage();
	lineList->updateLines();
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
	
	return NULL;
}

// ------------------------------------------------------------------------
void MainWindow::imageSelectionChanged()
{
	if(data == NULL)
	{
		std::cout << "Test" << std::endl;
		return;
	}

	data->setViewedImage(imageList->selectedIndex());
	updateSliders();
	data->setViewedTimeStep(tSlider->value());
	data->setViewedSlice(zSlider->value());
	updateAll(true);
}


// ------------------------------------------------------------------------
void MainWindow::sliderSelectionChanged()
{
	data->setViewedTimeStep(tSlider->value());
	data->setViewedSlice(zSlider->value());
	updateAll();
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
void MainWindow::setInfoPane()
{
	if(data == NULL) return;

	QString value = QString::fromStdString(data->getFolderName());
	folderValue->setText(value);

}


// ------------------------------------------------------------------------
void MainWindow::lockLine()
{
	if(!lineList->lineSelected())
		return; 

	this->data->lockLine(this->lineList->getSelectedLineType());
	updateAll();
}

// ------------------------------------------------------------------------
void MainWindow::updateAll(bool reset)
{
	lineChanged();
	imageViewer->updateImage(reset);
	lineList->updateLines();
}




// ------------------------------------------------------------------------
void MainWindow::closeEvent(QCloseEvent * event)
{
	delete data;
}

// ------------------------------------------------------------------------
void MainWindow::newActionPressed()
{
	// check that it is ok to overwrite
	if(!okToOverwrite())
		return;

	// get a filename
	QString fileName = getNewFolder();
	if(fileName.isEmpty())
		return;

	// delete old data 
	delete data;

	DataContainer * newData = new DataContainer;
	newData->LoadData(fileName.toStdString());
	initialiseFromData(newData);
}


// ------------------------------------------------------------------------
void MainWindow::loadActionPressed()
{
	// check that it is ok to overwrite
	if(!okToOverwrite())
		return;

	QString fileName = getLoadFile();

	loadData(fileName.toStdString());
}


// ------------------------------------------------------------------------
QString MainWindow::getLoadFile()
{
	QString fileName = QFileDialog::getOpenFileName(this, 
			tr("Select the file that you want to load"),
			"/Users/oliverferoze/Uni/Data/LineOutput",
			tr("Json Files (*.json *.txt)"));
	
	return fileName;
}


// ------------------------------------------------------------------------
QString MainWindow::getNewFolder()
{
	QString fileName = QFileDialog::getExistingDirectory(this, 
			tr("Set a folder to read from"),
			"/Users/oliverferoze/Uni/Data",
			QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
	return fileName;
}

// ------------------------------------------------------------------------
bool MainWindow::okToOverwrite()
{
	if(data == NULL)
		return true;

	int ret = QMessageBox::question(this, tr("Ok To Save"),
			tr("You are in the middle of a session\n do you want to save?"),
			QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);


	// check the return value from the box
	switch(ret)
	{
		case QMessageBox::Save :
			saveActionPressed();
			return true;
		case QMessageBox::Cancel :
			return false;
		case QMessageBox::Discard :
			return true;
		default:
			break;
	}

	return true;
}




// ------------------------------------------------------------------------
void MainWindow::initialiseFromData(DataContainer * data)
{
	this->data = data;
		
	QString folderName = QString::fromStdString(data->getFolderName());
	setInfoPane();
	
	imageViewer->setData(this->data);
	lineList->setData(this->data);
	imageList->setData(this->data);
	
	imageViewer->updateImage(true);

}

// ------------------------------------------------------------------------
MainWindow::~MainWindow()
{
	delete this->imageList;
	delete this->imageViewer;
}
