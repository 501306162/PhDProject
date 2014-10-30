#include "line_list.h"

#include <QtGui>


// ------------------------------------------------------------------------
LineList::LineList(DataContainer * data)
{
	setData(data);
	this->table = new QTableWidget;
	setUpTable();
}

// ------------------------------------------------------------------------
void LineList::setData(DataContainer * data)
{
	this->data = data;
}

// ------------------------------------------------------------------------
void LineList::setUpTable()
{
	QStringList headers;
	headers << "Type";
	headers << "Locked";

	this->table->setColumnCount(2);
	this->table->setHorizontalHeaderLabels(headers);
	this->table->setShowGrid(false);


	QHeaderView * view = this->table->horizontalHeader();
	view->setResizeMode(QHeaderView::Stretch);

	QHeaderView * rowView = this->table->verticalHeader();
	rowView->setDefaultSectionSize(rowView->fontMetrics().height()+6);


	this->table->setSelectionBehavior(QAbstractItemView::SelectRows);
	this->table->setSelectionMode(QAbstractItemView::SingleSelection);


}


// ------------------------------------------------------------------------
bool LineList::lineSelected()
{
	if(this->table->selectedItems().size() > 0)
		return true;
	else
		return false;
}


// ------------------------------------------------------------------------
unsigned int LineList::getSelectedLineIndex()
{
	if(!lineSelected())
		return -1;

	QTableWidgetItem * item = this->table->selectedItems().takeFirst();
	return item->row();
}

// ------------------------------------------------------------------------
void LineList::updateLines()
{
	// remove current entries
	removeEntrys();
	
	// set the new entries
	Line::Map lines = data->getLineData();
	table->setRowCount(lines.size());

	Line::Map::iterator it = lines.begin();
	int count =  0;
	while(it != lines.end())
	{
		Line * line = it->second;
		QString type = QString::fromStdString( Line::getTypeString(line->getType()) );

		QTableWidgetItem *item1 = new QTableWidgetItem(type);
		table->setItem(count,0,item1);

				

		QString lockedVal = "Unlocked";
		if(line->isLocked())
			lockedVal = "Locked";
	

		QTableWidgetItem *item2 = new QTableWidgetItem(lockedVal);
		table->setItem(count,1,item2);

		currentLines.push_back(line);

		++it; ++count;
	}

	if(currentIndex >= 0 && currentIndex < (int) currentLines.size())
	{
		QTableWidgetItem * item = table->item(currentIndex, 0);
		table->setItemSelected(item,true);
	}

}


// ------------------------------------------------------------------------
Line::Type LineList::getSelectedLineType()
{
	unsigned int index = getSelectedLineIndex();
	
	return currentLines[index]->getType();	
}


// ------------------------------------------------------------------------
void LineList::removeEntrys()
{
	while(table->rowCount() > 0)
	{
		table->removeRow(0);
	}

	// clear the line list
	currentLines.clear();
}
