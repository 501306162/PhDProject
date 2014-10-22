#include "line_list.h"

#include <QtGui>


// ------------------------------------------------------------------------
LineList::LineList(DataContainer * data)
{
	this->data = data;
	this->table = new QTableWidget;
	setUpTable();
}


// ------------------------------------------------------------------------
void LineList::setUpTable()
{
	QStringList headers;
	headers << "Num";
	headers << "Type";

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

	QTableWidgetItem * item = this->table->selectedItems().takeFirst();
	return item->row();
}

// ------------------------------------------------------------------------
void LineList::updateLines()
{
	// remove current entries
	removeEntrys();
	
	// set the new entries
	Line::List lines = data->getLineData();
	table->setRowCount(lines.size());


	for(unsigned int i = 0; i < lines.size(); i++)
	{
		Line * line = lines[i];
		QString num = QString::number(i+1);
		QString type = QString::fromStdString( Line::getTypeString(line->getType()) );

		QTableWidgetItem *item1 = new QTableWidgetItem(num);
		table->setItem(i,0,item1);
		QTableWidgetItem *item2 = new QTableWidgetItem(type);
		table->setItem(i,1,item2);
	}
}


// ------------------------------------------------------------------------
void LineList::removeEntrys()
{
	while(table->rowCount() > 0)
	{
		table->removeRow(0);
	}
}
