#ifndef LINE_LIST_H
#define LINE_LIST_H

#include "containers.h"
#include <QTableWidget>


class LineList
{
public:
	LineList(DataContainer * data);

	QTableWidget * getWidget() { return this->table; }

	void updateLines();
	
	void removeEntrys();
	void setUpTable();
	unsigned int getSelectedLineIndex();
	bool lineSelected();


private:
	DataContainer * data;
	QTableWidget * table;

};

#endif
