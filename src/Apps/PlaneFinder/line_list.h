#ifndef LINE_LIST_H
#define LINE_LIST_H

#include "containers.h"
#include <QTableWidget>


class LineList
{
public:
	LineList(DataContainer * data);

	QTableWidget * getWidget() { return this->table; }


	void setData(DataContainer * data);

	void updateLines();
	
	void removeEntrys();
	void setUpTable();
	unsigned int getSelectedLineIndex();
	bool lineSelected();

	Line::Type getSelectedLineType();

	void setCurrentIndex() { this->currentIndex = getSelectedLineIndex(); }



private:
	DataContainer * data;
	QTableWidget * table;

	std::vector<Line *> currentLines;
	int currentIndex;

};

#endif
