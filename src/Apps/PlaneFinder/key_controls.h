#ifndef KEY_CONTROLS_H
#define KEY_CONTROLS_H

#include <iostream>

#include <QObject>

class KeyControls : public QObject
{
	Q_OBJECT
public:
	KeyControls();

signals:
	void leftPressed();
	void rightPressed();
	void upPressed();
	void downPressed();
	void lockPressed();
	void propagatePressed();

protected:
	bool eventFilter(QObject * obj, QEvent * event);


};


#endif
