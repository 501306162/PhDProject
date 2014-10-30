#include "key_controls.h"

#include <QtGui>

// ------------------------------------------------------------------------
KeyControls::KeyControls()
{

}

// ------------------------------------------------------------------------
bool KeyControls::eventFilter(QObject * obj, QEvent * event)
{
	if(event->type() == QEvent::KeyPress)
	{
		QKeyEvent * keyEvent = static_cast<QKeyEvent*>(event);
		qDebug("Ate key press %d", keyEvent->key());

		switch(keyEvent->key())
		{
			case 16777236:
				emit rightPressed();
				break;
			case 16777234:
				emit leftPressed();
				break;
			case 16777235:
				emit upPressed();
				break;
			case 16777237:
				emit downPressed();
				break;
			case 76:
				emit lockPressed();
				break;
			case 80:
				emit propagatePressed();
				break;
			default:
				break;

		}

		return true;
	}
	return QObject::eventFilter(obj,event);

}
