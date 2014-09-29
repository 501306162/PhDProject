#include "PatientData.h"

#include <QFile>
#include <QTextStream>
#include <QStringList>

namespace utils
{
// ------------------------------------------------------------------------
PatientData::PatientData()
{

}

// ------------------------------------------------------------------------
PatientData * PatientData::Load(const std::string &filename)
{
	// try reading the file
	QFile file(QString::fromStdString(filename));
	if(!file.open(QIODevice::ReadOnly))
	{
		std::cout << filename << " is not a valid file"  << std::endl;
		return NULL;
	}
	
	// create a new patient data
	PatientData::Pointer data = PatientData::New();

	QTextStream in(&file);

	while(!in.atEnd())
	{
		QString line = in.readLine();
			
		// split the line
		QStringList parts = line.split(",");
		data->SetPatientId(parts[0].toInt());
	
		



	}

	file.close();




	return data;

}



} /*  utils */ 
