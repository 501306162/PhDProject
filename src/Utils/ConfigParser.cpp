#include "ConfigParser.h"

#include <QFile>
#include <qjson/parser.h>

namespace utils
{
// ------------------------------------------------------------------------
ConfigParser::ConfigParser(const std::string &configFilename)
{
	QFile file(configFilename.c_str());
	if(!file.open(QIODevice::Text | QIODevice::ReadOnly))
	{
		std::cout << "Couldn't open the config file: " << configFilename << std::endl;
		exit(1);
	}

	QJson::Parser parser;
	bool ok;
	m_Configs = parser.parse(&file, &ok).toMap();

	if(!ok)
	{
		std::cout << "Couldn't parse the config file: " << configFilename << std::endl;
		exit(1);
	}
}

// ------------------------------------------------------------------------
int ConfigParser::intValue(const std::string &name, int bk)
{
	if(hasValue(name)) return getValue(name).toInt();
	else return bk;
}

// ------------------------------------------------------------------------
std::string ConfigParser::strValue(const std::string &name, std::string bk)
{
	if(hasValue(name)) return getValue(name).toString().toStdString();
	else return bk;
}

// ------------------------------------------------------------------------
double ConfigParser::doubleValue(const std::string &name, double bk)
{
	if(hasValue(name)) return getValue(name).toDouble();
	else return bk;
}

// ------------------------------------------------------------------------
bool ConfigParser::hasValue(const std::string &value)
{
	if(m_Configs.contains(value.c_str()))
		return true;
	else
		return false;
}

QVariant ConfigParser::getValue(const std::string &type)
{
	return m_Configs[type.c_str()];
}


}
