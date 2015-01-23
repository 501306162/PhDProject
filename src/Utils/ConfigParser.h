#ifndef CONFIG_PARSER_H
#define CONFIG_PARSER_H

#include <QVariant>
#include <iostream>

namespace utils
{
class ConfigParser
{
public:
	ConfigParser(const std::string &configFilename);
	virtual ~ConfigParser() {};

	int intValue(const std::string &name, int bk=0);
	std::string strValue(const std::string &name, std::string bk="");
	double doubleValue(const std::string &name, double bk=0.0);


private:
	bool hasValue(const std::string &value);
	QVariant getValue(const std::string &value);

	QVariantMap m_Configs;

};



}



#endif
