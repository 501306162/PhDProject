#ifndef DICOM_PARSER_H
#define DICOM_PARSER_H

#include <iostream>
#include <vector>

class DicomParser
{
public:
	DicomParser();

	void SetInputFolder(const std::string &folder) { m_InputFolder = folder; }
	bool Parse();
	bool Check();



private:
	void Scan();

	std::string m_InputFolder;

	
	std::vector<std::string> m_Errors;

};

#endif
