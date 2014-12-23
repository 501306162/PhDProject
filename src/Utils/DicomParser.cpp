#include "DicomParser.h"


// ------------------------------------------------------------------------
DicomParser::DicomParser()
{
	m_InputFolder = "";

}


// ------------------------------------------------------------------------
bool DicomParser::Check()
{
	if(m_InputFolder.empty())
	{
		m_Errors.push_back("The input folder has not been set");
		return false;
	}


	return true;
}


// ------------------------------------------------------------------------
bool DicomParser::Parse()
{
	if(!Check())
		return false;

	


	return true;
}
