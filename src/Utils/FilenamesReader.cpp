#include "FilenamesReader.h"

#include <QFileInfo>
#include <QTextStream>

namespace utils
{
// ------------------------------------------------------------------------
FilenamesReader::FilenamesReader() 
{
	this->m_Filename = "";
	this->m_CheckFileContents = false;
}

// ------------------------------------------------------------------------
void FilenamesReader::SetInput(const FilenameType &filename)
{
	this->SetInputFilename(filename);
}

// ------------------------------------------------------------------------
void FilenamesReader::SetInputFilename(const FilenameType &filename)
{
	this->m_Filename = filename;
	this->Modified();
}

// ------------------------------------------------------------------------
void FilenamesReader::SetCheckFileContents(bool val)
{
	this->m_CheckFileContents = val;
}

// ------------------------------------------------------------------------
void FilenamesReader::Read() throw (itk::ExceptionObject)
{
	this->m_Ouptut.clear();


	// check that the input filename is set
	if(this->m_Filename.empty())
	{
		itkExceptionMacro(<< "Input Filename Not Set");
	}


	// check that the input filename is valid
	if(!this->IsValidFile(this->m_Filename))
	{
		itkExceptionMacro(<< "Input Filename Is Not A Valid File");
	}
	

	// now we read the file
	QFile file(QString::fromStdString(this->m_Filename));
	if(!file.open(QIODevice::ReadOnly))
	{
		itkExceptionMacro(<< "Error Loading File From Disk");
	}
	
	// read each line of the file one by one
	QTextStream ts(&file);
	while(!ts.atEnd())
	{
		std::string line = ts.readLine().toStdString();

		// check to see if the line is a valid file
		if(m_CheckFileContents 
				&& !this->IsValidFile(line))
		{
			itkExceptionMacro(<< "Error: " << line << " Is Not A Valid File");			
		}

		this->m_Ouptut.push_back(line);
	}

}


// ------------------------------------------------------------------------
FilenamesReader::FilenamesType FilenamesReader::GetOutput() const
{
	return m_Ouptut;
}


// ------------------------------------------------------------------------
bool FilenamesReader::IsValidFile(const std::string &filename)
{
	QFileInfo info(QString::fromStdString(filename));
	if(info.exists() && info.isFile()) return true;
	return false;
}


// ------------------------------------------------------------------------
FilenamesReader::FilenamesType FilenamesReader::Read(const std::string &input, bool checkFiles)
{
	FilenamesReader::Pointer reader = FilenamesReader::New();
	reader->SetCheckFileContents(checkFiles);
	reader->SetInputFilename(input);
	
	FilenamesType output;

	try
	{
		reader->Read();
	}
	catch(itk::ExceptionObject &e)
	{
		std::cerr << "Error in Filenames Reader: " << e << std::endl;
		return output;
	}

	output = reader->GetOutput();
	return output;
}



} /* utils */ 
