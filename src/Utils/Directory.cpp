#include "Directory.h"

#include <QDir>
#include <QString>
#include <QFileInfo>


namespace utils
{
// ------------------------------------------------------------------------
Directory::FilenamesType Directory::GetDirectories(const DirectoryNameType &dir)
{
	Directory::Pointer directory = Directory::New();
	directory->SetDirectory(dir);

	return directory->GetDirectories();
}

// ------------------------------------------------------------------------
Directory::FilenameType Directory::GetPath(const std::string &base, const std::string &ext)
{
	QDir dir(QString::fromStdString(base));
	return dir.absoluteFilePath(QString::fromStdString(ext)).toStdString();

}


// ------------------------------------------------------------------------
Directory::FilenamesType Directory::GetDirectories() const
{
	FilenamesType output;
	if(m_Directory.empty())
	{
		std::cout << "Directory::GetOutput - Directory Name Not Set" << std::endl;
		return output;
	}


	QDir dir(QString::fromStdString(m_Directory));

	QStringList qDirnames;
	qDirnames = dir.entryList(QDir::AllDirs | QDir::NoDotAndDotDot);

	for(int i = 0; i < qDirnames.size(); i++)
	{
		output.push_back(dir.absoluteFilePath(qDirnames[i]).toStdString());
	}

	return output;

}


// ------------------------------------------------------------------------
Directory::Directory()
{
	m_Directory = "";
	m_Extension = "";
}

// ------------------------------------------------------------------------
Directory::FilenamesType Directory::GetOutput() const
{
	FilenamesType output;
	if(m_Directory.empty())
	{
		std::cout << "Directory::GetOutput - Directory Name Not Set" << std::endl;
		return output;
	}


	QDir dir(QString::fromStdString(m_Directory));

	QStringList qFilenames;
	if(m_Extension.empty())
	{
		qFilenames = dir.entryList();
	}
	else
	{
		QStringList filters;
		filters << "*" + QString::fromStdString(m_Extension);
		qFilenames = dir.entryList(filters);
	}

	for(int i = 0; i < qFilenames.size(); i++)
	{
		output.push_back(dir.absoluteFilePath(qFilenames[i]).toStdString());
	}


	return output;
}


// ------------------------------------------------------------------------
Directory::FilenamesType Directory::GetFiles(const DirectoryNameType &dir, 
		const FilenameType & extension)
{
	Directory::Pointer directory = Directory::New();
	directory->SetDirectory(dir);
	directory->SetExtension(extension);
	return directory->GetOutput();
}



}
