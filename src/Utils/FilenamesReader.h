#ifndef FILENAMES_READER_H
#define FILENAMES_READER_H

#include <itkObject.h>
#include <itkImage.h>

#include <QString>

namespace utils
{
class FilenamesReader : public itk::Object
{
public:
	/** Typical itk typedefs */
	typedef FilenamesReader 				Self;
	typedef itk::Object     				Superclass;
	typedef itk::SmartPointer<Self> 		Pointer;
	typedef itk::SmartPointer<const Self> 	ConstPointer;

	itkTypeMacro(FilenamesReader, itk::Object);
	itkNewMacro(Self);


	typedef std::string 					FilenameType;
	typedef std::vector<FilenameType> 		FilenamesType;

	void SetInputFilename(const FilenameType & filename);
	void SetInput(const FilenameType & filename);
	void Read() throw (itk::ExceptionObject);
	void SetCheckFileContents(bool value);
	FilenamesType GetOutput() const;


	static FilenamesType Read(const std::string &input, bool checkFiles=false);

protected:
	FilenamesReader();	
	virtual ~FilenamesReader() {}


	bool IsValidFile(const std::string &filename);

private:
	FilenamesReader(const Self&);
	void operator=(const Self&);


	FilenameType m_Filename;
	FilenamesType m_Ouptut;
	bool m_CheckFileContents;
};



} /* utils */ 

#endif
