#ifndef IMAGE_IO_HPP
#define IMAGE_IO_HPP

#include "ImageIO.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

namespace utils
{

// ------------------------------------------------------------------------
template<typename TInputImage>
ImageIO<TInputImage>::ImageIO()
{
	m_ImageIO = NULL;
	m_Image   = NULL;	
}

// ------------------------------------------------------------------------
template<typename TInputImage>
ImageIO<TInputImage>::~ImageIO()
{
	
}

// ------------------------------------------------------------------------
template<typename TImageType>
void ImageIO<TImageType>::SetInput(const ImagePointer &image)
{
	this->m_Image = image;	
}

// ------------------------------------------------------------------------
template<typename TInputImage>
void ImageIO<TInputImage>::SetImageIO(const itk::ImageIOBase::Pointer &imageIO)
{
	this->m_ImageIO = imageIO;
}

// ------------------------------------------------------------------------
template<typename TInputImage>
typename ImageIO<TInputImage>::ImagePointer ImageIO<TInputImage>::GetOutput() const
{
	return m_Image;
}

// ------------------------------------------------------------------------
template<typename TInputImage>
void ImageIO<TInputImage>::SetFilename(const std::string &filename)
{
	m_Filename = filename;	
}

// ------------------------------------------------------------------------
template<typename TInputImage>
void ImageIO<TInputImage>::Read() throw (itk::ExceptionObject)
{
	if(m_Filename.empty())
	{
		itkExceptionMacro(<< "Filename Is Not Set");
	}	

	
	// create the reader object
	typedef itk::ImageFileReader<ImageType> ReaderType;
	typename ReaderType::Pointer reader = ReaderType::New();

	// set image io?
	if(m_ImageIO) reader->SetImageIO(m_ImageIO);
	reader->SetFileName(m_Filename);

	try
	{
		reader->Update();
	}
	catch(itk::ExceptionObject &e)
	{
		itkExceptionMacro(<< e);
	}


	m_Image = reader->GetOutput();

}

// ------------------------------------------------------------------------
template<typename TInputImage>
void ImageIO<TInputImage>::Write() throw (itk::ExceptionObject)
{
	if(m_Filename.empty())
	{
		itkExceptionMacro(<< "Filename Is Not Set");
	}	

	if(!m_Image)
	{
		itkExceptionMacro(<< "Image Is Not Set");
	}
	//TODO Implement the write function

	typedef itk::ImageFileWriter<ImageType> WriterType;
	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(m_Filename);
	writer->SetInput(m_Image);

	try
	{
		writer->Write();
	}
	catch(itk::ExceptionObject e)
	{
		itkExceptionMacro(<< e);
	}
}

// ------------------------------------------------------------------------
template<typename TInputImage>
typename ImageIO<TInputImage>::ImagePointer ImageIO<TInputImage>::Read(const std::string &filename,
		const ImageIOPointer &imageIO)
{
	Self::Pointer io = Self::New();
	io->SetFilename(filename);
	if(imageIO) io->SetImageIO(imageIO);
	
	try
	{
		io->Read();
	}
	catch(itk::ExceptionObject &e)
	{
		std::cerr << "Error In Static Load Image Function" << e << std::endl;
		exit(1);
	}

	return io->GetOutput();
}


// ------------------------------------------------------------------------
template<typename TInputImage>
void ImageIO<TInputImage>::Write(const std::string &filename,
		const ImagePointer &image,
		const ImageIOPointer &imageIO)
{
	Self::Pointer io = Self::New();
	io->SetFilename(filename);
	io->SetInput(image);
	if(imageIO) io->SetImageIO(imageIO);
	
	try
	{
		io->Write();
	}
	catch(itk::ExceptionObject &e)
	{
		std::cerr << "Error In Static Load Image Function" << e << std::endl;
		exit(1);
	}
}


} /* utils */ 

#endif
