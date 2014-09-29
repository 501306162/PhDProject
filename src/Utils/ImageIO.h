#ifndef IMAGE_IO_H
#define IMAGE_IO_H

#include <iostream>

#include <itkImage.h>
#include <itkObject.h>

#include <itkImageIOBase.h>

namespace utils
{
template <typename TImageType>
class ImageIO : public itk::Object
{
public:
	typedef ImageIO<TImageType> Self;
	typedef itk::Object    		Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;


	itkNewMacro(Self);
	itkTypeMacro(ImageIO, itk::Object);

	/** Image Typedefs */
	typedef TImageType ImageType;
	typedef typename ImageType::Pointer ImagePointer;
	typedef itk::ImageIOBase ImageIOBase;
	typedef ImageIOBase::Pointer ImageIOPointer;

	void SetInput(const ImagePointer &image);
	ImagePointer GetOutput() const;

	void SetFilename(const std::string &filename);
	void SetImageIO(const ImageIOPointer &imageIO);
	void Read() throw(itk::ExceptionObject);
	void Write() throw(itk::ExceptionObject);

	static ImagePointer Read(const std::string &filename, const ImageIOPointer &imageIO=NULL);
	static void Write(const std::string &filename, const ImagePointer &images,  const ImageIOPointer &imageIO=NULL);

protected:
	ImageIO();
	virtual ~ImageIO();

private:
	ImageIOPointer m_ImageIO;
	ImagePointer m_Image;
	std::string m_Filename;

};

} /* utils */ 


#include "ImageIO.hpp"


#endif
