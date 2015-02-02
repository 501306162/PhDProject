#ifndef CMR_FILE_NAME_EXTRACTOR_H
#define CMR_FILE_NAME_EXTRACTOR_H

#include <itkImage.h>

#include <QString>

namespace vt
{
class CMRFileExtractor : public itk::Object
{
public:
	typedef CMRFileExtractor Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(CMRFileExtractor, Object);
	itkNewMacro(Self);

	typedef itk::Image<unsigned short, 4> ImageSeriesType;
	typedef itk::Image<unsigned short, 3> ImageType;

	void SetFolderName(const std::string &foldername) { m_FolderName = foldername; }
	void Extract();

	void LoadImage(const QString &filename, ImageSeriesType::Pointer &image);
	
	bool Matches(const QString &filename, const QStringList &filters);

	ImageType::Pointer Get2CImage(const unsigned int index);
	ImageType::Pointer Get3CImage(const unsigned int index);
	ImageType::Pointer Get4CImage(const unsigned int index);
	ImageType::Pointer GetR2CImage(const unsigned int index);
	ImageType::Pointer GetR3CImage(const unsigned int index);
	ImageType::Pointer GetStackImage(const unsigned int index);

	bool Has2CImage() { return Exist(m_2CImage); }
	bool Has3CImage() { return Exist(m_3CImage); }
	bool Has4CImage() { return Exist(m_4CImage); }
	bool HasR2CImage() { return Exist(m_R2CImage); }
	bool HasStackImage() { return Exist(m_Stack); }

	itkSetMacro(Flip, bool);

	void FlipImage(const ImageType::Pointer &input, ImageType::Pointer &output);

	ImageType::Pointer ExtractImageTimeStep(const ImageSeriesType::Pointer &image, 
			const unsigned int index);

	void SetDebug(bool val) { m_Debug = val; }

protected:
	CMRFileExtractor() : m_Debug(false), m_Flip(false) {};
	virtual ~CMRFileExtractor() {}

private:
	CMRFileExtractor(const Self&);
	void operator=(const Self&);

	bool Exist(ImageSeriesType::Pointer &image);

	std::string m_FolderName;

	bool m_Flip;

	ImageSeriesType::Pointer m_Stack;
	ImageSeriesType::Pointer m_2CImage;
	ImageSeriesType::Pointer m_3CImage;
	ImageSeriesType::Pointer m_4CImage;
	ImageSeriesType::Pointer m_R2CImage;
	ImageSeriesType::Pointer m_R3CImage;


	bool m_Debug;
};


}

#endif
