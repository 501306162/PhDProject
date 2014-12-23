#ifndef CARDIAC_DICOM_EXTRACTOR_H
#define CARDIAC_DICOM_EXTRACTOR_H

#include <itkImage.h>
#include <itkObject.h>
#include <itkObjectFactory.h>


#include "SeriesExtractor.h"

namespace utils
{
class CardiacDicomExtractor : public itk::Object 
{
public:
	typedef CardiacDicomExtractor Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;



	

	itkTypeMacro(CardiacDicomExtractor, Object);
	itkNewMacro(Self);

	itkSetMacro(FolderName, std::string);
	itkGetMacro(FolderName, std::string);

	typedef std::string FlagType;
	typedef std::vector<FlagType> FlagListType;
	typedef std::map<FlagType, FlagListType> FlagContainerType;

	void Extract(); 

protected:
	CardiacDicomExtractor() {}
	virtual ~CardiacDicomExtractor() {}

private:

	void BuildFlagLists();
	void RunExtractor();

	CardiacDicomExtractor(const Self&);
	void operator=(const Self&);

	std::string m_FolderName;

	DicomSeriesList m_SeriesList;

	FlagContainerType m_Flags;



};





}


#endif
