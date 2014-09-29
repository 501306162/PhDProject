#ifndef PATIENT_DATA_H
#define PATIENT_DATA_H

#include <itkObject.h>
#include <itkObjectFactory.h>

namespace utils
{
class PatientData : public itk::Object
{
public:
	enum PhaseType { ED, ES };

	typedef PatientData Self;
	typedef itk::Object Superlcass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(PatientData, itk::Object);
	itkNewMacro(Self);

	typedef std::string FilenameType;
	typedef std::pair<double, double> PointType;
	typedef std::pair<PointType, PointType> LandmarksType;


	typedef struct slice_info_
	{
		unsigned int zNum;
		unsigned int sliceNum;
		FilenameType imageFilename;
		FilenameType labelFilename;
		typedef std::vector<slice_info_> List;

	} SliceInfo;

	static PatientData * Load(const FilenameType &filename);


	itkGetMacro(PatientId, unsigned int);
	itkSetMacro(PatientId, unsigned int);
	itkSetMacro(NumberOfSlices, unsigned int);
	itkGetMacro(NumberOfSlices, unsigned int);

	void SetSlice(unsigned int id, SliceInfo &info);


protected:
	PatientData();
	~PatientData() {}


private:
	unsigned int m_PatientId;
	unsigned int m_NumberOfSlices;
	PhaseType m_Phase;

	SliceInfo::List m_SliceList;



	PatientData(const Self&);
	void operator=(const Self&);
};
} /*  utils */ 

#endif
