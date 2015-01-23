#ifndef LENGTH_DATA_H
#define LENGTH_DATA_H

#include <itkObject.h>
#include <itkObjectFactory.h>
#include <MatrixCommon.h>
#include <itkListSample.h>
#include <itkGaussianMembershipFunction.h>
#include <itkCovarianceSampleFilter.h>

namespace vt
{
class LengthData : public itk::Object 
{
public:
	typedef LengthData Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(LengthData, Object);
	itkNewMacro(Self);

	static LengthData::Pointer Load(const std::string &folderName);
	void LoadData(const std::string &folderName);

	typedef itk::Vector<double, 1> MeasurementType;
	typedef itk::Statistics::ListSample<MeasurementType> SampleType;
	typedef itk::Statistics::GaussianMembershipFunction<MeasurementType> MembershipType;
	typedef std::vector<MembershipType::Pointer> MembershipListType;
	typedef std::map<std::string, MembershipListType> MembershipMapType;
	typedef utils::DoubleMatrixType MatrixType;
	typedef utils::IntMatrixType IntMatrixType;

	typedef std::pair<MatrixType, IntMatrixType> MatrixPair;
	typedef std::map<std::string, MatrixPair> MatrixMap;


	double Predict(const std::string &type, const unsigned int time, const double length);

	void Prepare(const unsigned int exclude);

protected:
	LengthData() {}
	virtual ~LengthData() {}

private:
	LengthData(const Self&);
	void operator=(const Self&);

	MatrixMap m_Data;
	MembershipMapType m_Functions;

};


}

#endif
