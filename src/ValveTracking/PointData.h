#ifndef POINT_DATA_H
#define POINT_DATA_H

#include <itkNormalVariateGenerator.h>
#include <itkObject.h>
#include <itkObjectFactory.h>
#include <MatrixCommon.h>
#include <itkPoint.h>
#include <ValveLine.h>
#include <itkCenteredAffineTransform.h>
#include <itkGaussianMembershipFunction.h>
#include <itkCovarianceSampleFilter.h>
#include <itkListSample.h>

namespace vt
{
class PointData : public itk::Object 
{
public:
	typedef PointData Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkNewMacro(Self);
	itkTypeMacro(PointData, Object);

	typedef utils::DoubleMatrixType MatrixType;
	typedef itk::Point<double, 3> PointType;
	typedef itk::Vector<double, 3> VectorType;
	typedef utils::IntMatrixType IntMatrixType;
	typedef std::pair<MatrixType, MatrixType> PointPairType;
	typedef std::vector<PointPairType> ValvePoints;
	typedef std::map<std::string, ValvePoints> PointMap;
	typedef itk::CenteredAffineTransform<double, 3> TransformType;

	typedef struct dist_
	{
		MatrixType rot;
		MatrixType scl;
		MatrixType mean;
	} Distribution;

	typedef std::pair<Distribution, Distribution> PointDistribution;
	typedef std::vector<PointDistribution> PointDistributionList;
	typedef std::map<std::string, PointDistributionList> DistributionData;
	typedef itk::Statistics::NormalVariateGenerator RNGType;


	typedef itk::Statistics::ListSample<VectorType> SampleType;
	typedef itk::Statistics::CovarianceSampleFilter<SampleType> CovFilterType;
	typedef itk::Statistics::GaussianMembershipFunction<VectorType> MembershipType;

	void TransformToOrigin(const PointType &p1, const PointType &p2, TransformType::Pointer &transform);
	static PointData::Pointer Load(const std::string &folderName, const unsigned int exclude, bool flat=false);
	void LoadData(const std::string &folderName, const unsigned int exclude);
	void LoadDataFlat(const std::string &folderName, const unsigned int exclude);

	void GetMembers(const std::string type, const unsigned int timeStep,
			const PointType &p1, const PointType &p2, MembershipType::Pointer &p1Member, 
			MembershipType::Pointer &p2Member);


	void GetUpdate(const std::string &type, const unsigned int timeStep,
			PointType &p1, PointType &p2);

	void GetSample(const Distribution & dist, PointType &sample);



protected:
	PointData() {
		m_Generator = RNGType::New();
		m_Generator->Initialize((int) time(NULL));
	}
	virtual ~PointData() {}

private:

	void LoadFolder(const std::string &folderName, ValvePoints &points);
	void ComputeDistributions(const ValvePoints &points, PointDistributionList &dists);
	void LoadFolderFlat(const std::string &folderName, ValvePoints &points);
	void GetMemberFunction(const MatrixType &points, const TransformType::Pointer &transform,
			MembershipType::Pointer &member);
	void GetDist(const MatrixType &points, Distribution &dist);
	std::string m_Exclude;
	
	PointData(const Self&);
	void operator=(const Self&);

	PointMap m_Data;
	DistributionData m_Distributions;
	RNGType::Pointer m_Generator;

};


}

#endif
