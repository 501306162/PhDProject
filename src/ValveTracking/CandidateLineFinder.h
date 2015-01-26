#ifndef CANDIDATE_LINE_FINDER_H
#define CANDIDATE_LINE_FINDER_H

#include "ValvePointCandidateFinder.h"
#include "LengthData.h"

#include <CommonDefinitions.h>

namespace vt 
{
class CandidateLineFinder : public itk::Object 
{
public:
	typedef CandidateLineFinder Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;	

	itkNewMacro(Self);

	typedef ValvePointCandidateFinder::PointListType PointListType;
	typedef ValvePointCandidateFinder::PointPairType PointPairType;
	typedef ValvePointCandidateFinder::PointType PointType;

	typedef std::pair<PointType, PointType> OutputPointPair;
	typedef std::pair<double, OutputPointPair> OutputPointCostPair;
	typedef std::vector<OutputPointCostPair> OutputListType;

	typedef Eigen::Vector3d VectorType;

	void SetP1List(const PointListType &list) { m_P1List = list; }
	void SetP2List(const PointListType &list) { m_P2List = list; }
	void SetLengthData(const LengthData::Pointer &lengths) { m_Lengths = lengths; }
	void SetOriginalLine(const PointType &p1, const PointType &p2) { m_OriginalP1 = p1; m_OriginalP2 = p2; }
	void SetType(const std::string &type) { m_Type = type; }
	void SetTimeStep(const unsigned int time) { m_TimeStep = time; }

	void Compute();
	OutputPointPair GetOutput(unsigned int index=0) { return m_OutputList[index].second; }
	void GetOutput(std::vector<OutputPointPair> &output, const unsigned int num=10);
	OutputPointPair GetWeightedOutput();


protected:
	CandidateLineFinder() {}
	virtual ~CandidateLineFinder() {}
	

private:
	CandidateLineFinder(const Self&);
	void operator=(const Self&);

	double LineCost(const PointPairType &p1, const PointPairType &p2);

	std::string m_Type;
	unsigned m_TimeStep;
	PointListType m_P1List;
	PointListType m_P2List;
	LengthData::Pointer m_Lengths;
	PointType m_OriginalP1;
	PointType m_OriginalP2;

	OutputListType m_OutputList;

};



}

#endif
