#ifndef VALVE_POINT_CANDIDATE_FINDER_H
#define VALVE_POINT_CANDIDATE_FINDER_H

#include <itkObject.h>
#include <itkObjectFactory.h>
#include <itkImage.h>


namespace vt
{
class ValvePointCandidateFinder : public itk::Object 
{
public: 
	typedef ValvePointCandidateFinder Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	typedef itk::Image<double, 3> RealType;
	typedef itk::Image<unsigned char, 3> LabelType;
	
	typedef RealType::PointType PointType;
	typedef std::pair<double, PointType> PointPairType;
	typedef std::vector<PointPairType> PointListType;

	itkNewMacro(Self);
	void SetInput(const RealType::Pointer &input) { m_Input = input; }
	void Compute();

	void GetOutput(PointListType &list, const unsigned int num=5);

protected:
	ValvePointCandidateFinder() {}
	virtual ~ValvePointCandidateFinder() {}


private:
	ValvePointCandidateFinder(const Self&);
	void operator=(const Self&);

	RealType::Pointer m_Input;
	PointListType m_Output;

};



}

#endif
