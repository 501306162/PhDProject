#ifndef VALVE_NORMALISER_H
#define VALVE_NORMALISER_H

#include <ValveLine.h>

namespace vt
{
class ValveNormaliser : public itk::Object 
{
public:
	typedef ValveNormaliser Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkNewMacro(Self);
	itkTypeMacro(ValveNormaliser, Object);

	typedef ValveLine<3> ValveType;
	typedef ValveType::ImageType ImageType;
	typedef ValveType::PointType PointType;
	typedef ValveType::ContIndexType ContIndexType;

	void SetInput(const ValveType::Pointer &input) { m_Valve = input; }
	void SetFlip(bool flip) { m_Flip = flip; }
	void Normalise();
	ValveType::Pointer GetOutput() const { return m_Output; }

	void FlipValve(const ValveType::Pointer &input, ValveType::Pointer &output);
	void AlignValve(const ValveType::Pointer &input, ValveType::Pointer &output);


protected:
	ValveNormaliser() : m_Flip(false) {}
	virtual ~ValveNormaliser() {}



private:
	ValveNormaliser(const Self&);
	void operator=(const Self&);

	
	ValveType::Pointer m_Valve;
	ValveType::Pointer m_Output;
	bool m_Flip;


};



}

#endif
