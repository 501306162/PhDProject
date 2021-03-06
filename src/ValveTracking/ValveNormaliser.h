#ifndef VALVE_NORMALISER_H
#define VALVE_NORMALISER_H

#include <ValveLine.h>
#include <itkCenteredAffineTransform.h>

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

	typedef itk::CenteredAffineTransform<double, 3> TransformType;
	void SetInput(const ValveType::Pointer &input) { m_Valve = input; }
	void SetFlip(bool flip) { m_Flip = flip; }
	void SetFlipPoints(bool flip) { m_FlipPoints = flip; }
	void Normalise();
	void UnNormalise();
	ValveType::Pointer GetOutput() const { return m_Output; }

	void FlipValve(const ValveType::Pointer &input, ValveType::Pointer &output);
	void AlignValve(const ValveType::Pointer &input, ValveType::Pointer &output);
	void FlipPoints(const ValveType::Pointer &input, ValveType::Pointer &output);


	TransformType::Pointer GetTransform() const { return m_Transform; }

	

protected:
	ValveNormaliser();
	virtual ~ValveNormaliser() {}



private:
	ValveNormaliser(const Self&);
	void operator=(const Self&);

	
	ValveType::Pointer m_Valve;
	ValveType::Pointer m_Output;
	bool m_Flip;
	bool m_FlipPoints;
	TransformType::Pointer m_Transform;


};



}

#endif
