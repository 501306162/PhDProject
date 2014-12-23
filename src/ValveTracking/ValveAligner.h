#ifndef VALVE_ALIGNER_H
#define VALVE_ALIGNER_H

#include <itkObject.h>
#include "ValveLine.h"

#include <itkSimilarity2DTransform.h>


namespace vt
{
class ValveAlignementTransformExtractor : public itk::Object
{
public:
	typedef ValveAlignementTransformExtractor Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	typedef itk::Similarity2DTransform<double> TransformType;
		
	typedef ValveLine<2> ValveType;

	itkTypeMacro(ValveAlignementTransformExtractor, Object);
	itkNewMacro(Self);

	void SetFixedValve(const ValveType::Pointer &valve) { m_FixedValve = valve; }
	void SetMovingValve(const ValveType::Pointer &valve) { m_MovingValve = valve; }

	void Extract();

	TransformType::Pointer GetTransform() const { return m_Transform; }

protected:
	ValveAlignementTransformExtractor() {}
	virtual ~ValveAlignementTransformExtractor() {}



private:

	ValveType::Pointer m_FixedValve;
	ValveType::Pointer m_MovingValve;
	TransformType::Pointer m_Transform;


	ValveAlignementTransformExtractor(const Self&);
	void operator=(const Self&);
};


class ValveLineResampler : public itk::Object
{
public:
	typedef ValveLineResampler Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	typedef itk::Similarity2DTransform<double> TransformType;
		
	typedef ValveLine<2> ValveType;
	typedef ValveType::ImageType ImageType;

	itkTypeMacro(ValveLineResampler, Object);
	itkNewMacro(Self);

	void SetInputValve(const ValveType::Pointer &valve) { m_InputValve = valve; }
	void SetReferenceValve(const ValveType::Pointer &valve) { m_ReferenceValve = valve; }
	void SetTransform(const TransformType::Pointer trans) { m_Transform = trans; }

	void Resample();
	ValveType::Pointer GetOutput() const { return m_Output; }

protected:
	ValveLineResampler() {}
	virtual ~ValveLineResampler() {}



private:

	ValveType::Pointer m_InputValve;
	ValveType::Pointer m_ReferenceValve;
	ValveType::Pointer m_Output;
	TransformType::Pointer m_Transform;


	ValveLineResampler(const Self&);
	void operator=(const Self&);

};


// ------------------------------------------------------------------------
class ValveSequenceAligner : public itk::Object
{
public:
	typedef ValveSequenceAligner Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	typedef itk::Similarity2DTransform<double> TransformType;
		
	typedef ValveLine<2> ValveType;
	typedef ValveType::ImageType ImageType;
	typedef ValveSequence<2> SequenceType;

	itkTypeMacro(ValveSequenceAligner, Object);
	itkNewMacro(Self);

	void SetFixedSequence(const SequenceType::Pointer &seq) { m_FixedSequence = seq; }
	void SetMovingSequence(const SequenceType::Pointer &seq) { m_MovingSequence = seq; }

	void Align();
	SequenceType::Pointer GetOutput() const { return m_Output; }

protected:
	ValveSequenceAligner() {}
	virtual ~ValveSequenceAligner() {}



private:

	SequenceType::Pointer m_FixedSequence;
	SequenceType::Pointer m_MovingSequence;
	SequenceType::Pointer m_Output;

	ValveSequenceAligner(const Self&);
	void operator=(const Self&);

};

}

#endif
