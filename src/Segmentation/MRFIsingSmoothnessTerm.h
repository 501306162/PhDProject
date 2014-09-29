#ifndef MRF_ISING_SMOOTHNESS_TERM_H
#define MRF_ISING_SMOOTHNESS_TERM_H

#include "MRFSmoothnessTermBase.h"

namespace segmentation
{
template<typename TInputValueType, typename TOutputValueType>
class MRFIsingSmoothnessTerm : public MRFSmoothnessTermBase<TInputValueType, TOutputValueType>
{
public:
	typedef MRFIsingSmoothnessTerm Self;
	typedef MRFSmoothnessTermBase<TInputValueType, TOutputValueType> Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(MRFIsingSmoothnessTerm, MRFSmoothnessTermBase);
	itkNewMacro(Self);


	typedef typename Superclass::InputType InputType;
	typedef TInputValueType InputValueType;
	typedef TOutputValueType OutputValueType;
	typedef typename Superclass::OutputType OutputType;


	virtual OutputType Evaluate(const InputType &input) const;

	itkSetMacro(Sigma, double);
	



protected:
	MRFIsingSmoothnessTerm(); 
	virtual ~MRFIsingSmoothnessTerm() {}

private:
	MRFIsingSmoothnessTerm(const Self&);
	void operator=(const Self &);


	double m_Sigma;



};
} /* segmentation */ 

#include "MRFIsingSmoothnessTerm.hpp"

#endif
