#ifndef MRF_ISING_SMOOTHNESS_TERM_HPP
#define MRF_ISING_SMOOTHNESS_TERM_HPP


#include "MRFIsingSmoothnessTerm.h"

namespace segmentation
{
// ------------------------------------------------------------------------
template<typename TInputValueType, typename TOutputValueType>
MRFIsingSmoothnessTerm<TInputValueType, TOutputValueType>::
MRFIsingSmoothnessTerm()
{
	m_Sigma = 1.0;
	Superclass::SetNumberOfParameters(3);
}


// ------------------------------------------------------------------------
template<typename TInputValueType, typename TOutputValueType>
typename MRFIsingSmoothnessTerm<TInputValueType, TOutputValueType>::OutputType
MRFIsingSmoothnessTerm<TInputValueType, TOutputValueType>::
Evaluate(const InputType &input) const
{
	if(input.Size() != Superclass::m_NumberOfParameters)
	{
		itkExceptionMacro(<< "Not the expected number of paramters");
	}

	InputValueType v1 = input(0);
	InputValueType v2 = input(1);
	InputValueType dist = input(2);

	// check if dist is unset
	if(dist <= 0) dist = 1.0;


	// compute the smoothness prior
	double numer = (v1-v2)*(v1-v2);
	double denom = 2*m_Sigma*m_Sigma;
	double expVal = exp(- (numer / denom)); 

	return static_cast<TOutputValueType>( expVal * (1.0/dist) );
}


} /* segmentation */ 

#endif
