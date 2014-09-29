#ifndef MRFGMM_DATA_TERM_HPP
#define MRFGMM_DATA_TERM_HPP

#include "MRFGMMDataTerm.h"


namespace segmentation
{

// ------------------------------------------------------------------------
template<typename TInputValueType, typename TOutputValueType, int TLabelNum>
MRFGMMDataTerm<TInputValueType, TOutputValueType, TLabelNum>::
MRFGMMDataTerm()
{
	m_Parameters.resize(TLabelNum);
	m_Distributions.resize(TLabelNum);
	Superclass::SetNumberOfParameters(1);
}

// ------------------------------------------------------------------------
template<typename TInputValueType, typename TOutputValueType, int TLabelNum>
typename MRFGMMDataTerm<TInputValueType, TOutputValueType, TLabelNum>::OutputType
MRFGMMDataTerm<TInputValueType, TOutputValueType, TLabelNum>::Evaluate(const InputType &input) const
{
	// create the output
	OutputType output;

	if(input.Size() != Superclass::m_NumberOfParameters)
	{
		itkExceptionMacro(<< "Number of paremeters is not as expected");
	}
	
	double sum = 0.0;

	// loop through the parameters array
	for(unsigned int i = 0; i < TLabelNum; i++)
	{

		unsigned int modes = m_Parameters[i].rows();
		for(unsigned int j = 0; j < modes; j++)
		{	
			MeasurementVectorType mv;
			mv[0] = input(0);
			double prob = m_Distributions[i][j]->Evaluate(mv) * m_Parameters[i](j,2);
			sum += prob;
			output[i] = prob;
		}		
	}

	for(unsigned int i = 0; i < TLabelNum; i++)
	{
		output[i] /= sum;
		output[i] = -log(output[i]);
	}

	return output;
}


// ------------------------------------------------------------------------
template<typename TInputValueType, typename TOutputValueType, int TLabelNum>
void
MRFGMMDataTerm<TInputValueType, TOutputValueType, TLabelNum>::
SetParameters(unsigned int idx, const GMMParametersType & parameters)
{
	if(idx < 0 || idx >= TLabelNum)
	{
		std::cout << "Index " << idx << " is out of range" << std::endl;
	}

	m_Parameters[idx] = parameters;
}


// ------------------------------------------------------------------------
template<typename TInputValueType, typename TOutputValueType, int TLabelNum>
void
MRFGMMDataTerm<TInputValueType, TOutputValueType, TLabelNum>::
Initialise()
{
	for(unsigned int i = 0; i < TLabelNum; i++)
	{
		if(m_Parameters[i].rows() == 0)
		{
			itkExceptionMacro(<< "Not all parameters have been set");
		}		
	}

	// set up the gaussian membership function
	for(unsigned int i = 0; i < TLabelNum; i++)
	{
		DisttributionListType distributions;
		unsigned int modes = m_Parameters[i].rows();
		for(unsigned int j = 0; j < modes; j++)
		{
			typename DistributionType::Pointer dist = DistributionType::New();
			
			MeanVectorType mean(1);
			mean[0] = m_Parameters[i](j,0);
			
			CovarianceMatrixType cov;
			cov.SetSize(1,1);
			cov[0][0] = m_Parameters[i](j,1);			
			dist->SetMean(mean);
			dist->SetCovariance(cov);

			distributions.push_back(dist);

		}

		m_Distributions[i] = distributions;
	}


}

} /* segmentation */ 

#endif
