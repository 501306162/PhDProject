#ifndef MRFGMM_DATA_TERM_H
#define MRFGMM_DATA_TERM_H

#include "MRFDataTermBase.h"
#include <MatrixCommon.h>
#include <itkGaussianDistribution.h>
#include <itkGaussianMembershipFunction.h>

namespace segmentation
{
template<typename TInputValueType, typename TOutputValueType, int TLabelNum>
class MRFGMMDataTerm : public MRFDataTermBase<TInputValueType, TOutputValueType, TLabelNum>
{
public:
	typedef MRFGMMDataTerm Self;
	typedef MRFDataTermBase<TInputValueType, TOutputValueType, TLabelNum> Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(MRFGMMDataTerm, MRFDataTermBase);
	itkNewMacro(Self);

	typedef utils::DoubleMatrixType GMMParametersType;
	typedef std::vector<GMMParametersType> GMMParametersListType;
	typedef itk::Vector<TInputValueType, 1> MeasurementVectorType;
	typedef itk::Statistics::GaussianMembershipFunction<MeasurementVectorType> DistributionType;
	typedef typename DistributionType::MeanVectorType MeanVectorType;
	typedef typename DistributionType::CovarianceMatrixType CovarianceMatrixType;
	typedef std::vector<typename DistributionType::Pointer> DisttributionListType;




	typedef typename Superclass::InputType InputType;
	typedef TOutputValueType OutputValueType;
	typedef typename Superclass::OutputType OutputType;

	virtual OutputType Evaluate(const InputType &input) const;

	void SetParameters(unsigned int idx, const GMMParametersType & parameters);

	void Initialise();
		
protected:
	MRFGMMDataTerm();
	~MRFGMMDataTerm() {}

private:
	MRFGMMDataTerm(const Self &);
	void operator=(const Self &);
	
	bool m_Initialised;


	GMMParametersListType m_Parameters;

	std::vector<DisttributionListType> m_Distributions;


};

} /* segmentation */ 

#include "MRFGMMDataTerm.hpp"

#endif
