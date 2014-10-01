#ifndef STRUCTURE_TENSOR_KEY_POINT_VALIDATOR_HPP
#define STRUCTURE_TENSOR_KEY_POINT_VALIDATOR_HPP

#include "StructureTensorKeyPointValidator.h"

#include <itkImageRegionConstIterator.h>
#include <MatrixCommon.h>

//#include <itkDiscreteHessianGaussianImageFunction.h>



namespace filter
{
// ------------------------------------------------------------------------
template<typename TImageType>
StructureTensorKeyPointValidator<TImageType>::
StructureTensorKeyPointValidator()
{
	this->SetNumberOfRequiredInputs(1);
	m_Ratio = 10.0;
	m_Valid = false;
}


// ------------------------------------------------------------------------
template<typename TImageType>
void 
StructureTensorKeyPointValidator<TImageType>::
SetKeyPoint(const KeyPointType &keyPoint)
{
	m_KeyPoint = keyPoint;
}

// ------------------------------------------------------------------------
template<typename TImageType>
void 
StructureTensorKeyPointValidator<TImageType>::
SetInput(const ImageType * input)
{
	this->SetNthInput(0,const_cast<ImageType*>(input));
}

// ------------------------------------------------------------------------
template<typename TImageType>
const typename StructureTensorKeyPointValidator<TImageType>::ImageType *
StructureTensorKeyPointValidator<TImageType>::
GetInput() const
{
	return itkDynamicCastInDebugMode<const ImageType *>(this->ProcessObject::GetInput(0));
}

// ------------------------------------------------------------------------
template<typename TImageType>
void 
StructureTensorKeyPointValidator<TImageType>::
Update()
{
	/**
	typedef itk::DiscreteHessianGaussianImageFunction<ImageType> HessianFunctionType;
	typename HessianFunctionType::Pointer hessianFunction = HessianFunctionType::New();
	hessianFunction->SetInputImage(GetInput());
	hessianFunction->SetSigma(m_KeyPoint.scale);
	hessianFunction->SetNormalizeAcrossScale(true);
	hessianFunction->SetUseImageSpacing(true);
	hessianFunction->SetMaximumKernelWidth(10000);
	hessianFunction->Initialize();

	typedef typename HessianFunctionType::OutputType HessianType;
	HessianType hessian = hessianFunction->Evaluate(m_KeyPoint.location);

	typedef typename HessianType::EigenValuesArrayType EigenValues;
	EigenValues values; 
	hessian.ComputeEigenValues(values);
	double alpha = values[1];
	double beta = values[0];


	double det   = alpha + beta;
	double trace = alpha * beta;

	if(det < 0.0)
	{
		m_Valid = false;
		return;
	}

	double lhs = std::pow(trace,2.0) / det;
	double rhs = std::pow(m_Ratio+1, 2.0) / m_Ratio;

	if(lhs < rhs)
	{
		m_Valid = true;
	}
	else
	{
		m_Valid = false;
	}

	*/
}


// ------------------------------------------------------------------------
template<typename TImageType>
bool 
StructureTensorKeyPointValidator<TImageType>::
IsValid()
{
	if(m_Valid)
		return true;
	else
		return false;
}



} /* filter */ 

#endif
