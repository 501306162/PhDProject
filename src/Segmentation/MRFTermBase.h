#ifndef MRF_TERM_BASE_H
#define MRF_TERM_BASE_H

#include <itkFunctionBase.h>
#include <itkVector.h>
#include <itkOptimizerParameters.h>

namespace segmentation
{
template<typename TInputValueType, typename TOutputValueType, int TLabelNumber >
class MRFTermBase : public itk::FunctionBase<itk::OptimizerParameters<TInputValueType>, itk::Vector<TOutputValueType, TLabelNumber> > 
{
public:
	typedef MRFTermBase 			Self;
	typedef itk::FunctionBase<itk::OptimizerParameters<TInputValueType>, itk::Vector<TOutputValueType, TLabelNumber> > Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(MRFTermBase, FunctionBase);

	typedef itk::OptimizerParameters<TInputValueType> InputType;
	typedef TOutputValueType OutputValueType;
	typedef typename Superclass::OutputType OutputType;


	itkSetMacro(NumberOfParameters, unsigned int);
	virtual OutputType Evaluate(const InputType &input) const = 0;

protected:
	MRFTermBase() {}
	virtual ~MRFTermBase(){}

	unsigned int m_NumberOfParameters;


private:
	MRFTermBase(const Self &);
	void operator=(const Self &);
};
} /* segmentation */ 

#endif
