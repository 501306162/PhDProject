#ifndef MRF_DATA_TERM_BASE_H
#define MRF_DATA_TERM_BASE_H


#include "MRFTermBase.h"




namespace segmentation
{
template<typename TInputValueType, typename TOutputValueType, int TLabelNum>
class MRFDataTermBase : public MRFTermBase<TInputValueType, TOutputValueType, TLabelNum>
{
public: 	
	typedef MRFDataTermBase Self;
	typedef MRFTermBase<TInputValueType, TOutputValueType, TLabelNum> Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(MRFDataTermBase, MRFTermBase);

	typedef typename Superclass::InputType InputType;
	typedef TOutputValueType OutputValueType;
	typedef typename Superclass::OutputType OutputType;

	virtual OutputType Evaluate(const InputType & input) const = 0;


protected:
	MRFDataTermBase() {}
	virtual ~MRFDataTermBase() {}

private:
	MRFDataTermBase(const Self &);
	void operator=(const Self &);
};
} /* segmentation */ 

#endif
