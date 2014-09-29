#ifndef MRF_SMOOTHNESS_TERM_BASE_H
#define MRF_SMOOTHNESS_TERM_BASE_H

#include "MRFTermBase.h"

namespace segmentation
{
template<typename TInput, typename TOutputValueType>
class MRFSmoothnessTermBase : public MRFTermBase<TInput, TOutputValueType, 1>
{
public:
	typedef MRFSmoothnessTermBase Self;
	typedef MRFTermBase<TInput, TOutputValueType, 1> Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(MRFSmoothnessTermBase, MRFTermBase);


	typedef typename Superclass::InputType InputType;
	typedef TOutputValueType OutputValueType;
	typedef typename Superclass::OutputType OutputType;


	virtual OutputType Evaluate(const InputType &input) const = 0;

protected:
	MRFSmoothnessTermBase() {}
	virtual ~MRFSmoothnessTermBase() {}

private:
	MRFSmoothnessTermBase(const Self&);
	void operator=(const Self &);

};
} /* segmentation */ 

#endif
