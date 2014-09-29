#ifndef SHAPE_GENERATOR_H
#define SHAPE_GENERATOR_H

#include <itkImageToImageFilter.h>

namespace manifold
{
template<typename TInputType, typename TOutputType>
class ShapeGenerator : public itk::ImageToImageFilter<TInputType, TOutputType>
{
public:
	typedef ShapeGenerator Self;
	typedef itk::ImageToImageFilter<TInputType, TOutputType> Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkTypeMacro(ShapeGenerator, ImageToImageFilter);
	itkNewMacro(Self);

	typedef TInputType InputType;
	typedef typename InputType::ValueType InputValueType;

	typedef TOutputType OutputType;
	typedef typename OutputType::RegionType OutputRegionType;
	typedef typename OutputType::ValueType OutputValueType;

	typedef std::vector<double> WeightListType;

	void SetWeightList(const WeightListType &list);

protected:
	ShapeGenerator();
	virtual ~ShapeGenerator() {}

	void BeforeThreadedGenerateData();
	void ThreadedGenerateData(const OutputRegionType &region,
			itk::ThreadIdType id);

private:
	ShapeGenerator(const Self&);
	void operator=(const Self&);

	WeightListType m_WeightList;
};
} /* manifold */ 


#include "ShapeGenerator.hpp"

#endif
