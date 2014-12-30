#ifndef MRF_SEGMENTATION_H
#define MRF_SEGMENTATION_H


#include <itkProcessObject.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkImageToImageFilter.h>
#include "MRFDataTermBase.h"
#include "MRFSmoothnessTermBase.h"
#include "MRFPriorTermBase.h"
#include "mrf.h"

namespace segmentation
{
template<typename TImageType, typename TPriorType, typename TOutputType, int TLabelNum=2, typename TCostType=double>
class MRFSegmentation : public itk::ImageToImageFilter<TImageType, TOutputType> 
{
public:
	typedef MRFSegmentation 					Self;
	typedef itk::ImageToImageFilter<TImageType, TOutputType> Superclass;
	typedef itk::SmartPointer<Self>				Pointer;
	typedef itk::SmartPointer<const Self>       ConstPointer;
	itkTypeMacro(MRFSegmentation, ImageToImageFilter); 
	itkNewMacro(Self);

	/** Standard typedefs */
	typedef TImageType ImageType;
	typedef typename ImageType::Pointer ImagePointer;
	typedef typename ImageType::ConstPointer ImageConstPointer;
	typedef typename ImageType::PixelType PixelType;
	typedef typename ImageType::SizeType SizeType;
	typedef typename ImageType::IndexType IndexType;
	typedef typename ImageType::SpacingType SpacingType;
	typedef typename ImageType::PointType PointType;

	typedef TOutputType OutputType;
	typedef typename OutputType::RegionType OutputRegionType;

	typedef TPriorType PriorImageType;
	typedef typename PriorImageType::Pointer PriorImagePointer;
	typedef typename PriorImageType::ConstPointer PriorImageConstPointer;

	typedef itk::LinearInterpolateImageFunction<PriorImageType, double> PriorInterpolatorType;
	typedef typename PriorInterpolatorType::Pointer PriorInterpolatorPointer;
	


	typedef TOutputType OutputImageType;
	typedef typename OutputImageType::Pointer OutputImagePointer;
	typedef typename OutputImageType::ConstPointer OutputImageConstPointer;


	/** Typdef the terms */
	
	typedef MRFDataTermBase<PixelType, TCostType, TLabelNum> DataTermType;
	typedef typename DataTermType::InputType TermInputType;
	typedef typename DataTermType::OutputType DataTermOutputType;

	typedef typename DataTermType::Pointer DataTermPointer;
	typedef MRFSmoothnessTermBase<PixelType, TCostType> SmoothnessTermType;
	typedef typename SmoothnessTermType::Pointer SmoothnessTermPointer;




	itkSetMacro(PriorWeight, double);
	itkSetMacro(SmoothnessWeight, double);
	itkSetMacro(DataWeight, double);

	/** typedef COnnections */
	typedef std::pair<unsigned int, unsigned int> EdgeType;
	typedef std::vector<EdgeType> EdgeListType;
	typedef typename SmoothnessTermType::InputType SmoothnessInputType;
	typedef std::vector<TCostType> SmoothnessPriorValuesType;

	void SetInput(const ImageType * image);
	void SetPrior(const PriorImageType * prior);
	
	void SetDataTerm(DataTermType * term);
	void SetSmoothnessTerm(SmoothnessTermType * term);

protected:
	/** Generate Data functions */
	virtual void BeforeThreadedGenerateData();
	virtual void AfterThreadedGenerateData();
	virtual void ThreadedGenerateData(const OutputRegionType & outputRegionForThread,
			itk::ThreadIdType threadId);


	void ComputeEdges();
	unsigned int GetLinearIndex(const IndexType &index);

	MRFSegmentation();
	virtual ~MRFSegmentation(){};


private:
	MRFSegmentation(const Self&);
	void operator=(const Self&);

	PriorImageConstPointer m_Prior;
	SmoothnessTermPointer m_SmoothnessTerm;
	DataTermPointer m_DataTerm;

	unsigned int m_VertexNumber;

	/** variables for holding the cost values */
	MRF::CostVal * m_DataCosts;
	MRF::CostVal * m_SmoothnessCosts;

	EdgeListType m_Edges;
	SmoothnessPriorValuesType m_SmoothnessPriorValues;


	PriorInterpolatorPointer m_PriorInterpolator;
	bool m_UsePrior;

	
	/** Data term weights */
	double m_PriorWeight;
	double m_SmoothnessWeight;
	double m_DataWeight;


	/**
	 * Input Image
	 * Input Mask
	 * Data Term
	 * Smoothness Term
	 * Prior Term
	 */

};	
} /* mrf */ 


#include "MRFSegmentation.hpp"

#endif
