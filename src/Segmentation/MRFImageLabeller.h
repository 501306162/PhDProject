#ifndef MRF_IMAGE_LABELLER_H
#define MRF_IMAGE_LABELLER_H

#include <itkImageSource.h>

namespace segmentation
{
template<typename TOutputImageType>
class MRFImageLabeller : public itk::ImageSource<TOutputImageType>
{
public:
	typedef MRFImageLabeller 				Self;
	typedef itk::ImageSource<TOutputImageType> Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	
	itkTypeMacro(MRFImageLabeller, ImageSource);
	itkNewMacro(Self);

	typedef TOutputImageType ImageType;
	typedef typename ImageType::Pointer ImagePointer;
	typedef typename ImageType::RegionType ImageRegionType;
	typedef typename ImageType::PointType PointType;
	typedef typename ImageType::SpacingType SpacingType;
	typedef typename ImageType::DirectionType DirectionType;
	typedef typename ImageType::SizeType SizeType;
	typedef typename ImageType::PixelType LabelType;
	typedef typename ImageType::IndexType IndexType;
	typedef std::vector<LabelType> LabelSetType;

	typedef itk::ImageBase<ImageType::ImageDimension> ImageBaseType;


	/** Set the output image parameters */
	itkSetMacro(Spacing, SpacingType);
	itkSetMacro(Direction, DirectionType);
	itkSetMacro(Origin, PointType);
	itkSetMacro(Size, SizeType);

	void UseImageParameters(const ImageBaseType * image);

	void SetLabelSet(const LabelSetType &labels);

protected:
	MRFImageLabeller();
	virtual ~MRFImageLabeller() {}

	virtual void GenerateOutputInformation();

	virtual void
  	ThreadedGenerateData(const ImageRegionType &
                       outputRegionForThread, itk::ThreadIdType threadId);


private:
	MRFImageLabeller(const Self&);
	void operator=(const Self&);

	SizeType      m_Size;      //size of the output image
	SpacingType   m_Spacing;   //spacing
	PointType     m_Origin;    //origin
	DirectionType m_Direction; //direction
	LabelSetType m_LabelSet;

};
} /* segmentation */ 


#include "MRFImageLabeller.hpp"

#endif
