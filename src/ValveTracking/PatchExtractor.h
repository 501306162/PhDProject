#ifndef PATCH_EXTRACTOR_H
#define PATCH_EXTRACTOR_H

#include <itkImage.h>
#include <ValveLine.h>

namespace vt 
{
template<typename PValue>
class PatchExtractor : public itk::Object 
{
public:
	typedef PatchExtractor Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(PatchExtractor, Object);
	itkNewMacro(Self);

	typedef ValveLine<3> ValveType;
	typedef itk::Image<PValue, 3> ImageType;
	typedef ValveType::PointType PointType;
	typedef ValveType::ContIndexType ContIndexType;
	typedef typename ImageType::SizeType SizeType;
	typedef itk::FixedArray<double, 2> LengthsType;
	typedef itk::Image<unsigned char, 3> MaskType;
	typedef MaskType::RegionType RegionType;
	typedef MaskType::ValueType ValueType;



	void SetImage(const typename ImageType::Pointer &image) { m_Image = image; }
	void SetPatchCenter(const ContIndexType & center) { m_Center = center; }
	void SetPatchSize(const SizeType & size) { m_Size = size; }

	typename ImageType::Pointer ExtractPatch();
	MaskType::Pointer ExtractMask();
	RegionType ExtractionRegion();

	itkSetMacro(MaskValue, ValueType);



protected:
	PatchExtractor() : m_MaskValue(1) {}
	virtual ~PatchExtractor() {}


private:

	typename ImageType::Pointer m_Image;
	SizeType m_Size;
	ContIndexType m_Center;
	ValueType m_MaskValue;

	PatchExtractor(const Self&);
	void operator=(const Self&);


};


typedef PatchExtractor<unsigned short> ImagePatchExtractor;
typedef PatchExtractor<unsigned char> MaskPatchExtractor;

}

#endif
