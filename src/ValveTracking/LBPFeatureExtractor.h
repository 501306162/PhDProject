#ifndef LBP_FEATURE_EXTRACTOR_H
#define LBP_FEATURE_EXTRACTOR_H

#include <FeatureExtractor.h>


namespace vt
{
class LBPFeatureExtractor : public FeatureExtractor<itk::Image<unsigned short, 3> >
{
public:
	typedef LBPFeatureExtractor Self;
	typedef FeatureExtractor<itk::Image<unsigned short, 3> > Superclass;	
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	typedef itk::Image<unsigned int, 3> ImageType;
	

	itkTypeMacro(LBPFeatureExtractor, FeatureExtractor);
	itkNewMacro(Self);

	typedef Superclass::FeatureType FeatureType;
	virtual void Extract(FeatureType &feature);
	itkSetMacro(NumNeighbours, unsigned int);
	itkSetMacro(Radius, double);
	itkSetMacro(GridX, unsigned int);
	itkSetMacro(GridY, unsigned int);

protected:
	LBPFeatureExtractor();
	virtual ~LBPFeatureExtractor() {}


private:
	LBPFeatureExtractor(const Self&);
	void operator=(const Self&);

	unsigned int m_NumNeighbours;
	unsigned int m_Radius;
	unsigned int m_GridX;
	unsigned int m_GridY;

};


}

#endif
