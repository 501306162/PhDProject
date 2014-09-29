#ifndef DO_G_FEATURE_POINT_EXTRACTOR_H
#define DO_G_FEATURE_POINT_EXTRACTOR_H

#include <itkProcessObject.h>
#include <itkImage.h>
#include <itkRecursiveMultiResolutionPyramidImageFilter.h>

namespace filter
{
template<typename TInputType>
class DoGKeyPointExtractor : public itk::ProcessObject
{
public:
	typedef DoGKeyPointExtractor Self;
	typedef itk::ProcessObject Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkTypeMacro(DoGKeyPointExtractor, ProcessObject);
	itkNewMacro(Self);

	typedef TInputType InputType;
	typedef typename InputType::SpacingType SpacingType;
	typedef typename InputType::SizeType SizeType;
	typedef typename InputType::SizeValueType SizeValueType;
	typedef typename InputType::PointType PointType;
	typedef itk::Image<double, InputType::ImageDimension> RealType;
	typedef typename RealType::Pointer RealPointer;
	typedef typename RealType::RegionType RealRegionType;
	typedef typename RealType::IndexType RealIndexType;
	typedef typename RealType::SizeType RealSizeType;

	typedef std::pair<double, std::vector<PointType> > KeyPointSet;
	typedef std::map<double, std::vector<PointType> > OutputType;
	
	typedef itk::CovariantVector<double, InputType::ImageDimension> GradientPixelType;
	typedef itk::Image<GradientPixelType, InputType::ImageDimension> GradientImageType;

	typedef itk::RecursiveMultiResolutionPyramidImageFilter<InputType, RealType> PyramidFilterType;
	typedef typename PyramidFilterType::Pointer PyramidFilterPointer;
	



	void SetInput(const InputType * input);
	const InputType * GetInput() const;

	void SetDistanceMap(const RealType * real);

	itkSetMacro(KeypointThreshold, double);
	itkSetMacro(StartingSigma, double);
	itkSetMacro(SplitsPerOctave, unsigned int);
	itkSetMacro(StoppingPercentage, double);
	itkSetMacro(DistanceThreshold, double);




	OutputType GetOutput() { return m_Output; }


	void Update();
protected:

	DoGKeyPointExtractor();
	virtual ~DoGKeyPointExtractor() {}

	void GetKeyPointsAtLevel(const std::vector<RealPointer> &dogImages,
			int level, OutputType & keypoints);

	void GetKeyPoints(const RealType * input, const RealType *rb, const RealType * ra,
			OutputType &keyPoints, double sigma);

	void PrepareScaleSpace();

private:
	DoGKeyPointExtractor(const Self&);
	void operator=(const Self&);

	typename RealType::Pointer m_DistanceMap;
	double m_DistanceThreshold;
	bool m_UseDistanceThreshold;

	double m_StartingSigma;
	double m_MaxSigma;
	unsigned int m_SplitsPerOctave;
	double m_StoppingPercentage;
	unsigned int m_Octaves;
	unsigned int m_Dimensions;
	std::vector<double> m_SigmaValues;

	PyramidFilterPointer m_Pyramid;
	unsigned int m_KeyPointSearchSize;
	double m_KeypointThreshold;
	OutputType m_Output;
	
	
};

} /* filter */ 

#include "DoGKeyPointExtractor.hpp"

#endif
