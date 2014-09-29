#ifndef FEATURE_POINT_ORIENTATOR_H
#define FEATURE_POINT_ORIENTATOR_H

#include <itkProcessObject.h>
#include <itkImage.h>
#include "FeatureCommon.h"

namespace filter
{
template<typename TInputType>
class KeyPointOrientator : public itk::ProcessObject 
{
public:
	typedef KeyPointOrientator Self;
	typedef itk::ProcessObject Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(KeyPointOrientator, ProcessObject);
	itkNewMacro(Self);

	typedef TInputType InputType;
	typedef typename InputType::PointType PointType;
	typedef std::pair<double, std::vector<PointType> > KeyPointSet;
	typedef std::map<double, std::vector<PointType> > KeyPointInputType;
	



	typedef itk::Image<double, InputType::ImageDimension> RealType;
	typedef itk::CovariantVector<double, InputType::ImageDimension> VectorType;
	typedef itk::Image<VectorType, InputType::ImageDimension> GradientType;


	typedef filter::OrientatedKeyPoint<InputType::ImageDimension> OrientatedKeyPoint;
	

	typedef std::vector<OrientatedKeyPoint> OrientatedKeyPointList;
	typedef std::pair<double, OrientatedKeyPointList> OrientedKeyPointPair;
	typedef std::map<double, OrientatedKeyPointList> OrientatedKeyPointMap; 
	typedef OrientatedKeyPointMap OutputType;

	OutputType GetOutput() const;



	void SetInput(const InputType * input);
	const InputType * GetInput() const;

	void SetKeyPoints(const KeyPointInputType &input);

	itkSetMacro(HistogramBins, unsigned int);
	itkSetMacro(SigmaScale, double);
	itkSetMacro(SampleRadius, unsigned int);
	itkSetMacro(SecondaryPeakThreshold, double);

	void Update();

protected:
	KeyPointOrientator();
	~KeyPointOrientator() {}

	void IdentifyPeaks(const std::vector<double> &histogram,
			std::vector<int> &peakBinNumbers);


	void ComputeOrientation(const GradientType * gradient, 
			double sigma, const PointType &point, 
			OrientatedKeyPointList &features);

	unsigned int GetHistogramIndex(double value);
	double GetFinalPeakPosition(const std::vector<double> &histogram, int index);



private:
	KeyPointOrientator(const Self&);
	void operator=(const Self&);

	unsigned int m_HistogramBins;
	double m_HistgramSpacing;
	double m_SigmaScale;
	unsigned int m_SampleRadius;
	unsigned int m_Dimensions;
	double m_SecondaryPeakThreshold;
	
	KeyPointInputType m_KeyPoints;
	
	OutputType m_Output;
};
} /* filter */ 


#include "KeyPointOrientator.hpp"

#endif
