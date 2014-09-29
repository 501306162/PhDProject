#ifndef FEATURE_POINT_ORIENTATOR_HPP
#define FEATURE_POINT_ORIENTATOR_HPP

#include "KeyPointOrientator.h"
#include "CommonDefinitions.h"

#include <itkDiscreteGaussianImageFilter.h>
#include <itkGradientRecursiveGaussianImageFilter.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkGradientImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include "ImagePatchExtractor.h"
#include <Eigen/Dense>

namespace filter
{
// ------------------------------------------------------------------------
template<typename TInputType>
KeyPointOrientator<TInputType>::
KeyPointOrientator()
{
	this->SetNumberOfRequiredInputs(1);
	m_HistogramBins = 32;
	m_SigmaScale = 2.0;
	m_SampleRadius = 5;
	m_SecondaryPeakThreshold = 0.2;

	m_Dimensions = InputType::ImageDimension;
}


// ------------------------------------------------------------------------
template<typename TInputType>
void
KeyPointOrientator<TInputType>::
SetKeyPoints(const KeyPointInputType &input)
{
	m_KeyPoints = input;
}

// ------------------------------------------------------------------------
template<typename TInputType>
void
KeyPointOrientator<TInputType>::
SetInput(const InputType * input)
{
	this->SetNthInput(0, const_cast<InputType*>(input));
}


// ------------------------------------------------------------------------
template<typename TInputType>
const typename KeyPointOrientator<TInputType>::InputType *
KeyPointOrientator<TInputType>::
GetInput() const
{
	return itkDynamicCastInDebugMode<const InputType*>(this->ProcessObject::GetInput(0));
}


// ------------------------------------------------------------------------
template<typename TInputType>
void
KeyPointOrientator<TInputType>::
Update()
{
	// some checks
	if(m_KeyPoints.size() == 0)
	{
		itkExceptionMacro(<< "KeyPoints have not been set");
	}

	// compute the histgram spacing
	m_HistgramSpacing = 360.0 / static_cast<double>(m_HistogramBins);

	// now we loop through the key points
	typename KeyPointInputType::iterator mapIt = m_KeyPoints.begin();
	std::cout << m_KeyPoints.size() << std::endl;
	while(mapIt != m_KeyPoints.end())
	{
		double sigma = mapIt->first;
		unsigned int numFeatures = mapIt->second.size();

		// if no features ignore
		if(numFeatures == 0)
		{
			++mapIt;
			continue;
		}

		OrientedKeyPointPair outputPair;
		outputPair.first = sigma;

		typedef itk::RecursiveGaussianImageFilter<InputType, RealType> Smoother;
		typename Smoother::Pointer smoother = Smoother::New();
		smoother->SetInput(this->GetInput());
		smoother->SetSigma(sigma);
		smoother->Update();

		typedef itk::GradientImageFilter<RealType, double, double> GradientFilterType;
		typename GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
		gradientFilter->SetInput(smoother->GetOutput());
		gradientFilter->SetUseImageSpacing(true);
		gradientFilter->Update();

		for(unsigned int i = 0; i < mapIt->second.size(); i++)
		{
			PointType point = mapIt->second[i];
			ComputeOrientation(gradientFilter->GetOutput(), sigma, point, outputPair.second);
		}

		m_Output.insert(outputPair);
		
		++mapIt;
	}
	
}


// ------------------------------------------------------------------------
template<typename TInputType>
void
KeyPointOrientator<TInputType>::
ComputeOrientation(const GradientType * gradient, double sigma, const PointType &point,
		OrientatedKeyPointList &features)
{
	// extract the patch
	typedef ImagePatchExtractor<GradientType> ExtractorType;
	typename ExtractorType::Pointer extractor = ExtractorType::New();
	extractor->SetInput(gradient);
	extractor->SetScale(sigma*m_SigmaScale);

	typename ExtractorType::SizeType patchSize;
	patchSize.Fill(m_SampleRadius*2);
	extractor->SetPatchSize(patchSize);
	extractor->SetInputPoint(point);
	extractor->Update();

	typedef itk::ImageRegionConstIterator<GradientType> IteratorType;
	IteratorType it(extractor->GetOutput(),
			extractor->GetOutput()->GetLargestPossibleRegion());


	// create the histogram
	std::vector<double> histogram(m_HistogramBins, 0.0);


	it.GoToBegin();
	while(!it.IsAtEnd())
	{
		// evalute the gradient at this point
		typename GradientType::PixelType grad = it.Get();

		// normalise the gradient
		double mag = utils::MagnitudeFrom2DVector(grad[0], grad[1]);
		double gx = grad[0] / mag;
		double gy = grad[1] / mag;


		double angle = utils::DegreeAngleFrom2DVector(gx, gy);

		unsigned int histBin = GetHistogramIndex(angle);
		histogram[histBin] += mag;

		++it;
	}


	std::vector<int> peaks;
	IdentifyPeaks(histogram, peaks);


	for(unsigned int p = 0; p < peaks.size(); p++)
	{
		double angle = GetFinalPeakPosition(histogram, peaks[p]);

		// add to the output
		OrientatedKeyPoint feature;
		feature.scale = sigma;
		feature.degrees = angle;

		// wrap and convert angle to radians
		angle = utils::WrapAngle(angle);

		feature.angle = utils::DegreesToRadians(angle);
		feature.location = point;

		features.push_back(feature);
	}

}

// ------------------------------------------------------------------------
template<typename TInputType>
typename KeyPointOrientator<TInputType>::OutputType
KeyPointOrientator<TInputType>::
GetOutput() const
{
	return m_Output;
}

// ------------------------------------------------------------------------
template<typename TInputType>
double
KeyPointOrientator<TInputType>::
GetFinalPeakPosition(const std::vector<double> &histogram, int index)
{
	// get the values of the index and the surrounding points
	double y2 = histogram[index];
	double y1=0.0,y3=0.0;

	double x2 = index * m_HistgramSpacing;
	double x1=0.0,x3=0.0;

	if(index == m_HistogramBins-1)
	{
		y3 = histogram[0];
		x3 = 0.0;
	}
	else
	{
		y3 = histogram[index+1];
		x3 = x2+m_HistgramSpacing;
	}

	if(index == 0)
	{
		y1 = histogram[m_HistogramBins-1];
		x1 = (m_HistogramBins-1)*m_HistgramSpacing;
	}
	else
	{
		y1 = histogram[index-1];
		x1 = x2-m_HistgramSpacing;
	}


	// build the matrix 
	Eigen::Matrix3f A;
	Eigen::Vector3f b;
	b << y1, y2, y3;
	A << x1*x1, x1, 1,   x2*x2, x2, 1,   x3*x3, x3, 1;
	Eigen::Vector3f x = A.colPivHouseholderQr().solve(b);


	// get the parabola vertex
	double xout = -(x(1) / (2*x(0)));
	return xout;
}

// ------------------------------------------------------------------------
template<typename TInputType>
void
KeyPointOrientator<TInputType>::
IdentifyPeaks(const std::vector<double> &histogram,
			std::vector<int> &peakBinNumbers)
{
	// get the max peak
	int maxId = 0;
	double maxVal = 0.0;
	for(unsigned int i = 0; i < histogram.size(); i++)
	{
		if(histogram[i] > maxVal)
		{
			maxVal = histogram[i];
			maxId = i;
		}		
	}

	peakBinNumbers.push_back(maxId);

	// go through again to find secondary peaks
	for(unsigned int i = 0; i < histogram.size(); i++)
	{
		if(abs(i-maxId) < 2) continue;

		double val = histogram[i];
		double diff = maxVal - val;
		if(diff < (maxVal * m_SecondaryPeakThreshold))
		{
			peakBinNumbers.push_back(i);
		}
		
	}
}


// ------------------------------------------------------------------------
template<typename TInputType>
unsigned int
KeyPointOrientator<TInputType>::
GetHistogramIndex(double value)
{
	if(value <= 0.0)
		return 0;

	if(value >= 360.0)
		return m_HistogramBins-1;

	double dev = value / m_HistgramSpacing;
	unsigned int bin = static_cast<unsigned int>(std::floor(dev));
	return bin;
}




} /* filter */ 

#endif
