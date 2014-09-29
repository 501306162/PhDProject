#ifndef DO_G_FEATURE_POINT_EXTRACTOR_HPP
#define DO_G_FEATURE_POINT_EXTRACTOR_HPP

#include "DoGKeyPointExtractor.h"
#include "ImageToDoGImageFilter.h"
#include "CommonDefinitions.h"

#include <itkDiscreteGaussianImageFilter.h>
#include <itkGradientImageFilter.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>

namespace filter
{
// ------------------------------------------------------------------------
template<typename TInputType>
DoGKeyPointExtractor<TInputType>::
DoGKeyPointExtractor()
{
	this->SetNumberOfRequiredInputs(1);
	m_StartingSigma = 1.6;
	m_SplitsPerOctave = 3.0;
	m_StoppingPercentage = 2.0;
	m_KeyPointSearchSize = 2;
	m_KeypointThreshold = 5.0;

	m_MaxSigma = 0.0;
	m_Octaves = 0;
	m_Dimensions = InputType::ImageDimension;

	m_DistanceMap = 0;
	m_DistanceThreshold = 0.0;
	m_UseDistanceThreshold = false;

}


// ------------------------------------------------------------------------
template<typename TInputType>
const typename DoGKeyPointExtractor<TInputType>::InputType *
DoGKeyPointExtractor<TInputType>::
GetInput() const
{
	return itkDynamicCastInDebugMode<const InputType*>(this->ProcessObject::GetInput(0));
}

// ------------------------------------------------------------------------
template<typename TInputType>
void
DoGKeyPointExtractor<TInputType>::
SetInput(const InputType * input)
{
	this->SetNthInput(0, 
			const_cast<InputType*>(input));
}


// ------------------------------------------------------------------------
template<typename TInputType>
void
DoGKeyPointExtractor<TInputType>::
SetDistanceMap(const RealType * map)
{
	m_DistanceMap = const_cast<RealType*>(map);
	m_UseDistanceThreshold = true;
}

// ------------------------------------------------------------------------
template<typename TInputType>
void
DoGKeyPointExtractor<TInputType>::
Update()
{
	PrepareScaleSpace();


	// prepare the pyramid
	m_Pyramid = PyramidFilterType::New();
	m_Pyramid->SetInput(this->GetInput());
	m_Pyramid->SetNumberOfLevels(m_Octaves);
	m_Pyramid->Update();


	

	OutputType keyPoints;
	int count = 0;
	for(int i = m_Octaves-1; i >= 0; i--)
	{
		// compute the dog images for each octave
		typename RealType::Pointer image = m_Pyramid->GetOutput(i);

		std::vector<RealPointer> dogImages;

		for(unsigned int j = 0; j < m_SigmaValues.size()-1; j++)
		{
			double s1 = m_SigmaValues[j];
			double s2 = m_SigmaValues[j+1];

			//std::cout << i << " " << s1 << " " << s2 << std::endl;

			typedef ImageToDoGImageFilter<RealType, RealType> DoGFilterType;
			typename DoGFilterType::Pointer dogFilter = DoGFilterType::New();
			dogFilter->SetInput(m_Pyramid->GetOutput(i));
			dogFilter->SetSigma1(s1);
			dogFilter->SetSigma2(s2);
			dogFilter->SetUseImageSpacing(true);
			dogFilter->Update();

			dogImages.push_back(dogFilter->GetOutput());
		}

		GetKeyPointsAtLevel(dogImages, count, keyPoints);
		count++;
	}


	m_Output = keyPoints;
}

// ------------------------------------------------------------------------
template<typename TInputType>
void
DoGKeyPointExtractor<TInputType>::
GetKeyPointsAtLevel(const std::vector<RealPointer> &dogImages,
			int level, OutputType & keyPoints)
{
	for(unsigned int i = 1; i < dogImages.size()-1; i++)
	{
		double octaveSigma = m_SigmaValues[i-1];
		double actualSigma = octaveSigma * pow(2.0, static_cast<double>(level));
		
		GetKeyPoints(dogImages[i], dogImages[i-1], dogImages[i+1], keyPoints, actualSigma);
	}

}


// ------------------------------------------------------------------------
template<typename TInputType>
void
DoGKeyPointExtractor<TInputType>::
GetKeyPoints(const RealType * input, const RealType *rb, const RealType * ra,
			OutputType &keyPoints, double sigma)
{

	// create the output pair
	KeyPointSet outputSet;
	outputSet.first = sigma;
	
	
	

	// create the neighbourhood iterator	
	typedef itk::ConstNeighborhoodIterator<RealType> IteratorType;
	typename IteratorType::SizeType radius;
	radius.Fill(m_KeyPointSearchSize);

	// set the search region so as not to worry about the image boundaries
	RealRegionType region = input->GetLargestPossibleRegion();
	RealSizeType size = region.GetSize();
	RealIndexType index = region.GetIndex();
	
	for(unsigned int i = 0; i < m_Dimensions; i++)
	{
		size[i] -= m_KeyPointSearchSize;
		index[i] += m_KeyPointSearchSize;				
	}

	region.SetSize(size);
	region.SetIndex(index);
	

	// create the iterators
	IteratorType it[3];
	it[0] = IteratorType(radius, input, region);
	it[1] = IteratorType(radius, rb, region);
	it[2] = IteratorType(radius, ra, region);


	const unsigned int regionSize = pow(1+m_KeyPointSearchSize*2, m_Dimensions);
	const unsigned int centerIndex  = regionSize / 2;





	while(!it[0].IsAtEnd())
	{
		if(m_UseDistanceThreshold)
		{
			PointType point;
			input->TransformIndexToPhysicalPoint(it[0].GetIndex(), point);
			RealIndexType index;
			m_DistanceMap->TransformPhysicalPointToIndex(point, index);

			double distanceValue = m_DistanceMap->GetPixel(index);
			if(distanceValue > m_DistanceThreshold)
			{
				++(it[0]); ++(it[1]); ++(it[2]);
				continue;
			}
		}


		bool isMin = true;
		bool isMax = true;
		double testValue = it[0].GetCenterPixel();

		if(fabs(testValue) > m_KeypointThreshold)
		{

			// iterate through all the neighbourhoods
			for(unsigned int i = 0; i < 3; i++)
			{
				for(unsigned int j = 0; j < regionSize;	j++)
				{
					// check for the center pixels
					if(i==0 && j == centerIndex) continue;

					double value = it[i].GetPixel(j);
					if(value < testValue)
					{
						isMin = false;					
					}

					if(value > testValue)
					{
						isMax = false;
					}

					if(!isMax && !isMin) break;
				}
				if(!isMax && !isMin) break;
			}

			if(isMax || isMin)
			{		
				RealIndexType index = it[0].GetIndex();

				// add the point to the keypoint list
				PointType point;
				input->TransformIndexToPhysicalPoint(index,point);
				outputSet.second.push_back(point);
			}

		}

		++(it[0]); ++(it[1]); ++(it[2]);
	}


	keyPoints.insert(outputSet);
}


// ------------------------------------------------------------------------
template<typename TInputType>
void
DoGKeyPointExtractor<TInputType>::
PrepareScaleSpace()
{
	// get the smallest dimension of the input
	const InputType * input = GetInput();
	SizeType inputSize = input->GetLargestPossibleRegion().GetSize();
	SizeValueType min = itk::NumericTraits<SizeValueType>::max();
	unsigned int minIndex = itk::NumericTraits<unsigned int>::max();
	for(unsigned int i = 0; i < m_Dimensions; i++)
	{
		if(inputSize[i] < min)
		{
			min = inputSize[i];
			minIndex = i;
		}
	}

	// set the stopping value as half the smallest dimension
	SpacingType inputSpacing = input->GetSpacing();
	double largestDistance = min * inputSpacing[minIndex];
	m_MaxSigma = largestDistance / m_StoppingPercentage;


	// compute the sigma increments per octave
	m_SigmaValues.clear();
	double currentSigma = m_StartingSigma;
	for(unsigned int i = 0; i <= m_SplitsPerOctave; i++)
	{
		m_SigmaValues.push_back(currentSigma);
		
		double k = pow(2, 1.0 / static_cast<double>(m_SplitsPerOctave));
		currentSigma *= k;
		
	}

	// compute the number of ocaves
	currentSigma = m_StartingSigma;
	
	int iterationCount = 0;
	m_Octaves = 0;
	while(currentSigma < m_MaxSigma)
	{
		if(iterationCount == m_SplitsPerOctave)
		{		
			m_Octaves++;
			iterationCount = 0;
		}
		
		double k = pow(2, 1.0 / static_cast<double>(m_SplitsPerOctave));
		currentSigma *= k;
		iterationCount++;

	}
}

} /* filter */ 

#endif
