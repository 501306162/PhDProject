#ifndef MRF_MEAN_SQUARE_UNARY_POTENTIALS_METRIC_HPP
#define MRF_MEAN_SQUARE_UNARY_POTENTIALS_METRIC_HPP


#include "MRFMeanSquaredUnaryPotentialsMetric.h"


namespace myitk
{
using namespace itk;
template<typename TFixedImage, typename TMovingImage>
void
MRFMeanSquaredUnaryPotentialsMetric<TFixedImage, TMovingImage>::
ComputeUnaryPotentials(const SampleRegion & region,  int threadId)
{
	int index = region.threadStartIndex;
	int endIndex = region.threadRange + index;

	//std::cout << index << " " << endIndex << " " << threadId << std::endl;

	for(; index < endIndex; index++)
	{
		FixedImageSamples samples = this->GetFixedImageSamplesList()[index];




		double msdSum = 0.0;
		int msdCount = 0;
		int miss = 0;

		m_OutputValue[index] = 0.0;

		for(unsigned int i = 0; i < samples.size(); i++)
		{
			FixedImageIndexType fixedIndex = samples[i].first;
			FixedImagePointType fixedPoint, transformedPoint;

			// transform the index to a point then transform it
			this->GetFixedImage()->TransformIndexToPhysicalPoint(fixedIndex, fixedPoint);
			transformedPoint = this->GetTransform()->TransformPoint(fixedPoint);

			//std::cout << fixedPoint << " " << transformedPoint << " ";


			// check that the transformed point is within the image bounds
			if(this->GetInterpolator()->IsInsideBuffer(transformedPoint))
			{
				MovingImageValueType movingValue = this->GetInterpolator()->Evaluate(transformedPoint);
				FixedImageValueType  fixedValue  = this->GetFixedImage()->GetPixel(fixedIndex);


				if( movingValue == fixedValue )
				{
					miss++;
				}

				//std::cout << movingValue << " " << fixedValue << std::endl;

				double diff = (double) fixedValue - (double) movingValue;
				double value = diff * diff;


				if(this->GetUseDistanceWeighting())
				{
					//value *= samples[i].second;
				}

				msdCount++;
				msdSum += value;
			}
		}


		// finalise the sum
		if(msdCount == 0)
		{
			//std::cout << "ERRRRORRRR" << std::endl;
		}
		else
		{
			//m_OutputValue[index] = msdSum / (double) msdCount;
			m_OutputValue[index] = msdCount - miss;
		}
	}
}


}


#endif
