#include "BinaryPatchFeatureExtractor.h"
#include <itkImageRegionConstIterator.h>


namespace vt
{
// ------------------------------------------------------------------------
void BinaryPatchFeatureExtractor::Extract(FeatureType & feature)
{
	unsigned int numFeatures = m_Input->GetLargestPossibleRegion().GetNumberOfPixels();
	feature = Superclass::FeatureType::Zero(1, numFeatures);

	itk::ImageRegionConstIterator<ImageType> it(m_Input, m_Input->GetLargestPossibleRegion());
	unsigned int count = 0;
	while(!it.IsAtEnd())
	{
		feature(count) = it.Get();
		++it; ++count;
	}
}

}
