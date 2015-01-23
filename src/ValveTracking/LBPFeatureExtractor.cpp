#include "LBPFeatureExtractor.h"

#include "lbp.h"
#include "lbp_hist.h"




namespace vt
{

// ------------------------------------------------------------------------
LBPFeatureExtractor::LBPFeatureExtractor()
{
	m_NumNeighbours = 8;
	m_Radius = 1;
	m_GridX = 2;
	m_GridY = 2;
}

// ------------------------------------------------------------------------
void LBPFeatureExtractor::Extract(FeatureType &feature)
{
	// conver the patch to an opencv image 
	unsigned int rows = m_Input->GetLargestPossibleRegion().GetSize()[1];
	unsigned int cols = m_Input->GetLargestPossibleRegion().GetSize()[0];

	cv::Mat cvPatch(rows, cols, CV_16SC1);

	for(unsigned int i = 0; i < rows; i++)
	{
		for(unsigned int j = 0; j < cols; j++)
		{
			ImageType::IndexType index;
			index[0] = j;
			index[1] = i;
			index[2] = 0;
			cvPatch.at<unsigned short>(i,j) = m_Input->GetPixel(index);
		}
	}

	cv::GaussianBlur(cvPatch, cvPatch, cv::Size(3,3), 1);
	cv::Mat lbp = lbp::ELBP(cvPatch, m_Radius, m_NumNeighbours);

	unsigned int numPatterns = static_cast<unsigned int>(std::pow(2.0, static_cast<double>(m_NumNeighbours)));
	cv::Mat hist = lbp::spatial_histogram(lbp, numPatterns, m_GridX, m_GridY, 0);

	feature = FeatureType(1, hist.cols);
	for(int i = 0; i < hist.cols; i++)
	{
		feature(0,i) = (double) hist.at<int>(0,i);			
	}
}




}
