#include "ComputeProbImage.h"
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <LBPFeatureExtractor.h>


namespace vt
{
// ------------------------------------------------------------------------
ComputeProbImage::ComputeProbImage()
{
	m_HasMask = false;
}



// ------------------------------------------------------------------------
void ComputeProbImage::ThreadedGenerateData(const RegionType & outputRegionForThread, ThreadIdType threadId)
{
	const ImageType * input = this->GetInput();
	OutputType * output = this->GetOutput();

	itk::ImageRegionConstIterator<ImageType> inIt(input, outputRegionForThread);
	itk::ImageRegionIterator<OutputType> outIt(output, outputRegionForThread);
	itk::ImageRegionConstIterator<MaskType> maskIt(m_Mask, outputRegionForThread);


	while(!inIt.IsAtEnd())
	{
		IndexType index = inIt.GetIndex();
		if(maskIt.Get() == 255)
		{
			PointType point;
			input->TransformIndexToPhysicalPoint(index, point);	

			MatrixType feature;
			ExtractLBPFeature(point, feature);

			MatrixType probs;
			IntMatrixType classes;
			m_Classifier->PredictProbability(feature, classes, probs);
			
			outIt.Set(probs(0,1));
		}
		else
		{
			outIt.Set(0);
		}


		++inIt;
		++outIt;
		++maskIt;
	}

}

// ------------------------------------------------------------------------
void ComputeProbImage::ExtractLBPFeature(const PointType &point, MatrixType &feature)
{
	ExtractorType::Pointer extractor  = ExtractorType::New();
	extractor->SetInput(this->GetInput());
	extractor->SetSize(m_PatchSize);
	extractor->SetLine(m_Line);
	extractor->SetDistance(m_PatchDistance);
	extractor->SetCenter(point);

	extractor->Update();

	LBPFeatureExtractor::Pointer featureExtractor = LBPFeatureExtractor::New();
	featureExtractor->SetRadius(m_Radius);
	featureExtractor->SetGridX(m_GridSize);
	featureExtractor->SetGridY(m_GridSize);
	featureExtractor->SetNumNeighbours(m_Neighbors);
	featureExtractor->SetInput(extractor->GetOutput());
	featureExtractor->Extract(feature);
}





}
