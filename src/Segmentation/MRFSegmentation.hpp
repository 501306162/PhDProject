#ifndef MRF_SEGMENTATION_HPP
#define MRF_SEGMENTATION_HPP

#include "MRFSegmentation.h"
#include "GCoptimization.h"

#include <itkImageRegionConstIterator.h>
#include "MRFImageLabeller.h"

#include <itkLinearInterpolateImageFunction.h>


namespace segmentation
{
// ------------------------------------------------------------------------
template<typename TImageType, typename TPriorType, typename TOutputType, int TLabelNum, typename TCostType>
MRFSegmentation<TImageType, TPriorType, TOutputType, TLabelNum, TCostType>::
MRFSegmentation()
{
	this->ProcessObject::SetNumberOfRequiredInputs(1);	
	m_Prior = 0;
	m_SmoothnessTerm = 0;
	m_DataTerm = 0;
	m_UsePrior = false;

	m_DataWeight = 1.0;
	m_SmoothnessWeight = 1.0;
	m_PriorWeight = 1.0;
}


// ------------------------------------------------------------------------
template<typename TImageType, typename TPriorType, typename TOutputType, int TLabelNum, typename TCostType>
void
MRFSegmentation<TImageType, TPriorType, TOutputType, TLabelNum, TCostType>::
SetInput(const ImageType * image)
{
	// add the images one by one
	this->ProcessObject::SetNthInput(0, 
			const_cast<ImageType *>(image));	
}


// ------------------------------------------------------------------------
template<typename TImageType, typename TPriorType, typename TOutputType, int TLabelNum, typename TCostType>
void
MRFSegmentation<TImageType, TPriorType, TOutputType, TLabelNum, TCostType>::
SetPrior(const PriorImageType * prior)
{
	m_Prior = prior;
	m_UsePrior = true;
	m_PriorInterpolator = PriorInterpolatorType::New();
	m_PriorInterpolator->SetInputImage(m_Prior);
}


// ------------------------------------------------------------------------
template<typename TImageType, typename TPriorType, typename TOutputType, int TLabelNum, typename TCostType>
void
MRFSegmentation<TImageType, TPriorType, TOutputType, TLabelNum, TCostType>::
SetSmoothnessTerm(SmoothnessTermType * term)
{
	m_SmoothnessTerm = term;
}



// ------------------------------------------------------------------------
template<typename TImageType, typename TPriorType, typename TOutputType, int TLabelNum, typename TCostType>
void
MRFSegmentation<TImageType, TPriorType, TOutputType, TLabelNum, TCostType>::
SetDataTerm(DataTermType  * term)
{
	m_DataTerm = term;
}



// ------------------------------------------------------------------------
template<typename TImageType, typename TPriorType, typename TOutputType, int TLabelNum, typename TCostType>
void
MRFSegmentation<TImageType, TPriorType, TOutputType, TLabelNum, TCostType>::
BeforeThreadedGenerateData()
{
	// do the checking
	//
	//
	//
	if(!m_DataTerm)
	{
		itkExceptionMacro(<< "Data Term Not Set");
	}

	if(!m_SmoothnessTerm)
	{
		itkExceptionMacro(<< "Smoothness Term Not Set");
	}
	
	
	// set up some of the variables
	const ImageType * image = this->GetInput(0);
	m_VertexNumber = image->GetLargestPossibleRegion().GetNumberOfPixels();


	// initialise the cost arrays
	m_DataCosts = new MRF::CostVal[m_VertexNumber * TLabelNum];
	m_SmoothnessCosts = new MRF::CostVal[TLabelNum * TLabelNum];

	// assign the smoothness cost
	for(unsigned int i = 0; i < TLabelNum; i++)
	{
		for(unsigned int j = 0; j <TLabelNum; j++)
		{
			if(i == j)
				m_SmoothnessCosts[i*TLabelNum+j] = 0.0;
			else
				m_SmoothnessCosts[i*TLabelNum+j] = m_SmoothnessWeight;

		}
	}



	// compute the connections
	this->ComputeEdges();

	
}

// ------------------------------------------------------------------------
template<typename TImageType, typename TPriorType, typename TOutputType, int TLabelNum, typename TCostType>
void
MRFSegmentation<TImageType, TPriorType, TOutputType, TLabelNum, TCostType>::
ComputeEdges()
{
	const ImageType * input = this->GetInput(0);		
	SizeType imageSize = input->GetLargestPossibleRegion().GetSize();
	SpacingType spacing = input->GetSpacing();

	typedef itk::ImageRegionConstIterator<ImageType> InputIteratorType;
	InputIteratorType it(input,input->GetLargestPossibleRegion());
	
	while(!it.IsAtEnd())
	{
		IndexType index = it.GetIndex();
		unsigned int idx = GetLinearIndex(index);



		for(unsigned int i = 0; i < ImageType::ImageDimension; i++)
		{
			if(index[i] < imageSize[i]-1)
			{
				// create the connection
				int idx2 = 0;
				if(i == 0)
					idx2 = idx+1;
				else
					idx2 = idx+imageSize[i-1];


				EdgeType con(idx, idx2);
				m_Edges.push_back(con);


				// compute the smoothness weight
				IndexType index2 = index;
				index2[i]+=1;
				
				PixelType v1 = it.Get();
				PixelType v2 = input->GetPixel(index2);

			
				SmoothnessInputType params(3);
				params(0) = v1;
				params(1) = v2;
				params(2) = spacing[i];

				double w = m_SmoothnessTerm->Evaluate(params)[0];
				m_SmoothnessPriorValues.push_back(w);
			}
		}

	
		++it;
	}
	


}


// ------------------------------------------------------------------------
template<typename TImageType, typename TPriorType, typename TOutputType, int TLabelNum, typename TCostType>
unsigned int
MRFSegmentation<TImageType, TPriorType, TOutputType, TLabelNum, TCostType>::
GetLinearIndex(const IndexType &index)
{
	SizeType size = this->GetInput(0)->GetLargestPossibleRegion().GetSize();
	unsigned int idx = 0;
	if(ImageType::ImageDimension == 2)
	{
		idx = (index[1]*size[0])+index[0];
	}
	else if(ImageType::ImageDimension == 3)
	{
		idx = (index[2]*size[1]*size[0])+(index[1]*size[0])+index[0];
	}

	return idx;
}

// ------------------------------------------------------------------------
template<typename TImageType, typename TPriorType, typename TOutputType, int TLabelNum, typename TCostType>
void
MRFSegmentation<TImageType, TPriorType, TOutputType, TLabelNum, TCostType>::
ThreadedGenerateData(const OutputRegionType & outputRegionForThread,
			itk::ThreadIdType threadId)
{
	// iterate through the image region and compute the 
	typedef itk::ImageRegionConstIterator<ImageType> InputIteratorType;
	InputIteratorType it(this->GetInput(0), outputRegionForThread);


	while(!it.IsAtEnd())
	{
		// get the pixel index
		unsigned int idx = GetLinearIndex(it.GetIndex());

		// get the data term value
		TermInputType input(1);
		input(0) = it.Get();

		DataTermOutputType output = m_DataTerm->Evaluate(input);


		// assign the cost values
		for(unsigned int j = 0; j < TLabelNum; j++)
		{
			m_DataCosts[idx*TLabelNum+j] = output[j] * m_DataWeight;
		}


		// do the prior values? Not this will only work with a two label problem (at the moment)
		if(m_UsePrior)
		{
			PointType point;
			this->GetInput(0)->TransformIndexToPhysicalPoint(it.GetIndex(), point);

			// set very small value
			double priorValue = 0.00000001;
			if(m_PriorInterpolator->IsInsideBuffer(point))
			{
				priorValue = m_PriorInterpolator->Evaluate(point);
				//if(priorValue > 0.9)
					//std::cout << priorValue << std::endl;
			}

			if(priorValue < 0.00000001)
				priorValue = 0.00000001;


					
			m_DataCosts[idx*TLabelNum+0] += -log(1-priorValue) * m_PriorWeight;
			m_DataCosts[idx*TLabelNum+1] += -log(priorValue)   * m_PriorWeight;




		}


		++it;
	}
}


// ------------------------------------------------------------------------
template<typename TImageType, typename TPriorType, typename TOutputType, int TLabelNum, typename TCostType>
void
MRFSegmentation<TImageType, TPriorType, TOutputType, TLabelNum, TCostType>::
AfterThreadedGenerateData()
{
	DataCost * data = new DataCost(m_DataCosts);
	SmoothnessCost * smooth = new SmoothnessCost(m_SmoothnessCosts);
	EnergyFunction * eng = new EnergyFunction(data,smooth);


	// here we do the actual segmentation
	MRF * mrf = new Expansion(m_VertexNumber, TLabelNum, eng);

	for(unsigned int i = 0; i < m_Edges.size(); i++)
	{
		mrf->setNeighbors(m_Edges[i].first, m_Edges[i].second, m_SmoothnessPriorValues[i]);
		//mrf->setNeighbors(m_Edges[i].first, m_Edges[i].second, 1.0);
	}

	mrf->initialize();
	mrf->clearAnswer();
	float t;
	mrf->optimize(5,t);

	MRF::EnergyVal E_smooth = mrf->smoothnessEnergy();
    MRF::EnergyVal E_data   = mrf->dataEnergy();

	std::cout << E_smooth << " " << E_data << std::endl;

	// do the labelling operation
	typedef MRFImageLabeller<OutputType> LabellerType;
	typedef typename LabellerType::LabelSetType LabelSetType;
	typedef typename LabellerType::LabelType LabelType;
	LabelSetType labels;

	for(unsigned int i = 0; i < m_VertexNumber; i++)
	{
		LabelType label = static_cast<LabelType>(mrf->getLabel(i));
		labels.push_back(label);
	}


	typename LabellerType::Pointer labeller = LabellerType::New();
	labeller->SetLabelSet(labels);
	labeller->UseImageParameters(this->GetInput(0));
	labeller->Update();


	this->GraftOutput(labeller->GetOutput());

	delete mrf;

}



} /* segmentation */ 

#endif
