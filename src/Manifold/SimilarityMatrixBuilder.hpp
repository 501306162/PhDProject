#ifndef SIMILARITY_MATRIX_BUILDER_HPP
#define SIMILARITY_MATRIX_BUILDER_HPP

#include "SimilarityMatrixBuilder.h"

namespace manifold
{
// ------------------------------------------------------------------------
template<typename TInputImage>
SimilarityMatrixBuilder<TInputImage>::SimilarityMatrixBuilder()
{
	m_NumberOfImages = 0;
	m_Measure = 0;
	m_IsSymmetrical = false;
}

// ------------------------------------------------------------------------
template<typename TInputImage>
void SimilarityMatrixBuilder<TInputImage>::SetImages(const ImageListType &images)
{
	// add the images one by one
	for(unsigned int i = 0; i < images.size(); i++)
	{
		this->ProcessObject::PushBackInput( images[i] );
			//const_cast<ImageType *>(images[i]));	
		m_NumberOfImages++;
	}	
}

// ------------------------------------------------------------------------
template<typename TInputImage>
void SimilarityMatrixBuilder<TInputImage>::AddImage(const ImageType *image)
{
	// add the images one by one
	this->ProcessObject::PushBackInput( 
			const_cast<ImageType *>(image));	
	m_NumberOfImages++;
}


// ------------------------------------------------------------------------
template<typename TInputImage>
const typename SimilarityMatrixBuilder<TInputImage>::ImageType *
SimilarityMatrixBuilder<TInputImage>::GetInput(unsigned int index) const
{
	return itkDynamicCastInDebugMode<const ImageType*>(
			this->ProcessObject::GetInput(index));
}

// ------------------------------------------------------------------------
template<typename TInputImage>
void SimilarityMatrixBuilder<TInputImage>::GenerateData()
{
	// create the output matrix
	m_Output = MatrixType::Zero(m_NumberOfImages, m_NumberOfImages);
	
	// start looping through the data
	for(unsigned int i = 0; i < m_NumberOfImages; i++)
	{
		for(unsigned int j = 0; j < i; j++)
		{
			if(i == j) 
			{
				m_Output(i,j) = 0.0;
			}	
			else if((j < i) && m_IsSymmetrical)
			{
				m_Output(i,j) = m_Output(j,i);
			}
			else
			{
				
				// compute the measure
				ImageConstPointer im1 = this->GetInput(i);
				ImageConstPointer im2 = this->GetInput(j);
				m_Measure->SetInput1(im1);
				m_Measure->SetInput2(im2);
				m_Output(i,j) = m_Measure->GetValue();
			}
		}
	}
}


// ------------------------------------------------------------------------
template<typename TInputImage>
void SimilarityMatrixBuilder<TInputImage>::SetIsSymmetrical(bool val)
{
	m_IsSymmetrical = val;
}


// ------------------------------------------------------------------------
template<typename TInputImage>
void SimilarityMatrixBuilder<TInputImage>::Update()
{
	this->Initialise();
	this->GenerateData();
}

// ------------------------------------------------------------------------
template<typename TInputImage>
void SimilarityMatrixBuilder<TInputImage>::Initialise()
throw(itk::ExceptionObject)
{
	if(m_NumberOfImages == 0)
	{
		itkExceptionMacro( << "No Images Have Been Set");
	}

	if(!m_Measure)
	{
		itkExceptionMacro( << "No Distance Measure Has Been Set");
	}
}

// ------------------------------------------------------------------------
template<typename TInputImage>
SimilarityMatrixBuilder<TInputImage>::~SimilarityMatrixBuilder()
{
}

// ------------------------------------------------------------------------
template<typename TInputImage>
const typename SimilarityMatrixBuilder<TInputImage>::OutputType
SimilarityMatrixBuilder<TInputImage>::GetOutput() const
{
	return m_Output;
}	


// ------------------------------------------------------------------------
template<typename TInputImage>
void SimilarityMatrixBuilder<TInputImage>::SetDistanceMeasure(DistanceMeasureType * measure)
{
	m_Measure = measure;
}




} /* manifold */ 

#endif
