#ifndef MYITKMRFIMAGEREGISTRATIONMETHODCONNECTIVITYFILTER_HPP
#define MYITKMRFIMAGEREGISTRATIONMETHODCONNECTIVITYFILTER_HPP

#include "MRFImageRegistrationMethodConnectivityFilter.h"
#include <itkImageRegionConstIterator.h>

namespace myitk
{
// ---------------------------------------------------------------------------------
template<typename TGridImage >
MRFImageRegistrationMethodConnectivityFilter<TGridImage>
::MRFImageRegistrationMethodConnectivityFilter() : m_GridImage(0)
{
	m_Dimensions = ImageDimension;
}

// ---------------------------------------------------------------------------------
template<typename TGridImage >
const typename MRFImageRegistrationMethodConnectivityFilter<TGridImage>::ConnectionListType
MRFImageRegistrationMethodConnectivityFilter<TGridImage>
::ComputeConnections()
{

	if( !m_GridImage )
	{
		itkExceptionMacro( << "Grid Image has not been set");
	}

	// to compute the connectivity we just iterate over the image. We visit
	// each grid point once, and log the forward connections
	typedef itk::ImageRegionConstIterator< GridImageType > IteratorType;
	IteratorType it( m_GridImage, m_GridImage->GetLargestPossibleRegion() );
	it.GoToBegin();

	ConnectionListType connectivityList;


	// get the size to accomodate the bounds of the grid
	typedef typename GridImageType::SizeType SizeType;
	typedef typename GridImageType::IndexType IndexType;
	SizeType gridSize = m_GridImage->GetLargestPossibleRegion().GetSize();

	while( !it.IsAtEnd() )
	{
		IndexType index = it.GetIndex();
		for( unsigned int dim = 0; dim < m_Dimensions; dim++ )
		{

			// check that the neighbour is not out of bounds
			if( index[dim] < gridSize[dim]-1 )
			{
				IndexType newIndex = index;
				newIndex[dim]+=1;				

				//compute the linear indices
				unsigned int p1 = ComputeLinearIndex( index, gridSize );
				unsigned int p2 = ComputeLinearIndex( newIndex, gridSize );

				//std::cout << index << " " << newIndex << std::endl;
				//std::cout << p1 << " " << p2  << std::endl;

				// create a new pair
				ConnectionType con(p1,p2);
				connectivityList.push_back( con );
			}
		}

		++it;
	}

	return connectivityList;
}


// ---------------------------------------------------------------------------------
template<typename TGridImage >
unsigned int 
MRFImageRegistrationMethodConnectivityFilter<TGridImage>
::ComputeLinearIndex( IndexType &index, SizeType &size )
{
	unsigned int outIndex = 0;

	for( unsigned int dim = 0; dim < m_Dimensions; dim++ )
	{
		int cval = index[dim];
		int mult = dim;
		while( mult > 0 )
		{
			cval *= size[mult];
			mult--;
		}
		outIndex+=cval;
	}

	return outIndex;
}



// ---------------------------------------------------------------------------------
template<typename TGridImage >
void 
MRFImageRegistrationMethodConnectivityFilter<TGridImage>
::SetGridImage( const TGridImage * grid )
{
	m_GridImage = grid;
}
}


#endif /* End of MYITKMRFIMAGEREGISTRATIONMETHODCONNECTIVITYFILTER_HPP */
