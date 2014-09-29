#ifndef MYITKMRFIMAGEREGISTRATIONMETHODCONNECTIVITYFILTER_H
#define MYITKMRFIMAGEREGISTRATIONMETHODCONNECTIVITYFILTER_H

#include <itkObject.h>
#include <itkImage.h>

namespace myitk
{
/**
 * @brief Class to compute the connectivity of the grid used in the MRF
 * registration class
 */
template<typename TGridImage >
class MRFImageRegistrationMethodConnectivityFilter : public itk::Object
{
public:
	typedef MRFImageRegistrationMethodConnectivityFilter	Self;
	typedef itk::Object										Superclass;
	typedef itk::SmartPointer< Self >						Pointer;
	typedef itk::SmartPointer< const Self >					ConstPointer;

	/** factory initialisation macro */
	itkNewMacro( Self );
	itkTypeMacro( MRFImageRegistrationMethodConnectivityFilter, itk::Object );

	/** typedefs for the grid image */
	typedef				TGridImage					GridImageType;
	typedef typename	GridImageType::ConstPointer GridImageConstPointer;
	typedef typename	GridImageType::Pointer		GridImagePointer;
	typedef typename	GridImageType::SizeType		SizeType;
	typedef typename	GridImageType::IndexType	IndexType;


	/** typedef the connectivity type */
	typedef std::pair<int, int>				ConnectionType;
	typedef std::vector< ConnectionType >	ConnectionListType;

	/** Get the image dimension */
	itkStaticConstMacro(ImageDimension,
                      unsigned int,
                      TGridImage::ImageDimension);

	/** Set / Get for the grid Image */
	void SetGridImage( const GridImageType * gridImage );
	itkGetConstObjectMacro(GridImage, GridImageType);

	const ConnectionListType ComputeConnections();

protected:
	MRFImageRegistrationMethodConnectivityFilter();
	~MRFImageRegistrationMethodConnectivityFilter() {}

	unsigned int ComputeLinearIndex( IndexType &index, SizeType &size );


private:
	MRFImageRegistrationMethodConnectivityFilter( const Self & );
	void operator=( const Self & );
	GridImageConstPointer m_GridImage;
	int m_Dimensions;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "MRFImageRegistrationMethodConnectivityFilter.hpp"
#endif

#endif /* End of MYITKMRFIMAGEREGISTRATIONMETHODCONNECTIVITYFILTER_H */
