#ifndef SIMILARITY_MATRIX_BUILDER_H
#define SIMILARITY_MATRIX_BUILDER_H

#include <itkProcessObject.h>
#include <itkVariableSizeMatrix.h>
#include <Eigen/Dense>
#include <vector>

#include "ImageToImageDistanceMeasure.h"



namespace manifold
{
template<typename TInputImage>
class SimilarityMatrixBuilder : public itk::ProcessObject
{
public:
	typedef SimilarityMatrixBuilder 		Self;
	typedef itk::ProcessObject				Superclass;
	typedef itk::SmartPointer<Self> 		Pointer;
	typedef itk::SmartPointer<const Self> 	ConstPointer;
	itkTypeMacro(SimilarityMatrixBuilder, itk::ProcessObject);
	itkNewMacro(Self);



	/** Image Typedefs */
	typedef TInputImage ImageType;
	typedef typename ImageType::Pointer ImagePointer;
	typedef typename ImageType::ConstPointer ImageConstPointer;
	typedef std::vector<ImagePointer> ImageListType;

	/** Output Typedef */
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
	typedef MatrixType OutputType;
	
	/** Distance Measure Typedef */
	typedef ImageToImageDistanceMeasure<ImageType> DistanceMeasureType;
	typedef typename DistanceMeasureType::Pointer DistanceMeasurePointer;
	typedef typename DistanceMeasureType::ConstPointer DistanceMeasureConstPointer;



	void AddImage(const ImageType * image);
	void SetImages(const ImageListType &imageList);
	void SetDistanceMeasure(DistanceMeasureType * measure);
	void SetIsSymmetrical(bool val);
	const ImageType * GetInput(unsigned int index) const;
	

	const OutputType GetOutput() const;
	void Update();

	virtual ~SimilarityMatrixBuilder();


protected:

	void GenerateData();
	void Initialise() throw(itk::ExceptionObject);

	SimilarityMatrixBuilder();

private:
	DistanceMeasurePointer m_Measure;

	unsigned int m_NumberOfImages;
	bool m_IsSymmetrical;
	OutputType m_Output;


};


} /* manifold */ 

#include "SimilarityMatrixBuilder.hpp"

#endif
