#ifndef MANIFOLD_SHAPE_MODEL_H
#define MANIFOLD_SHAPE_MODEL_H

#include <itkObject.h>
#include <itkSimilarity3DTransform.h>
#include <itkObjectFactory.h>

#include <vector>

#include "Manifold.h"
#include "ShapeEmbedder.h"
#include "ImageToImageDistanceMeasure.h"

namespace manifold
{
template<int TDimension, typename TTransformType>
class ManifoldShapeModel : public itk::Object
{
public:
	typedef ManifoldShapeModel 			Self;
	typedef itk::Object 				Superclass;
	typedef itk::SmartPointer<Self>		Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkTypeMacro(ManifoldShapeModel, itk::Object);
	itkNewMacro(Self);

	typedef itk::Image<unsigned char, TDimension> LabelType;
	typedef typename LabelType::Pointer LabelPointer;
	typedef std::vector<LabelPointer> LabelList;

	typedef Manifold::Pointer ManifoldPointer;
	typedef Manifold::MatrixType MatrixType;
	

	typedef itk::Image<unsigned short, TDimension> ImageType;
	typedef itk::Image<double, TDimension> RealType;
	typedef typename RealType::Pointer RealPointer;
	typedef std::vector<RealPointer> RealList;

	typedef ImageToImageDistanceMeasure<RealType> DistanceMeasureType;
	typedef typename DistanceMeasureType::Pointer DistanceMeasurePointer;
	
	typedef TTransformType TransformType;
	typedef typename TransformType::Pointer TransformPointer;
	typedef typename TransformType::ParametersType TransformParametersType;


	typedef ShapeEmbedder<TDimension> ShapeEmbedderType;
	typedef typename ShapeEmbedderType::IndexDistancePair IndexDistancePair;
	typedef typename ShapeEmbedderType::IndexDistanceList IndexDistanceList;
	typedef typename ShapeEmbedderType::IndexListType IndexListType;
	typedef typename ShapeEmbedderType::WeightListType WeightListType;






	/**
	 * Set the images that will constitute the manifold
	 */
	void SetLabelData(const LabelList &images);

	void SetDistanceMaps(const RealList &maps);
	/**
	 * Set the mean image for the shape model
	 */
	void SetMeanLabel(const LabelType * mean);


	void ComputeWeightsFromDistances(const IndexDistanceList &input,
			std::vector<double> &weights);



	void Optimise(const LabelType * segmentationEstimate, 
			const ImageType * targetImage, 
			TransformPointer &transformParams,
			LabelPointer &output);


	void SetDistanceMeasure(const DistanceMeasureType * measure);




	itkSetMacro(Knn, unsigned int);



	/**
	 * Function to set the manifold building method
	 */
	void SetManifoldBuilder(Manifold * manifoldBuilder);

	itkSetMacro(EntropyWeight, double);
	itkSetMacro(OptimiserIterations, unsigned int);
		
	typedef itk::ImageBase<TDimension> BaseType;


protected:
	ManifoldShapeModel();

	void ApplyTransform(const LabelType * input,
			const TransformType *transform,
			const BaseType * ref,
			LabelPointer &output);

	void ApplyFinalTransform(const LabelType * input,
			const TransformType *transform,
			const BaseType * ref,
			LabelPointer &output);


	void ApplyTransform(const RealType * input,
			const TransformType *transform,
			const BaseType * ref,
			RealPointer &output);

	double ComputeCost(
			const LabelType * estimation, 
			const ImageType * target,
			const TransformType *transform,
			const IndexListType &indices,
			const WeightListType &theta);


	void ParseX(const std::vector<double> &x, WeightListType &weights,
			TransformParametersType &params);

	void ParseX(const WeightListType &weights, const TransformParametersType &params,
			std::vector<double> &x);

	
	void CreateShape(WeightListType &weights, IndexListType &indices,
			RealPointer &shape);

	void CreateLabel(WeightListType &weights, IndexListType &indices,
			LabelPointer &label);



	double ComputeConditionOneCost(const LabelType * estimation,
			const RealType * current);

	// function to compute the distance map
	void ComputeDistanceMap(const LabelType * input, 
			RealPointer &output);

	void ThresholdDistance(const RealType * input, LabelPointer &output);


	static double func(const std::vector<double> &x, std::vector<double> &grad, void *f_data);


private:
   	const ImageType *m_Target;
	const LabelType * m_Estimate;
	IndexListType m_Indices;
	TransformType * m_Transform;


	unsigned int m_OptimiserIterations;
	LabelList m_LabelData;
	RealList m_DistanceMaps;
	LabelPointer m_MeanLabel;
	DistanceMeasurePointer m_DistanceMeasure;
	ManifoldPointer m_Manifold;
	unsigned int m_Knn;
	double m_EntropyWeight;

};



} /* manifold */ 

#include "ManifoldShapeModel.hpp"
#endif
