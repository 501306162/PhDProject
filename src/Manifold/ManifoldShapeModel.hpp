#ifndef MANIFOLD_SHAPE_MODEL_HPP
#define MANIFOLD_SHAPE_MODEL_HPP

#include "ManifoldShapeModel.h"

#include <itkResampleImageFilter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBinaryThresholdImageFilter.h>

#include "ShapeGenerator.h"
#include "ShapeEmbedder.h"
#include "ShapeEntropyCalculator.h"

#include "CommonDefinitions.h"

#include <nlopt.hpp>


namespace manifold
{
// ------------------------------------------------------------------------
template<int TDimension, typename TTransform> 
ManifoldShapeModel<TDimension, TTransform>::
ManifoldShapeModel()
{

	m_EntropyWeight = 1.0;
	m_OptimiserIterations = 100;
	m_Knn = 4;
}

// ------------------------------------------------------------------------
template<int TDimension, typename TTransform> 
void
ManifoldShapeModel<TDimension, TTransform>::
SetLabelData(const LabelList &labels)
{
	m_LabelData = labels;
}


// ------------------------------------------------------------------------
template<int TDimension, typename TTransform> 
void
ManifoldShapeModel<TDimension, TTransform>::
SetDistanceMaps(const RealList &maps)
{
	m_DistanceMaps = maps;
}

// ------------------------------------------------------------------------
template<int TDimension, typename TTransform> 
void
ManifoldShapeModel<TDimension, TTransform>::
SetMeanLabel(const LabelType *image)
{
	m_MeanLabel = const_cast<LabelType*>(image);
}

// ------------------------------------------------------------------------
template<int TDimension, typename TTransform>
void  ManifoldShapeModel<TDimension, TTransform>::
SetManifoldBuilder(Manifold * manifoldBuilder)
{
	m_Manifold = manifoldBuilder;
}


// ------------------------------------------------------------------------
template<int TDimension, typename TTransform>
void  
ManifoldShapeModel<TDimension, TTransform>::
Optimise(const LabelType * segmentationEstimate, 
			const ImageType * targetImage, 
			TransformPointer &transform,
			LabelPointer &output)
{

	typename LabelType::Pointer se = LabelType::New();
	se = const_cast<LabelType*>(segmentationEstimate);
	typename ImageType::Pointer re = ImageType::New();
	re = const_cast<ImageType*>(targetImage);


	// transform the segmentation estimation into the space of the model
	//
	//
	//
	TransformPointer imageToModelTransform = TransformType::New(); 
	imageToModelTransform->SetCenter(transform->GetCenter());
	transform->GetInverse(imageToModelTransform);

	LabelPointer estimationTransformed = LabelType::New();
	ApplyTransform(segmentationEstimate, imageToModelTransform, m_MeanLabel, estimationTransformed);


	typename ShapeEmbedderType::Pointer embedder = ShapeEmbedderType::New();
	embedder->SetInput(estimationTransformed);
	embedder->SetDistanceMeasure(m_DistanceMeasure);
	embedder->SetManifold(m_Manifold);
	embedder->SetDistanceMaps(m_DistanceMaps);
	embedder->SetKnn(m_Knn);
	embedder->Update();

	WeightListType weights = embedder->GetWeights();
	IndexListType indices = embedder->GetIndices();
	

	for(unsigned int i = 0; i < indices.size(); i++)
	{
		std::cout << indices[i] << " ";
	}

	std::cout << "" << std::endl;

	// parse out the parameters
	unsigned int transformParamNum = transform->GetNumberOfParameters();
	unsigned int weightsNum = weights.size();
	unsigned int totalParamNum = transformParamNum + weightsNum;
	std::vector<double> x(totalParamNum);

	for(unsigned int i = 0; i < weightsNum; i++)
	{
		x[i] = weights[i];
	}

	for(unsigned int i = weightsNum; i < totalParamNum; i++)
	{
		x[i] = transform->GetParameters()(i-weightsNum);
	}

	for(unsigned int i = 0; i < x.size(); i++)
	{
		std::cout << x[i] << std::endl;
	}

	// set up the optimiser
	nlopt::opt optimiser(nlopt::LN_BOBYQA, totalParamNum);
	optimiser.set_maxeval(m_OptimiserIterations);
	optimiser.set_xtol_rel(1e-4);
	optimiser.set_min_objective(func, (void*) this);

	// set up the bounds
	std::vector<double> lb(totalParamNum);
	std::vector<double> ub(totalParamNum);

	// set the weight bounds
	for(unsigned int i = 0; i < weightsNum; i++)
	{
		lb[i] = 0.0;
		ub[i] = 1.0;
	}


	double rotationTolerance = 0.5;
	double translationTolerance = 10.0;
	double scaleTolerance = 0.2;


	if(TDimension == 3)
	{
		lb[weightsNum] = x[weightsNum]-rotationTolerance;
		lb[weightsNum+1] = x[weightsNum+1]-rotationTolerance;
		lb[weightsNum+2] = x[weightsNum+2]-rotationTolerance;
		lb[weightsNum+3] = x[weightsNum+3]-translationTolerance;
		lb[weightsNum+4] = x[weightsNum+4]-translationTolerance;
		lb[weightsNum+5] = x[weightsNum+5]-translationTolerance;
		lb[weightsNum+6] = x[weightsNum+6]-scaleTolerance;

		ub[weightsNum] = x[weightsNum]+rotationTolerance;
		ub[weightsNum+1] = x[weightsNum+1]+rotationTolerance;
		ub[weightsNum+2] = x[weightsNum+2]+rotationTolerance;
		ub[weightsNum+3] = x[weightsNum+3]+translationTolerance;
		ub[weightsNum+4] = x[weightsNum+4]+translationTolerance;
		ub[weightsNum+5] = x[weightsNum+5]+translationTolerance;
		ub[weightsNum+6] = x[weightsNum+6]+scaleTolerance;
	}
	
	optimiser.set_lower_bounds(lb);
	optimiser.set_upper_bounds(ub);

	m_Target = targetImage;
	m_Indices = indices;
	m_Estimate = estimationTransformed;
	m_Transform = transform;


	double minf;
	try
	{
		nlopt::result res = optimiser.optimize(x, minf);
		std::cout << res << std::endl;
	}
	catch(std::exception e)
	{
		std::cout << e.what() << std::endl;
		return;
	}


	// parse out the optimised values
	WeightListType finalWeights;
	TransformParametersType finalTransformParameters;
	ParseX(x, finalWeights, finalTransformParameters);

	// generate a shape
	LabelPointer outLabel = LabelType::New();
	CreateLabel(finalWeights, indices, outLabel);

	transform->SetParameters(finalTransformParameters);
	ApplyFinalTransform(outLabel, transform, se, output);
}


// ------------------------------------------------------------------------
template<int TDimension, typename TTransform>
void  
ManifoldShapeModel<TDimension, TTransform>::
ParseX(const std::vector<double> &x, WeightListType &weights,
			TransformParametersType &params)
{
	unsigned int weightNum = m_Knn;
	for(unsigned int i = 0; i < weightNum; i++)
	{
		weights.push_back(x[i]);
	}

	unsigned int tparamNum = m_Transform->GetNumberOfParameters();
	params = TransformParametersType(tparamNum);
	for(unsigned int i = weightNum; i < tparamNum+weightNum; i++)
	{
		params(i-weightNum) = x[i];	
	}
}



// ------------------------------------------------------------------------
template<int TDimension, typename TTransform>
void  
ManifoldShapeModel<TDimension, TTransform>::
CreateShape(WeightListType &weights, IndexListType &indices,
			RealPointer &shape)
{
	// generate the estimation of the shape
	typedef ShapeGenerator<RealType, RealType> GeneratorType;
	typename GeneratorType::Pointer generator = GeneratorType::New();
	
	for(unsigned int i = 0; i < indices.size(); i++)
	{
		generator->SetInput(i, m_DistanceMaps[indices[i]]);
	}
	generator->SetWeightList(weights);
	generator->Update();

	shape = generator->GetOutput();
}


// ------------------------------------------------------------------------
template<int TDimension, typename TTransform>
void  
ManifoldShapeModel<TDimension, TTransform>::
CreateLabel(WeightListType &weights, IndexListType &indices,
			LabelPointer &shape)
{
	RealPointer real = RealType::New();
	CreateShape(weights, indices, real);

	ThresholdDistance(real, shape);

}

// ------------------------------------------------------------------------
template<int TDimension, typename TTransform>
double  
ManifoldShapeModel<TDimension, TTransform>::
func(const std::vector<double> &x, std::vector<double> &grad, void *f_data)
{
	// parse out the parameters
	Self * data = (Self*) f_data;
	grad.empty();

	
	//set up the wieght
	WeightListType weights;
	for(unsigned int i = 0; i < data->m_Knn; i++)
	{
		weights.push_back(x[i]);
	}

	unsigned int paramNum = data->m_Transform->GetNumberOfParameters();
	TransformParametersType transparams(paramNum);
	for(unsigned int i = data->m_Knn; i < paramNum + data->m_Knn; i++)
	{
		transparams(i-data->m_Knn) = x[i];
	}


	
	

	data->m_Transform->SetParameters(transparams);

	double cost = data->ComputeCost(data->m_Estimate, data->m_Target, data->m_Transform,
			data->m_Indices, weights);

	return cost;

	
}



// ------------------------------------------------------------------------
template<int TDimension, typename TTransform>
double  
ManifoldShapeModel<TDimension, TTransform>::
ComputeCost(const LabelType * estimation, 
			const ImageType * target,
			const TransformType  * transformToImage,
			const IndexListType &indices,
			const WeightListType &theta)
{
	// generate the estimation of the shape
	typedef ShapeGenerator<RealType, RealType> GeneratorType;
	typename GeneratorType::Pointer generator = GeneratorType::New();
	
	for(unsigned int i = 0; i < indices.size(); i++)
	{
		generator->SetInput(i, m_DistanceMaps[indices[i]]);
	}
	generator->SetWeightList(theta);
	generator->Update();
	
	


	// compute the divergence from the original shape
	double conditionOneCost = ComputeConditionOneCost(estimation, 
			generator->GetOutput());


	typename RealType::Pointer transformedShape = RealType::New();
	ApplyTransform(generator->GetOutput(), transformToImage, target, transformedShape);

	typename LabelType::Pointer transformedLabel = LabelType::New();
	ThresholdDistance(transformedShape, transformedLabel);
	

	typedef ShapeEntropyCalculator<TDimension> EntropyCalculatorType;
	typename EntropyCalculatorType::Pointer entropyCalculator = EntropyCalculatorType::New();
	entropyCalculator->SetInputImage(target);
	entropyCalculator->SetInputShape(transformedLabel);
	double entropy = entropyCalculator->GetValue();
	

	return conditionOneCost + (m_EntropyWeight * entropy);

}

// ------------------------------------------------------------------------
template<int TDimension, typename TTransform>
double  
ManifoldShapeModel<TDimension, TTransform>::
ComputeConditionOneCost(const LabelType * estimation,
		const RealType * current)
{
	// compute the distance map
	typename RealType::Pointer map = RealType::New();
	ComputeDistanceMap(estimation, map);

	// apply the transform to the estimation	
	m_DistanceMeasure->SetInput1(map);
	m_DistanceMeasure->SetInput2(current);

	return m_DistanceMeasure->GetValue();

}

// ------------------------------------------------------------------------
template<int TDimension, typename TTransform>
void
ManifoldShapeModel<TDimension, TTransform>::
ComputeDistanceMap(const LabelType * input, RealPointer &output)
{
	// compute the distance map
	typedef itk::SignedMaurerDistanceMapImageFilter<LabelType, RealType> DistanceFilterType;
	typename DistanceFilterType::Pointer distanceFilter = DistanceFilterType::New();
	distanceFilter->SetInput(input);
	distanceFilter->SetUseImageSpacing(true);
	distanceFilter->SetSquaredDistance(false);
	distanceFilter->SetBackgroundValue(0);
	distanceFilter->Update();

	output = distanceFilter->GetOutput();

}


// ------------------------------------------------------------------------
template<int TDimension, typename TTransform>
void  
ManifoldShapeModel<TDimension, TTransform>::
ThresholdDistance(const RealType *input, LabelPointer &output)
{
	typedef itk::BinaryThresholdImageFilter<RealType, LabelType> ThresholderType;
	typename ThresholderType::Pointer thresholder = ThresholderType::New();
	thresholder->SetInput(input);
	thresholder->SetUpperThreshold(0.0);
	thresholder->SetOutsideValue(0.0);
	thresholder->SetInsideValue(1.0);

	thresholder->Update();
	output = thresholder->GetOutput();


}

// ------------------------------------------------------------------------
template<int TDimension, typename TTransform>
void  
ManifoldShapeModel<TDimension, TTransform>::
SetDistanceMeasure(const DistanceMeasureType *measure)
{
	m_DistanceMeasure = const_cast<DistanceMeasureType*>(measure);
}


// ------------------------------------------------------------------------
template<int TDimension, typename TTransform>
void  
ManifoldShapeModel<TDimension, TTransform>::
ApplyTransform(const RealType * input,
			const TransformType * transform,
			const BaseType * ref,
			RealPointer &output)
{
	typedef itk::ResampleImageFilter<RealType, RealType> ResamplerType;
	typedef itk::NearestNeighborInterpolateImageFunction<RealType, double> InterpolatorType;

	typename ResamplerType::Pointer resampler = ResamplerType::New();
	typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

	resampler->SetTransform( transform );
	resampler->SetInterpolator( interpolator );
	resampler->SetInput( input );
	resampler->SetOutputSpacing( ref->GetSpacing() );
	resampler->SetOutputDirection( ref->GetDirection() );
	resampler->SetOutputOrigin( ref->GetOrigin() );
	resampler->SetSize( ref->GetLargestPossibleRegion().GetSize() );
	resampler->SetDefaultPixelValue(std::numeric_limits<double>::max());
	resampler->Update();

	output = resampler->GetOutput();
}

// ------------------------------------------------------------------------
template<int TDimension, typename TTransform>
void  
ManifoldShapeModel<TDimension, TTransform>::
ApplyFinalTransform(const LabelType * input,
			const TransformType * transform,
			const BaseType * ref,
			LabelPointer &output)
{
	typedef itk::ResampleImageFilter<LabelType, LabelType> ResamplerType;
	typedef itk::NearestNeighborInterpolateImageFunction<LabelType, double> InterpolatorType;

	typename ResamplerType::Pointer resampler = ResamplerType::New();
	typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

	resampler->SetTransform( transform );
	resampler->SetInterpolator( interpolator );
	resampler->SetInput( input );
	resampler->SetOutputSpacing( input->GetSpacing() );
	resampler->SetOutputDirection( ref->GetDirection() );
	resampler->SetOutputOrigin( ref->GetOrigin() );
	resampler->SetSize( input->GetLargestPossibleRegion().GetSize() );
	resampler->SetDefaultPixelValue(0);
	resampler->Update();

	output = resampler->GetOutput();
}


// ------------------------------------------------------------------------
template<int TDimension, typename TTransform>
void  
ManifoldShapeModel<TDimension, TTransform>::
ApplyTransform(const LabelType * input,
			const TransformType * transform,
			const BaseType * ref,
			LabelPointer &output)
{
	typedef itk::ResampleImageFilter<LabelType, LabelType> ResamplerType;
	typedef itk::NearestNeighborInterpolateImageFunction<LabelType, double> InterpolatorType;

	typename ResamplerType::Pointer resampler = ResamplerType::New();
	typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

	resampler->SetTransform( transform );
	resampler->SetInterpolator( interpolator );
	resampler->SetInput( input );
	resampler->SetOutputSpacing( ref->GetSpacing() );
	resampler->SetOutputDirection( ref->GetDirection() );
	resampler->SetOutputOrigin( ref->GetOrigin() );
	resampler->SetSize( ref->GetLargestPossibleRegion().GetSize() );
	resampler->SetDefaultPixelValue(0);
	resampler->Update();

	output = resampler->GetOutput();
}

} /* manifold */ 
#endif
