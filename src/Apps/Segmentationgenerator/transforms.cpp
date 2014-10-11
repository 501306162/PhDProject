#include "transforms.h"

#include <itkLinearInterpolateImageFunction.h>
#include <itkImageRegionIterator.h>

#include <vnl/vnl_inverse.h>


// ------------------------------------------------------------------------
void createLabelImage(const SeriesTransform &series,
	   const ImageType::Pointer &reference,
	   const LevelSetType::Pointer &levelSet,
	   const std::vector<int> &roiOffset,
	   ImageType::Pointer &label)
{
	// initialise the label image
	createOutput(series.images[0], label);


	// create the level set interpolator
	typedef itk::LinearInterpolateImageFunction<LevelSetType> InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	interpolator->SetInputImage(levelSet);


	// iterate through the output 
	itk::ImageRegionIterator<ImageType> it(label, label->GetLargestPossibleRegion());
	it.GoToBegin();

	while(!it.IsAtEnd())
	{
		ImageType::IndexType index = it.GetIndex();
		ImageType::PointType p1, p2;

		// get the point in the level set space
		label->TransformIndexToPhysicalPoint(index, p1);
		transformPointToPatient(reference, series, roiOffset, p1, p2);

		// interpolate the level set value
		if(interpolator->IsInsideBuffer(p2))
		{
			float val = interpolator->Evaluate(p2);
			if(val >= 0) it.Set(1);
		}

		++it;
	}
}

// ------------------------------------------------------------------------
void createOutput(const ImageType::Pointer &input, ImageType::Pointer &output)
{
	output->SetDirection(input->GetDirection());
	output->SetSpacing(input->GetSpacing());
	output->SetRegions(input->GetLargestPossibleRegion());
	output->SetOrigin(input->GetOrigin());
	output->Allocate();
	output->FillBuffer(0);
}



// ------------------------------------------------------------------------
void transformPointToPatient(const ImageType::Pointer &reference,
		const SeriesTransform &series, const std::vector<int> &roiOffset,
		const ImageType::PointType &input, ImageType::PointType &output)
{
	// rotation is easy
	MatrixType rotTemp = reference->GetDirection().GetVnlMatrix();
	MatrixType rotation;
	rotation.set_size(3,3);
	for(unsigned int i = 0; i < 3; i++)
	{
		rotation(i,0) = rotTemp(i,1);
		rotation(i,1) = rotTemp(i,0);
		rotation(i,2) = rotTemp(i,2);
	}

	VectorType refOrig = reference->GetOrigin().GetVnlVector();
	VectorType offset  = vnl_inverse(rotation) * refOrig;

	vnl_vector<double> vnl_output;

	// translate the point
	vnl_output = vnl_inverse(rotation) * input.GetVnlVector();
	vnl_output = vnl_output - offset;
	vnl_output /= reference->GetSpacing()[0];
	for(unsigned int i = 0; i < 3; i++)
	{
		vnl_output[i] -= roiOffset[i];
		vnl_output[i] -= series.translation[i];
	}

	output[0] = vnl_output[0];
	output[1] = vnl_output[1];
	output[2] = vnl_output[2];
}

