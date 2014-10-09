#include "transforms.h"

#include <vnl/vnl_inverse.h>

// ------------------------------------------------------------------------
void transformPointToPatient(const ImageType::Pointer &reference,
		const SeriesTransform &series, const OptionsData &options,
		VectorType &input, VectorType &output)
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

	// translate the point
	output = vnl_inverse(rotation) * input;
	output = output - offset;
	output /= reference->GetSpacing()[0];
	for(unsigned int i = 0; i < 3; i++)
	{
		output[i] -= options.roiOffset[i];
		output[i] -= series.translation[i];
	}



	/*
	for(unsigned int i = 0; i < 3; i++)
	{
		output[i] = input[i];
		
		// first negate the registration offsets
		output[i] += series.translation[i];

		// now scale
		output[i] *= reference->GetSpacing()[0];

		// now offset
		output[i] += offset[i];
	}

	// now apply the rotation to the point
	output = rotation * output;
	*/


	
}
