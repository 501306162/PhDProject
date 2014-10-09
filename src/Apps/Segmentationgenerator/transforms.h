#ifndef TRANSFORMS_H
#define TRANSFORMS_H

#include "common.h"

/**
 * Function to get the transform from the reference space to the 
 * patient space
 */
void transformPointToPatient(const ImageType::Pointer &reference,
		const SeriesTransform &series, const OptionsData &options,
		VectorType &input, VectorType &output);


#endif
