#ifndef TRANSFORMS_H
#define TRANSFORMS_H

#include "common.h"

/**
 * Function to get the transform from the reference space to the 
 * patient space
 */
void transformPointToPatient(const ImageType::Pointer &reference,
		const SeriesTransform &series, const std::vector<int> &roiOffset,
		const ImageType::PointType &input, ImageType::PointType &output);

/**
 * Funciton to build the label image from the level set and the transforms
 */ 
void createLabelImage(const SeriesTransform &series,
	   const ImageType::Pointer &reference,
	   const LevelSetType::Pointer &levelSet,
	   const std::vector<int> &roiOffset,
	   ImageType::Pointer &label);
		
/**
 * Function to build the output from the input
 */
void createOutput(const ImageType::Pointer &input, ImageType::Pointer &output);

#endif
