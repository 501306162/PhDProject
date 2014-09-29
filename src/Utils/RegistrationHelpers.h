#ifndef REGISTRATION_HELPERS_H
#define REGISTRATION_HELPERS_H

#include "CommonDefinitions.h"

#include <itkSimilarity3DTransform.h>
#include <itkSimilarity2DTransform.h>

namespace utils
{

typedef itk::Similarity3DTransform<double> Rigid3DTransformType;
typedef itk::Similarity2DTransform<double> Rigid2DTransformType;

void RegisterLabelVolumeImages(
		LabelVolume::Pointer &fixed, 
		LabelVolume::Pointer &moving,
		Rigid3DTransformType::Pointer &transform);

void RegisterLabelSliceImages(
		LabelSlice::Pointer &fixed, 
		LabelSlice::Pointer &moving,
		Rigid2DTransformType::Pointer &transform);

void Apply2DRigidTransformToLabelSlice(
		LabelSlice::Pointer &input,
		LabelSlice::Pointer &reference,
		Rigid2DTransformType::Pointer &transform,
		LabelSlice::Pointer &output);


void Apply3DRigidTransformTOLabelVolume(
		LabelVolume::Pointer &input,
		LabelVolume::Pointer &reference,
		Rigid3DTransformType::Pointer &transform,
		LabelVolume::Pointer &output);


void Apply3DRigidTransformToImageVolume(
		ImageVolume::Pointer &input,
		ImageVolume::Pointer &reference,
		Rigid3DTransformType::Pointer &transform,
		ImageVolume::Pointer &output);


} /* utils */ 

#endif
