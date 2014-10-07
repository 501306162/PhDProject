#ifndef INITIAL_TRANSFORM_EXTRACTOR_H
#define INITIAL_TRANSFORM_EXTRACTOR_H

#include <itkSimilarity3DTransform.h>
#include <gdcmScanner.h>
#include <gdcmDirectory.h>
#include <gdcmTag.h>

#include <itkImage.h>
#include <itkAffineTransform.h>

#include "common.h"

class InitialTransformExtractor
{
public:

	typedef gdcm::Directory::FilenamesType FilenamesType;
	typedef itk::Image<unsigned short, 3> ImageType;
	typedef itk::AffineTransform<double, 3> TransformType;
	typedef TransformType::TranslationType TranslationType;

	InitialTransformExtractor();

	void SetOptions(const OptionsData &options);
	void SetDicomDir(const std::string &dir);
	void Compute();

	ImageType::PointType GetTranslation() { return m_Translation; }
	ImageType::DirectionType GetRotation() { return m_Rotation; }
	ImageType::Pointer GetReference() { return m_Reference; }
	

private:
	
	std::string GetReferenceImageFilename();
	void BuildMapping(gdcm::Scanner::MappingType &mapping);
	void GetFilenames(FilenamesType &filenames);
	TranslationType ComputeTranslation(const ImageType::Pointer &image, const OptionsData &options);

	ImageType::Pointer m_Reference;
	
	ImageType::PointType m_Translation;
	ImageType::DirectionType m_Rotation;
	std::string m_DicomDir;
	OptionsData m_Options;
	gdcm::Scanner::MappingType m_Mapping;
	FilenamesType m_Filenames;

	gdcm::Directory m_Dir;
	gdcm::Scanner m_Scanner;

	static gdcm::Tag seriesDescription;
	static gdcm::Tag seriesNumber;
	static gdcm::Tag instanceNumber;
	static gdcm::Tag sliceThickness;
	static gdcm::Tag triggerTime;
	static gdcm::Tag numberOfImages;
	static gdcm::Tag slicePosition;
	static gdcm::Tag studyId;
	static gdcm::Tag imagePosition;
	static gdcm::Tag imageOrientation;
	static gdcm::Tag sliceLocation;
	static gdcm::Tag rows;
	static gdcm::Tag cols;
	static gdcm::Tag pixelSpacing;


};



#endif
