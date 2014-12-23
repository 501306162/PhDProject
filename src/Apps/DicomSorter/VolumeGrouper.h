#ifndef VOLUME_GROUPER_H
#define VOLUME_GROUPER_H


#include "common.h"
#include <itkImage.h>
#include <map>

class VolumeGrouper
{
public:

	typedef std::map<std::string, DicomSeriesList> NamedSeriesMap;
	typedef itk::Image<unsigned short, 3> ImageType;

	VolumeGrouper();

	void SetInput(const DicomSeriesList &series) { this->m_Series = series; }

	void Group();

private:

	void GroupByName(NamedSeriesMap &seriesMap);

	void CheckPositionAndOrientation(const DicomSeriesList &series, NamedSeriesMap &outputSeries);

	double ComputeImageDistance(ImageType::Pointer &im1, ImageType::Pointer &im2);

	DicomSeriesList m_Series;
	double m_DistanceTolerance;
};

#endif
