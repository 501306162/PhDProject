#include "VolumeGrouper.h"


#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkGDCMImageIO.h>

// ------------------------------------------------------------------------
VolumeGrouper::VolumeGrouper()
{
	m_DistanceTolerance = 0.5;

}


// ------------------------------------------------------------------------
void VolumeGrouper::Group()
{
	NamedSeriesMap seriesMap;
	GroupByName(seriesMap);

	// loop through the map and process groups that have more than
	// one series
	NamedSeriesMap::iterator mapIt = seriesMap.begin();
	while(mapIt != seriesMap.end())
	{
		if(mapIt->second.size() > 1)
		{
			NamedSeriesMap newSeries;
			CheckPositionAndOrientation(mapIt->second, newSeries);
		}


		++mapIt;
		
	}
}


// ------------------------------------------------------------------------
void VolumeGrouper::CheckPositionAndOrientation(const DicomSeriesList &series, 
		NamedSeriesMap &outputSeries)
{
	// first lets load all the images 
	typedef std::vector<ImageType::Pointer> ImageList;
	typedef itk::ImageFileReader<ImageType> ReaderType;

	ImageList images;

	for(unsigned int i = 0; i < series.size(); i++)
	{
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(series[i].images.front().filename);
		reader->SetImageIO(itk::GDCMImageIO::New());
		reader->Update();
		images.push_back(reader->GetOutput());		
	}

	// first check the origins are not equal within a tolerance
	for(unsigned int i = 0; i < images.size(); i++)
	{
		ImageType::Pointer im1 = images[i];
		for(unsigned int j = i+1; j < images.size(); j++)
		{
			ImageType::Pointer im2 = images[j];
			double distance = ComputeImageDistance(im1, im2);
			if(distance < m_DistanceTolerance)
			{
				std::cout << series[i].description << std::endl;
			}
			
		}
	}


}


// ------------------------------------------------------------------------
double VolumeGrouper::ComputeImageDistance(ImageType::Pointer &im1, 
		ImageType::Pointer &im2)
{
	ImageType::PointType p1 = im1->GetOrigin();
	ImageType::PointType p2 = im2->GetOrigin();

	return p1.EuclideanDistanceTo(p2);
}



// ------------------------------------------------------------------------
void VolumeGrouper::GroupByName(NamedSeriesMap &seriesMap)
{
	for(unsigned int i = 0; i < m_Series.size(); i++)
	{
		DicomSeries & series = m_Series[i];
		std::string name = series.description;

		seriesMap[name].push_back(series);
	}
}
