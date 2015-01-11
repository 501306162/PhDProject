#include "CMRFileExtractor.h"

#include <QDir>
#include <QFileInfo>
#include <itkOrientImageFilter.h>


#include <itkImageFileReader.h>
#include <itkNrrdImageIO.h>
#include <itkPermuteAxesImageFilter.h>
#include <itkFlipImageFilter.h>

#include <itkExtractImageFilter.h>

namespace vt
{
// ------------------------------------------------------------------------
void CMRFileExtractor::Extract()
{
	// parse all the files in the folder
	QDir dir(QString::fromStdString(m_FolderName));
	QStringList filters;
	filters << "*.nrrd";
	QStringList files = dir.entryList(filters);

	// set up the filters 
	QStringList c2Filters;
	c2Filters << "_2C";
	c2Filters << "_2c";
	
	QStringList c3Filters;
 	c3Filters << "_3C";
 	c3Filters << "_3c";

	QStringList c4Filters;
	c4Filters << "_4C";
	c4Filters << "_4c";

	QStringList stackFilters;
	stackFilters << "stack";

	QStringList r2cFilters;
	r2cFilters << "_R2C";
	r2cFilters << "_RT_2C";
	r2cFilters << "_r2c";
	r2cFilters << "_RV2C";

	QStringList r3cFilters;
	r3cFilters << "_R3C";
	r3cFilters << "_RT_3C";
	r3cFilters << "_r3c";
	r3cFilters << "_RV3C";

	for (int i = 0; i < files.size(); ++i)
	{

		if(Matches(files[i], r2cFilters))
		{
			if(m_Debug) std::cout << "R2C : " << files[i].toStdString() << std::endl;
			LoadImage(dir.absoluteFilePath(files[i]), m_R2CImage);
			continue;
		}			

		if(Matches(files[i], r3cFilters))
		{
			if(m_Debug) std::cout << "R3C : " << files[i].toStdString() << std::endl;
			LoadImage(dir.absoluteFilePath(files[i]), m_R3CImage);
			continue;
		}			


		if(Matches(files[i], c2Filters))
		{
			if(m_Debug) std::cout << "2C : " << files[i].toStdString() << std::endl;
			LoadImage(dir.absoluteFilePath(files[i]), m_2CImage);
			continue;
		}			

		if(Matches(files[i], c3Filters))
		{
			if(m_Debug) std::cout << "3C : " << files[i].toStdString() << std::endl;
			LoadImage(dir.absoluteFilePath(files[i]), m_3CImage);
			continue;
		}			

		if(Matches(files[i], c4Filters))
		{
			if(m_Debug) std::cout << "4C : " << files[i].toStdString() << std::endl;
			LoadImage(dir.absoluteFilePath(files[i]), m_4CImage);
			continue;
		}			

		if(Matches(files[i], stackFilters))
		{
			if(m_Debug) std::cout << "Stack : " << files[i].toStdString() << std::endl;
			LoadImage(dir.absoluteFilePath(files[i]), m_Stack);
			continue;
		}			

	
		//std::cout << "No Match : " << files[i].toStdString() << std::endl;
	}


	if(!m_2CImage)
	{
		std::cout << "No 2C Image Found" << std::endl;
	}


	if(!m_3CImage)
	{
		std::cout << "No 3C Image Found" << std::endl;
	}


	//if(!m_4CImage)
	//{
		//std::cout << "No 4C Image Found" << std::endl;
	//}


	//if(!m_R2CImage)
	//{
		//std::cout << "No R2C Image Found" << std::endl;
	//}


	//if(!m_R3CImage)
	//{
		//std::cout << "No R3C Image Found" << std::endl;
	//}

	if(!m_Stack)
	{
		std::cout << "No Stack Image Found" << std::endl;
	}
}

// ------------------------------------------------------------------------
bool CMRFileExtractor::Matches(const QString &filename, const QStringList &filters)
{
	for (int i = 0; i < filters.size(); ++i)
	{
		if(filename.contains(filters[i]))
			return true;
	}

	return false;

}


// ------------------------------------------------------------------------
void CMRFileExtractor::FlipImage(const ImageType::Pointer &input, ImageType::Pointer &output)
{
	typedef itk::OrientImageFilter<ImageType, ImageType> OreinterType;
	OreinterType::Pointer orienter = OreinterType::New();
	orienter->SetInput(input);
	orienter->SetGivenCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSA);
	orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSA);
	orienter->Update();

	output = orienter->GetOutput();

	/*
	ImageSeriesType::DirectionType dir = input->GetDirection();
	typedef itk::PermuteAxesImageFilter<ImageSeriesType> PermType;
	PermType::Pointer permer = PermType::New();
	permer->SetInput(input);

	PermType::PermuteOrderArrayType axes;
	axes[0] = 1;
	axes[1] = 0;
	axes[2] = 2;
	axes[3] = 3;

	permer->SetOrder(axes);
	permer->Update();

	output = permer->GetOutput();


	typedef itk::FlipImageFilter<ImageSeriesType> FlipperType;
	FlipperType::Pointer flipper = FlipperType::New();
	FlipperType::FlipAxesArrayType flipAxes;
	flipAxes.Fill(false);

	*/


}


// ------------------------------------------------------------------------
void CMRFileExtractor::LoadImage(const QString &filename, ImageSeriesType::Pointer &image)
{
	typedef itk::ImageFileReader<ImageSeriesType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename.toStdString());
	reader->SetImageIO(itk::NrrdImageIO::New());
	reader->Update();

	image = reader->GetOutput();

}

// ------------------------------------------------------------------------
CMRFileExtractor::ImageType::Pointer 
CMRFileExtractor::ExtractImageTimeStep(const ImageSeriesType::Pointer &image, 
		const unsigned int index)
{
	typedef itk::ExtractImageFilter<ImageSeriesType, ImageType> ExtractorType;
	ExtractorType::Pointer extractor = ExtractorType::New();
	extractor->SetInput(image);

	ImageSeriesType::RegionType exRegion = image->GetLargestPossibleRegion();
	ImageSeriesType::SizeType exSize = exRegion.GetSize();
	ImageSeriesType::IndexType exIndex = exRegion.GetIndex();

	exSize[3] = 0;
	exIndex[3] = index;
	exRegion.SetSize(exSize);
	exRegion.SetIndex(exIndex);

	extractor->SetExtractionRegion(exRegion);
	extractor->SetDirectionCollapseToSubmatrix();
	extractor->Update();

	ImageType::Pointer output = ImageType::New();
	if(m_Flip)
	{
		FlipImage(extractor->GetOutput(), output);
		return output;
	}

	return extractor->GetOutput();

	

}


CMRFileExtractor::ImageType::Pointer  CMRFileExtractor::Get2CImage(const unsigned int index)
{
	return ExtractImageTimeStep(m_2CImage, index);
}

CMRFileExtractor::ImageType::Pointer CMRFileExtractor::Get3CImage(const unsigned int index)
{
	return ExtractImageTimeStep(m_3CImage, index);
}

CMRFileExtractor::ImageType::Pointer CMRFileExtractor::Get4CImage(const unsigned int index)
{
	return ExtractImageTimeStep(m_4CImage, index);
}

CMRFileExtractor::ImageType::Pointer CMRFileExtractor::GetR2CImage(const unsigned int index)
{
	return ExtractImageTimeStep(m_R2CImage, index);
}

CMRFileExtractor::ImageType::Pointer CMRFileExtractor::GetR3CImage(const unsigned int index)
{
	return ExtractImageTimeStep(m_R3CImage, index);
}

CMRFileExtractor::ImageType::Pointer CMRFileExtractor::GetStackImage(const unsigned int index)
{
	return ExtractImageTimeStep(m_Stack, index);
}

}
