#include <iostream>

#include <PatchExtractor2.h>

#include <ValveIO.h>
#include <itkVector.h>
#include <itkImageRegionIterator.h>
#include <itkLinearInterpolateImageFunction.h>
#include <CommonDefinitions.h>
#include <FlipChecker.h>
#include <itkPNGImageIO.h>
#include <itkImageFileWriter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkNrrdImageIO.h>

#include <Directory.h>


using namespace vt;

typedef ValveSequence<3> SequenceType;
typedef ValveLine<3> ValveType;
typedef ValveType::ImageType ImageType;
typedef ValveType::PointType PointType;
typedef itk::Vector<double, 3> VectorType;


int main(int argc, char **argv)
{


	utils::Directory::FilenamesType filenames = utils::Directory::GetFiles(argv[1], ".txt");
	for(unsigned int k = 0; k < filenames.size(); k++)
	{
		// get the id 
		const std::string filename = filenames[k];
		std::string fname = utils::Directory::GetFileName(filename);
		unsigned int id = QString::fromStdString(fname).replace("d","").replace(".txt","").toInt();


		FlipChecker::Pointer checker = FlipChecker::New();
		bool flipImage = checker->FlipImage(argv[2], id);
		bool flipPoints = checker->FlipPoints(argv[2], id);

		flipPoints = false;
		flipImage = false;


		ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
		reader->SetFileName(filename);
		SequenceType::Pointer sequence = reader->GetOutput();

		ValveType::Pointer valve = sequence->GetValveLine(0);

		PointType p1 = valve->GetP1();
		PointType p2 = valve->GetP2();

		std::cout << id << " " << flipImage << " " << flipPoints << std::endl;

		if(flipPoints)
		{
			p1 = valve->GetP2();
			p2 = valve->GetP1();
		}

		VectorType xDir = (p2-p1);
		typedef PatchExtractor2<ImageType> PatchExtractorType;
		PatchExtractorType::Pointer extractor = PatchExtractorType::New();

		PatchExtractorType::SizeType size;
		size.Fill(200);
		size[2] = 1;

		PatchExtractorType::DistanceType distance;
		distance.Fill(200);
		distance[2] = 0.0;

		extractor->SetSize(size);
		extractor->SetDistance(distance);
		extractor->SetInput(valve->GetImage());
		extractor->SetCenter(p1);
		extractor->SetLine(xDir);
		extractor->Update();

		ImageType::Pointer output = extractor->GetOutput();
		
		bool check = false;

		ImageType::DirectionType tdir = valve->GetImage()->GetDirection();
		double sum = 0.0;
		for(unsigned int i = 0; i < 3; i++)
		{
			sum += tdir(i,i);			
		}


		if(sum < 0.0)
			check = true;


		//std::cout << sum << " " << flipImage  << std::endl;

		//if(flipImage != check)
		//{
			//std::cout << id << " " << flipImage << " " << check << std::endl;
		//}

	


		/*



		ImageType::Pointer image = valve->GetImage();
		typedef ImageType::DirectionType DirectionType;
		DirectionType imageDir = image->GetDirection();


		const unsigned int patchSize = 200;
		const unsigned int patchLength = 200;

		ImageType::SpacingType spacing;
		for(unsigned int i = 0; i < 2; i++)
		{
			spacing[i] = static_cast<double>(patchLength) / static_cast<double>(patchSize);
		}

		spacing[2] = image->GetSpacing()[2];


		
		xDir.Normalize();

		VectorType zDir;
		for(unsigned int i = 0; i < 3; i++)
		{
			zDir[i] = imageDir(i,2);
		}

		VectorType yDir = itk::CrossProduct(xDir, zDir);
		if(flipImage)
			yDir = -yDir;

		DirectionType newDir;
		for(unsigned int i = 0; i < 3; i++)
		{
			newDir(i,0) = xDir[i];
			newDir(i,1) = yDir[i];
			newDir(i,2) = zDir[i];
		}

		PointType np;
		for(unsigned int i = 0; i < 3; i++)
		{
			np[i] = p1[i] - (0.5*patchLength*yDir[i]);
			np[i] +=  -(0.5*patchLength*xDir[i]);
		}

		std::cout << newDir << std::endl;

		ImageType::Pointer output = ImageType::New();
		output->SetDirection(newDir);
		output->SetSpacing(spacing);
		output->SetOrigin(np);

		ImageType::RegionType region;
		ImageType::SizeType size;
		size.Fill(patchSize);
		size[2] = 1;
		ImageType::IndexType index;
		index.Fill(0);

		region.SetSize(size);
		region.SetIndex(index);

		output->SetRegions(region);
		output->Allocate();


		typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		interpolator->SetInputImage(image);

		itk::ImageRegionIterator<ImageType> it(output, output->GetLargestPossibleRegion());
		while(!it.IsAtEnd())
		{
			ImageType::IndexType index = it.GetIndex();
			PointType loc;
			output->TransformIndexToPhysicalPoint(index, loc);

			if(interpolator->IsInsideBuffer(loc))
				it.Set(interpolator->Evaluate(loc));
			else
				it.Set(0);

			++it;

		}
		*/

		typedef itk::Image<unsigned char, 3> OutputType;
		typedef itk::RescaleIntensityImageFilter<ImageType, OutputType> RescalerType;
		RescalerType::Pointer rescaler = RescalerType::New();
		rescaler->SetInput(output);
		rescaler->SetOutputMaximum(255);
		rescaler->SetOutputMinimum(0);

		typedef itk::ImageFileWriter<OutputType> WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetInput(rescaler->GetOutput());
		writer->SetImageIO(itk::PNGImageIO::New());

		std::stringstream ss;
		ss << id << ".png";
		writer->SetFileName(ss.str());
		writer->Update();


		
	}



}
