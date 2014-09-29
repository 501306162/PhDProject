#include <iostream>

#include <CommonDefinitions.h>
#include <DoGKeyPointExtractor.h>
#include <KeyPointOrientator.h>

#include <vtkRegularPolygonSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkAppendPolyData.h>

#include <itkDiscreteGaussianImageFilter.h>
#include <itkGradientRecursiveGaussianImageFilter.h>
#include <itkNrrdImageIO.h>
#include <itkImageFileWriter.h>
#include <ImagePatchExtractor.h>

#include <itkCenteredRigid2DTransform.h>

#include <vtkLine.h>
#include <HistogramOfGradientsFeatureExtractor.h>
#include <vtkLineSource.h>

int main(int argc, char *argv[])
{

	utils::ImageSlice::Pointer image = utils::ImageSliceIO::Read(argv[1]);

	utils::ImageSliceIO::Write("image.nrrd", image);

	typedef filter::DoGKeyPointExtractor<utils::ImageSlice> ExtractorType;
	ExtractorType::Pointer extractor = ExtractorType::New();
	extractor->SetInput(image);
	extractor->SetKeypointThreshold(atof(argv[3]));
	extractor->SetSplitsPerOctave(atoi(argv[2]));
	extractor->SetStartingSigma(1.6);
	extractor->Update();




	// draw the output
	ExtractorType::OutputType output = extractor->GetOutput();




	typedef filter::KeyPointOrientator<utils::ImageSlice> Orientator;
	Orientator::Pointer orientator  = Orientator::New();
	orientator->SetInput(image);
	orientator->SetKeyPoints(output);
	orientator->SetHistogramBins(32);
	orientator->SetSigmaScale(2);
	orientator->SetSampleRadius(5);
	orientator->Update();


	Orientator::OrientatedKeyPointMap fout = orientator->GetOutput();




	

	// mean histogram over all the features
	typedef Orientator::GradientType GradientType;
	typedef filter::HistogramOfGradeintsFeatureExtractor<GradientType> FeatureBuilderType;
	typedef FeatureBuilderType::HistogramType FeatureType;
	FeatureType meanFeature;
	std::vector<FeatureType> featureList;
	

	unsigned int fcount = 0;

	Orientator::OrientatedKeyPointMap::iterator it = fout.begin();
	while(it != fout.end())
	{

		vtkSmartPointer<vtkPolyData> poutput =
			vtkSmartPointer<vtkPolyData>::New();

		// create the smoothed image
		typedef itk::DiscreteGaussianImageFilter<utils::ImageSlice, utils::RealSlice> Smoother;
		Smoother::Pointer smoother = Smoother::New();
		smoother->SetInput(image);
		smoother->SetVariance(it->first*it->first);
		smoother->SetUseImageSpacingOff();
		smoother->Update();



		for(unsigned int i = 0; i < it->second.size(); i++)
		{



			Orientator::OrientatedKeyPointList features = it->second;

			double p1 = features[i].location[0];
			double p2 = features[i].location[1];


			typedef filter::ImagePatchExtractor<utils::ImageSlice> ImageExtractor;
			ImageExtractor::Pointer imExtractor = ImageExtractor::New();
			imExtractor->SetInput(image);
			ImageExtractor::SizeType pSize;
			pSize.Fill(10);
			imExtractor->SetPatchSize(pSize);
			imExtractor->SetScale(features[i].scale);
			imExtractor->SetInputPoint(features[i].location);
			imExtractor->Update();

			utils::ImageSliceIO::Write("imp.nrrd", imExtractor->GetOutput());




			// draw the line in as well
			double x = cos(features[i].angle);
			double y = sin(features[i].angle);

			double np1 = p1+(features[i].scale*x);
			double np2 = p2+(features[i].scale*y);

			double origin[3];
			origin[0] = p1;
			origin[1] = p2;
			origin[2] = 0.0;

			double line[3];
			line[0] = np1;
			line[1] = np2;
			line[2] = 0.0;

			vtkSmartPointer<vtkLineSource> lineSource =
				vtkSmartPointer<vtkLineSource>::New();
			lineSource->SetPoint1(origin);
			lineSource->SetPoint2(line);
			lineSource->Update();


			// Create a circle
			vtkSmartPointer<vtkRegularPolygonSource> polygonSource =
				vtkSmartPointer<vtkRegularPolygonSource>::New();

			//polygonSource->GeneratePolygonOff();
			polygonSource->SetNumberOfSides(50);
			polygonSource->SetRadius(features[i].scale);
			polygonSource->SetCenter(p1,p2,0);
			polygonSource->GeneratePolygonOff();
			polygonSource->Update();	


			vtkSmartPointer<vtkAppendPolyData> appender = 
				vtkSmartPointer<vtkAppendPolyData>::New();
			appender->AddInputData(polygonSource->GetOutput());
			appender->AddInputData(poutput);
			appender->AddInputData(lineSource->GetOutput());
			appender->Update();

			poutput = appender->GetOutput();


			typedef itk::CenteredRigid2DTransform<double> TransformType;
			TransformType::Pointer transform = TransformType::New();
			transform->SetCenter(features[i].location);
			transform->SetAngleInDegrees(360-features[i].degrees);


			std::cout << features[i].degrees << std::endl;



			typedef itk::GradientRecursiveGaussianImageFilter<utils::RealSlice, GradientType> GradientFilterType;
			GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
			gradientFilter->SetInput(smoother->GetOutput());
		


			typedef filter::ImagePatchExtractor<GradientType> PatchExtractorType;
			PatchExtractorType::Pointer patchExtractor = PatchExtractorType::New();
			patchExtractor->SetInput(gradientFilter->GetOutput());
			patchExtractor->SetTransform(transform);
			patchExtractor->SetScale(features[i].scale*3);
			
			PatchExtractorType::SizeType patchSize;
			patchSize.Fill(10);
			patchExtractor->SetPatchSize(patchSize);
			patchExtractor->SetInputPoint(features[i].location);
			patchExtractor->Update();



			FeatureBuilderType::Pointer builder = FeatureBuilderType::New();
			builder->SetInput(patchExtractor->GetOutput());
			builder->Update();

			typedef FeatureBuilderType::FeatureType HoGFeatureType;
			HoGFeatureType hist = builder->GetOutput();
			fcount++;
		}


		

		std::stringstream ss;
		ss << it->first;

		std::string polyName = "poly" + ss.str() + ".vtk";
		std::cout << polyName << std::endl;

		vtkSmartPointer<vtkPolyDataWriter> writer = 
			vtkSmartPointer<vtkPolyDataWriter>::New();
		writer->SetInputData(poutput);
		writer->SetFileName(polyName.c_str());
		writer->Write();
		writer->Update();

		/*

		typedef Orientator::GradientType GradientType;
		gradientFilter->S
*/


		
		utils::RealSliceIO::Write("sm" + ss.str() + ".nrrd", smoother->GetOutput());

		it++;
	}


	for(unsigned int a = 0; a < meanFeature.size(); a++)
	{
		for(unsigned int b = 0; b < meanFeature.front().size(); b++)
		{
			meanFeature[a][b] /= static_cast<double>(fcount);
		}
	}


	// rank order the 


	return 0;
}
