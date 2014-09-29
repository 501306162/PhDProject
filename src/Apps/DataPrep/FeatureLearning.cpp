#include "FeatureLearning.h"
#include "FeatureCommon.h"


#include <FeaturePointGrouper.h>
#include <FeatureClusterLearner.h>
#include <FeatureClusterCuller.h>
#include <StructureTensorKeyPointValidator.h>
#include <ClusterDistributionGenerator.h>


using namespace utils;
int main(int argc, char *argv[])
{
	std::string inputFilenamesFilename = argv[1];
	double keyPointIntensityThreshold = atof(argv[2]);
	double dogSplitsPerOctave = atof(argv[3]);
	double statingScale = atof(argv[4]);
	double eLocation = atof(argv[5]);
	double eScale = std::log(atof(argv[6]));
	double eOrientation = atof(argv[7]);
	double gammaValue = atof(argv[8]);

	std::string distanceMapFilenamesFilename = argv[9];
	double extractionDistanceThreshold = atof(argv[10]);




	// load up the set of aligned images
	FilenamesReader::FilenamesType inputFilenames = FilenamesReader::Read(inputFilenamesFilename);
	FilenamesReader::FilenamesType distanceMapFilenames   = FilenamesReader::Read(distanceMapFilenamesFilename);
	ImageVolumeList images;
	RealVolumeList distanceMaps;
	for(unsigned int i = 0; i < inputFilenames.size(); i++)
	{
		ImageVolume::Pointer image = ImageVolumeIO::Read(inputFilenames[i]);
		images.push_back(image);
		RealVolume::Pointer distMap = RealVolumeIO::Read(distanceMapFilenames[i]);
		distanceMaps.push_back(distMap);
	}


	unsigned int sliceToTest = 7;

	// for each slice we want to learn the features
	const unsigned int sliceNum = images.front()->GetLargestPossibleRegion().GetSize()[2];
	for(unsigned int slice = sliceToTest; slice < sliceNum; slice++)
	{

		// get the set of slices that have some image data in them
		ImageSliceList validImages;
		RealSliceList validDistanceMaps;
		
		for(unsigned int im = 0; im < images.size(); im++)
		{
			ImageSlice::Pointer extractedSlice = ImageSlice::New();
			RealSlice::Pointer distanceSlice = RealSlice::New();
			ExtractSlice<ImageVolume, ImageSlice>(images[im], slice, extractedSlice);
			ExtractSlice<RealVolume, RealSlice>(distanceMaps[im], slice, distanceSlice);

			if(ImageContainsData(extractedSlice))
			{
				validDistanceMaps.push_back(distanceSlice);
				validImages.push_back(extractedSlice);
			}
		}

		/*
		if(validImages.size() < 3)
			continue;
		*/

		std::cout << "Slice Num: " << slice << " Image Number: " << validImages.size() << std::endl;



		typedef itk::Vector<double, 2> VectorType;
		typedef itk::Image<VectorType, 2> GradientType;
		typedef filter::HistogramOfGradeintsFeatureExtractor<GradientType> FeatureBuilderType;
		typedef FeatureBuilderType::FeatureType HoGFeatureType;
		std::vector<HoGFeatureType> allFeatures;
		std::vector<HoGFeatureType> allFeatures1;
		std::vector<HoGFeatureType> allFeatures2;

		

		unsigned int featureCount = 0;
		for(unsigned int im = 0; im < validImages.size(); im++)
		{
			ImageSlice::Pointer extractedSlice = validImages[im];

			// first we extract all of the keypoints points
			typedef filter::DoGKeyPointExtractor<utils::ImageSlice> ExtractorType;
			ExtractorType::Pointer extractor = ExtractorType::New();
			extractor->SetInput(extractedSlice);
			extractor->SetKeypointThreshold(keyPointIntensityThreshold);
			extractor->SetSplitsPerOctave(dogSplitsPerOctave);
			extractor->SetStartingSigma(statingScale);
			extractor->SetDistanceMap(validDistanceMaps[im]);
			extractor->SetDistanceThreshold(extractionDistanceThreshold);
			extractor->Update();

			// orientate the feature points
			typedef filter::KeyPointOrientator<utils::ImageSlice> Orientator;
			Orientator::Pointer orientator  = Orientator::New();
			orientator->SetInput(extractedSlice);
			orientator->SetKeyPoints(extractor->GetOutput());
			orientator->SetHistogramBins(32);
			orientator->SetSigmaScale(2);
			orientator->SetSampleRadius(5);
			orientator->Update();

			Orientator::OrientatedKeyPointMap orientateKeyPoints = orientator->GetOutput();


			

			// now we go through the features and compute the HOG descriptors
			Orientator::OrientatedKeyPointMap::iterator keyPointIt = orientateKeyPoints.begin();
			std::cout << orientateKeyPoints.size() << std::endl;
			while(keyPointIt != orientateKeyPoints.end())
			{
				double sigma = keyPointIt->first;
				Orientator::OrientatedKeyPointList keyPoints = keyPointIt->second;

				// smooth the image to the sigma level
				typedef itk::DiscreteGaussianImageFilter<utils::ImageSlice, utils::RealSlice> Smoother;
				Smoother::Pointer smoother = Smoother::New();
				smoother->SetInput(extractedSlice);
				smoother->SetVariance(sigma*sigma);
				smoother->SetUseImageSpacingOn();

				typedef itk::GradientRecursiveGaussianImageFilter<RealSlice, GradientType> GradientFilterType;
				GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
				gradientFilter->SetInput(smoother->GetOutput());
				//gradientFilter->SetSigma(sigma);
				gradientFilter->Update();



		
				


				std::cout << "Doing Sigma " << sigma << " Key Point Number: " << keyPoints.size() << std::endl;



				for(unsigned int fnum = 0; fnum < keyPoints.size(); fnum++)
				{
					Orientator::OrientatedKeyPoint keyPoint = keyPoints[fnum];

					// build the tranform
					typedef itk::CenteredRigid2DTransform<double> TransformType;
					TransformType::Pointer transform = TransformType::New();
					transform->SetCenter(keyPoint.location);
					transform->SetAngleInDegrees(360-keyPoint.degrees);

					// extract the patch from the gradient image
					typedef filter::ImagePatchExtractor<GradientType> PatchExtractorType;
					PatchExtractorType::Pointer patchExtractor = PatchExtractorType::New();
					patchExtractor->SetInput(gradientFilter->GetOutput());
					patchExtractor->SetTransform(transform);
					patchExtractor->SetScale(keyPoint.scale*2);

					PatchExtractorType::SizeType patchSize;
					patchSize.Fill(10);
					patchExtractor->SetPatchSize(patchSize);
					patchExtractor->SetInputPoint(keyPoint.location);
					patchExtractor->Update();


					/*
					// validate the keypoint
					typedef filter::StructureTensorKeyPointValidator<utils::ImageSlice> ValidatorType;
					ValidatorType::Pointer validator = ValidatorType::New();
					validator->SetInput(extractedSlice);	
					validator->SetKeyPoint(keyPoint);
					validator->SetRatio(validatorBeta);
					validator->Update();


					bool valid = validator->IsValid();
					*/


					// create the descriptor
					FeatureBuilderType::Pointer builder = FeatureBuilderType::New();
					builder->SetInput(patchExtractor->GetOutput());
					builder->SetOrientationBins(8);
					builder->SetKeyPoint(keyPoint);
					builder->Update();


					// add the feature to the list
					FeatureBuilderType::FeatureType feature = builder->GetOutput();
					feature.featureId = featureCount;
					allFeatures.push_back(feature);

					featureCount++;

				}



				
				keyPointIt++;
			}
		}

		




		std::cout << "Computing Distance Matrix" << std::endl;
		// compute the distance matrix
		typedef utils::DoubleMatrixType MatrixType;
		MatrixType distanceMatrix = MatrixType::Zero(allFeatures.size(), allFeatures.size());
		ComputeDifferenceMatrix(allFeatures, distanceMatrix);


		std::cout << "Grouping Features" << std::endl;
		// now we group the features by thier geometry
		typedef filter::FeaturePointGrouper<2> GrouperType;
		GrouperType::Pointer grouper = GrouperType::New();
		grouper->SetInput(allFeatures);
		grouper->SetAngleThreshold(eOrientation);
		grouper->SetLocationThreshold(eLocation);
		grouper->SetScaleThreshold(eScale);
		grouper->Update();





		std::cout << "Creating Clusters" << std::endl;
		GrouperType::FeatureGroupList clusters = grouper->GetOutput();
		std::sort(clusters.begin(), clusters.end());

		GrouperType::FeatureGroupList newClusters;

		for(unsigned int i = 0; i < clusters.size(); i++)
		{
			typedef filter::FeatureClusterLearner<2> ClusterLearnerType;
			ClusterLearnerType::Pointer learner = ClusterLearnerType::New();
			learner->SetInput(clusters[i]);
			learner->SetFeatures(allFeatures);
			learner->SetDistanceMatrix(distanceMatrix);
			learner->SetGamma(gammaValue);
			learner->Update();
			
			ClusterLearnerType::ClusterType newCluster = learner->GetOutput();
			newClusters.push_back(newCluster);



		}

		std::cout << "Culling Clusters" << std::endl;
		
		typedef filter::FeatureClusterCuller<2> Culler;
		Culler::Pointer culler = Culler::New();
		culler->SetInput(newClusters);
		culler->Update();


		Culler::ClusterList culledClusters = culler->GetOutput();
		std::sort(culledClusters.begin(), culledClusters.end());
		for(unsigned int i = 0; i < culledClusters.size(); i++)
		{
			typedef filter::ClusterDistributionGenerator<2> DistributionGenerator;
			DistributionGenerator::Pointer generator = DistributionGenerator::New();
			generator->SetInput(culledClusters[i]);
			generator->Update();

			exit(1);
		}



		/*


		ImageSlice::Pointer extractedSlice = ImageSlice::New();
		ExtractSlice<ImageVolume, ImageSlice>(images[0], sliceToTest, extractedSlice);

		std::vector<std::pair<int, ImageSlice::PointType> > testOut;
		for(unsigned int i = 0; i < culledClusters.size(); i++)
		{
			for(unsigned int j = 0; j < culledClusters[i].clusterMembers.size(); j++)
			{
			
				std::pair<int, ImageSlice::PointType> p(i, culledClusters[i].clusterMembers[j].keyPoint.location);
				testOut.push_back(p);				
			}			
		}




		utils::DoubleMatrixType pOut = utils::DoubleMatrixType::Zero(testOut.size(), 2);
		utils::IntMatrixType iOut = utils::IntMatrixType::Zero(testOut.size(),1);
		for(unsigned int i = 0; i < testOut.size(); i++)
		{
			itk::ContinuousIndex<double, 2> contIndex;
			extractedSlice->TransformPhysicalPointToContinuousIndex(testOut[i].second, contIndex);

			pOut(i,0) = contIndex[0];			
			pOut(i,1) = contIndex[1];

			iOut(i,0) = testOut[i].first;			
		}


	


		utils::MatrixDataSet::Pointer dout = utils::MatrixDataSet::New();
		dout->AddData("locations", pOut);
		dout->AddData("index", iOut);
		utils::MatrixWriter::Pointer writer = utils::MatrixWriter::New();
		writer->SetInput(dout);
		writer->SetFilename("data.hdf5");
		writer->Update();



		exit(1);
		*/

		


		/*
		// compute the affinity matrix between all of the features
		unsigned int numFeatures = allFeatures.size();
		double sum = 0.0;
		int count = 0;
		int max = 0;
		int maxId = 0;
		std::vector<int> counts;

		std::vector<Cluster> allGroups;

		// groupd all the features that have a similar location / scale / orientation
		for(unsigned int i = 0; i < numFeatures; i++)
		{
			Cluster cluster;
			cluster.featureIndex = i;
			cluster.feature = allFeatures[i];
			cluster.e = 0.0;
			cluster.members.push_back(allFeatures[i]);
			
			for(unsigned int j = 0; j < numFeatures; j++)
			{
				if(i == j) continue;
				if(AreSimilar(allFeatures[i], allFeatures[j],
							eLocation, eOrientation, eScale))
				{
					cluster.members.push_back(allFeatures[j]);
				}
			}
			
			counts.push_back(cluster.members.size());
			if(cluster.members.size() > max)
			{
				max = cluster.members.size();
				maxId = i;
			}

			allGroups.push_back(cluster);
			sum += cluster.members.size();
			std::sort(cluster.members.begin(), cluster.members.end(), member_sort);
			count++;
		}



		std::sort(counts.begin(), counts.end());
		for(unsigned int i = 0; i < counts.size(); i++)
		{
			std::cout << counts[i] << std::endl;
		}


		// compute the difference matrix
		utils::DoubleMatrixType diffMat;
		ComputeDifferenceMatrix(allFeatures, diffMat);

		// loop through the groups to form the clusters
		std::vector<Cluster> allClusters;
		for(unsigned int i = 0; i < allGroups.size(); i++)
		{
			Cluster cluster;
			CreateCluster(allGroups[i], allFeatures, diffMat, cluster);
			allClusters.push_back(cluster);
		}

		// remove duplicates
		std::vector<int> toRemove;
		for(unsigned int i = 0; i < allClusters.size(); i++)
		{
			bool duplicate = false;
			for(unsigned int j = i; j < allClusters.size(); j++)
			{
				if(i == j) continue;
				if(allClusters[i].members.size() != allClusters[j].members.size())
					continue;

				bool sameMembers = true;
				for(unsigned int k = 0; k < allClusters[i].members.size(); k++)
				{
					if(allClusters[i].members[k].index != allClusters[j].members[k].index)
					{
						sameMembers = false;				
						break;
					}
				}

				if(sameMembers)
				{
					duplicate = true;
				}

				if(duplicate) break;
			}
			if(duplicate) toRemove.push_back(i);

		}

		
		std::cout << allClusters.size() << std::endl;
		for(unsigned int i = 0; i < toRemove.size(); i++)
		{
	//		allClusters.erase(allGroups.begin()+(toRemove[i]-i));
		}
		std::cout << allClusters.size() << std::endl;

		// trim the clusters
		std::vector<Cluster> trimmedClusters;
		TrimClusters(allClusters, trimmedClusters);
		std::cout << trimmedClusters.size() << std::endl;

		std::vector<std::pair<int, ImageSlice::PointType> > testOut;


		for(unsigned int i = 0; i < trimmedClusters.size(); i++)
		{
			for(unsigned int j = 0; j < trimmedClusters[i].members.size(); j++)
			{
				std::pair<int, ImageSlice::PointType> p(i, trimmedClusters[i].members[j].location);
				testOut.push_back(p);				
			}			
			std::cout << trimmedClusters[i].members.size() << std::endl;
		}

		ImageSlice::Pointer extractedSlice = ImageSlice::New();
		ExtractSlice<ImageVolume, ImageSlice>(images[0], sliceToTest, extractedSlice);


		utils::DoubleMatrixType pOut = utils::DoubleMatrixType::Zero(testOut.size(), 2);
		utils::IntMatrixType iOut = utils::IntMatrixType::Zero(testOut.size(),1);
		for(unsigned int i = 0; i < testOut.size(); i++)
		{
			itk::ContinuousIndex<double, 2> contIndex;
			extractedSlice->TransformPhysicalPointToContinuousIndex(testOut[i].second, contIndex);

			pOut(i,0) = contIndex[0];			
			pOut(i,1) = contIndex[1];
			
			// compute the image indexes


			
	
			iOut(i,0) = testOut[i].first;			
		}


	


		utils::MatrixDataSet::Pointer dout = utils::MatrixDataSet::New();
		dout->AddData("locations", pOut);
		dout->AddData("index", iOut);
		utils::MatrixWriter::Pointer writer = utils::MatrixWriter::New();
		writer->SetInput(dout);
		writer->SetFilename("data.hdf5");
		writer->Update();
		exit(1);
		*/
	}

	return 0;
}

