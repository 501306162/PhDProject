
#include <iostream>
#include <CommonDefinitions.h>
#include <FilenamesReader.h>

#include <Helpers.h>

#include <DoGKeyPointExtractor.h>
#include <KeyPointOrientator.h>
#include <itkStatisticsImageFilter.h>
#include <itkGradientRecursiveGaussianImageFilter.h>
#include <HistogramOfGradientsFeatureExtractor.h>
#include <itkCenteredRigid2DTransform.h>

#include <MatrixWriter.h>

#include <Eigen/Dense>


#include <vtkRegularPolygonSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkAppendPolyData.h>
#include <vtkLineSource.h>

#include <FeatureCommon.h>

typedef struct whole_feature_
{
	unsigned int imageNum;
	double scale;
	double angle;
	std::vector<std::vector<double> > histogram;
	utils::ImageSlice::PointType location;
	unsigned int index;
} WholeFeature;

typedef struct cluster_
{
	int featureIndex;
	WholeFeature feature;
	std::vector<WholeFeature> members;
	double e;
} Cluster;


bool member_sort(const WholeFeature &f1, const WholeFeature &f2)
{
	return f1.index < f2.index;
}




double FeatureDistance(const filter::HoGFeature<2> &f1, const filter::HoGFeature<2> &f2);
void GroupNonClusterFeatures(const Cluster &cluster, const std::vector<WholeFeature> &allFeatures, std::vector<WholeFeature> &output);
void CreateCluster(const Cluster &cluster, const std::vector<WholeFeature> &allFeatures,
	   const utils::DoubleMatrixType &diffMat,	Cluster &output);

// ------------------------------------------------------------------------
void ComputeDifferenceMatrix(const std::vector<filter::HoGFeature<2> > &input, utils::DoubleMatrixType &diffMat)
{
	diffMat = utils::DoubleMatrixType::Zero(input.size(), input.size());
	for(unsigned int i = 0; i < input.size(); i++)
	{
		for(unsigned int j = 0; j < input.size(); j++)
		{
			diffMat(i,j) = FeatureDistance(input[i], input[j]);			
		}
	}
}



// ------------------------------------------------------------------------
void GroupNonClusterFeatures(const Cluster &cluster, const std::vector<WholeFeature> &allFeatures, std::vector<WholeFeature> &output)
{
	// get the set of features that dont belong to the cluster
	for(unsigned int i = 0; i < allFeatures.size(); i++)
	{
		bool member = false;
		for(unsigned int j = 0; j < cluster.members.size(); j++)
		{
			if(i == cluster.members[j].index)
			{
				member = true;
				break;
			}
		}
		if(!member)
			output.push_back(allFeatures[i]);				
	}
}

// ------------------------------------------------------------------------
double ComputeGamma(const Cluster &cluster, const utils::DoubleMatrixType &mat, const double e)
{
	unsigned int featureNumber = mat.cols();
	unsigned int clusterNumber = cluster.featureIndex;	
	
	unsigned int fCount = 0;
	unsigned int bCount = 0;
	for(unsigned int i = 0; i < featureNumber; i++)
	{
		double diff = mat(clusterNumber,i);		


		if(diff > e) continue;

		bool member = false;
		for(unsigned int j = 0; j < cluster.members.size(); j++)
		{
			if(i == cluster.members[j].index)
			{
				member = true;
				break;
			}
		}

		if(member)
			fCount++;
		else
			bCount++;

		
	}
	
	
	double output = static_cast<double>(fCount) / static_cast<double>(bCount);



	return output;
}

// ------------------------------------------------------------------------
void CreateCluster(const Cluster &cluster, const std::vector<WholeFeature> &allFeatures,
	   const utils::DoubleMatrixType &diffMat,	Cluster &output)
{
	std::vector<WholeFeature> nonClusterFeatures;
	GroupNonClusterFeatures(cluster,allFeatures,nonClusterFeatures);

	double eInc = 1.0;

	double e = 0.0;
	double gammaThreshold = 0.05;
	double gamma = std::numeric_limits<double>::max();
	while(gamma >= gammaThreshold)
	{
		double newGamma = ComputeGamma(cluster, diffMat, e);

		if(newGamma >= gammaThreshold)
			gamma = newGamma;
		else
			break;
		
		e+=eInc;		
	}



	// create the new cluster
	output.e = e;
	output.feature = cluster.feature;
	output.featureIndex = cluster.featureIndex;
	

	// get the cluster members
	for(unsigned int i = 0; i < cluster.members.size(); i++)
	{
		if(diffMat(cluster.featureIndex, cluster.members[i].index) <= e)		
			output.members.push_back(cluster.members[i]);
	}

	std::sort(output.members.begin(), output.members.end(), member_sort);
}


// ------------------------------------------------------------------------
bool ClustersAreTheSame(const Cluster &ca, const Cluster &cb)
{
	if(ca.members.size() != cb.members.size())
		return false;

	for(unsigned int i = 0; i < cb.members.size(); i++)
	{
		if(ca.members[i].index != cb.members[i].index)
			return false;
	}


	return true;

}


// ------------------------------------------------------------------------
void TrimClusters(const std::vector<Cluster> &input, std::vector<Cluster> &output)
{
	for(unsigned int i = 0; i < input.size(); i++)
	{
		Cluster c = input[i];
		if(c.members.size() <= 2)
			continue;

		unsigned int size = c.members.size();
		bool inLargerCluster = false;
		bool duplicate = false;
		for(unsigned int j = 0; j < input.size(); j++)
		{
			if(i== j) continue;
			
			if(j > i)	
				duplicate = ClustersAreTheSame(c, input[j]);
			
			

			for(unsigned int k = 0; k < input[j].members.size(); k++)
			{
				if(input[j].members[k].index == c.featureIndex)
				{
					if(input[j].members.size() > size)
					{
						inLargerCluster = true;
						break;
					}					
				}				
			}

			if(inLargerCluster) break;
			if(duplicate) break;

		}	

		if(!inLargerCluster && !duplicate)
			output.push_back(c);
	}
}



// ------------------------------------------------------------------------
bool AreSimilar(const WholeFeature &f1, const WholeFeature &f2, double eLocation, double eOrientation, double eScale)
{
	// compute the locations distance
	utils::ImageSlice::PointType locationDiff = f1.location - f2.location;
	double locDiff = std::sqrt((locationDiff[0]*locationDiff[0]) + (locationDiff[1] *locationDiff[1]));


	locDiff /= f1.scale;


	// compute the angular difference
	double angleDiff = std::atan2(sin(f1.angle-f2.angle), cos(f1.angle-f2.angle));
	angleDiff = utils::RadiansToDegrees(utils::UnwrapAngle(angleDiff));
	if(angleDiff > 180)
		angleDiff = 360 - angleDiff;

	angleDiff = utils::DegreesToRadians(angleDiff);


	// scale diff
	double scaleDiff = std::fabs(std::log(f2.scale / f1.scale));



	if(scaleDiff <= eScale && locDiff <= eLocation && angleDiff < eOrientation)
	{
		return true;
	}
	return false;
}



double FeatureDistance(const filter::HoGFeature<2> &f1, const filter::HoGFeature<2> &f2)
{
	// compute the norms of the histograms
	unsigned int spatialBins = f1.histogram.size();
	unsigned int orientationBins = f1.histogram.front().size();
	double sum = 0.0;
	for(unsigned int i = 0; i < spatialBins; i++)
	{
		for(unsigned int j = 0; j < orientationBins; j++)
		{
			double diff = f1.histogram[i][j]-f2.histogram[i][j];
			sum += diff*diff;		
		}
	}

	double histDiff = std::sqrt(sum);

	return histDiff;
}

double ImageContainsData(const utils::ImageSlice * image)
{
	typedef itk::StatisticsImageFilter<utils::ImageSlice> StatsFilter;
	StatsFilter::Pointer stats = StatsFilter::New();
	stats->SetInput(image);
	stats->Update();

	if(stats->GetMean() > 0)
		return true;
	else
		return false;
		
}




		/*
		std::cout << max << " " << maxId << std::endl;

		vtkSmartPointer<vtkPolyData> poutput =
			vtkSmartPointer<vtkPolyData>::New();


		for(unsigned int i = 0; i < allGroups[maxId].size(); i++)
		{
			WholeFeature feat = allGroups[maxId][i];

			// write the polydata
			// draw the line in as well
			double x = cos(feat.angle);
			double y = sin(feat.angle);

			double p1 = feat.location[0];
			double p2 = feat.location[1];




			double np1 = p1+(feat.scale*x);
			double np2 = p2+(feat.scale*y);

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
			polygonSource->SetRadius(feat.scale);
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





		}

		vtkSmartPointer<vtkPolyDataWriter> writer = 
			vtkSmartPointer<vtkPolyDataWriter>::New();
		writer->SetInputData(poutput);
		writer->SetFileName("polyGroup.vtk");
		writer->Write();
		writer->Update();
		ImageSlice::Pointer extractedSlice = ImageSlice::New();
		ExtractSlice<ImageVolume, ImageSlice>(images[0], sliceToTest, extractedSlice);


		ImageSliceIO::Write("image.nrrd", extractedSlice);
		*/

