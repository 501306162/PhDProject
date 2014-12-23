#include <iostream>

#include <OpenCVValve.h>
#include <ValveIO.h>
#include <Directory.h>

#include <opencv2/features2d/features2d.hpp>
#include <opencv2/nonfree/features2d.hpp>
#include <opencv2/ml/ml.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace vt;

typedef std::vector<cv::KeyPoint> KeyPointsType;
void computeDescriptors(KeyPointsType &keyPoints, cv::Mat &image, std::vector<cv::Mat> &descriptors);
void computePositiveKeyPoints(const OpenCVValve::Pointer &valve, KeyPointsType &keyPoint, float size);
void computeNegativeKeyPoints(const OpenCVValve::Pointer &valve, KeyPointsType &keyPoint, float size, unsigned int numKeyPoints);
void computeTestingDescriptors(const std::string &filename, cv::Mat &descritors, cv::Mat &indices, float size);
void computeProbImage(const std::string &filename, CvRTrees * tree, cv::Mat &descriptors, cv::Mat &indices);



int main(int argc, char ** argv)
{




	utils::Directory::FilenamesType fnames = utils::Directory::GetFiles(argv[1], ".txt");

	float size = 10;
	std::vector<cv::Mat> desc;
	std::vector<cv::Mat> positiveDescriptors, negativeDescriptors;
	for(unsigned int i = 0; i < fnames.size()-1; i++)
	{
		std::cout  << fnames[i] << std::endl;
		ValveSequenceReader<2>::Pointer reader = ValveSequenceReader<2>::New();
		reader->SetFileName(fnames[i]);
		ValveSequence<2>::Pointer valves = reader->GetOutput();

		OpenCVValve::Pointer cvValve = OpenCVValve::New();
		cvValve->InitialiseFromValve(valves->GetValveLine(0));

		KeyPointsType positiveKeyPoints;
		computePositiveKeyPoints(cvValve, positiveKeyPoints, size);
		computeDescriptors(positiveKeyPoints, *cvValve->GetImage(), positiveDescriptors);

		KeyPointsType negativeKeyPoints;
		computeNegativeKeyPoints(cvValve, negativeKeyPoints, size, 1);
		computeDescriptors(negativeKeyPoints, *cvValve->GetImage(), negativeDescriptors);

	}


	// build the training matrix 
	cv::Mat classifications(positiveDescriptors.size() + negativeDescriptors.size(), 1, CV_32F);
	cv::Mat trainingData(positiveDescriptors.size() + negativeDescriptors.size(), positiveDescriptors.front().cols, CV_32F);

	unsigned int rowCount = 0;
	for(unsigned int i = 0; i < positiveDescriptors.size(); i++)
	{
		classifications.at<float>(rowCount,0) = 1;
		trainingData.row(rowCount) = positiveDescriptors[i].row(0);
		rowCount++;
	}

	for(unsigned int i = 0; i < negativeDescriptors.size(); i++)
	{
		classifications.at<float>(rowCount,0) = 0;
		trainingData.row(rowCount) = negativeDescriptors[i].row(0);
		rowCount++;
	}



	// compute the testing descriptors
	cv::Mat testIndices;
	cv::Mat testDescriptors;

	computeTestingDescriptors(fnames.front(), testDescriptors, testIndices, size);



	cv::FileStorage fs("test.yml", cv::FileStorage::WRITE);
	fs << "Training" << trainingData;
	fs << "Labels" << classifications;
	fs << "TestD" << testDescriptors;
	fs << "TestI" << testIndices;
	fs.release();



	cv::FileStorage fsIn("test.yml", cv::FileStorage::READ);
	cv::Mat tdata, ldata, testD, testI;
	fsIn["Training"] >> tdata;
	fsIn["Labels"] >> ldata;
	fsIn["TestD"] >> testD;
	fsIn["TestI"] >> testI;
	fsIn.release();
	cv::Mat var_type = cv::Mat(tdata.cols + 1, 1, CV_8U );
	var_type.setTo(cv::Scalar(CV_VAR_NUMERICAL) ); // all inputs are numerical

	// this is a classification problem (i.e. predict a discrete number of class
	// outputs) so reset the last (+1) output var_type element to CV_VAR_CATEGORICAL

	var_type.at<uchar>(tdata.cols, 0) = CV_VAR_CATEGORICAL;


	float w1 = (float) positiveDescriptors.size() / (float) tdata.rows;
	float w2 = (float) negativeDescriptors.size() / (float) tdata.rows;


	// start building the tree
	float priors[] = {1.0,1.0};
	std::cout << w1 << " " << w2 << std::endl;
	CvRTParams params = CvRTParams(5, // max depth
			2, // min sample count
			0.0, // regression accuracy: N/A here
			true, // compute surrogate split, no missing data
			2, // max number of categories (use sub-optimal algorithm for larger numbers)
			priors, // the array of priors
			false,  // calculate variable importance
			4,       // number of variables randomly selected at node and used to find the best split(s).
			1000,	 // max number of trees in the forest
			0.0001f,				// forrest accuracy
			CV_TERMCRIT_ITER |	CV_TERMCRIT_EPS // termination cirteria
			);

	CvRTrees * rtree = new CvRTrees;

	std::cout << "training" << std::endl;
	rtree->train(tdata, CV_ROW_SAMPLE, ldata, cv::Mat(), cv::Mat(), var_type, cv::Mat(), params);

	std::cout << "testing" <<  std::endl;

	computeProbImage(fnames.front(), rtree, testD, testI);

	


	return 0;

}

// ------------------------------------------------------------------------
void computeProbImage(const std::string &filename, CvRTrees * tree, cv::Mat &descriptors, cv::Mat &indices)
{
	ValveSequenceReader<2>::Pointer reader = ValveSequenceReader<2>::New();
	reader->SetFileName(filename);
	ValveSequence<2>::Pointer sequence = reader->GetOutput();

	OpenCVValve::Pointer cvValve = OpenCVValve::New();
	cvValve->InitialiseFromValve(sequence->GetValveLine(0));

	OpenCVValve::MatPtr image = cvValve->GetImage();
	cv::Mat probIm(image->rows, image->cols, CV_8UC1);

	for(int i = 0; i < descriptors.rows; i++)
	{
		double prob = tree->predict_prob(descriptors.row(i), cv::Mat());
		unsigned char val = (unsigned char) (prob * 255.0);
		std::cout << (int) val << " " << prob << std::endl;
		probIm.at<unsigned char>(indices.at<unsigned char>(i,0), indices.at<unsigned char>(i,1)) = val;
	}


	cv::imshow("test", probIm);
	cv::waitKey();

}



// ------------------------------------------------------------------------
void computeTestingDescriptors(const std::string &filename, cv::Mat &descriptors, cv::Mat &indices, float size)
{
	ValveSequenceReader<2>::Pointer reader = ValveSequenceReader<2>::New();
	reader->SetFileName(filename);
	ValveSequence<2>::Pointer sequence = reader->GetOutput();

	OpenCVValve::Pointer cvValve = OpenCVValve::New();
	cvValve->InitialiseFromValve(sequence->GetValveLine(0));

	OpenCVValve::MatPtr image = cvValve->GetImage();
	KeyPointsType keyPoints;

	indices = cv::Mat(image->rows*image->cols, 2, CV_8UC1);

	int pcount = 0;
	for(int i = 0; i < image->rows; i++)
	{
		for(int j = 0; j < image->cols; j++)
		{
			cv::KeyPoint kp;
			kp.pt.x = j;
			kp.pt.y = i;
			kp.size = size;

			keyPoints.push_back(kp);			
			indices.at<unsigned char>(pcount,0) = i;
			indices.at<unsigned char>(pcount,1) = j;

			pcount++;
		}
	}

	// compute the descriptors
	cv::SiftDescriptorExtractor extractor;
	extractor.compute(*image, keyPoints, descriptors);
}



// ------------------------------------------------------------------------
void computePositiveKeyPoints(const OpenCVValve::Pointer &valve, KeyPointsType &keyPoint, float size)
{
	cv::KeyPoint key;
	key.pt = valve->GetP1();
	key.size = size;
	keyPoint.push_back(key);
}


// ------------------------------------------------------------------------
void computeNegativeKeyPoints(const OpenCVValve::Pointer &valve, KeyPointsType &keyPoints, float size, unsigned int numKeyPoints)
{
	OpenCVValve::MatPtr mat = valve->GetImage();
	unsigned int rows = mat->rows;
	unsigned int cols = mat->cols;
	cv::Point2f positive = valve->GetP1();

	unsigned int count = 0;
	while(count < numKeyPoints)
	{
		// get a random position
		float rand1 = ((float) std::rand() / (float) RAND_MAX);
		float rand2 = ((float) std::rand() / (float) RAND_MAX);

		float x = rand1 * (float) cols;
		float y = rand2 * (float) rows;

		cv::Point2f newPoint(x,y);
		

		// compute the distance between the two points
		cv::Point2f diff = positive - newPoint;
		double dist = sqrt((diff.x * diff.x) + (diff.y * diff.y)); 

		if(dist > size)
		{
			cv::KeyPoint kp;
			kp.pt = newPoint;
			kp.size = size;
			keyPoints.push_back(kp);
			count++;
		}
	}
}


// ------------------------------------------------------------------------
void computeDescriptors(KeyPointsType &keyPoints, cv::Mat &image, std::vector<cv::Mat> &descriptors)
{
	cv::SiftDescriptorExtractor extractor;
	cv::Mat desc;
	extractor.compute(image, keyPoints, desc);
	for(int i = 0; i < desc.rows; i++)
	{
		descriptors.push_back(desc.row(i));		
	}
}


