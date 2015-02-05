#include <iostream>
#include <CommonDefinitions.h>
#include <MatrixCommon.h>
#include <MatrixReader.h>
#include <random>
#include <SVMClassifier.h>
#include <ValveIO.h>
#include <opencv2/ml/ml.hpp>
#include <opencv2/core/core.hpp>

#include "ExtractTrainingFeatures.h"

#include <itkImageRegionIterator.h>
#include <ParameterHelpers.h>

typedef std::pair<int, double> IndexDistPair;
bool index_sort(const IndexDistPair &a, const IndexDistPair &b)
{
	return (a.second < b.second);
}

typedef utils::DoubleMatrixType MatrixType;
typedef utils::IntMatrixType IntMatrixType;
void combine(const MatrixType &pos, const MatrixType &neg, MatrixType &all, IntMatrixType &labels);
cv::Mat convert(const MatrixType &mat);
cv::Mat convert(const IntMatrixType &mat);

using namespace vt;
int main(int argc, char ** argv)
{
	// load the input matrix
	utils::MatrixReader::Pointer reader = utils::MatrixReader::New();
	reader->SetFilename("/home/om0000/ValveTracking/TrainingData/conf6/conf6-MV-2C-0.hdf5");
	reader->Read();
	

	Params pparams("/home/om0000/ValveTracking/TrainingData/conf6/extract_features.json");

	MatrixType pos = reader->GetOutput()->doubleData["p1-pos"];
	MatrixType neg = reader->GetOutput()->doubleData["p1-neg"];
	

	for(unsigned int tNum = 0; tNum < pos.rows(); tNum++)
	{
		MatrixType test = pos.row(tNum);		
		MatrixType others(pos.rows()-1, pos.cols());
		unsigned int count = 0;
		for(unsigned int i = 0; i < pos.rows(); i++)
		{
			if(i != tNum)
			{
				others.row(count) = pos.row(i);			
				count++;
			}
		}

		unsigned int T = others.rows();
		unsigned int Knn = 5;
		unsigned int newCount = 0;
		unsigned int NN = 50;
		srand (time(NULL));

		MatrixType smotePos(NN*T, pos.cols());

		for(unsigned int i = 0; i < T; i++)
		{
			MatrixType s1 = others.row(i);
			MatrixType dists(T,1);
			std::vector<IndexDistPair> sortList;
			for(unsigned int j = 0; j < T; j++)
			{
				if(i == j) continue;
				std::pair<int, double> p(j, (others.row(j)-s1).norm());
				sortList.push_back(p);
			}

			std::sort(sortList.begin(), sortList.end(), index_sort);

			unsigned int N = NN;
			while(N != 0)
			{
				int nn = sortList[rand() % Knn].first;
				MatrixType newF(1, others.cols());
				for(unsigned int j = 0; j < others.cols(); j++)
				{
					double diff = others(nn,j) - s1(0,j);
					double gap = (double) rand() / RAND_MAX;
					smotePos(newCount, j) = round(s1(0,j) + (gap*diff));
				}
				newCount++;
				N--;
			}
		}

		// combine the pos and neg 
		MatrixType newAll, all;
		IntMatrixType newLabels, labels;
		combine(others, neg,  all, labels);
		combine(smotePos, neg, newAll, newLabels);

		/*
		cv::Mat vType = cv::Mat(pos.cols()+1, 1, CV_8U);
		vType.setTo(cv::Scalar(CV_VAR_NUMERICAL));
		vType.at<uchar>(pos.cols(), 0) = CV_VAR_CATEGORICAL;

		cv::Mat X = convert(newAll);
		cv::Mat y = convert(newLabels);
		cv::Mat cvtest = convert(test);

		cv::BoostParams params;
		params.boost_type = CvBoost::GENTLE;
		cv::Boost cls;
		std::cout << "Training" << std::endl;
		//cls.train(X, CV_ROW_SAMPLE, y,
				//cv::Mat(), cv::Mat(), vType, cv::Mat(), params);

		cls.load("cls.cls", "hello");

		std::cout << "Testing" << std::endl;
		std::cout << cls.predict(cvtest, cv::Mat(), cv::Range::all(), false, true) << std::endl;

		/*
		cv::RandomTreeParams params;
		params.max_depth = 25;
		params.regression_accuracy = 0.00001;
		

		cv::NormalBayesClassifier cls;
		cls.train(X, y);

		cv::Mat out;
		cls.predict(cvtest, &out);
		std::cout << out << std::endl;
		*/

		//cv::RandomTrees trees;
		//trees.train(X, CV_ROW_SAMPLE, y,
				//cv::Mat(), cv::Mat(), vType, cv::Mat(), params);

		//std::cout << trees.predict(cvtest) << std::endl;
		/*
		MatrixType probs1, probs2;
		IntMatrixType cls1, cls2;

		*/
		vt::SVMClassifier::Pointer svm1 = vt::SVMClassifier::New();
		svm1->Train(all, labels);


		vt::SVMClassifier::Pointer svm2 = vt::SVMClassifier::New();
		svm2->Train(newAll, newLabels);




		// load the image 
		ValveSequenceReader<3>::Pointer reader = ValveSequenceReader<3>::New();
		reader->SetFileName(argv[1]);
		ValveSequence<3>::Pointer sequence = reader->GetOutput();

		typedef ValveLine<3> ValveType;
		typedef ValveType::ImageType ImageType;
		ValveType::Pointer valve = sequence->GetValveLine(0);
		ImageType::Pointer image = valve->GetImage();
		
		typedef itk::Image<double, 3> RealImageType;
		RealImageType::Pointer output = RealImageType::New();
		output->SetOrigin(image->GetOrigin());
		output->SetSpacing(image->GetSpacing());
		output->SetRegions(image->GetLargestPossibleRegion());
		output->SetDirection(image->GetDirection());
		output->Allocate();

		RealImageType::Pointer output2 = RealImageType::New();
		output2->SetOrigin(image->GetOrigin());
		output2->SetSpacing(image->GetSpacing());
		output2->SetRegions(image->GetLargestPossibleRegion());
		output2->SetDirection(image->GetDirection());
		output2->Allocate();
	
		
		itk::ImageRegionIterator<RealImageType> it(output, output->GetLargestPossibleRegion());
		itk::ImageRegionIterator<RealImageType> it2(output2, output2->GetLargestPossibleRegion());
		while(!it.IsAtEnd())
		{
			ImageType::IndexType index = it.GetIndex();
			ImageType::PointType pt;
			output->TransformIndexToPhysicalPoint(index, pt);

			MatrixType feature;
			extractLBPFeature(pparams, valve, pt, feature);

			MatrixType p1,p2;
			IntMatrixType c1, c2;

			svm1->PredictProbability(feature, c1, p1);
			svm2->PredictProbability(feature, c2, p2);

			it.Set(p1(0,0));
			it2.Set(p2(0,0));


			++it2;
			++it;
		}

		utils::RealVolumeIO::Write("votes.nrrd", output);
		utils::RealVolumeIO::Write("votes2.nrrd", output);
		utils::ImageVolumeIO::Write("image.nrrd", image);
		return 0;


	}







	
	return 0;



}

// ------------------------------------------------------------------------
void combine(const MatrixType &pos, const MatrixType &neg, MatrixType &all, IntMatrixType &labels)
{
	unsigned int total = pos.rows() + neg.rows();
	all = MatrixType::Zero(total, pos.cols());
	labels = IntMatrixType::Zero(total, 1);

	unsigned int count = 0;
	for(unsigned int i = 0; i < pos.rows(); i++)
	{
		all.row(count) = pos.row(i);
		labels(count, 0) = 1;
		count++;
	}

	for(unsigned int i = 0; i < neg.rows(); i++)
	{
		all.row(count) = neg.row(i);
		labels(count, 0) = 0;
		count++;
	}

}

// ------------------------------------------------------------------------
cv::Mat convert(const IntMatrixType &mat)
{
	cv::Mat output(mat.rows(), mat.cols(), CV_32FC1);
	for(unsigned int i = 0; i < mat.rows(); i++)
	{
		for(unsigned int j = 0; j < mat.cols(); j++)
		{
			output.at<float>(i,j) = mat(i,j);
		}
	}

	return output;
}

// ------------------------------------------------------------------------
cv::Mat convert(const MatrixType &mat)
{
	cv::Mat output(mat.rows(), mat.cols(), CV_32FC1);
	for(unsigned int i = 0; i < mat.rows(); i++)
	{
		for(unsigned int j = 0; j < mat.cols(); j++)
		{
			output.at<float>(i,j) = mat(i,j);
		}
	}


	return output;
}


