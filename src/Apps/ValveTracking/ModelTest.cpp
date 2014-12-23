#include <iostream>

#include <OpenCVValve.h>
#include <ValveIO.h>
#include <Directory.h>

#include <opencv2/features2d/features2d.hpp>
#include <opencv2/nonfree/features2d.hpp>
#include <opencv2/ml/ml.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>


using namespace vt;


void extractPatch(cv::Mat &input, cv::Point2f location, unsigned int size, cv::Mat &patch);


int main(int argc, char ** argv)
{
	utils::Directory::FilenamesType fnames = utils::Directory::GetFiles(argv[1], ".txt");

	std::vector<cv::Mat> patches;
	std::vector<cv::Mat> descriptors;
	for(unsigned int i = 0; i < fnames.size()-1; i++)
	{
		std::cout  << fnames[i] << std::endl;
		ValveSequenceReader<2>::Pointer reader = ValveSequenceReader<2>::New();
		reader->SetFileName(fnames[i]);
		ValveSequence<2>::Pointer valves = reader->GetOutput();

		OpenCVValve::Pointer cvValve = OpenCVValve::New();
		cvValve->InitialiseFromValve(valves->GetValveLine(0));

		cv::Mat desc;
		cv::KeyPoint kp;
		kp.pt = cvValve->GetP1();
		std::vector<cv::KeyPoint> kps;
		kps.push_back(kp);
		cv::SurfDescriptorExtractor extractor;
		extractor.compute(*cvValve->GetImage(), kps, desc);
		descriptors.push_back(desc);


		cv::Mat patch;
		extractPatch(*cvValve->GetImage(), cvValve->GetP1(), 20, patch);
		cv::Mat normPatch;
		cv::normalize(patch ,normPatch,0,255, cv::NORM_MINMAX, -1, cv::Mat()); 

		cv::Mat thresh;
		cv::adaptiveThreshold(patch, thresh, 255, CV_ADAPTIVE_THRESH_MEAN_C, CV_THRESH_BINARY, 3, 1);

		//cv::imshow("patch", normPatch);
		//cv::waitKey();
		patches.push_back(patch);
	}


	ValveSequenceReader<2>::Pointer reader = ValveSequenceReader<2>::New();
	reader->SetFileName(fnames.back());
	ValveSequence<2>::Pointer valves = reader->GetOutput();

	OpenCVValve::Pointer cvValve = OpenCVValve::New();
	cvValve->InitialiseFromValve(valves->GetValveLine(0));

	cv::Mat test = *(cvValve->GetImage());
	cv::Mat finalOut;
	for(unsigned int i = 0; i < patches.size(); i++)
	{
		cv::Mat output;
		cv::matchTemplate(test, patches[i], output,CV_TM_CCOEFF_NORMED);
		if(i == 0)
			finalOut = output;
		else
			finalOut.mul(output, (double) patches.size());


		// get the min and max
		double min, max;
		cv::Point minLoc, maxLoc;	
		cv::minMaxLoc(output,&min,&max,&minLoc, &maxLoc);
		std::cout << i << " " << max << std::endl;
		std::cout << maxLoc << std::endl;

		//if(max > 0.9)
			//cv::rectangle(test, cv::Point(maxLoc.x, maxLoc.y), cv::Point(maxLoc.x+15, maxLoc.y+15), cv::Scalar::all(0),1,8,0);
	}


	double min, max;
	cv::Point minLoc, maxLoc;	
	cv::minMaxLoc(finalOut,&min,&max,&minLoc, &maxLoc);
	cv::rectangle(test, cv::Point(maxLoc.x, maxLoc.y), cv::Point(maxLoc.x+15, maxLoc.y+15), cv::Scalar::all(0),1,8,0);

	//finalOut /= (float) patches.size();
	cv::Mat normed;
	cv::normalize(finalOut , normed, 0, 255, cv::NORM_MINMAX, CV_8UC1, cv::Mat()); 
	

	cv::imshow("out", normed);
	cv::imshow("orig", test);
	cv::waitKey();



	return 0;
}

// ------------------------------------------------------------------------
void extractPatch(cv::Mat &input, cv::Point2f location, unsigned int size, cv::Mat &patch)
{
	// get the indices 
	int x = (int) round(location.x);
	int y = (int) round(location.y);


	patch = cv::Mat(size, size, CV_8UC1);
	for(unsigned int i = 0; i < size; i++)
	{
		int row = (y - ((size/2)-1)) + i;
		for(unsigned int j = 0; j < size; j++)
		{
			int col = (x - ((size/2)-1)) + j;
			patch.at<uchar>(i,j) = input.at<uchar>(row,col);
			
		}
	}
}
