#ifndef OPEN_CV_VALVE_H
#define OPEN_CV_VALVE_H

#include <opencv2/core/core.hpp>

#include <itkObject.h>
#include <itkObjectFactory.h>

#include "ValveLine.h"

namespace vt
{
class OpenCVValve : public itk::Object 
{	
public:
	typedef OpenCVValve Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(OpenCVValve, Object);
	itkNewMacro(Self);

	typedef cv::Mat MatType;
	typedef cv::Ptr<MatType> MatPtr;
	typedef ValveLine<2> ValveType;
	typedef ValveType::ImageType ImageType;

	itkSetMacro(P1, cv::Point2f);
	itkGetMacro(P1, cv::Point2f);
	itkSetMacro(P2, cv::Point2f);
	itkGetMacro(P2, cv::Point2f);

	itkGetMacro(Image, MatPtr);
	itkSetMacro(Image, MatPtr);

	void InitialiseFromValve(const ValveLine<2>::Pointer &input);

protected:
	OpenCVValve() {}
	virtual ~OpenCVValve() {}

	void ConvertImage(const ImageType::Pointer &input, MatPtr &mat);

private:

	cv::Ptr<cv::Mat> m_Image;
	cv::Point2f m_P1;
	cv::Point2f m_P2;


	OpenCVValve(const Self&);
	void operator=(const Self&);

};
}


#endif
