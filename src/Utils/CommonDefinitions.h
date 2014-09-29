#ifndef COMMON_DEFINITIONS_H
#define COMMON_DEFINITIONS_H

#include <itkImage.h>

#include "ImageIO.h"

namespace utils
{
typedef itk::Image<unsigned short, 2> ImageSlice;
typedef itk::Image<unsigned short, 3> ImageVolume;
typedef itk::Image<unsigned char, 2> LabelSlice;
typedef itk::Image<unsigned char, 3> LabelVolume;
typedef itk::Image<double, 2> RealSlice;
typedef itk::Image<double, 3> RealVolume;

typedef std::vector<LabelVolume::Pointer> LabelVolumeList;
typedef std::vector<RealVolume::Pointer> RealVolumeList;
typedef std::vector<RealSlice::Pointer> RealSliceList;
typedef std::vector<ImageVolume::Pointer> ImageVolumeList;
typedef std::vector<ImageSlice::Pointer> ImageSliceList;


typedef ImageIO<ImageSlice> ImageSliceIO;
typedef ImageIO<ImageVolume> ImageVolumeIO;
typedef ImageIO<LabelSlice> LabelSliceIO;
typedef ImageIO<LabelVolume> LabelVolumeIO;
typedef ImageIO<RealVolume> RealVolumeIO;
typedef ImageIO<RealSlice> RealSliceIO;


/* some useful functions */
double DegreesToRadians(double degree)
{
	double mult = M_PI / 180.0;
	return degree * mult;
}

double RadiansToDegrees(double rad)
{
	double mult = 180.0 / M_PI;
	return rad * mult;

}

double UnwrapAngle(double angle)
{
	if(angle >= 0) return angle;
	return 2*M_PI + angle;
}

double WrapAngle(double angle)
{
	angle = fmod(angle + 180, 360);
	if(angle < 0)
		angle+=360;
	return angle-180;	
}


double RadianAngleFrom2DVector(double x, double y, bool unwrap=false)
{
	if(!unwrap)
		return std::atan2(y,x);
	else 
		return UnwrapAngle(std::atan2(x,y));
}

double DegreeAngleFrom2DVector(double x, double y)
{
	double angle = std::atan2(y,x);
	angle = (angle > 0 ? angle : (2*M_PI+angle));
   	return angle * (360.0 / (2*M_PI));
}

double MagnitudeFrom2DVector(double x, double y)
{
	return std::sqrt((x*x)+(y*y));
}


} /* utils */ 


#endif
