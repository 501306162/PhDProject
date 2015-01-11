#include "ValveOriginFinder.h"

#include <itkRegionOfInterestImageFilter.h>
#include <itkPermuteAxesImageFilter.h>
#include <itkFlipImageFilter.h>

#include <Eigen/Dense>

namespace vt
{
// ------------------------------------------------------------------------
void ValveOriginFinder::Compute()
{
	// extract the SA central image
	ImageType::Pointer saSlice = ImageType::New();
	ExtractSASlice(m_Stack, saSlice);

	Plane saPlane;
	ComputePlane(saSlice, saPlane);

	// compute the plane equations of the 2c and 3c images
	Plane c2Plane;
	ComputePlane(m_2CImage, c2Plane);

	Plane c3Plane;
	ComputePlane(m_3CImage, c3Plane);


	// build the matrixes
	Eigen::Matrix3f A;
	Eigen::Vector3f b;

	for(unsigned int i = 0; i < 3; i++)
	{
		A(0,i) = saPlane.normal[i]; 
		A(1,i) = c2Plane.normal[i]; 
		A(2,i) = c3Plane.normal[i]; 
	}

	b(0) = -saPlane.d;
	b(1) = -c2Plane.d;
	b(2) = -c3Plane.d;

	Eigen::Vector3f x = A.fullPivLu().solve(b);

	m_Origin[0] = x(0,0);
	m_Origin[1] = x(1,0);
	m_Origin[2] = x(2,0);

	m_Translation[0] = x(0,0);
	m_Translation[1] = x(1,0);
	m_Translation[2] = x(2,0);

	// the Y axis is the stack normal
	m_YAxis = saPlane.normal;
	m_Rotation = m_Stack->GetDirection();
}

// ------------------------------------------------------------------------
void ValveOriginFinder::ComputeNormalise()
{

	// extract the SA central image
	ImageType::Pointer saSlice = ImageType::New();
	ExtractSASlice(m_Stack, saSlice);

	Plane saPlane;
	ComputePlane(saSlice, saPlane);

	// compute the plane equations of the 2c and 3c images
	Plane c2Plane;
	ComputePlane(m_2CImage, c2Plane);

	Plane c3Plane;
	ComputePlane(m_3CImage, c3Plane);

	// build the matrixes
	Eigen::Matrix3f A;
	Eigen::Vector3f b;

	for(unsigned int i = 0; i < 3; i++)
	{
		A(0,i) = saPlane.normal[i]; 
		A(1,i) = c2Plane.normal[i]; 
		A(2,i) = c3Plane.normal[i]; 
	}

	b(0) = -saPlane.d;
	b(1) = -c2Plane.d;
	b(2) = -c3Plane.d;

	Eigen::Vector3f x = A.fullPivLu().solve(b);

	m_Origin[0] = x(0,0);
	m_Origin[1] = x(1,0);
	m_Origin[2] = x(2,0);


	VectorType perp = itk::CrossProduct(saPlane.normal, c2Plane.normal);
	perp.Normalize();

	VectorType p1;
   	p1[0] = m_Origin[0];
   	p1[1] = m_Origin[1];
   	p1[2] = m_Origin[2];



	m_Center = p1;


	VectorType upVec;
	VectorType rotVec;

	upVec[0] = m_2CImage->GetDirection()(0,1);
	upVec[1] = m_2CImage->GetDirection()(1,1);
	upVec[2] = m_2CImage->GetDirection()(2,1);

	rotVec[0] = -m_2CImage->GetDirection()(0,2);
	rotVec[1] = -m_2CImage->GetDirection()(1,2);
	rotVec[2] = -m_2CImage->GetDirection()(2,2);

	upVec.Normalize();
	m_Angle = acos(perp*upVec);
	m_Axis = rotVec;



}

// ------------------------------------------------------------------------
void ValveOriginFinder::FlipImage(const ImageType::Pointer &input, ImageType::Pointer &output)
{
	ImageType::DirectionType dir = input->GetDirection();
	typedef itk::PermuteAxesImageFilter<ImageType> PermType;
	PermType::Pointer permer = PermType::New();
	permer->SetInput(input);

	PermType::PermuteOrderArrayType axes;
	axes[0] = 1;
	axes[1] = 0;
	axes[2] = 2;

	permer->SetOrder(axes);
	permer->Update();

	std::cout << input->GetOrigin() << std::endl;
	output = permer->GetOutput();
	output->SetDirection(dir);
	std::cout << output->GetOrigin() << std::endl;
}



// ------------------------------------------------------------------------
void ValveOriginFinder::ComputePlane(const ImageType::Pointer &image, Plane &plane)
{
	PointListType points;
	ExtractPlanePoints(image, points);
	

	// get vec AB
	itk::Vector<double, 3> AB, AC;
	AB = points[1]-points[0];
	AC = points[2]-points[0];

	// compute the normal
	itk::Vector<double, 3> normal;
   	normal[0] = image->GetDirection()(0,2);
   	normal[1] = image->GetDirection()(1,2);
   	normal[2] = image->GetDirection()(2,2);

	// compute d
	double sum = 0;
	for(unsigned int i = 0; i < 3; i++)
	{
		sum += normal[i] * points[0][i];
	}

	plane.normal = normal;
	plane.d = -sum;
}

// ------------------------------------------------------------------------
void ValveOriginFinder::ExtractPlanePoints(const ImageType::Pointer &image, PointListType &points)
{
	points.push_back(image->GetOrigin());

	ImageType::SizeType imageSize = image->GetLargestPossibleRegion().GetSize();
	ImageType::IndexType ind1, ind2;
	ind1.Fill(0); ind2.Fill(0);
	ind1[1] = imageSize[0]-1;
	ind2[2] = imageSize[1]-1;

	PointType p1, p2;
	image->TransformIndexToPhysicalPoint(ind1,p1);
	image->TransformIndexToPhysicalPoint(ind2,p2);

	points.push_back(p1);
	points.push_back(p2);
}

// ------------------------------------------------------------------------
void ValveOriginFinder::ExtractSASlice(const ImageType::Pointer &image, ImageType::Pointer &slice)
{
	typedef itk::RegionOfInterestImageFilter<ImageType, ImageType> ROIType;
	ROIType::Pointer roiFilter = ROIType::New();
	roiFilter->SetInput(image);

	// create the region 
	ImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
	double middle = (double) size[2] / 2.0;
	int index = (int) floor(middle);

	ImageType::IndexType start;
	start.Fill(0);
	start[2] = index;

	size[2] = 1;
	ImageType::RegionType region;
	region.SetSize(size);
	region.SetIndex(start);

	roiFilter->SetRegionOfInterest(region);
	roiFilter->Update();
	slice = roiFilter->GetOutput();
}


}	
