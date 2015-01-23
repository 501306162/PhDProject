#include "LineIntersectionFinder.h"


namespace vt
{
// ------------------------------------------------------------------------
void LineIntersectionFinder::Compute()
{
	// get the image normal and point 
	VectorType imageNormal, imagePoint;
	imageNormal(0) = m_Image->GetDirection()(0,2);
	imageNormal(1) = m_Image->GetDirection()(1,2);
	imageNormal(2) = m_Image->GetDirection()(2,2);

	imagePoint(0) = m_Image->GetOrigin()[0];
	imagePoint(1) = m_Image->GetOrigin()[1];
	imagePoint(2) = m_Image->GetOrigin()[2];


	// find the line at the intersection of the planes
	PlaneType plane1(m_PlaneNormal, m_PlanePoint);
	PlaneType plane2(imageNormal, imagePoint);
	LineType intersectingLine;
	FindIntersectingLine(plane1, plane2, intersectingLine);


	// get the planes that bound the bounding box
	double minT, maxT;
	GetBoundingPlaneIntersections(intersectingLine, m_BoundingBox, minT, maxT);
	
	// set the output points
	VectorType start = intersectingLine.pointAt(minT);
	VectorType end = intersectingLine.pointAt(maxT);

	double trace = 0.0;
	for(unsigned int i = 0; i < 3; i++)
	{
		trace += m_Image->GetDirection()(i,i);
	}

	for(unsigned int i = 0; i < 3; i++)
	{
		if(trace > 0.0)
		{
			m_Output.p1[i] = end(i);
			m_Output.p2[i] = start(i);
		}
		else
		{
			m_Output.p2[i] = end(i);
			m_Output.p1[i] = start(i);
		}
	}

	m_Output.direction = -intersectingLine.direction();
}

// ------------------------------------------------------------------------
void LineIntersectionFinder::GetBoundingPlaneIntersections(const LineType &line, 
		vtkBoundingBox &box, double &minT, double &maxT)
{
	// bounding box is aligned with standard coords
	VectorType xNormal, yNormal, zNormal;
	xNormal[0] = 1.0;
	xNormal[1] = 0.0;
	xNormal[2] = 0.0;
	yNormal[0] = 0.0;
	yNormal[1] = 1.0;
	yNormal[2] = 0.0;
	zNormal[0] = 0.0;
	zNormal[1] = 0.0;
	zNormal[2] = 1.0;

	const double * maxPoint = box.GetMaxPoint();
	const double * minPoint = box.GetMinPoint();

	VectorType minp;
	VectorType maxp;

	for(unsigned int i = 0; i < 3; i++)
	{
		minp(i) = minPoint[i];
		maxp(i) = maxPoint[i];
	}
	
	PlaneType minXPlane(xNormal, minp);
	PlaneType minYPlane(yNormal, minp);
	PlaneType minZPlane(zNormal, minp);
	PlaneType maxXPlane(xNormal, maxp);
	PlaneType maxYPlane(yNormal, maxp);
	PlaneType maxZPlane(zNormal, maxp);



	// get the intersections points with the line 
	VectorType i1,i2,i3,i4,i5,i6;
	i1 = line.intersectionPoint(minXPlane);
	i2 = line.intersectionPoint(minYPlane);
	i3 = line.intersectionPoint(minZPlane);
	i4 = line.intersectionPoint(maxXPlane);
	i5 = line.intersectionPoint(maxYPlane);
	i6 = line.intersectionPoint(maxZPlane);

	double ip1[3],ip2[3],ip3[3],ip4[3],ip5[3],ip6[3];

	for(unsigned int i = 0; i < 3; i++)
	{
		ip1[i] = i1(i);		
		ip2[i] = i2(i);		
		ip3[i] = i3(i);		
		ip4[i] = i4(i);		
		ip5[i] = i5(i);		
		ip6[i] = i6(i);		
	}

	box.Inflate(0.1);
	minT = std::numeric_limits<float>::max();
	maxT = std::numeric_limits<float>::min();

	if(box.ContainsPoint(ip1))
	{
		minT = std::min(line.intersectionParameter(minXPlane), minT);
		maxT = std::max(line.intersectionParameter(minXPlane), maxT);
	}

	if(box.ContainsPoint(ip2))
	{
		minT = std::min(line.intersectionParameter(minYPlane), minT);
		maxT = std::max(line.intersectionParameter(minYPlane), maxT);
	}
	
	if(box.ContainsPoint(ip3))
	{
		minT = std::min(line.intersectionParameter(minZPlane), minT);
		maxT = std::max(line.intersectionParameter(minZPlane), maxT);
	}

	if(box.ContainsPoint(ip4))
	{
		minT = std::min(line.intersectionParameter(maxXPlane), minT);
		maxT = std::max(line.intersectionParameter(maxXPlane), maxT);
	}

	if(box.ContainsPoint(ip5))
	{
		minT = std::min(line.intersectionParameter(maxYPlane), minT);
		maxT = std::max(line.intersectionParameter(maxYPlane), maxT);
	}

	if(box.ContainsPoint(ip6))
	{
		minT = std::min(line.intersectionParameter(maxZPlane), minT);
		maxT = std::max(line.intersectionParameter(maxZPlane), maxT);
	}


}


// ------------------------------------------------------------------------
void LineIntersectionFinder::FindIntersectingLine(const PlaneType &p1, const PlaneType &p2, LineType &line)
{
	Eigen::Matrix<double,2,3> A;
	Eigen::Matrix<double,2,1> b;
	b(0) = -p1.offset();
	b(1) = -p2.offset();

	for(unsigned int i = 0; i < 3; i++)
	{
		A(0,i) = p1.normal()(i);
		A(1,i) = p2.normal()(i);
	}

	VectorType x = A.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b);
	VectorType outputNormal = p1.normal().cross(p2.normal());
	line = LineType(x, outputNormal);
}

}
