#include "BoundingBox.h"
#include <qjson/parser.h>
#include <QFile>
#include <QVariant>

#include <itkImageRegionIterator.h>
#include <vtkCubeSource.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkMatrix4x4.h>
#include <vtkPointData.h>
#include <vtkSelectEnclosedPoints.h>

namespace vt
{
// ------------------------------------------------------------------------
void BoundingBox::Load(const std::string &filename, const int exclude)
{
	//  open the file 
	QFile file(QString::fromStdString(filename));
	file.open(QIODevice::ReadOnly | QIODevice::Text);


	// parse the input file
	QJson::Parser parser;
	bool ok;
	QVariantMap pointMap = parser.parse(&file, &ok).toMap();
	file.close();

	std::cout << pointMap.size() << std::endl;


	// get the exclude key. If it is in the map, remove it
	QString excludeKey =  "d" + QString::number(exclude);
	if(pointMap.contains(excludeKey))
	{
		pointMap.remove(excludeKey);
	}

	// start reading the points from the list
	QVariantMap::iterator mapIt = pointMap.begin();
	PointContainerType::Pointer containerP1 = PointContainerType::New();
	PointContainerType::Pointer containerP2 = PointContainerType::New();
	
	vtkSmartPointer<vtkPoints> points1 = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> points2 = vtkSmartPointer<vtkPoints>::New();

	while(mapIt != pointMap.end())
	{
		QVariantMap points = mapIt.value().toMap();
		QVariantList p1List = points["p1"].toList();
		QVariantList p2List = points["p2"].toList();

		PointType p1, p2;

		for(unsigned int i = 0; i < 3; i++)
		{
			p1[i] = p1List[i].toDouble();			
			p2[i] = p2List[i].toDouble();			
		}
		points1->InsertNextPoint(p1.GetDataPointer());
		points2->InsertNextPoint(p2.GetDataPointer());

		containerP1->push_back(p1);
		containerP2->push_back(p2);

		++mapIt;
	}

	SetBoundingBox(points1, m_BoundingBoxP1);
	SetBoundingBox(points2, m_BoundingBoxP2);
}

// ------------------------------------------------------------------------
void BoundingBox::SetBoundingBox(vtkSmartPointer<vtkPoints> &points, vtkSmartPointer<vtkPolyData> &box)
{
	// initialise if needed
	vtkBoundingBox * bounding = new vtkBoundingBox;
	for(unsigned int i = 0; i < points->GetNumberOfPoints(); i++)
	{
		bounding->AddPoint(points->GetPoint(i));
	}

	bounding->Inflate(m_Inflation);
	double bounds[6];
	bounding->GetBounds(bounds);

	vtkSmartPointer<vtkCubeSource> cubeSource = vtkSmartPointer<vtkCubeSource>::New();
	cubeSource->SetBounds(bounds);
	cubeSource->Update();

	if(!box) box = vtkSmartPointer<vtkPolyData>::New();
	box->DeepCopy(cubeSource->GetOutput());
}

// ------------------------------------------------------------------------
void BoundingBox::TransformBoundingBox(const TransformType::Pointer &transform)
{
	vtkSmartPointer<vtkPolyData> newP1;
	vtkSmartPointer<vtkPolyData> newP2;

	ApplyTransform(m_BoundingBoxP1, transform, newP1);
	ApplyTransform(m_BoundingBoxP2, transform, newP2);

	m_BoundingBoxP1 = newP1;
	m_BoundingBoxP2 = newP2;
}

// ------------------------------------------------------------------------
void BoundingBox::ApplyTransform(vtkSmartPointer<vtkPolyData> &input, 
		const TransformType::Pointer &transform, 
		vtkSmartPointer<vtkPolyData> &output)
{
	// get the point containers
	vtkSmartPointer<vtkTransform> trans = vtkSmartPointer<vtkTransform>::New();
	vtkSmartPointer<vtkMatrix4x4> transMat = vtkSmartPointer<vtkMatrix4x4>::New();
	

	// create the vtk transformation matrix
	TransformType::MatrixType rot = transform->GetMatrix();
	TransformType::OutputVectorType translation = transform->GetTranslation();
	
	for(unsigned int i = 0; i < 3; i++)
	{
		for(unsigned int j = 0; j < 3; j++)
		{
			transMat->SetElement(i,j, rot(i,j));
		}
		transMat->SetElement(i,3, translation[i]);
	}
	transMat->SetElement(3,3,1);
	trans->SetMatrix(transMat);

	
	// perform the transform
	vtkSmartPointer<vtkTransformPolyDataFilter> transformer = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformer->SetTransform(trans);
	transformer->SetInputData(input);
	transformer->Update();

	if(!output) output = vtkSmartPointer<vtkPolyData>::New();
	output->DeepCopy(transformer->GetOutput());
}

// -----------------------------------------------------------------------
void BoundingBox::ApplyInflation(BoundingBoxType::Pointer &box, double &inflation)
{
	PointType max = box->GetMaximum();
	PointType min = box->GetMinimum();

	for(unsigned int i = 0; i < 3; i++)
	{
		max[i] += inflation;
		min[i] -= inflation;
	}
	
	box->ConsiderPoint(max);
	box->ConsiderPoint(min);

}

// ------------------------------------------------------------------------
void BoundingBox::ComputeImageMask(const ImageType::Pointer &image, unsigned int pointType, MaskType::Pointer &mask)
{
	mask->SetDirection(image->GetDirection());
	mask->SetOrigin(image->GetOrigin());
	mask->SetSpacing(image->GetSpacing());
	mask->SetRegions(image->GetLargestPossibleRegion());
	mask->Allocate();
	mask->FillBuffer(0);

	// now we iterate through the mask pixels and check which ones fall inside the 
	// bounding boxes
	itk::ImageRegionIterator<MaskType> maskIt(mask, mask->GetLargestPossibleRegion());
	while(!maskIt.IsAtEnd())
	{
		MaskType::IndexType index = maskIt.GetIndex();
		PointType point;

		mask->TransformIndexToPhysicalPoint(index, point);

		bool inside1 = IsInside(point, m_BoundingBoxP1);
		bool inside2 = IsInside(point, m_BoundingBoxP2);

		if(pointType == 0)
		{
			if(inside1) maskIt.Set(1);
			if(inside2) maskIt.Set(1);
			if(inside1 && inside2) maskIt.Set(1);
		}
		else if(pointType == 1)
		{
			if(inside1) maskIt.Set(1);
		}
		else if(pointType == 2)
		{
			if(inside2) maskIt.Set(1);
		}
		
		++maskIt;
	}

	mask->SetDirection(image->GetDirection());
}

// ------------------------------------------------------------------------
bool BoundingBox::IsInside(const PointType &point, vtkSmartPointer<vtkPolyData> &box)
{
	vtkSmartPointer<vtkPolyData> points = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> p = vtkSmartPointer<vtkPoints>::New();
	p->InsertNextPoint(point.GetDataPointer());
	points->SetPoints(p);

	vtkSmartPointer<vtkSelectEnclosedPoints> selector = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
	selector->SetSurfaceData(box);
	selector->SetInputData(points);
	selector->Update();

	if(selector->IsInside(0))
		return true;
	return false;

}


// ------------------------------------------------------------------------
void BoundingBox::ExtractMaskMinMaxIndex(const ImageType::Pointer &image, 
		const BoundingBoxType::Pointer &box,
		ImageType::RegionType &region)
{
	PointType max = box->GetMaximum();
	PointType min = box->GetMinimum();

	IndexType maxIndex, minIndex;
	image->TransformPhysicalPointToIndex(max, maxIndex);
	image->TransformPhysicalPointToIndex(min, minIndex);

	// zero the z index
	minIndex[2] = 0;
	maxIndex[2] = 0;

	// get the size 
	ImageType::SizeType size;
	size.Fill(0);

	for(unsigned int i = 0; i < 3; i++)
	{
		size[i] = maxIndex[i]-minIndex[i];
	}

	std::cout << min << " " << max  << std::endl;
	std::cout << size << " " << minIndex << " " << maxIndex <<  std::endl;

	region.SetIndex(minIndex);
	region.SetSize(size);

}



}
