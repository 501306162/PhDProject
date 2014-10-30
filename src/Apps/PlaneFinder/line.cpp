#include "line.h"

#include <vtkLineSource.h>
#include <vtkUnsignedCharArray.h>
#include <vtkCellData.h>
#include <vtkProperty.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkVertex.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkLine.h>

// ------------------------------------------------------------------------
Line::Line()
{
	this->poly = vtkSmartPointer<vtkPolyData>::New();
	this->mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	this->actor = vtkSmartPointer<vtkActor>::New();
}


// ------------------------------------------------------------------------
Line * Line::copy()
{
	Line * line = new Line();
	line->x = x;
	line->x2 = x2;
	line->y = y;
	line->y2 = y2;
	line->z = z;
	line->z2 = z2;

	line->poly->DeepCopy(poly);

	line->mapper->SetInputData(line->poly);
	line->actor->SetMapper(line->mapper);
	line->actor->GetProperty()->SetLineWidth(2.0);
	line->actor->GetProperty()->SetPointSize(5.0);
	line->actor->GetProperty()->SetRepresentationToWireframe();

	double col[3];
	line->getColour(col);
	line->actor->GetProperty()->SetColor(col);


	line->type = type;

	return line;
}

// ------------------------------------------------------------------------
Line * Line::NewLine(Type type, double * p1, double * p2)
{
	Line * line = new Line();
	line->type = type;

	line->x = p1[0];
	line->x2 = p2[0];
	line->y = p1[1];
	line->y2 = p2[1];
	line->z = p1[2];
	line->z2 = p2[2];


	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->InsertNextPoint(p1[0],p1[1],p1[2]);
	points->InsertNextPoint(p2[0],p2[1],p2[2]);

	vtkSmartPointer<vtkVertex> v1 = vtkSmartPointer<vtkVertex>::New();
	vtkSmartPointer<vtkVertex> v2 = vtkSmartPointer<vtkVertex>::New();
	v1->GetPointIds()->SetId(0,0);
	v2->GetPointIds()->SetId(0,1);

	vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
	vertices->InsertNextCell(v1);
	vertices->InsertNextCell(v2);

	vtkSmartPointer<vtkLine> pline = vtkSmartPointer<vtkLine>::New();
	pline->GetPointIds()->SetId(0,0);
	pline->GetPointIds()->SetId(1,1);

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	lines->InsertNextCell(pline);
	


	// create the polydata 
	line->poly->SetPoints(points);
	line->poly->SetLines(lines);
	line->poly->SetVerts(vertices);
	line->poly->Modified();

	line->mapper->SetInputData(line->poly);
	

	line->actor->SetMapper(line->mapper);
	line->actor->GetProperty()->SetLineWidth(2.0);
	line->actor->GetProperty()->SetPointSize(5.0);
	line->actor->GetProperty()->SetRepresentationToWireframe();

	double col[3];
	line->getColour(col);
	line->actor->GetProperty()->SetColor(col);

	return line;
}


// ------------------------------------------------------------------------
Line * Line::NewLine(vtkImageData * image, Type type)
{
	Line * line = new Line();
	line->type = type;


	// get the image bounds 
	double bounds[6];
   	image->GetBounds(bounds);

	// set the initial points
	line->x = (bounds[1]-bounds[0])/2.0 + image->GetOrigin()[0];
	line->x2 = line->x;
	line->y = ((bounds[3]-bounds[2]) / 2.0) - 50 + image->GetOrigin()[1];
	line->y2 = ((bounds[3]-bounds[2]) / 2.0) + 50 + image->GetOrigin()[1];
	line->z = image->GetOrigin()[2];
	line->z2 = image->GetOrigin()[2];


	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->InsertNextPoint(line->x, line->y, line->z);
	points->InsertNextPoint(line->x2, line->y2, line->z2);

	vtkSmartPointer<vtkVertex> v1 = vtkSmartPointer<vtkVertex>::New();
	vtkSmartPointer<vtkVertex> v2 = vtkSmartPointer<vtkVertex>::New();
	v1->GetPointIds()->SetId(0,0);
	v2->GetPointIds()->SetId(0,1);

	vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
	vertices->InsertNextCell(v1);
	vertices->InsertNextCell(v2);

	vtkSmartPointer<vtkLine> pline = vtkSmartPointer<vtkLine>::New();
	pline->GetPointIds()->SetId(0,0);
	pline->GetPointIds()->SetId(1,1);

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	lines->InsertNextCell(pline);
	


	// create the polydata 
	line->poly->SetPoints(points);
	line->poly->SetLines(lines);
	line->poly->SetVerts(vertices);

	line->mapper->SetInputData(line->poly);
	

	line->actor->SetMapper(line->mapper);
	line->actor->GetProperty()->SetLineWidth(2.0);
	line->actor->GetProperty()->SetPointSize(5.0);
	line->actor->GetProperty()->SetRepresentationToWireframe();

	double col[3];
	line->getColour(col);
	line->actor->GetProperty()->SetColor(col);

	return line;

}

// ------------------------------------------------------------------------
Line::Type Line::getTypeEnum(unsigned int btnIndex)
{
	switch(btnIndex)
	{
		case 0:
			return Line::MV;
			break;
		case 1:
			return Line::TP;
			break;
		case 2:
			return Line::AV;
			break;
		default:
			return Line::MV;
			break;
	}
}

// ------------------------------------------------------------------------
std::string Line::getTypeString(Type type)
{
	switch(type)
	{
		case MV:
			return "MV";
			break;
		case TP:
			return "TP";
			break;
		case AV:
			return "AV";
			break;
		default:
			return "";
			break;
	}
}

// ------------------------------------------------------------------------
void Line::getColour(double * col)
{
	col[0] = 0; col[1] = 0; col[2] = 0;

	switch(this->type)
	{
		case MV:
			col[0] = 1.0;
			break;
		case TP:
			col[1] = 1.0;
			break;
		case AV:
			col[2] = 1.0;
			break;
		default:
			break;
	}			
}
