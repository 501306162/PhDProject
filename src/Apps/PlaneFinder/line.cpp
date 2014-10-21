#include "line.h"

#include <vtkLineSource.h>
#include <vtkUnsignedCharArray.h>
#include <vtkCellData.h>
#include <vtkProperty.h>
#include <vtkVertexGlyphFilter.h>

// ------------------------------------------------------------------------
Line::Line()
{
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
	line->z = bounds[4] + image->GetOrigin()[2];
	line->z = bounds[5] + image->GetOrigin()[2];



	// create the polydata 
	vtkSmartPointer<vtkLineSource> lineSource = 
		vtkSmartPointer<vtkLineSource>::New();
	lineSource->SetPoint1(line->x, line->y, line->z);
	lineSource->SetPoint2(line->x2, line->y2, line->z2);
	lineSource->Update();

	vtkSmartPointer<vtkVertexGlyphFilter> glypher = 
		vtkSmartPointer<vtkVertexGlyphFilter>::New();
	glypher->SetInputData(lineSource->GetOutput());
	glypher->Update();

	line->poly = vtkSmartPointer<vtkPolyData>::New();
	line->poly->DeepCopy(glypher->GetOutput());

	line->poly->Print(std::cout);
	line->mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	line->mapper->SetInputData(line->poly);
	

	line->actor = vtkSmartPointer<vtkActor>::New();
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
Line::Type Line::getType(unsigned int btnIndex)
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
