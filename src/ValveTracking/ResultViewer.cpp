#include "ResultViewer.h"

#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkPNGWriter.h>
#include <vtkWindowToImageFilter.h>
#include <itkImageToVTKImageFilter.h>
#include <itkFlipImageFilter.h>
#include <vtkImageData.h>
#include <vtkImageSlice.h>
#include <vtkImageSliceMapper.h>
#include <vtkRenderer.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkPlaneSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkMatrix4x4.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkVolumeProperty.h>
#include <LineIntersectionFinder.h>
#include <vtkLineSource.h>
#include <vtkProperty.h>
#include <vtkImageProperty.h>
#include <vtkCamera.h>
#include <itkResampleImageFilter.h>
#include <vtkVertexGlyphFilter.h>

namespace vt
{
// ------------------------------------------------------------------------
void ResultViewer::SetUp()
{
	m_RenderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	m_Interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

	unsigned int numViewports = m_Images.size();
	double xSplit = 1.0 / static_cast<double>(numViewports);
	m_Interactor->SetRenderWindow(m_RenderWindow);

	m_RenderWindow->SetSize(400*m_Images.size(), 500);
	
	//vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
	//m_Interactor->SetInteractorStyle(style);

	for(unsigned int i = 0; i < numViewports; i++)
	{
		{
			vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
			renderer->SetViewport(xSplit*i,0,xSplit*i+xSplit,0.5);		





			vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
			GetImage(m_Images[i], image);

			////GetLine(m_Images[i], m_Point, m_Normal, image,  line);

			//vtkSmartPointer<vtkPolyData> startLine = vtkSmartPointer<vtkPolyData>::New();
			////GetLine(m_Images[i], m_StartPoint, m_StartNormal, image,  startLine);
			//GetPlane(m_Point, m_Normal, m_Images[i], line);
			////GetPlane(m_StartPoint, m_StartNormal, m_Images[i], startLine);

			vtkSmartPointer<vtkImageSliceMapper> imapper = vtkSmartPointer<vtkImageSliceMapper>::New();
			imapper->SetInputData(image);
			vtkSmartPointer<vtkImageSlice> iactor = vtkSmartPointer<vtkImageSlice>::New();
			iactor->SetMapper(imapper);

			iactor->GetProperty()->SetColorLevel(177.48);
			iactor->GetProperty()->SetColorWindow(512.04);


			vtkSmartPointer<vtkPolyData> line = vtkSmartPointer<vtkPolyData>::New();
			vtkSmartPointer<vtkPoints> inP = vtkSmartPointer<vtkPoints>::New();
			for(unsigned int pn = 0; pn < m_Points[i].size(); pn++)
			{
				inP->InsertNextPoint(m_Points[i][pn].first.GetDataPointer());
				inP->InsertNextPoint(m_Points[i][pn].second.GetDataPointer());
			}

			line->SetPoints(inP);

			vtkSmartPointer<vtkVertexGlyphFilter> filter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
			filter->SetInputData(line);
			filter->Update();


			vtkSmartPointer<vtkPolyDataMapper> pmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			vtkSmartPointer<vtkActor> pactor = vtkSmartPointer<vtkActor>::New();
			pmapper->SetInputData(filter->GetOutput());
			pactor->SetMapper(pmapper);
			pactor->GetProperty()->SetLineWidth(2.0);
			pactor->GetProperty()->SetPointSize(3.0);
			pactor->GetProperty()->SetColor(0,1,0);


			//vtkSmartPointer<vtkPolyDataMapper> pmapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
			//vtkSmartPointer<vtkActor> pactor2 = vtkSmartPointer<vtkActor>::New();
			//pmapper2->SetInputData(startLine);
			//pactor2->SetMapper(pmapper2);
			//pactor2->GetProperty()->SetLineWidth(2.0);
			//pactor2->GetProperty()->SetColor(1,0,0);


			// apply transform to the image actor
			ImageType::DirectionType d=m_Images[i]->GetDirection();
			vtkMatrix4x4 *mat=vtkMatrix4x4::New(); //start with identity matrix
			for (int j=0; j<3; j++)
				for (int k=0; k<3; k++)
					mat->SetElement(j,k, d(j,k));

			//counteract the built-in translation by origin
			ImageType::PointType origin=m_Images[i]->GetOrigin();
			iactor->SetPosition(-origin[0], -origin[1], -origin[2]);

			//add translation to the user matrix
			for (int j=0; j<3; j++)
				mat->SetElement(j,3, origin[j]);
			iactor->SetUserMatrix(mat);	

			vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();

			ImageType::Pointer im = m_Images[i];
			itk::ContinuousIndex<double, 3> contIndex;
			for(unsigned int i = 0; i < 2; i++)
			{
				contIndex[i] = ((double) im->GetLargestPossibleRegion().GetSize()[i]) / 2.0;
			}
			contIndex[2] = 0.0;


			ImageType::PointType p;
			im->TransformContinuousIndexToPhysicalPoint(contIndex, p);


			camera->SetFocalPoint(p.GetDataPointer());

			double pos[3];
			double up[3];
			double * linepos = pactor->GetPosition();
			for(unsigned int k = 0; k < 3; k++)
			{
				pos[k] = p[k] + (-400.0*mat->GetElement(k,2));			
				up[k] = -mat->GetElement(k,1);
				linepos[k] += -10.0 * mat->GetElement(k,2);
			}


			pactor->SetPosition(linepos);
			camera->SetPosition(pos);
			camera->SetViewUp(up);

			//renderer->AddActor(pactor2);
			renderer->AddActor(iactor);
			renderer->AddActor(pactor);
			renderer->SetActiveCamera(camera);
			m_RenderWindow->AddRenderer(renderer);
		}


		{

			vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
			renderer->SetViewport(xSplit*i,0.5,xSplit*i+xSplit,1.0);		





			vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
			GetImage(m_Images[i], image);

			////GetLine(m_Images[i], m_Point, m_Normal, image,  line);

			//vtkSmartPointer<vtkPolyData> startLine = vtkSmartPointer<vtkPolyData>::New();
			////GetLine(m_Images[i], m_StartPoint, m_StartNormal, image,  startLine);
			//GetPlane(m_Point, m_Normal, m_Images[i], line);
			////GetPlane(m_StartPoint, m_StartNormal, m_Images[i], startLine);

			vtkSmartPointer<vtkImageSliceMapper> imapper = vtkSmartPointer<vtkImageSliceMapper>::New();
			imapper->SetInputData(image);
			vtkSmartPointer<vtkImageSlice> iactor = vtkSmartPointer<vtkImageSlice>::New();
			iactor->SetMapper(imapper);

			iactor->GetProperty()->SetColorLevel(177.48);
			iactor->GetProperty()->SetColorWindow(512.04);


			vtkSmartPointer<vtkPolyData> line = vtkSmartPointer<vtkPolyData>::New();
			vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
			vtkSmartPointer<vtkPoints> inP = vtkSmartPointer<vtkPoints>::New();

			//PointType p1Mean;
			//p1Mean.Fill(0.0);
			//PointType p2Mean;
			//p2Mean.Fill(0.0);
			

			//double W = 0.0;
			//for(unsigned int wn = 0; wn < m_Weights[i].size(); wn++)
			//{
				//W += m_Weights[i][wn];				
			//}

			//for(unsigned int pn = 0; pn < m_Points[i].size(); pn++)
			//{
				//for(unsigned int n = 0; n < 3; n++)
				//{
					//double w = m_Weights[i][pn] / W;
					//p1Mean[n] += (m_Points[i][pn].first[n] * w);
					//p2Mean[n] += (m_Points[i][pn].second[n] * w);
				//}
			//}

			lineSource->SetPoint1(m_Res[i].first.GetDataPointer());
			lineSource->SetPoint2(m_Res[i].second.GetDataPointer());
			lineSource->Update();



			vtkSmartPointer<vtkPolyDataMapper> pmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			vtkSmartPointer<vtkActor> pactor = vtkSmartPointer<vtkActor>::New();
			pmapper->SetInputData(lineSource->GetOutput());
			pactor->SetMapper(pmapper);
			pactor->GetProperty()->SetLineWidth(2.0);
			pactor->GetProperty()->SetPointSize(3.0);
			pactor->GetProperty()->SetColor(0,0,1);



			vtkSmartPointer<vtkLineSource> gt = vtkSmartPointer<vtkLineSource>::New();
			gt->SetPoint1(m_GT[i].first.GetDataPointer());
			gt->SetPoint2(m_GT[i].second.GetDataPointer());
			gt->Update();

			vtkSmartPointer<vtkPolyDataMapper> pmapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
			vtkSmartPointer<vtkActor> pactor2 = vtkSmartPointer<vtkActor>::New();
			pmapper2->SetInputData(gt->GetOutput());
			pactor2->SetMapper(pmapper2);
			pactor2->GetProperty()->SetLineWidth(2.0);
			pactor2->GetProperty()->SetColor(1,0,0);


			// apply transform to the image actor
			ImageType::DirectionType d=m_Images[i]->GetDirection();
			vtkMatrix4x4 *mat=vtkMatrix4x4::New(); //start with identity matrix
			for (int j=0; j<3; j++)
				for (int k=0; k<3; k++)
					mat->SetElement(j,k, d(j,k));

			//counteract the built-in translation by origin
			ImageType::PointType origin=m_Images[i]->GetOrigin();
			iactor->SetPosition(-origin[0], -origin[1], -origin[2]);

			//add translation to the user matrix
			for (int j=0; j<3; j++)
				mat->SetElement(j,3, origin[j]);
			iactor->SetUserMatrix(mat);	

			vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();

			ImageType::Pointer im = m_Images[i];
			itk::ContinuousIndex<double, 3> contIndex;
			for(unsigned int i = 0; i < 2; i++)
			{
				contIndex[i] = ((double) im->GetLargestPossibleRegion().GetSize()[i]) / 2.0;
			}
			contIndex[2] = 0.0;


			ImageType::PointType p;
			im->TransformContinuousIndexToPhysicalPoint(contIndex, p);


			camera->SetFocalPoint(p.GetDataPointer());

			double pos[3];
			double up[3];
			double * linepos = pactor->GetPosition();
			for(unsigned int k = 0; k < 3; k++)
			{
				pos[k] = p[k] + (-400.0*mat->GetElement(k,2));			
				up[k] = -mat->GetElement(k,1);
				linepos[k] += -10.0 * mat->GetElement(k,2);
			}


			pactor->SetPosition(linepos);
			camera->SetPosition(pos);
			camera->SetViewUp(up);

			renderer->AddActor(pactor2);
			renderer->AddActor(iactor);
			renderer->AddActor(pactor);
			renderer->SetActiveCamera(camera);
			m_RenderWindow->AddRenderer(renderer);

		}

	}
}


// ------------------------------------------------------------------------
void ResultViewer::View()
{
	m_RenderWindow->Render();
	m_Interactor->Start();
}

// ------------------------------------------------------------------------
void ResultViewer::Save(const std::string &filename)
{

	m_RenderWindow->OffScreenRenderingOn();
	m_RenderWindow->Render();
	vtkSmartPointer<vtkWindowToImageFilter> saver = vtkSmartPointer<vtkWindowToImageFilter>::New();
	saver->SetInput(m_RenderWindow);
	saver->SetInputBufferTypeToRGBA();
	saver->SetReadFrontBuffer(false);
	saver->Update();

	vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputConnection(saver->GetOutputPort());
	writer->Write();
	

}

// ------------------------------------------------------------------------
void ResultViewer::GetLine(const ImageType::Pointer &image, VectorType &point, VectorType &normal, 
			vtkSmartPointer<vtkImageData> &vtkImage,
		   	vtkSmartPointer<vtkPolyData> &line)
{
	LineIntersectionFinder::Pointer finder = LineIntersectionFinder::New();

	itk::Point<double, 3> tp;
	itk::CovariantVector<double, 3> tn;

	for(unsigned int i = 0; i < 3; i++)
	{
		tp[i] = point(i);
		tn[i] = normal(i);
	}

	tp = m_Transform->GetInverseTransform()->TransformPoint(tp);
	tn = m_Transform->GetInverseTransform()->TransformCovariantVector(tn);

	VectorType np, nn;
	for(unsigned int i = 0; i < 3; i++)
	{
		np(i) = tp[i];
		nn(i) = tn[i];
	}

	finder->SetImage(image);
	finder->SetPlane(nn, np);

	vtkBoundingBox newBox;


	finder->SetBoundingBox(m_BoundingBox);
	finder->Compute2();	


	LineIntersectionFinder::OutputLineType tmpLine = finder->GetOutput();

	tmpLine.p1 = m_Transform->TransformPoint(tmpLine.p1);
	tmpLine.p2 = m_Transform->TransformPoint(tmpLine.p2);

	vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
	lineSource->SetPoint1(tmpLine.p1.GetDataPointer());
	lineSource->SetPoint2(tmpLine.p2.GetDataPointer());
	lineSource->Update();
	line->DeepCopy(lineSource->GetOutput());

}


// ------------------------------------------------------------------------
void ResultViewer::GetImage(const ImageType::Pointer &in,vtkSmartPointer<vtkImageData> &image)
{

	typedef itk::ImageToVTKImageFilter<ImageType> ConverterType;
	ConverterType::Pointer converter = ConverterType::New();
	converter->SetInput(in);
	converter->Update();

	image->DeepCopy(converter->GetOutput());
		

}


// ------------------------------------------------------------------------
void ResultViewer::GetPlane(VectorType &point, VectorType &normal, 
		ImageType::Pointer &image, vtkSmartPointer<vtkPolyData> &outpoly)
{
	vtkSmartPointer<vtkPlaneSource> source = vtkSmartPointer<vtkPlaneSource>::New();
	source->SetCenter(point(0), point(1), point(2));
	source->SetNormal(normal(0), normal(1), normal(2));
	source->SetResolution(1,1);
	source->Update();

	vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
	poly->DeepCopy(source->GetOutput());

	vtkSmartPointer<vtkPoints> p = poly->GetPoints();
	double cent[3];
	cent[0] = 0.0; cent[1] = 0.0; cent[2] = 0.0;

	for(unsigned int j = 0; j < p->GetNumberOfPoints(); j++)
	{
		cent[0] += p->GetPoint(j)[0] / (double) p->GetNumberOfPoints();
		cent[1] += p->GetPoint(j)[1] / (double) p->GetNumberOfPoints();
		cent[2] += p->GetPoint(j)[2] / (double) p->GetNumberOfPoints();
	}

	for(unsigned int j = 0; j < p->GetNumberOfPoints(); j++)
	{
		double pin[3];
		p->GetPoint(j,pin);		

		for(unsigned int k = 0; k < 3; k++)
		{
			pin[k] -= cent[k];
			pin[k] *= 100;
			pin[k] += cent[k];
		}

		p->SetPoint(j,pin);

	}


	vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
	cutter->SetInputData(poly);

	vtkSmartPointer<vtkPlane> cutPlane = vtkSmartPointer<vtkPlane>::New();
	cutPlane->SetNormal(image->GetDirection()(0,2), image->GetDirection()(1,2), image->GetDirection()(2,2));
	cutPlane->SetOrigin(image->GetOrigin()[0], image->GetOrigin()[1], image->GetOrigin()[2]);
	cutter->SetCutFunction(cutPlane);

	cutter->Update();
	outpoly= vtkSmartPointer<vtkPolyData>::New();
	outpoly->DeepCopy(cutter->GetOutput());

}

}
