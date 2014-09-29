#include <iostream>

#include <CommonDefinitions.h>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkRenderer.h>
#include <vtkImageActor.h>
#include <vtkImageResliceMapper.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageSlice.h>
#include <vtkImageShiftScale.h>
#include <vtkLookupTable.h>
#include <vtkImageMapToColors.h>

#include <itkImageToVTKImageFilter.h>
#include <itkFlipImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>


using namespace utils;
int main(int argc, char *argv[])
{
	// load an image and a label
	ImageVolume::Pointer image = ImageVolumeIO::Read(argv[1]);
	LabelVolume::Pointer label = LabelVolumeIO::Read(argv[2]);

	// first we need to flip the image and rescale it
	typedef itk::RescaleIntensityImageFilter< ImageVolume, ImageVolume > RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetInput( image );
	rescaler->SetOutputMaximum(255);
	rescaler->SetOutputMinimum(0);
	rescaler->Update();


	typedef itk::FlipImageFilter<ImageVolume> ImageFlipperType;
	typedef ImageFlipperType::FlipAxesArrayType ImageFlipAxes;
	ImageFlipperType::Pointer imageFlipper = ImageFlipperType::New();
	imageFlipper->SetInput(image);

	ImageFlipAxes imageAxes;
	imageAxes.Fill(false);
	imageAxes[1] = true;
	
	imageFlipper->SetFlipAxes(imageAxes);



	typedef itk::FlipImageFilter<LabelVolume> LabelFlipperType;
	LabelFlipperType::Pointer labelFlipper = LabelFlipperType::New();
	labelFlipper->SetInput(label);
	labelFlipper->SetFlipAxes(imageAxes);



	// convert the label and the image to vtk representation
	typedef itk::ImageToVTKImageFilter<ImageVolume> ImageConverter;
	ImageConverter::Pointer imageConverter = ImageConverter::New();
	imageConverter->SetInput(imageFlipper->GetOutput());
	imageConverter->Update();

	
	typedef itk::ImageToVTKImageFilter<LabelVolume> LabelConverter;
	LabelConverter::Pointer labelConverter = LabelConverter::New();
	labelConverter->SetInput(labelFlipper->GetOutput());
	labelConverter->Update();


	vtkSmartPointer<vtkImageResliceMapper> mapper = 
		vtkSmartPointer<vtkImageResliceMapper>::New();
	mapper->SetInputData(imageConverter->GetOutput());
	mapper->SliceFacesCameraOn();
	mapper->SliceAtFocalPointOn();

	vtkSmartPointer<vtkImageSlice> actor = 
		vtkSmartPointer<vtkImageSlice>::New();
	actor->SetMapper(mapper);


	// 
	// Set up the opactiy
	//
	vtkSmartPointer< vtkLookupTable > lookupTable =
		vtkSmartPointer< vtkLookupTable >::New();
	lookupTable->SetNumberOfTableValues(2);
	lookupTable->SetRange(0,1);
	lookupTable->SetTableValue(0,0.0,0.0,0.0,0.0);
	lookupTable->SetTableValue(1,1,0,0,0.5);
	lookupTable->Build();

	vtkSmartPointer< vtkImageMapToColors > imageMapToColors =
		vtkSmartPointer< vtkImageMapToColors >::New();
	imageMapToColors->SetLookupTable( lookupTable );
	imageMapToColors->PassAlphaToOutputOn();
	imageMapToColors->SetInputData( labelConverter->GetOutput() );
	imageMapToColors->Update();


	vtkSmartPointer<vtkImageResliceMapper> mapper2 = 
		vtkSmartPointer<vtkImageResliceMapper>::New();
	mapper2->SetInputData(imageMapToColors->GetOutput());
	mapper2->SliceFacesCameraOn();
	mapper2->SliceAtFocalPointOn();

	vtkSmartPointer<vtkImageSlice> actor2 = 
		vtkSmartPointer<vtkImageSlice>::New();
	actor2->SetMapper(mapper2);


	// Setup renderers
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->AddViewProp(actor);
	//renderer->AddViewProp(actor2);
	renderer->ResetCamera();

	// Setup render window
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetSize(300, 300);
	renderWindow->AddRenderer(renderer);

	// Setup render window interactor
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();

	vtkSmartPointer<vtkInteractorStyleImage> style =
		vtkSmartPointer<vtkInteractorStyleImage>::New();
	style->SetInteractionModeToImageSlicing();

	renderWindowInteractor->SetInteractorStyle(style);

	// Render and start interaction
	renderWindowInteractor->SetRenderWindow(renderWindow);
	renderWindowInteractor->Initialize();

	renderWindowInteractor->Start();
	return 0;
}




