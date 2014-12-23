#include "ValveHelpers.h"

#include <itkExtractImageFilter.h>

namespace vt
{
// ------------------------------------------------------------------------
void FlattenValve(const ValveLine<3>::Pointer &input, ValveLine<2>::Pointer &output)
{
	// make sure everyhting is up to date
	input->UpdateIndexs();
	
	// extract the 2d image
	typedef itk::ExtractImageFilter<ValveLine<3>::ImageType, ValveLine<2>::ImageType> ExtractorType;
	ExtractorType::Pointer extractor = ExtractorType::New();

	ValveLine<3>::ImageType::RegionType exRegion = input->GetImage()->GetLargestPossibleRegion();
	ValveLine<3>::ImageType::SizeType exSize = exRegion.GetSize();
	ValveLine<3>::ImageType::IndexType exIndex = exRegion.GetIndex();

	exSize[2] = 0;
	exIndex[2] = 0;
	exRegion.SetSize(exSize);
	exRegion.SetIndex(exIndex);

	extractor->SetExtractionRegion(exRegion);
	extractor->SetInput(input->GetImage());
	extractor->SetDirectionCollapseToIdentity();

	extractor->Update();
	output->SetImage(extractor->GetOutput());


	// get the indexes
	ValveLine<2>::ContIndexType ind1, ind2;
	for(unsigned int i = 0; i < 2; i++)
	{
		ind1[i] = input->GetInd1()[i];
		ind2[i] = input->GetInd2()[i];
	}

	output->SetInd1(ind1);
	output->SetInd2(ind2);
	output->UpdatePoints();
}

// ------------------------------------------------------------------------
void FlattenValveSequence(const ValveSequence<3>::Pointer &input, 
		ValveSequence<2>::Pointer &output)
{
	for(unsigned int i = 0; i < input->GetNumberOfLines(); i++)
	{
		ValveLine<2>::Pointer newValve = ValveLine<2>::New();
		FlattenValve(input->GetValveLine(i), newValve);		
		output->AddValveLine(newValve);
	}
}


}
