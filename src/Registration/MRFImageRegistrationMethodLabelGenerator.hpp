#ifndef MYITKMRFIMAGEREGISTRATIONMETHODLABELGENERATOR_HPP
#define MYITKMRFIMAGEREGISTRATIONMETHODLABELGENERATOR_HPP

#include "MRFImageRegistrationMethodLabelGenerator.h"

namespace myitk
{
// ---------------------------------------------------------------------------------
template<int VDimensions>
MRFImageRegistrationMethodLabelGenerator<VDimensions>
::MRFImageRegistrationMethodLabelGenerator()
{
	// set some defualts
	m_MaxDisplacement = 0.0;
	m_LabelSteps      = 0;
	m_DenseSampling   = false;
}

// ---------------------------------------------------------------------------------
template<int VDimensions>
const typename MRFImageRegistrationMethodLabelGenerator<VDimensions>::LabelListType 
MRFImageRegistrationMethodLabelGenerator<VDimensions>
::GenerateLabels()
{
	// make sure that the dimensions make sense
	if( VDimensions < 2 || VDimensions > 4 )
	{
		itkExceptionMacro( << "Dimensions are strange, should be 2, 3 or 4");
	}

	// make sure the other parameters are set
	if( m_MaxDisplacement <= 0.0 )
	{
		itkExceptionMacro( << "Max Displacement should be larger than 0");
	}

	if( m_LabelSteps <= 0 )
	{
		itkExceptionMacro( << "Label Steps Should be greater than 0");
	}


	// create the output object
	LabelListType labelList;



	// first we generate as if for non dense sampling
	if( m_DenseSampling == false )
	{
		// generate the full range of numbers
		int totalSteps = m_LabelSteps * 2;

		std::vector< std::vector< double > > valueRange(VDimensions);
		for( unsigned int dim = 0; dim < VDimensions; dim++ )
		{
			// find the actual extent of the labels given the spacing of the
			// grid
			double actualExtent = m_MaxDisplacement * m_Spacing[dim];
			double stepSize = actualExtent / static_cast<double>( m_LabelSteps );
			double currentLocation = - actualExtent;
			for( unsigned int i = 0; i < totalSteps+1; i++ )
			{
				if( i != m_LabelSteps )
				{
					valueRange[dim].push_back(currentLocation);
				}
				currentLocation += stepSize;
			}
		}


		unsigned int lcount = 0;

		// find out the total label size, and initailise the list
		int totalLabels = (VDimensions * totalSteps) + 1;
		labelList.resize(totalLabels);

		// the first value will always be the zero displacement
		LabelType zeroLabel(VDimensions);
		zeroLabel.Fill(0.0);
		labelList[lcount] = zeroLabel;
		lcount++;

		// now we build the labels
		for( unsigned int dim = 0; dim < VDimensions; dim++ )
		{
			for( unsigned int l = 0; l < totalSteps; l++ )
			{
				LabelType newLabel(VDimensions);
				newLabel.Fill(0.0);
				newLabel(dim) = valueRange[dim][l];
				labelList[lcount] = newLabel;
				lcount++;
			}
		}
		
	}
	else
	{

		// generate the full range of numbers
		double stepSize = m_MaxDisplacement / static_cast<double>( m_LabelSteps );
		int totalSteps = m_LabelSteps * 2 + 1;


		std::vector< std::vector< double > > valueRange(VDimensions);
		for( unsigned int dim = 0; dim < VDimensions; dim++ )
		{
			// find the actual extent of the labels given the spacing of the
			// grid
			double actualExtent = m_MaxDisplacement * m_Spacing[dim];
			double stepSize = actualExtent / static_cast<double>( m_LabelSteps );
			double currentLocation = - actualExtent;
			for( unsigned int i = 0; i < totalSteps+1; i++ )
			{
				valueRange[dim].push_back(currentLocation);
				currentLocation += stepSize;
			}
		}


		// find out the total size
		int totalLabels = pow( totalSteps, VDimensions );
		labelList.resize( totalLabels );

		unsigned int lcount = 0;





		// fill up the rest of the list with zero labels
		for( unsigned int l = lcount; l < totalLabels; l++ )
		{
			LabelType label(VDimensions);
			label.Fill(0.0);
			labelList[l] = label;
		}

		// now we put in the values
		for( unsigned int dim = 0; dim < VDimensions; dim++ )
		{
			for( unsigned int l = lcount; l < totalLabels; l++ )
			{
				// given the dimension and the label find out what to put where

				// compute the devisor
				double exponent = static_cast<double>( VDimensions-dim-1 );
				int devisor = static_cast<int>( pow( totalSteps, exponent ));
				int ind = (l / devisor) % totalSteps;
				labelList[l](dim) = valueRange[dim][ind];
			}
		}
	}

	return labelList;
}

}

#endif /* End of MYITKMRFIMAGEREGISTRATIONMETHODLABELGENERATOR_HPP */
