#ifndef MYITKMRFIMAGEREGISTRATIONMETHODLABELGENERATOR_H
#define MYITKMRFIMAGEREGISTRATIONMETHODLABELGENERATOR_H

#include <itkImage.h>
#include <itkTranslationTransform.h>


namespace myitk
{
/**
 * @brief This is a helper class for the MRFImageRegistrationMethod that
 * facilitates the creation of the label set. It takes as an argument the
 * Bspline Grid and generates the labels based on parameters passed in
 */
template <int VDimensions >
class MRFImageRegistrationMethodLabelGenerator : public itk::Object
{
public:
	typedef MRFImageRegistrationMethodLabelGenerator	Self;
	typedef itk::Object									Superclass;
	typedef itk::SmartPointer< Self >					Pointer;
	typedef itk::SmartPointer< const Self >				ConstPointer;



	/** Method for creation through object factory */
	itkNewMacro( Self );

	/** Runtime type info */
	itkTypeMacro( MRFImageRegistrationMethodLabelGenerator, itk::Object );

	/** typedefs for the labels */
	typedef typename itk::TranslationTransform<double, VDimensions>::ParametersType LabelType;
	typedef std::vector< LabelType > LabelListType;
	typedef itk::Vector<double, VDimensions > SpacingType;

	/** Set / Get for MaxDisplacement */
	itkSetMacro( MaxDisplacement, double );
	itkGetConstMacro( MaxDisplacement, double );

	/** Set / Get for LabelSteps */
	itkSetMacro( LabelSteps, int );
	itkGetConstMacro( LabelSteps, int );

	/** Set / Get for Spacing */
	itkSetMacro( Spacing, SpacingType );
	itkGetConstMacro( Spacing, SpacingType );

	/** Set / Get for DenseSampling */
	itkSetMacro( DenseSampling, bool );
	itkGetConstMacro( DenseSampling, bool );

	/** Method that actually generates the labels */
	const LabelListType  GenerateLabels();


protected:
	MRFImageRegistrationMethodLabelGenerator();
	~MRFImageRegistrationMethodLabelGenerator() {}
	

private:
	MRFImageRegistrationMethodLabelGenerator( const Self & );
	void operator=( const Self & );



	SpacingType m_Spacing;
	double m_MaxDisplacement;
	unsigned int m_LabelSteps;
	bool m_DenseSampling;

};
} // end namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "MRFImageRegistrationMethodLabelGenerator.hpp"
#endif

#endif /* End of MYITKMRFIMAGEREGISTRATIONMETHODLABELGENERATOR_H */
