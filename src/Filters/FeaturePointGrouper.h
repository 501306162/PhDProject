#ifndef FEATURE_POINT_GROUPER_H
#define FEATURE_POINT_GROUPER_H

#include <itkProcessObject.h>
#include "FeatureCommon.h"

namespace filter
{
template<int TDimension>
class FeaturePointGrouper : public itk::ProcessObject
{
public:
	typedef FeaturePointGrouper Self;
	typedef itk::ProcessObject Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	itkTypeMacro(FeaturePointGrouper, ProcessObject);
	itkNewMacro(Self);

	typedef HoGFeature<TDimension> FeatureType;
	typedef std::vector<FeatureType> FeatureListType;
	typedef HoGFeatureCluster<TDimension> FeatureGroup;
	typedef std::vector<FeatureGroup> FeatureGroupList;

	void SetInput(const FeatureListType &input);


	itkSetMacro(LocationThreshold, double);
	itkSetMacro(AngleThreshold, double);
	itkSetMacro(ScaleThreshold, double);
	itkSetMacro(AdjustLocationByScale, bool);
	itkSetMacro(GroupSizeThreshold, unsigned int);
	FeatureGroupList GetOutput() const;

	void Update();

protected:
	FeaturePointGrouper();
	~FeaturePointGrouper(){}


	bool AreSimilarEnough(const FeatureType &f1, const FeatureType &f2);

private:
	double m_LocationThreshold;
	double m_AngleThreshold;
	double m_ScaleThreshold;
	bool m_AdjustLocationByScale;
	unsigned int m_GroupSizeThreshold;

	FeatureListType m_Input;
	FeatureGroupList m_Output;

	FeaturePointGrouper(const Self&);
	void operator=(const Self&);
	

};
} /* filter */ 

#include "FeaturePointGrouper.hpp"

#endif
