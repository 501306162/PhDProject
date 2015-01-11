#ifndef VALVE_LINE_H
#define VALVE_LINE_H

#include <iostream>
#include <itkObject.h>
#include <itkObjectFactory.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

#include <itkImage.h>

namespace vt
{
template<int VDimensions>
class ValveLine : public itk::Object
{
public:
	typedef ValveLine Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	typedef itk::ContinuousIndex<double, VDimensions> ContIndexType;

	itkTypeMacro(ValveLine, Object);
	itkNewMacro(Self);

	typedef itk::Image<unsigned short, VDimensions> ImageType;
	typedef typename ImageType::PointType PointType;

	void SetImage(const typename ImageType::Pointer &image) { m_Image = image; }
	typename ImageType::Pointer GetImage() const { return m_Image; }

	PointType GetPoint(const unsigned int index) const
	{
		if(index == 1)
			return m_P1;
		else
			return m_P2;	
	}

	itkSetMacro(P1, PointType);
	itkGetMacro(P1, PointType);
	itkSetMacro(P2, PointType);
	itkGetMacro(P2, PointType);

	ContIndexType GetIndex(const unsigned int index) const 
	{
		if(index == 1)
			return m_Ind1;
		else
			return m_Ind2;	
	}
	
	itkSetMacro(Ind1, ContIndexType);
	itkGetMacro(Ind1, ContIndexType);
	itkSetMacro(Ind2, ContIndexType);
	itkGetMacro(Ind2, ContIndexType);


	void UpdateIndexs();
	void UpdatePoints();

	vtkSmartPointer<vtkPolyData> GetPolyData() const;

protected:

	ValveLine() {};
	virtual ~ValveLine() {}

private:	

	ValveLine(const Self&);
	void operator=(const Self&);

	std::string m_ImageType;
	std::string m_LineType;
	PointType m_P1;
	PointType m_P2;

	ContIndexType m_Ind1;
	ContIndexType m_Ind2;

	typename ImageType::Pointer m_Image;

};


template<int VDimensions>
class ValveSequence : public itk::Object 
{
public:
	typedef ValveSequence Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(ValveSequence, Object);
	itkNewMacro(Self);

	typedef ValveLine<VDimensions> ValveLineType;
	typedef typename ValveLineType::ImageType ImageType;

	void AddValveLine(const typename ValveLineType::Pointer &valve) { m_ValveLines.push_back(valve); }
	void SetValveLine(const typename ValveLineType::Pointer &valve, 
			const unsigned int &index) { m_ValveLines[index] = valve; }
	typename ValveLineType::Pointer GetValveLine(const unsigned int &index) const { return m_ValveLines[index]; }
	unsigned int GetNumberOfLines() const { return m_ValveLines.size(); }



protected:
	ValveSequence() {}
	virtual ~ValveSequence() {}

	

private:

	ValveSequence(const Self&);
	void operator=(const Self&);

	std::vector<typename ValveLineType::Pointer> m_ValveLines;

};




}


#endif
