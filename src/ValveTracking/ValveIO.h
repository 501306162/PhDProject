#ifndef VALVE_IO_H
#define VALVE_IO_H

#include "ValveLine.h"

#include <itkImage.h>
#include <itkObject.h>
#include <itkObjectFactory.h>


namespace vt
{
template<int VDimensions>
class ValveWriter : public itk::Object
{
public:
	typedef ValveWriter Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(ValveWriter, Object);
	itkNewMacro(Self);

	typedef ValveLine<VDimensions> ValveLineType;

	void SetFileName(const std::string &filename) { m_Filename = filename; }
	void Write();
	void SetInput(const typename ValveLineType::Pointer &valve) { m_Valve = valve; }


protected:
	ValveWriter() {}
	virtual ~ValveWriter() {}


private:

	std::string m_Filename;
	typename ValveLineType::Pointer m_Valve;

	ValveWriter(const Self&);
	void operator=(const Self&);

};	

// ------------------------------------------------------------------------
template<int VDimensions>
class ValveReader : public itk::Object
{
public:
	typedef ValveReader Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(ValveReader, Object);
	itkNewMacro(Self);


	typedef ValveLine<VDimensions> ValveLineType;

	void SetFileName(const std::string &filename) { m_Filename = filename; }
	void Read();
	typename ValveLineType::Pointer GetOutput() const;
protected:
	ValveReader() {}
	virtual ~ValveReader() {}

private:
	std::string m_Filename;
	typename ValveLineType::Pointer m_Valve;

	ValveReader(const Self&);
	void operator=(const Self&);
};	


// ------------------------------------------------------------------------
template<int VDimensions>
class ValveSequenceWriter : public itk::Object
{
public:
	typedef ValveSequenceWriter Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(ValveSequenceWriter, Object);
	itkNewMacro(Self);

	typedef ValveSequence<VDimensions> ValveSequenceType;

	void Write();
	void SetFileName(const std::string &filename) { m_Filename = filename; }
	void SetInput(const typename ValveSequenceType::Pointer &input) { m_Sequence = input; }

protected:
	ValveSequenceWriter() {}
	virtual ~ValveSequenceWriter() {}

private:

	typename ValveSequenceType::Pointer m_Sequence;

	std::string m_Filename;

	ValveSequenceWriter(const Self&);
	void operator=(const Self&);

};


// ------------------------------------------------------------------------
template<int VDimensions>
class ValveSequenceReader : public itk::Object
{
public:
	typedef ValveSequenceReader Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(ValveSequenceReader, Object);
	itkNewMacro(Self);

	typedef ValveSequence<VDimensions> ValveSequenceType;

	void Write();
	void SetFileName(const std::string &filename) { m_Filename = filename; }
	typename ValveSequenceType::Pointer GetOutput() const;

protected:
	ValveSequenceReader() {}
	virtual ~ValveSequenceReader() {}

private:

	typename ValveSequenceType::Pointer m_Sequence;

	std::string m_Filename;

	ValveSequenceReader(const Self&);
	void operator=(const Self&);

};



}

#endif
