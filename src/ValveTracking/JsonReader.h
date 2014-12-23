#ifndef JSON_READER_H
#define JSON_READER_H

#include <iostream>
#include <itkObject.h>
#include <itkObjectFactory.h>


namespace vt
{
class JsonReader : public itk::Object
{
public:
	typedef JsonReader Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	typedef struct line_
	{
		int t;
		double p1[3];
		double p2[3];
		std::string imageName;
	} Line;


	typedef struct line_sequence_
	{
		std::string lineType;
		std::string imageType;
		std::vector<Line> lines;
		std::string imageFilename;
	} LineSequence;

	typedef std::map<std::string, LineSequence> LineSet;


	itkTypeMacro(JsonReader, Object);
	itkNewMacro(Self);

	itkSetMacro(Filename, std::string);
	itkGetMacro(Filename, std::string);


	LineSet Read() const;

	

protected:
	JsonReader();
	virtual ~JsonReader() {}

private:

	std::string m_Filename;

	JsonReader(const Self&);
	void operator=(const Self&);




};


}


#endif
