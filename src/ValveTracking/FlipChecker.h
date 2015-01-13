#ifndef FLIP_CHECKER_H
#define FLIP_CHECKER_H

#include <itkObject.h>
#include <itkObjectFactory.h>
#include <QVariant>

namespace vt
{
class FlipChecker : public itk::Object
{
public:
	typedef FlipChecker Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	typedef std::vector<int> FlipKeys;
	typedef std::map<std::string, FlipKeys> FlipMap;

	itkTypeMacro(FlipChecker, Object);
	itkNewMacro(Self);

	bool FlipImage(const std::string &type, const int number);
	bool FlipImageFromFileName(const std::string type, const std::string &filename);
	bool FlipPoints(const std::string &type, const int number);


protected:
	FlipChecker();
	virtual ~FlipChecker() {}


private:

	void BuildMap(QVariantMap input, FlipMap &map);


	FlipChecker(const Self&);
	void operator=(const Self&);
	FlipMap m_FlipMap;
	FlipMap m_PointFlipMap;
};


}


#endif



