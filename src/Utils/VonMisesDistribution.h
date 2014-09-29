#ifndef VON_MISES_DISTRIBUTION_H
#define VON_MISES_DISTRIBUTION_H

#include <itkListSample.h>

namespace utils
{
class VonMisesDistribution
{
public:
	VonMisesDistribution();
	~VonMisesDistribution() {}

	typedef std::vector<double> DataArray;

	void SetData(const DataArray &data);

	void ComputeMean();
	void ComputeK();

	double ComputeR();

	double Evaluate(double x);

protected: 

private:
	DataArray m_Data;
	double m_Mean;
	double m_K;




};
} /* utils */ 


#endif
