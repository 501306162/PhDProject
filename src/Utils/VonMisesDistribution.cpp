#include "VonMisesDistribution.h"

#include <boost/math/special_functions/bessel.hpp>

namespace utils
{
// ------------------------------------------------------------------------
VonMisesDistribution::VonMisesDistribution()
{
}

// ------------------------------------------------------------------------
void VonMisesDistribution::SetData(const DataArray &data)
{
	m_Data = data;	
}

// ------------------------------------------------------------------------
void VonMisesDistribution::ComputeMean()
{
	unsigned int numPoints = m_Data.size();
	std::complex<double> angleSum(0,0);
	for(unsigned int i = 0; i < numPoints; i++)
	{
		double val = m_Data[i];
		std::complex<double> a(0,1);
		a *= val;
		angleSum+= std::exp(a);
	}

	m_Mean = std::arg(angleSum);
}

// ------------------------------------------------------------------------
double VonMisesDistribution::Evaluate(double x)
{
	double c = 1.0 / (2*M_PI*boost::math::cyl_bessel_i(0, m_K));
	return c * std::exp(m_K*std::cos(x-m_Mean));
}

// ------------------------------------------------------------------------
void VonMisesDistribution::ComputeK()
{
	double r = ComputeR();

	if(r < 0.53)
		m_K = (2*r) + std::pow(r,3.0) + (5*std::pow(r,5.0) / 6);
	else if(r >= 0.53 && r < 0.85)
		m_K = -0.4 + 1.39*r + 0.43/(1-r);
	else
		m_K = 1/(std::pow(r,3.0) - 4*std::pow(r,2.0) + 3*r);

	if(m_Data.size() < 15 && m_Data.size() > 1)
	{
		if(m_K < 2)
		{	
			m_K = std::max(m_K-2*std::pow((m_Data.size()*m_K),-1.0),0.0);
		}
		else
		{
			m_K = std::pow((m_Data.size()-1),3.0) * m_K / (std::pow( (double) m_Data.size(), 3.0) + m_Data.size());
		}
	}
	std::cout << m_K << std::endl;
}

// ------------------------------------------------------------------------
double VonMisesDistribution::ComputeR()
{
	unsigned int numPoints = m_Data.size();
	std::complex<double> angleSum(0,0);
	for(unsigned int i = 0; i < numPoints; i++)
	{
		double val = m_Data[i];
		std::complex<double> a(0,1);
		a *= val;
		angleSum+= std::exp(a);
	}


	return std::abs(angleSum) / static_cast<double>(numPoints);

}

} /* utils */ 
