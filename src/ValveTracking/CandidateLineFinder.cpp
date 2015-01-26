#include "CandidateLineFinder.h"

namespace vt
{
bool cost_sort(const CandidateLineFinder::OutputPointCostPair &c1, const CandidateLineFinder::OutputPointCostPair &c2)
{
	return (c1.first<c2.first);
}


// ------------------------------------------------------------------------
void CandidateLineFinder::Compute()
{
	// iterate through the lsit of points and find out the cost for each 
	for(unsigned int i = 0; i < m_P1List.size(); i++)
	{
		for(unsigned int j = 0; j < m_P2List.size(); j++)
		{
			double cost = LineCost(m_P1List[i], m_P2List[j]);
			OutputPointPair points;
			points.first = m_P1List[i].second;
			points.second = m_P2List[j].second;

			OutputPointCostPair costs;
			costs.first = cost;
			costs.second = points;

			m_OutputList.push_back(costs);
		}
	}

	std::sort(m_OutputList.begin(), m_OutputList.end(), cost_sort);
}


// ------------------------------------------------------------------------
double CandidateLineFinder::LineCost(const PointPairType &p1, const PointPairType &p2)
{
	double p1Prob = p1.first;
	double p2Prob = p2.first;
	double lineDist = p1.second.EuclideanDistanceTo(p2.second);
	double lengthProb = m_Lengths->Predict(m_Type, m_TimeStep, lineDist);

	// compute the angle difference between the input line and the new line 
	VectorType v1,v2,ov1,ov2;
	for(unsigned int i = 0; i < 3; i++)
	{
		v1(i) = p1.second[i];
		v2(i) = p2.second[i];
		ov1(i) = m_OriginalP1[i];
		ov2(i) = m_OriginalP2[i];
	}

	
	// get the direction vectors 
	VectorType d1 = (v2-v1);
	VectorType d2 = (ov2-ov1);
	d1.normalize();
	d2.normalize();

	double angle = std::acos((d1.dot(d2)) / (d1.norm() * d2.norm()));

	double vecDiff = (d2-d1).norm();
	
	double cost = -log(p1Prob) - log(p2Prob) - log(lengthProb) + vecDiff;
	
	return cost;

}

// ------------------------------------------------------------------------
void CandidateLineFinder::GetOutput(std::vector<OutputPointPair> &output, const unsigned int num)
{
	for(unsigned int i = 0; i < num; i++)
	{
		output.push_back(m_OutputList[i].second);
	}
}

// ------------------------------------------------------------------------
CandidateLineFinder::OutputPointPair CandidateLineFinder::GetWeightedOutput()
{
	double maxCost = 0.0;
	double totalCost = 0.0;
	for(unsigned int i = 0; i < m_OutputList.size(); i++)
	{
		double cost = m_OutputList[i].first;
		if(cost > maxCost)
			maxCost = cost;

		totalCost += cost;
	}


	OutputPointPair output;
	output.first.Fill(0);
	output.second.Fill(0);

	for(unsigned int i = 0; i < m_OutputList.size(); i++)
	{
		double cost = m_OutputList[i].first;
		double weight = (maxCost - cost) / totalCost;

		for(unsigned int j = 0; j < 3; j++)
		{
			output.first[j] += weight * m_OutputList[i].second.first[j];
			output.second[j] += weight * m_OutputList[i].second.second[j];
		}
	}


	return output;

}


}
