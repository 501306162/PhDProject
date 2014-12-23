#include "CardiacDicomExtractor.h"


namespace utils
{

// ------------------------------------------------------------------------
void CardiacDicomExtractor::Extract()
{
	BuildFlagLists();
	RunExtractor();

	DicomSeriesList::iterator it = m_SeriesList.begin();
	while(it != m_SeriesList.end())
	{
		std::cout << it->description << std::endl;
		++it;
	}

}


// ------------------------------------------------------------------------
void CardiacDicomExtractor::RunExtractor()
{
	SeriesExtractor extractor;
	extractor.SetDirectory(m_FolderName);
	extractor.ExtractSeries();
	m_SeriesList = extractor.GetOutput();
}


// ------------------------------------------------------------------------
void CardiacDicomExtractor::BuildFlagLists()
{
	FlagListType v2cFlags, v3cFlags, v4cFlags, vr2cFlags, vsaFlags;

	v2cFlags.push_back("Cine_ipat 2C");
	v2cFlags.push_back("Cine_2C");
	v2cFlags.push_back("cine_2C");

	m_Flags["2C"] = v2cFlags;

	v3cFlags.push_back("Cine_ipat 3C");
	v3cFlags.push_back("Cine_3C");
	v3cFlags.push_back("cine_3C");

	m_Flags["3C"] = v3cFlags;

	v4cFlags.push_back("Cine_ipat 4C");
	v4cFlags.push_back("Cine_4C");
	v4cFlags.push_back("cine_4C");

	m_Flags["4C"] = v4cFlags;

	vr2cFlags.push_back("Cine_R2C");
	vr2cFlags.push_back("Cine_RV2C");
	vr2cFlags.push_back("cine_RV2C");

	m_Flags["R2C"] = vr2cFlags;

	vsaFlags.push_back("Cine_ipat SA");
	vsaFlags.push_back("Cine_SA");
	vsaFlags.push_back("Cine_SA 1");
	vsaFlags.push_back("cine_SA");
	vsaFlags.push_back("Cine_SA RPT");
	
	m_Flags["SA"] = vsaFlags;

}


}	
