#include "ParameterHelpers.h"
#include <Directory.h>
#include <cmath>

namespace vt
{

// -----------------------------------------------------------------------
void PatchExtractorParameters::print()
{
	std::cout << std::endl;
	std::cout << "==================================================" << std::endl;
	std::cout << "Name: " << this->name << std::endl;
	std::cout << "Input Directory: " << this->inputDirectory << std::endl;
	std::cout << "Output Directory: " << this->outputDirectory << std::endl;
	std::cout << "Valve Type: " << this->valveType << std::endl;
	
	std::cout << "Valve Sub Directories: "  << std::endl;
	for(unsigned int i = 0; i < this->valveSubDirectories.size(); i++)
	{
		std::cout << "\t" << this->valveSubDirectories[i] << std::endl;		
	}


	std::cout << "Feature Type: " << this->featureType << std::endl;

	if(this->featureType == "LBP")
	{
		std::cout << "\tLBP grid Size: " << this->lbpGridSize << std::endl;
		std::cout << "\tLBP neighbors: " << this->lbpNeighbors << std::endl;
		std::cout << "\tLBP radius: " << this->lbpRadius << std::endl;
	}

	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "Segmentation Patch Size: " << this->segmentationPatchSize << std::endl;
	std::cout << "Negative Feature Region Size: " << this->negativeFeatureRegionSize << std::endl;
	std::cout << "Feature Patch Size: " << this->featurePatchSize << std::endl;
	std::cout << "Feature Patch Distance: " << this->featurePatchDistance << std::endl;
	std::cout << "Negative Ignore Distance: " << this->negativeIgnoreDistance << std::endl;

	std::cout << "Output Template: " << this->outputFilenameTemplate << std::endl;

	std::cout << "==================================================" << std::endl;
	std::cout << std::endl;
}




// ------------------------------------------------------------------------
PatchExtractorParameters::PatchExtractorParameters(const std::string &filename, bool all)
{

	utils::ConfigParser configs = utils::ConfigParser(filename);
	this->inputDirectory = configs.strValue("input_directory");
	this->outputDirectory = configs.strValue("output_directory");
	this->valveType = configs.strValue("valve_type");
	this->name = configs.strValue("name");
	
	if(all)
	{
		this->valveSubTypes.push_back("TP-R2C");
		this->valveSubTypes.push_back("TP-4C");
	}

	// create the list of valve directories
	if(this->valveType == "MV")
	{		
		this->valveSubTypes.push_back("MV-2C");
		this->valveSubTypes.push_back("MV-3C");
		this->valveSubTypes.push_back("MV-4C");
	}
	else
	{	
		std::cout << "Valve Type: " << this->valveType << " not recognised" << std::endl;
		exit(1);
	}


	for(unsigned int i = 0; i < this->valveSubTypes.size(); i++)
	{
		this->valveSubDirectories.push_back(utils::Directory::GetPath(
					this->inputDirectory, this->valveSubTypes[i]));
	}



	this->featureType = configs.strValue("feature_type");
	if(this->featureType == "LBP")
	{
		this->lbpGridSize = configs.intValue("lbp_grid_size");
		this->lbpNeighbors = configs.intValue("lbp_neighbors");
		this->lbpRadius = configs.doubleValue("lbp_radius");
		this->featureSize = static_cast<int>(std::pow(2.0, static_cast<double>(this->lbpNeighbors)));
		this->featureSize *= std::pow(static_cast<double>(this->lbpGridSize), 2.0);
	}



	this->segmentationPatchSize = configs.intValue("segmentation_patch_size");
	this->negativeFeatureRegionSize = configs.intValue("negative_feature_region_size");
	this->featurePatchSize = configs.intValue("feature_patch_size");
	this->featurePatchDistance = configs.doubleValue("feature_patch_distance", 30.0);
	this->negativeIgnoreDistance = configs.doubleValue("ignore_distance",
		   	static_cast<double>(this->featurePatchDistance / 2.0));
	
	QString templateName = "{Name}-{ValveType}-{ValveSubType}-{PointNumber}.hdf5";
	templateName = templateName.replace("{Name}", this->name.c_str());
	this->outputFilenameTemplate = templateName.toStdString();
	this->numberOfTimeSteps = configs.intValue("num_timesteps");

}

// ------------------------------------------------------------------------
std::string PatchExtractorParameters::getOutputFilename(const unsigned int time, const std::string &type)
{
	std::stringstream sFilename;
	sFilename << this->name << "-" << type << "-" << time << ".hdf5";
	std::string output = utils::Directory::GetPath(this->outputDirectory, sFilename.str()); 
	return output;
}




}
