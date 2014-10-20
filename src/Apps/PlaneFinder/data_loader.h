#ifndef DATA_LOADER_H
#define DATA_LOADER_H

#include <iostream>
#include <vector>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>

#include <QStringList>

#include "common.h"


class DataLoader
{
public:
	
	typedef std::string FilenameType;
	typedef std::vector<FilenameType> FilenamesType;


	DataLoader() {};
	virtual ~DataLoader() {}

	/**
	 * Set the folder that the data is stored in
	 */
	void SetDataFolder(const std::string &folder)
	{
		this->folderName = folder;
	}

	/**
	 * Load the data from the folder
	 */
	bool LoadData(ImageDataList &imageData);

	/**
	 * Add a filter term to the filename search
	 */
	void AddFilenameFilter(const std::string &filter);

protected: 
	/**
	 * Get the list of image filenames in  the folder
	 */
	void GetFilenames(FilenamesType &filenames);

	void LoadImage(const std::string &fname, ImageData &image);

	std::string GetFileBasename(const std::string &filename);

private:


	QStringList filters;

	std::string folderName;
};

#endif
