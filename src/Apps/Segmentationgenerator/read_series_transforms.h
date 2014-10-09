#ifndef READ_SERIES_TRANSFORMS_H
#define READ_SERIES_TRANSFORMS_H
#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <QString>
#include <QStringList>


#include "common.h"

typedef std::map<int, std::map<int, int> > LookupMap;
typedef std::vector<double> TranslationType;
typedef std::vector<double> RotationType;
typedef std::pair<TranslationType, RotationType> TransformationType;
typedef std::map<int, TransformationType> TransformationMap;


void readSeriesTransforms(const std::string &registrationFilename, 
		const std::string &lookupFilename, 
		SeriesTransform::Map &transforms);
void readLookUpFile(const std::string &filename, LookupMap &map);
void sortRegistrationFile(const std::string &filename, 
		LookupMap &lookup, TransformationMap &transforms);
void getImageIds(const QString &line, int &series, int &slice);
void getTranslation(const QString &line, TranslationType &trans);
void getRotation(const QString &line, RotationType &rot);

#endif
