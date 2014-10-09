#include "VolumeData.h"

using namespace std;
using namespace cv;

VolumeData::VolumeData(XMLreader& _xml_reader): xml_reader(_xml_reader)
{
	cout << "constructor VolumeData" << endl;

	volume = NULL;

	timeOfFrame = NULL;

	translationVect = NULL;
	rotationAngles = NULL;
	pixelCoordsInVolume = NULL;
	voxelCoordsInImageVolume = NULL;

	mat_serieRef = MAT_NULL;
	mat_serieRef_inv = MAT_NULL;
	offset_serieRef = VECT_NULL;

	imagePresence = NULL;
}

void VolumeData::loadData()
{
	repertory = xml_reader.get_repertory();

	periodic_motion = xml_reader.is_motion_periodic();

	try
	{
		listSeries();
		printListedSeries();
		selectSeries();
	}
	catch(string message)
	{
		cout << message << endl;
		exit(-1);
	}

	/// read series ///

	set<int>::iterator lit(selectedSeries.begin()), lend(selectedSeries.end()), lit_tmp;

	for(;lit != lend;)
	{
		try
		{
			readSerie(*lit);
			++lit;
		}
		catch(string message)
		{
			cout << message << endl;

			lit_tmp = lit;
			++lit;

			cout << "Discarding serie " << *lit_tmp << endl;
			selectedSeries.erase(lit_tmp);

			if (selectedSeries.empty())
				throw string("VolumeData() error: no readable serie");
		}
	}

	try
	{
		set<int>::iterator lit(selectedSeries.begin()), lend(selectedSeries.end());

		for(; lit != lend; )
		{
			int current_serieNumber = *lit;
			check_serie(lit);

			if (current_serieNumber == *lit) ++lit; //else, the serie has been erased and lit incremented
		}

		choose_serie_ref();

		/// build list of used images, for convenience

		for(lit = selectedSeries.begin(); lit != lend; ++lit)
		{
			cout << *lit << endl;
			serieFirstImageIndex[*lit] = selectedImages.size();

			for (int s=0; s<seriesNbSlices[*lit]; s++)
				selectedImages.push_back( make_pair(*lit, s) );
		}

		nbImages = selectedImages.size();

		/// place series in our coordinate system ///

		computeVolumeCoords();

		/// build and fill volume

		buildListTimeframes();
		buildVolume();

		loadRegistration();

		fillVolume();
		display3D(0);
	}
	catch(string message)
	{
		cout << message << endl;
		exit(-1);
	}
}

VolumeData::~VolumeData()
{
	if (imagePresence != NULL) delete[] imagePresence;

	map< pair<int, int>, Point3<float>* >::iterator mit(seriesDataVolumeCoords.begin()), mend(seriesDataVolumeCoords.end());
	for (; mit != mend; ++mit) delete[] mit->second;

	if ( volume != NULL )
	{
		for (int f=0; f<nbFrames; f++) delete[] volume[f];
		delete[] volume;
	}

	if (translationVect != NULL) delete[] translationVect;
	if (rotationAngles != NULL) delete[] rotationAngles;

	if (pixelCoordsInVolume != NULL) delete[] pixelCoordsInVolume;
	if (voxelCoordsInImageVolume != NULL) delete[] voxelCoordsInImageVolume;

	if (timeOfFrame != NULL)
	{
		for (int i=0; i<nbImages; i++) delete[] timeOfFrame[i];
		delete[] timeOfFrame;
	}
}

void VolumeData::listFilenames(DIR *rep, list<string>& filenames, string path)
{
	dirent *file;

	while ((file = readdir(rep)))
	{
		if (file->d_name[0] == '.') continue;

		string name = file->d_name;

		if (file->d_type == 8 && name.size() >= 4)
		{
			if ( ! name.compare(name.size()-4, 4, ".png") ) continue;
			if ( ! name.compare(name.size()-4, 4, ".jpg") ) continue;
			if ( ! name.compare(name.size()-4, 4, ".txt") ) continue;
			if ( ! name.compare(name.size()-4, 4, ".EXE") ) continue;
			if ( ! name.compare(name.size()-4, 4, ".HTM") ) continue;
		}
		if (file->d_type == 8 && name.size() >= 5)
			if ( ! name.compare(name.size()-5, 5, ".html") ) continue;
		if (file->d_type == 8 && name.size() >= 12)
			if ( ! name.compare(0, 12, "distribution") ) continue;
		if (file->d_type == 8 && name.size() >= 7)
			if ( ! name.compare(0, 7, "DIRFILE") ) continue;

		if (file->d_type == 4)
		{
			//enter the directory and list DICOM files
			//cout << "ouverture repertoire " << chemin + name << endl;

			DIR *rep2 = opendir ((path + name).c_str());
			if (rep2 != NULL)
			{
				listFilenames(rep2, filenames, path + name + "/");

				closedir (rep2);
				rep2 = NULL;
			}
			else
				throw(string("Error: cannot open the directory"));
		}
		else if (file->d_type == 8)
		{
			filenames.push_back(path + name);
		}
		else
		{
			cout << name << " " << file->d_type << endl;
			throw string("VolumeData::listFilenames() Error : unknown type");
		}
	}
}

void VolumeData::printListedSeries()
{
	cout << "listed series:" << endl;

	cout << left;
	cout << setw(15) << "Serie number" << setw(25) << "Serie description" << endl;
	cout << setw(15) << "Nb of images" << setw(25) << "Images size (rows x cols)" << "Pixel spacing (rows x cols)" << endl;
	cout << endl;

	int listeSeriesSize = seriesDescription.size();
	for (int n=0; n<listeSeriesSize; n++)
	{
		cout << setw(15) << n << setw(25) << seriesDescription[n] << endl;
		cout << setw(15) << seriesNbSlices[n] * nbFrames << setw(25) << seriesPixelSpacing[n].first << "x" << seriesPixelSpacing[n].second << endl;
		cout << endl;
	}
}

void VolumeData::queryChoiceSeries()
{
	cout << "Please, choose series. Enter -1 when done." << endl;

	int choice;
	bool stop = false;

	while (!stop)
	{
		cin >> choice;

		if (choice >= 0 && choice < (int)seriesImageSizes.size())
			selectedSeries.insert(choice);
		else if (choice == -1)
			stop = true;
		else
			cout << "Please, choose a serie between 0 and " << seriesImageSizes.size() - 1 << " or enter -1 to finish" << endl;
	}
}

void VolumeData::selectSeries()
{
	//cout << "VolumeData::selectSeries()" << endl;

	if (selectedSeries.empty()) queryChoiceSeries();

	printSelectedSeries();

	cout << "Enter v to validate this choice, or m to modify the selection" << endl;

	char choice;
	bool stop = false;
	while (!stop)
	{
		cin >> choice;

		if (choice == 'v')
		{
			cout << "Choice validated" << endl;
			stop = true;
		}
		else if (choice == 'm')
		{
			selectedSeries.clear();

			queryChoiceSeries();
		}
		else
		{
			cout << "Please, choose between v and m" << endl;
		}
	}
}

void VolumeData::printSelectedSeries()
{
	cout << "selected series:" << endl;

	cout << left;
	cout << setw(15) << "Serie number" << setw(25) << "Serie description" << endl;
	cout << setw(15) << "Nb of images" << setw(25) << "Images size (rows x cols)" << "Pixel spacing (rows x cols)" << endl;
	cout << endl;

	set<int>::const_iterator lit(selectedSeries.begin()), lend(selectedSeries.end());

	for (; lit != lend; ++lit)
	{
		cout << setw(15) << *lit << setw(25) << seriesDescription[*lit] << endl;
		cout << setw(15) << seriesNbSlices[*lit] * nbFrames << setw(25) << seriesPixelSpacing[*lit].first << " x " << seriesPixelSpacing[*lit].second << endl;
		cout << endl;
	}
}

/// iterator of the selectedSeries set
void VolumeData::check_serie(set<int>::iterator& it_serie)
{
	/// if 2 series have same name/description, image size, pixel size and image orientation, check if they should be fused

	set<int>::iterator lit, lend(selectedSeries.end());
	lit = it_serie;
	++lit;

	int serieNumber = *it_serie;

	string desc = seriesDescription[*it_serie];
	pair<int, int> dim = seriesImageSizes[*it_serie];
	vector<double> orientation = seriesImageOrientation[*it_serie];
	pair<double, double> pixelSp = seriesPixelSpacing[*it_serie];

	for (; lit != lend; )
	{
		if ( seriesDescription[*lit].compare(desc) ||
				seriesPixelSpacing[*lit] != pixelSp ||
				seriesImageSizes[*lit] != dim )
		{
			++lit;
			continue;
		}

		bool similar_orientation = true;

		for (int i=0; i<6; i++)
		{
			if ( fabs(seriesImageOrientation[*lit][i] - orientation[i]) > 1e-5 )
			{
				similar_orientation = false;
				break;
			}
		}

		if (!similar_orientation)
		{
			++lit;
			continue;
		}

		handle_similar_series(it_serie, lit);

		if (*it_serie != serieNumber) break; // the current serie has been erased, so it does not make sense to continue the comparisons
	}
}

/// iterators of the selectedSeries set
void VolumeData::handle_similar_series(set<int>::iterator& it_serie, set<int>::iterator& it_serie_add)
{
	cout << endl;
	cout << "Two similar series detected" << endl;
	cout << "series " << *it_serie << " and " << *it_serie_add << endl;

	/// ask the user what to do (discard one serie, fuse series or keep them?)

	int slice_mid_1 = seriesNbSlices[*it_serie]/2;

	myargs userData1;
	userData1.ptr = this;
	userData1.serie = *it_serie;
	userData1.slice = slice_mid_1;
	userData1.frame = 0;
	userData1.win_name = string("Serie #1");

	float maxVal = 0;
	int Rows = seriesImageSizes[*it_serie].first, Columns = seriesImageSizes[*it_serie].second;
	for(int r=0; r<Rows; r++)
		for(int c=0; c<Columns; c++)
		{
			if (seriesDataValues[*it_serie][slice_mid_1][0][c + Columns * r] > maxVal)
				maxVal = seriesDataValues[*it_serie][slice_mid_1][0][c + Columns * r];
		}

	int valIni = (int)maxVal;
	printImage(*it_serie, 0, 0, "Serie #1", valIni);
	cv::createTrackbar("range", "Serie #1", &valIni, valIni, updateRange, &userData1);

	int slice_mid_2 = seriesNbSlices[*it_serie_add]/2;

	myargs userData2;
	userData2.ptr = this;
	userData2.serie = *it_serie_add;
	userData2.slice = slice_mid_2;
	userData2.frame = 0;
	userData2.win_name = string("Serie #2");

	Rows = seriesImageSizes[*it_serie_add].first;
	Columns = seriesImageSizes[*it_serie_add].second;
	for(int r=0; r<Rows; r++)
		for(int c=0; c<Columns; c++)
		{
			if (seriesDataValues[*it_serie_add][slice_mid_2][0][c + Columns * r] > maxVal)
				maxVal = seriesDataValues[*it_serie_add][slice_mid_2][0][c + Columns * r];
		}

	int valIni2 = (int)maxVal;
	printImage(*it_serie_add, 0, 0, "Serie #2", valIni2);
	cv::createTrackbar("range", "Serie #2", &valIni2, valIni2, updateRange, &userData2);

	cv::waitKey(100);

	char choice;
	int choice_serie;
	bool stop;

	while(1)
	{
		cout << endl;
		cout << "Should we discard one serie (d), fuse them (f) or do nothing (n)? d/f/n/s(change display settings)" << endl;
		cin >> choice;

		switch(choice)
		{
			case 'n':
				{
					cvDestroyWindow("Choice #1");
					cvDestroyWindow("Choice #2");
					cv::waitKey(50);

					/// the series should not be discarded or fused, so do nothing
					++it_serie_add;
					return;
				}

			case 'd':
				{
					set<int>::iterator it_tmp;

					cout << "Which serie should be discarded? 1 or 2" << endl;
					stop = false;
					while (!stop)
					{
						cin >> choice_serie;

						switch(choice_serie)
						{
							case 1:
								{
									it_tmp = it_serie;
									++it_serie;
									selectedSeries.erase(it_tmp);
									stop = true;
									break;
								}

							case 2:
								{
									it_tmp = it_serie_add;
									++it_serie_add;
									selectedSeries.erase(it_tmp);
									stop = true;
									break;
								}

							default:
								cout << "Please choose between 1 and 2" << endl;
						}
					}

					cvDestroyWindow("Choice #1");
					cvDestroyWindow("Choice #2");
					cv::waitKey(50);

					return;
				}

			case 'f':
				{
					cout << "Fusing series..." << endl;

					cvDestroyWindow("Choice #1");
					cvDestroyWindow("Choice #2");
					cv::waitKey(50);

					fuse_series(it_serie, it_serie_add);
					return;
				}

			case 's':
				cout << "Set the position and grey value range of the images, then press a key (while images selected) to continue" << endl;
				cv::waitKey();
				cout << "display setting done" << endl;
				break;

			default:
				cout << "Please, choose between 'd', 'f', 'n' or reset with 's'" << endl;
				break;
		}
	}
}


void VolumeData::replace_slice(int serieNumber, int slice, int add_serie, int add_slice, double add_distance)
{
	cout << "replacing serie " << serieNumber << " slice " << slice << endl;
	seriesImagePosition[serieNumber][slice] = seriesImagePosition[add_serie][add_slice];
	seriesDataValues[serieNumber][slice] = seriesDataValues[add_serie][add_slice];
	listDistances[serieNumber].insert(add_distance);
}


void VolumeData_dicom::replace_slice(int serieNumber, int slice, int add_serie, int add_slice, double add_distance)
{
	VolumeData::replace_slice(serieNumber, slice, add_serie, add_slice, add_distance);

	vector< string>::iterator vit(seriesImagesFilenamesSorted[add_serie][add_slice].begin()), vend(seriesImagesFilenamesSorted[add_serie][add_slice].end());

	for (; vit != vend; ++vit) seriesImagesFilenames[serieNumber].remove(*vit);

	seriesImagesFilenamesSorted[serieNumber][slice] = seriesImagesFilenamesSorted[add_serie][add_slice];

	seriesImagesFilenames[serieNumber].insert( seriesImagesFilenames[serieNumber].end(), seriesImagesFilenamesSorted[add_serie][add_slice].begin(), seriesImagesFilenamesSorted[add_serie][add_slice].end() );

	seriesTriggerTimes[serieNumber][slice] = seriesTriggerTimes[add_serie][add_slice];
}

void VolumeData_nifti::replace_slice(int serieNumber, int slice, int add_serie, int add_slice, double add_distance)
{
	VolumeData::replace_slice(serieNumber, slice, add_serie, add_slice, add_distance);

	seriesImagesFilenames[serieNumber].push_back( seriesImagesFilenames[add_serie].front() );
}

void VolumeData_nrrd::replace_slice(int serieNumber, int slice, int add_serie, int add_slice, double add_distance)
{
	VolumeData::replace_slice(serieNumber, slice, add_serie, add_slice, add_distance);

	seriesImagesFilenames[serieNumber].push_back( seriesImagesFilenames[add_serie].front() );
}


/// iterators of the selectedSeries set
void VolumeData::fuse_series(set<int>::iterator& it_serie, set<int>::iterator& it_serie_add)
{
	/// first, check if both series have slice at the same position, and let the user resolve conflicts

	set< pair<double, pair<int, int> > > dists; // distance of the slice to the origine of the coordinates system along the normal axis: used to sort the slices. pair<dist, pair<serie, slice> >

	set<double>::const_iterator sit(listDistances[*it_serie].begin()), send(listDistances[*it_serie].end());

	// list the distances of the slices of the existing serie
	for (int s=0; sit != send; ++sit, s++)
	{
		dists.insert( make_pair(*sit, make_pair(*it_serie, s) ) );
	}

	// check each slice of the serie to add for conflicts with the slices of the existing serie
	set< pair<double, pair<int, int> > >::const_iterator sitdists, senddists;

	int conflictingSerie, conflictingSlice, choice;
	bool found, stop;

	sit = listDistances[*it_serie_add].begin();
	send = listDistances[*it_serie_add].end();

	for (int s=0; sit != send; ++sit, s++)
	{
		//is the slice distance already known?
		found = false;

		sitdists = dists.begin();
		senddists = dists.end();

		for (; sitdists!=senddists; ++sitdists)
		{
			cout << fabs(sitdists->first - *sit) << endl;
			if ( fabs(sitdists->first - *sit) < 1e-4 )
			{
				found = true;
				conflictingSerie = sitdists->second.first;
				conflictingSlice = sitdists->second.second;
				break;
			}
		}

		// if not create a new entry in dists
		if (!found) dists.insert( make_pair(*sit, make_pair(*it_serie_add, s) ) );
		else // else, let the user choose which slice to use
		{
			cout << endl;
			cout << "distance = " << *sit << endl;

			cout << left;
			cout << setw(15) << "Serie number" << setw(25) << "Serie description" << endl;
			cout << setw(15) << "Nb of images" << setw(25) << "Images size (rows x cols)" << "Pixel spacing (rows x cols)" << endl;
			cout << "Images orientation" << endl << endl;

			cout << setw(15) << *it_serie_add << setw(25) << seriesDescription[*it_serie_add] << endl;
			cout << setw(15) << seriesNbSlices[*it_serie_add]*nbFrames << setw(25) << seriesImageSizes[*it_serie_add].first << " x " << seriesImageSizes[*it_serie_add].second << seriesPixelSpacing[*it_serie_add].first << " x " << seriesPixelSpacing[*it_serie_add].second << endl;
			for (int i=0; i<6; i++) cout << seriesImageOrientation[*it_serie_add][i] << " ";
			cout << endl << endl;

			cout << "vs" << endl << endl;

			cout << setw(15) << conflictingSerie << setw(25) << seriesDescription[conflictingSerie] << endl;
			cout << setw(15) << seriesNbSlices[conflictingSerie]*nbFrames << setw(25) << seriesImageSizes[conflictingSerie].first << " x " << seriesImageSizes[conflictingSerie].second << seriesPixelSpacing[conflictingSerie].first << " x " << seriesPixelSpacing[conflictingSerie].second << endl;
			for (int i=0; i<6; i++) cout << seriesImageOrientation[conflictingSerie][i] << " ";
			cout << endl << endl;

			myargs userData1;
			userData1.ptr = this;
			userData1.serie = *it_serie_add;
			userData1.slice = s;
			userData1.frame = 0;
			userData1.win_name = string("Choice #1");

			float maxVal = 0;
			int Rows = seriesImageSizes[*it_serie_add].first, Columns = seriesImageSizes[*it_serie_add].second;
			for(int r=0; r<Rows; r++)
				for(int c=0; c<Columns; c++)
				{
					if (seriesDataValues[*it_serie_add][s][0][c + Columns * r] > maxVal)
						maxVal = seriesDataValues[*it_serie_add][s][0][c + Columns * r];
				}

			int valIni = (int)maxVal;
			printImage(*it_serie_add, s, 0, "Choice #1", valIni);
			cv::createTrackbar("trackbar1", "Choice #1", &valIni, valIni, updateRange, &userData1);

			myargs userData2;
			userData2.ptr = this;
			userData2.serie = conflictingSerie;
			userData2.slice = conflictingSlice;
			userData2.frame = 0;
			userData2.win_name = string("Choice #2");

			Rows = seriesImageSizes[conflictingSerie].first;
			Columns = seriesImageSizes[conflictingSerie].second;
			for(int r=0; r<Rows; r++)
				for(int c=0; c<Columns; c++)
				{
					if (seriesDataValues[conflictingSerie][conflictingSlice][0][c + Columns * r] > maxVal)
						maxVal = seriesDataValues[conflictingSerie][conflictingSlice][0][c + Columns * r];
				}

			int valIni2 = (int)maxVal;
			printImage(conflictingSerie, conflictingSlice, 0, "Choice #2", valIni2);
			cv::createTrackbar("trackbar2", "Choice #2", &valIni2, valIni2, updateRange, &userData2);
			cv::waitKey(100);

			stop = false;

			while(!stop)
			{
				cout << endl;
				cout << "Which image should be used?" << endl;
				cout << "1: use image #1" << endl;
				cout << "2: use image #2" << endl;
				cout << "0: cancel fusion" << endl;
				cout << "9: show display settings" << endl;

				cin >> choice;

				switch(choice)
				{
					case 1:
						cout << "Image #1 will be used" << endl;
						cout << "Discarding image #2" << endl;

						replace_slice(conflictingSerie, conflictingSlice, *it_serie_add, s, *sit);

						stop = true;
						break;					
					case 2:
						cout << "Image #2 will be used" << endl;
						cout << "Discarding image #1" << endl;

						stop = true;
						break;

					case 0:
						cout << "Cancelling fusion" << endl;
						++it_serie_add;
						return;

					case 9:
						cout << "Set the position and grey value range of the images, then press a key (while images selected) to continue" << endl;
						cv::waitKey();
						cout << "display setting done" << endl;
						break;

					default:
						cout << "Please, choose between #1 and #2 or reset display #9" << endl;
						break;
				} //switch(choice)
			} //while user

			cvDestroyWindow("Choice #1");
			cvDestroyWindow("Choice #2");
			cv::waitKey(100);
		} //if conflict
	} //foreach slice

	/// fill the arrays of the new fused serie

	add_slices_to_serie(*it_serie, dists);

	set<int>::iterator it_tmp = it_serie_add;
	++it_serie_add;
	selectedSeries.erase(it_tmp);
}

void VolumeData::add_slices_to_serie(int serieNumber, set< pair<double, pair<int, int> > >& dists)
{
	seriesNbSlices[serieNumber] = dists.size();
	cout << "The new fused serie has now " << seriesNbSlices[serieNumber] << " slices" << endl;

	set< pair<double, pair<int, int> > >::const_iterator sitdists(dists.begin()), senddists(dists.end());
	int serie, slice;

	for (int s=0; sitdists!=senddists; ++sitdists, s++)
	{
		serie = sitdists->second.first;
		slice = sitdists->second.second;

		if (serie == serieNumber) continue;

		cout << "adding slice " << s << " : serie " << serie << " slice " << slice << endl;

		seriesImagePosition[serieNumber].insert(seriesImagePosition[serieNumber].begin()+s, seriesImagePosition[serie][slice]); //serie, slice, xyz
		seriesDataValues[serieNumber].insert(seriesDataValues[serieNumber].begin()+s, seriesDataValues[serie][slice]); //serie, slice, timeframe, col + Columns * row

		listDistances[serieNumber].insert( sitdists->first );
	}
}

void VolumeData::choose_serie_ref()
{
	printSelectedSeries();

	cout << "Please, choose the serie that defines the orientation and resolution of the volume.." << endl;

	int choice;
	bool stop = false;

	while (!stop)
	{
		cin >> choice;

		if ( selectedSeries.count(choice) )
		{
			serie_ref = choice;
			stop = true;
		}
		else
			cout << "Please, choose a serie between 0 and " << seriesImageSizes.size() - 1 << endl;
	}
}

void VolumeData::computeVolumeCoords()
{
	cout << "VolumeData::computeVolumeCoords() " << endl;

	// orientation of serie ref
	Vector3 colVector_ref(seriesImageOrientation[serie_ref][0],
			seriesImageOrientation[serie_ref][1],
			seriesImageOrientation[serie_ref][2]);

	Vector3 rowVector_ref(seriesImageOrientation[serie_ref][3],
			seriesImageOrientation[serie_ref][4],
			seriesImageOrientation[serie_ref][5]);

	Vector3 normalVector_ref = vectProduct(colVector_ref, rowVector_ref);

	mat_serieRef = Matrice4x4(rowVector_ref[0], colVector_ref[0], normalVector_ref[0], 0,
			rowVector_ref[1], colVector_ref[1], normalVector_ref[1], 0,
			rowVector_ref[2], colVector_ref[2], normalVector_ref[2], 0,
			0, 0, 0, 1);

	mat_serieRef_inv = mat_serieRef.inv();

	offset_serieRef = mat_serieRef_inv * seriesImagePosition[serie_ref][0];
	cout << seriesImagePosition[serie_ref][0] << " " << offset_serieRef << endl;

	set<int>::const_iterator sit(selectedSeries.begin()), send(selectedSeries.end());

	for ( ; sit != send; ++sit)
	{
		// orientation of serie
		Vector3 colVector(seriesImageOrientation[*sit][0],
				seriesImageOrientation[*sit][1],
				seriesImageOrientation[*sit][2]);

		Vector3 rowVector(seriesImageOrientation[*sit][3],
				seriesImageOrientation[*sit][4],
				seriesImageOrientation[*sit][5]);

		Vector3 normalVector = vectProduct(colVector, rowVector);

		int Rows = seriesImageSizes[*sit].first, Columns = seriesImageSizes[*sit].second;
		int nbPixImg = Rows * Columns;

		Vector3 coordPatient, coord_vol;

		for (int s=0; s<seriesNbSlices[*sit]; s++)
		{
			seriesDataVolumeCoords[ make_pair(*sit, s) ] = new Point3<float>[nbPixImg];

			Vector3 ImagePositionPatient = seriesImagePosition[*sit][s];

			Matrice4x4 mat(rowVector[0], colVector[0], normalVector[0], ImagePositionPatient[0],
					rowVector[1], colVector[1], normalVector[1], ImagePositionPatient[1],
					rowVector[2], colVector[2], normalVector[2], ImagePositionPatient[2],
					0, 0, 0, 1);

			//store the inverse matrix (used later for computing the image coordinates of a voxel)
			mat_img_inv[ make_pair(*sit, s) ] = mat.inv();

			for (int r=0; r<Rows; r++)
				for (int c=0; c<Columns; c++)
				{
					/// compute image to patient coord

					Vector3 coordImage(r * seriesPixelSpacing[*sit].first,
							c * seriesPixelSpacing[*sit].second,
							0 );

					coordPatient = mat * coordImage;

					/// compute patient to volume coord

					coord_vol = mat_serieRef_inv * coordPatient;

					coord_vol[0] = (coord_vol[0] - offset_serieRef[0]) / seriesPixelSpacing[serie_ref].first;
					coord_vol[1] = (coord_vol[1] - offset_serieRef[1]) / seriesPixelSpacing[serie_ref].second;
					coord_vol[2] = (coord_vol[2] - offset_serieRef[2]) / seriesPixelSpacing[serie_ref].first;

					seriesDataVolumeCoords[ make_pair(*sit, s) ][ computeIndex(seriesImageSizes[*sit], r, c) ] = coord_vol;
				}
		}
	}
}

void VolumeData::buildListTimeframes()
{
	cout << "VolumeData::buildListTimeframes()" << endl;

	if ( xml_reader.use_first_timeframe_only() ) nbFrames = 1;

	if (nbFrames != 1)
		cout << "VolumeData::buildListTimeframes() Warning: this implementation assumes that frames of different slices are synchronised (e.g. by EMG). If this is not the case, you will need to modify the code in VolumeData::buildListTimeframes() to assign proper time values to the individual frames. If alignment is required in the time domain, IReSD will need to be modified as well (feature not currently implemented)." << endl;

	timeOfFrame = new float*[nbImages];

	int image;

	for (image = 0; image < nbImages; image++)
		timeOfFrame[image] = new float[nbFrames];

	if (nbFrames <= 1)
	{
#pragma omp parallel for default(shared) private(image) schedule(dynamic) num_threads(NB_THREADS)
		for (image=0; image<nbImages; image++)
			timeOfFrame[image][0] = 0;

		size[3] = 1;

		return;
	}

	/// fill array ///
	//time = frame number, we don't care about the triggerTime because we assume that the frames are synchronised (e.g. by ECG), so all the slices of a timeframe should show the same moment in time

	int f;

#ifdef _OPENMP
#if GCC_VERSION >= 40400 && _OPENMP >= 200805

#pragma omp parallel for default(shared) private(image, f) schedule(dynamic) collapse(2) num_threads(NB_THREADS)
	for (image=0; image<nbImages; image++)
		for (f=0; f<nbFrames; f++)
			timeOfFrame[image][f] = f;

#else

#pragma omp parallel for default(shared) private(image, f) schedule(dynamic) num_threads(NB_THREADS)
	for (image=0; image<nbImages; image++)
		for (f=0; f<nbFrames; f++)
			timeOfFrame[image][f] = f;

#endif
#else
	for (image=0; image<nbImages; image++)
		for (f=0; f<nbFrames; f++)
			timeOfFrame[image][f] = f;
#endif

	size[3] = nbFrames;
}

void VolumeData::buildVolume()
{
	cout << "VolumeData::buildVolume" << endl;

	if (xml_reader.size_volume_known())
	{
		Vector3 size_volume = xml_reader.get_volume_size();

		size[0] = (int)size_volume[0];
		size[1] = (int)size_volume[1];
		size[2] = (int)size_volume[2];

		offset_ROI = xml_reader.get_ROI_offset();

		xMin = 0;
		yMin = 0;
		zMin = 0;
		xMax = size[0]-1;
		yMax = size[1]-1;
		zMax = size[2]-1;
	}
	else
	{
		/// find volume's limits

		xMin = VAL_INVALID;
		xMax = -VAL_INVALID;
		yMin = VAL_INVALID;
		yMax = -VAL_INVALID;
		zMin = VAL_INVALID;
		zMax = -VAL_INVALID;

		list<int> list_series;

		if ( xml_reader.serie_ref_defines_volume_size() ) list_series.push_back(serie_ref);
		else list_series.assign(selectedSeries.begin(), selectedSeries.end());

		list<int>::const_iterator lit(list_series.begin()), lend(list_series.end());
		float x, y, z;
		int nbPixImg, s, indImg;

		for (; lit != lend; ++lit)
		{
			nbPixImg = seriesImageSizes[*lit].first * seriesImageSizes[*lit].second;

			for (s=0; s<seriesNbSlices[*lit]; s++)
			{
				for (indImg = 0; indImg < nbPixImg; indImg++)
				{
					x = seriesDataVolumeCoords[ make_pair(*lit, s) ][indImg].x();
					if (x == VAL_INVALID) continue;

					y = seriesDataVolumeCoords[ make_pair(*lit, s) ][indImg].y();
					z = seriesDataVolumeCoords[ make_pair(*lit, s) ][indImg].z();

					if (x > xMax) xMax = x;
					if (x < xMin) xMin = x;

					if (y > yMax) yMax = y;
					if (y < yMin) yMin = y;

					if (z > zMax) zMax = z;
					if (z < zMin) zMin = z;
				}
			}
		}

		cout << "initial volume limits: x " << xMin << " " << xMax << " y " << yMin << " " << yMax << " z " << zMin << " " << zMax << endl;

		/// ROI

		if ( xml_reader.ROI_selection() ) ROIselection();

		cout << "ROI volume limits: x " << xMin << " " << xMax << " y " << yMin << " " << yMax << " z " << zMin << " " << zMax << endl;

		size[0] = (int)floor(xMax - xMin)+1;
		size[1] = (int)floor(yMax - yMin)+1;
		size[2] = (int)floor(zMax - zMin)+1;

		offset_ROI[0] = xMin;
		offset_ROI[1] = yMin;
		offset_ROI[2] = zMin;

		xMin -= offset_ROI[0];
		yMin -= offset_ROI[1];
		zMin -= offset_ROI[2];

		xMax -= offset_ROI[0];
		yMax -= offset_ROI[1];
		zMax -= offset_ROI[2];

		cout << "new volume limits: x " << xMin << " " << xMax << " y " << yMin << " " << yMax << " z " << zMin << " " << zMax << endl;
	}

	cout << "volume size : " << size[0] << " " << size[1] << " " << size[2] << " " << size[3] << endl;

	nbPixVol = size[0]*size[1]*size[2];
	nbPix = nbPixVol * size[3];


	/// offset volume coords (because of ROI selection)

	vector< pair<int, int> >::const_iterator lit(selectedImages.begin()), lend(selectedImages.end());
	int nbPixImg, indImg;

	for (; lit != lend; ++lit)
	{
		nbPixImg = seriesImageSizes[lit->first].first * seriesImageSizes[lit->first].second;

		for (indImg = 0; indImg < nbPixImg; indImg++)
		{
			if (seriesDataVolumeCoords[*lit][indImg].x() == VAL_INVALID) continue;

			seriesDataVolumeCoords[*lit][indImg] -= offset_ROI;
		}
	}

	/// create volume

	if ( volume != NULL )
	{
		for (int f=0; f<nbFrames; f++) delete[] volume[f];
		delete[] volume;
	}

	volume = new float*[size[3]];

	int i, f;

	for (f=0; f<size[3]; f++)
	{
		volume[f] = new float[nbPixVol];

#pragma omp parallel for default(shared) private(i) schedule(dynamic) num_threads(NB_THREADS)
		for (i=0; i<nbPixVol; i++)
			volume[f][i] = VAL_INVALID;
	}

	imagePresence = new set<int>[nbPixVol];
}

void VolumeData::ROIselection()
{
	cout << "VolumeData::ROIselection()" << endl;

	vector< pair<int, int> >::const_iterator lit(selectedImages.begin()), lend(selectedImages.end());

	int indImg;
	int cMin, cMax, rMin, rMax;
	float xMin_tmp, xMax_tmp, yMin_tmp, yMax_tmp, zMin_tmp, zMax_tmp;
	Point3<float> point;

	xMin_tmp = VAL_INVALID;
	yMin_tmp = VAL_INVALID;
	zMin_tmp = VAL_INVALID;

	xMax_tmp = -VAL_INVALID;
	yMax_tmp = -VAL_INVALID;
	zMax_tmp = -VAL_INVALID;

	cv::namedWindow("ROI selection", CV_WINDOW_AUTOSIZE);

	mousse_args_box args;
	args.actif = false;
	cvSetMouseCallback("ROI selection", mouse_callback, (void*)&args);

	for (; lit != lend ; ++lit)
	{
		int Row = seriesImageSizes[lit->first].first;
		int Col = seriesImageSizes[lit->first].second;

		// look for an ROI file
		string filename = seriesImagesFilenamesSorted[lit->first][lit->second][0];
		int pos = filename.find_last_of(".");
		string descrFile = filename.substr(0, pos) + "_ROI.txt";

		ifstream fichier(descrFile.c_str(), ios::in);

		if (fichier)
		{
			fichier >> cMin;
			fichier >> cMax;
			fichier >> rMin;
			fichier >> rMax;

			cout << "ROI read" << endl;

			fichier.close();
		}
		else
		{
			// draw a rectangle to define an ROI

			cv::Mat imgDisp(Row, Col, CV_64FC1);
			float max = -VAL_INVALID;
			float min = VAL_INVALID;

			for (int i=0; i<Row; i++)
				for (int j=0; j<Col; j++)
				{
					float val = seriesDataValues[lit->first][lit->second][0][j + Col * i];

					if (val != VAL_INVALID)
					{
						if (val > max) max = val;
						if (val < min) min = val;
					}

					imgDisp.at<double>(i, j) = (double)val;
				}

			if (min < -1000) min = -1000;

			for (int i=0; i<Row; i++)
				for (int j=0; j<Col; j++)
				{
					if (imgDisp.at<double>(i, j) != VAL_INVALID)
						imgDisp.at<double>(i, j) = (imgDisp.at<double>(i, j) - min) / (max - min) * 0.99 + 0.01;
					else
						imgDisp.at<double>(i, j) = 0;
				}

			cv::imshow("ROI selection", imgDisp);

			args.drawing_box = false;
			args.img = &imgDisp;
			args.validated = false;
			args.cancelled = false;
			args.actif = true;

			cout << "Draw a rectangle, then press 'v' to validate, or 'a' to cancel and go to the next image" << endl;

			char reponse;

			while (!args.validated && !args.cancelled)
			{
				reponse = cv::waitKey(15);

				if (reponse == 'a')
				{
					args.cancelled = true;
					args.actif = false;
				}
			}

			if (args.validated)
			{
				cMin = args.box.x;
				cMax = args.box.x + args.box.width;
				rMin = args.box.y;
				rMax = args.box.y + args.box.height;
			}
			else
			{
				cMin = VAL_INVALID;
			}
		} //no ROI file

		if (cMin != VAL_INVALID)
		{
			/// select min and max coordinates in the ROI

			for (int r=rMin; r<=rMax; r++)
				for (int c=cMin; c<=cMax; c++)
				{
					indImg = c + Col * r;
					point = seriesDataVolumeCoords[*lit][indImg];

					if (point.x() == VAL_INVALID) continue;

					if (point.x() < xMin_tmp) xMin_tmp = point.x();
					if (point.y() < yMin_tmp) yMin_tmp = point.y();
					if (point.z() < zMin_tmp) zMin_tmp = point.z();

					if (point.x() > xMax_tmp) xMax_tmp = point.x();
					if (point.y() > yMax_tmp) yMax_tmp = point.y();
					if (point.z() > zMax_tmp) zMax_tmp = point.z();
				}
		}
	}

	/// compare with current min/max

	if (xMin_tmp > xMin) xMin = xMin_tmp;
	if (yMin_tmp > yMin) yMin = yMin_tmp;
	if (zMin_tmp > zMin) zMin = zMin_tmp;

	if (xMax_tmp < xMax) xMax = xMax_tmp;
	if (yMax_tmp < yMax) yMax = yMax_tmp;
	if (zMax_tmp < zMax) zMax = zMax_tmp;

	cvDestroyWindow("ROI selection");
	cv::waitKey(100);
}

void VolumeData::loadRegistration()
{
	/// initialise arrays

	int nb_images = selectedImages.size();

	translationVect = new Vector3[nb_images];
	rotationAngles = new Vector3[nb_images];

	pixelCoordsInVolume = new map<int, Point3<float> >[nb_images];
	voxelCoordsInImageVolume = new map<int, Point3<float> >[nb_images];


	/// read registration from xml file and fill translationVect and rotationAngles

	vector< pair<int, int> >::const_iterator itS(selectedImages.begin()), endS(selectedImages.end());

	for (int i=0; itS != endS; itS++, i++)
	{
		translationVect[i] = xml_reader.get_registration_translation(seriesDescription[itS->first], itS->second);
		rotationAngles[i] = xml_reader.get_registration_rotation(seriesDescription[itS->first], itS->second);
	}
}

void VolumeData::fillVolume()
{
	cout << "VolumeData::fillVolume()" << endl;
	if ( xml_reader.accurate_volume_filling() ) return fillVolume_interpolation();
	else return fillVolume_closestPoint();
}

void VolumeData::fillVolume_closestPoint()
{
	//cout << "VolumeData::fillVolume_closestPoint()" << endl;

	int i, t;

#pragma omp parallel for default(shared) private(i, t) schedule(dynamic) num_threads(NB_THREADS)
	for (i=0; i<nbPixVol; i++)
	{
		imagePresence[i].clear();

		for (t=0; t<size[3]; t++) volume[t][i] = VAL_INVALID;
	}

	int image, serie, slice;
	int nbPixImg, indImg;
	Point3<float> point;
	int indV;
	int fM, fP;
	float distM, distP;
	float res;

#pragma omp parallel for default(shared) private(image, serie, slice, nbPixImg, indImg, point, indV, t, fM, fP, distM, distP, res) schedule(dynamic) num_threads(NB_THREADS)
	for (image=0; image < nbImages; image++)
	{
		serie = selectedImages[image].first;
		slice = selectedImages[image].second;

		nbPixImg = seriesImageSizes[serie].first * seriesImageSizes[serie].second;

		for (indImg=0; indImg<nbPixImg; indImg++)
		{
			if ( seriesDataValues[serie][slice][0][indImg] == VAL_INVALID ) continue;

			point = getImagePointCoordsInVolume(image, indImg);

			point.round_coords();
			if (! isPointInVolume(point) ) continue;

			indV = computeIndex(size, (int)point.x(), (int)point.y(), (int)point.z());

			for (t=0; t<size[3]; t++)
			{
				getClosestFrames(t, image, fM, distM, fP, distP);

				if (distM > 1 && distP > 1) continue; // too far from any usable frame

				res = 0;

				if (distM <= 1)
				{
					if (seriesDataValues[serie][slice][fM][indImg] != VAL_INVALID)
					{
						res += seriesDataValues[serie][slice][fM][indImg] * (1. - distM);
					}
				}

				if (distP <= 1)
				{
					if (seriesDataValues[serie][slice][fP][indImg] != VAL_INVALID)
					{
						res += seriesDataValues[serie][slice][fP][indImg] * (1. - distP);
					}
				}

#pragma omp critical
				volume[t][indV] = res;
			}

#pragma omp critical
			imagePresence[indV].insert(image);
		}
	}
}

void VolumeData::fillVolume_interpolation()
{
	cout << "VolumeData::fillVolumes_interpolation()" << endl;

	int i, t;

#pragma omp parallel for default(shared) private(i, t) schedule(dynamic) num_threads(NB_THREADS)
	for (i=0; i<nbPixVol; i++)
	{
		imagePresence[i].clear();

		for (t=0; t<size[3]; t++) volume[t][i] = VAL_INVALID;
	}

	int nbPixImg, indImg;
	Point3<float> point;
	int xM, yM, zM, xP, yP, zP;
	set< Point3<int> > waitingList_set;
	Point3<float> imageCoord;
	int indV;
	int fM, fP;
	float distM, distP, valtM, valtP, ax, ay, az;
	int indMM, indMP, indPM, indPP;
	float valMM, valMP, valPM, valPP, valCM, valCP, res;

	vector< pair<int, int> >::const_iterator vit(selectedImages.begin()), vend(selectedImages.end());

	for (int image=0; vit != vend; ++vit, image++)
	{
		waitingList_set.clear();

		int serie = vit->first;
		int slice = vit->second;

		const pair<int, int>& imgSize = seriesImageSizes[serie];
		nbPixImg = imgSize.first * imgSize.second;

#pragma omp parallel for default(shared) private(indImg, point, xM, yM, zM, xP, yP, zP) num_threads(NB_THREADS)
		for (indImg=0; indImg<nbPixImg; indImg++)
		{
			if ( seriesDataValues[serie][slice][0][indImg] == VAL_INVALID ) continue;

			point = getImagePointCoordsInVolume(image, indImg);

			xM = (int)floor(point.x());
			yM = (int)floor(point.y());
			zM = (int)floor(point.z());

			xP = (int)ceil(point.x());
			yP = (int)ceil(point.y());
			zP = (int)ceil(point.z());

#pragma omp critical
			waitingList_set.insert( Point3<int>(xM, yM, zM) );
#pragma omp critical
			waitingList_set.insert( Point3<int>(xP, yM, zM) );
#pragma omp critical
			waitingList_set.insert( Point3<int>(xM, yP, zM) );
#pragma omp critical
			waitingList_set.insert( Point3<int>(xP, yP, zM) );
#pragma omp critical
			waitingList_set.insert( Point3<int>(xM, yM, zP) );
#pragma omp critical
			waitingList_set.insert( Point3<int>(xP, yM, zP) );
#pragma omp critical
			waitingList_set.insert( Point3<int>(xM, yP, zP) );
#pragma omp critical
			waitingList_set.insert( Point3<int>(xP, yP, zP) );
		}

		//convert waitingList into a vector for openMP
		vector< Point3<int> > waitingList(waitingList_set.begin(), waitingList_set.end());
		waitingList_set.clear();

		int indList, list_size = waitingList.size();

#pragma omp parallel for default(shared) private(indList, indV, imageCoord, xM, yM, zM, xP, yP, zP, t, fM, fP, distM, distP, valtM, valtP, ax, ay, az, indMM, indMP, indPM, indPP, valMM, valMP, valPM, valPP, valCM, valCP, res) schedule(dynamic) num_threads(NB_THREADS)
		for (indList=0; indList < list_size; indList++)
		{
			if ( !isPointInVolume(waitingList[indList]) ) continue;

			indV = computeIndex(size, waitingList[indList].x(), waitingList[indList].y(), waitingList[indList].z());

#pragma omp critical
			if ( !imagePresence[indV].size() )
			{
				for (t=0; t<size[3]; t++)
				{
					volume[t][indV] =  0;
				}
			}

			imageCoord = getVoxelCoordsInImageVolume(indV, image);

			zM = (int)floor(imageCoord.z());
			zP = (int)ceil(imageCoord.z());

			if (zM != 0 && zP != 0) continue; // too far from the plane of the image

			/// tri-linear interpolation of the pixel intensity at imageCoord

			xM = (int)floor(imageCoord.x());
			xP = (int)ceil(imageCoord.x());

			yM = (int)floor(imageCoord.y());
			yP = (int)ceil(imageCoord.y());

			if (xM < 0) xM = 0;
			if (xP < 0) xP = 0;
			if (xM >= imgSize.first) xM = imgSize.first - 1;
			if (xP >= imgSize.first) xP = imgSize.first - 1;
			if (yM < 0) yM = 0;
			if (yP < 0) yP = 0;
			if (yM >= imgSize.second) yM = imgSize.second - 1;
			if (yP >= imgSize.second) yP = imgSize.second - 1;

			indMM = computeIndex(imgSize, xM, yM);
			indMP = computeIndex(imgSize, xM, yP);
			indPM = computeIndex(imgSize, xP, yM);
			indPP = computeIndex(imgSize, xP, yP);

			ax = imageCoord.x() - xM;
			ay = imageCoord.y() - yM;
			az = fabs(imageCoord.z());

			for (t=0; t<size[3]; t++)
			{
				getClosestFrames(t, image, fM, distM, fP, distP);

				if (distM > 1 && distP > 1) continue; // too far from any usable frame

				if (ax == 0 && ay == 0) // the projection of the point on the image is a pixel, so no interpolation needed
				{
					if (distM > 1) valtM = 0;
					else valtM = seriesDataValues[serie][slice][fM][indMM];

					if (distP > 1) valtP = 0;
					else valtP = seriesDataValues[serie][slice][fP][indMM];
				}
				else /// bi-linear interpolation in the image plane (z=0)
				{
					if (distM > 1) valtM = 0;
					else
					{
						valMM = seriesDataValues[serie][slice][fM][indMM];
						valMP = seriesDataValues[serie][slice][fM][indMP];
						valPM = seriesDataValues[serie][slice][fM][indPM];
						valPP = seriesDataValues[serie][slice][fM][indPP];

						if (valMM == VAL_INVALID || valMP == VAL_INVALID || valPM == VAL_INVALID || valPP == VAL_INVALID) continue;

						valCM = ax * valPM + (1. - ax) * valMM;
						valCP = ax * valPP + (1. - ax) * valMP;

						valtM = ay * valCP + (1. - ay) * valCM;
					}

					if (distP > 1) valtP = 0;
					else
					{
						valMM = seriesDataValues[serie][slice][fP][indMM];
						valMP = seriesDataValues[serie][slice][fP][indMP];
						valPM = seriesDataValues[serie][slice][fP][indPM];
						valPP = seriesDataValues[serie][slice][fP][indPP];

						if (valMM == VAL_INVALID || valMP == VAL_INVALID || valPM == VAL_INVALID || valPP == VAL_INVALID) continue;

						valCM = ax * valPM + (1. - ax) * valMM;
						valCP = ax * valPP + (1. - ax) * valMP;

						valtP = ay * valCP + (1. - ay) * valCM;
					}
				}

				/// linear interpolation in the t direction at z=0

				res = (1. - distM) * valtM + (1. - distP) * valtP;

				/// linear interpolation in the z direction, using the value at z=0 and at z=+-1 (=0)

				res = (1. - az) * res;

				/// fill volume

#pragma omp critical
				volume[t][indV] += res;
			}

#pragma omp critical
			imagePresence[indV].insert(image);
		}
	}

	/// normalisation

	int nbImages_tmp;

#pragma omp parallel for default(shared) private(i, nbImages_tmp, t) schedule(dynamic) num_threads(NB_THREADS)
	for (i=0; i<nbPixVol; i++)
	{
		nbImages_tmp = imagePresence[i].size();

		if ( nbImages_tmp <= 1 ) continue;

		for (t=0; t<size[3]; t++)
		{
#pragma omp critical
			volume[t][i] /= (float)nbImages_tmp;
		}
	}
}

const set<int>& VolumeData::get_selectedSeries()
{
	return selectedSeries;
}

const set<int>& VolumeData::get_images_at_voxel(int voxel)
{
	return imagePresence[voxel];
}

Point3<float> VolumeData::getImagePointCoordsInVolume(int image, int pixel)
{
	if ( rotationAngles[image] != VECT_NULL || translationVect[image] != VECT_NULL )
	{
		/// first check if the coordinates have already been computed for this point

		map<int, Point3<float> >::const_iterator mit;
		Point3<float> result(VECT_INVALID);

#pragma omp critical
		{
			mit = pixelCoordsInVolume[image].find(pixel);
			if ( mit != pixelCoordsInVolume[image].end() ) result = mit->second;
		}

		if (result != VECT_INVALID) return result;
	}

	//original voume coordinates, before registration
	Point3<float> point = seriesDataVolumeCoords[ selectedImages[image] ][pixel];

	if (rotationAngles[image] == VECT_NULL && translationVect[image] == VECT_NULL ) return point;


	/// apply registration transformation

	if (point.x() != VAL_INVALID)
	{
		/// translation

		point -= translationVect[image];

		/// rotations around all three axes

		// warning : the centre of rotation in at the centre of the volume
		int cx = size[0]/2, cy = size[1]/2, cz = size[2]/2;

		Vector3 vect(point.x() - cx, point.y() - cy, point.z() - cz);

		if ( rotationAngles[image][2] )
		{
			Matrice3x3 matz = Matrice3x3( cos( -rotationAngles[image][2] ), -sin( -rotationAngles[image][2] ), 0,
					sin( -rotationAngles[image][2] ), cos( -rotationAngles[image][2] ), 0,
					0, 0, 1 );
			vect = matz * vect;
		}

		if ( rotationAngles[image][1] )
		{
			Matrice3x3 maty = Matrice3x3( cos( -rotationAngles[image][1] ), 0, sin( -rotationAngles[image][1] ),
					0, 1, 0,
					-sin( -rotationAngles[image][1] ), 0, cos( -rotationAngles[image][1] ) );
			vect = maty * vect;
		}

		if ( rotationAngles[image][0] )
		{
			Matrice3x3 matx = Matrice3x3( 1, 0, 0,
					0, cos( -rotationAngles[image][0] ), -sin( -rotationAngles[image][0] ),
					0, sin( -rotationAngles[image][0] ), cos( -rotationAngles[image][0] ) );
			vect = matx * vect;
		}

		point = Point3<float>(vect + Vector3(cx, cy, cz));
	}

	/// save this coordinate to avoid repeated computations

#pragma omp critical
	pixelCoordsInVolume[image][pixel] = point;

	return point;
}

Point3<float> VolumeData::getVoxelCoordsInImageVolume(int voxel, int image)
{
	/// first check if the coordinates have already been computed for this point

	Point3<float> result(VECT_INVALID);

#pragma omp critical
	{
		map<int, Point3<float> >::const_iterator mit = voxelCoordsInImageVolume[image].find(voxel);
		if ( mit != voxelCoordsInImageVolume[image].end() ) result = mit->second;
	}

	if (result != VECT_INVALID) return result;

	int x, y, z;
	computeCoords(voxel, size, x, y, z);

	Vector3 coord_volume(x, y, z);

	/// apply inverse registration transformation

	if ( rotationAngles[image] != VECT_NULL || translationVect[image] != VECT_NULL )
	{
		// rotations
		// warning : the centre of rotation in at the centre of the volume
		int cx = size[0]/2, cy = size[1]/2, cz = size[2]/2;

		Vector3 vect = coord_volume - Vector3(cx, cy, cz);

		if ( rotationAngles[image][0] )
		{
			Matrice3x3 matx = Matrice3x3( 1, 0, 0,
					0, cos( rotationAngles[image][0] ), -sin( rotationAngles[image][0] ),
					0, sin( rotationAngles[image][0] ), cos( rotationAngles[image][0] ) );
			vect = matx * vect;
		}

		if ( rotationAngles[image][1] )
		{
			Matrice3x3 maty = Matrice3x3( cos( rotationAngles[image][1] ), 0, sin( rotationAngles[image][1] ),
					0, 1, 0,
					-sin( rotationAngles[image][1] ), 0, cos( rotationAngles[image][1] ) );
			vect = maty * vect;
		}

		if ( rotationAngles[image][2] != 0 )
		{
			Matrice3x3 matz = Matrice3x3( cos( rotationAngles[image][2] ), -sin( rotationAngles[image][2] ), 0,
					sin( rotationAngles[image][2] ), cos( rotationAngles[image][2] ), 0,
					0, 0, 1 );
			vect = matz * vect;
		}

		coord_volume = vect + Vector3(cx, cy, cz);

		// translation
		coord_volume += translationVect[image];
	}

	/// offset used to select the ROI
	coord_volume += offset_ROI;

	/// offset and scaling of serie_ref

	coord_volume[0] = coord_volume[0] * seriesPixelSpacing[serie_ref].first + offset_serieRef[0];
	coord_volume[1] = coord_volume[1] * seriesPixelSpacing[serie_ref].second + offset_serieRef[1];
	coord_volume[2] = coord_volume[2] * seriesPixelSpacing[serie_ref].first + offset_serieRef[2];

	/// convert to patient coordinates

	Vector3 patient_coord = mat_serieRef * coord_volume;

	/// convert to image volume coordinates

	Vector3 image_coord = mat_img_inv[ selectedImages[image] ] * patient_coord;

	image_coord[0] /= seriesPixelSpacing[ selectedImages[image].first ].first;
	image_coord[1] /= seriesPixelSpacing[ selectedImages[image].first ].second;
	image_coord[2] /= seriesPixelSpacing[ selectedImages[image].first ].first;

	/// save for later

#pragma omp critical
	voxelCoordsInImageVolume[image][voxel] = image_coord;

	return image_coord;
}

bool VolumeData::isImagePointInVolume(int image, int pixel)
{
	return isPointInVolume( getImagePointCoordsInVolume(image, pixel) );
}

bool VolumeData::isPointInVolume(Point3<float> point)
{
	if (point.x() < 0)
		return false;

	if (point.x() >= size[0])
		return false;

	if(point.y() < 0)
		return false;

	if(point.y() >= size[1])
		return false;

	if(point.z() < 0)
		return false;

	if(point.z() >= size[2])
		return false;

	return true;
}

bool VolumeData::isPointInVolume(Point3<int> point)
{
	if (point.x() < 0)
		return false;

	if (point.x() >= size[0])
		return false;

	if(point.y() < 0)
		return false;

	if(point.y() >= size[1])
		return false;

	if(point.z() < 0)
		return false;

	if(point.z() >= size[2])
		return false;

	return true;
}

bool VolumeData::isPointInVolume(Point3<int> point, int t)
{
	if(t < 0)
		return false;

	if(t >= size[3])
		return false;

	return isPointInVolume(point);
}

void VolumeData::getClosestFrames(int t, int image, int& frameM, float& distM, int& frameP, float& distP)
{
	distM = VAL_INVALID;
	distP = VAL_INVALID;

	if (nbFrames == 1)
	{
		frameM = 0;
		distM = 0;
		return;
	}

	int frame = 0;

	if (timeOfFrame[image][frame] == t)
	{
		frameM = frame;
		distM = 0;

		return;
	}

	if (timeOfFrame[image][frame] > t)
	{
		frameP = frame;
		distP = timeOfFrame[image][frame] - t;

		return;
	}

	while (timeOfFrame[image][frame] < t && frame < nbFrames)
	{
		frame++;
	}

	if (frame == nbFrames)
	{
		frameM = frame - 1;
		distM = t - timeOfFrame[image][frame-1];

		if (periodic_motion)
		{
			frameP = 0;

			// distance to the second occurence of the first frame has to be estimated, based on the average time between frames
			float av_time_between_frames = (timeOfFrame[image][nbFrames-1] - timeOfFrame[image][0]) / (float)(nbFrames - 1);
			distP = timeOfFrame[image][nbFrames-1] + av_time_between_frames - t;

			if (distP < 0)
				throw string("VolumeData::getClosestFrames() Error: time above superior limit");
		}

		return;
	}

	if (timeOfFrame[image][frame] == t)
	{
		frameM = frame;
		distM = 0;
		return;
	}

	frameM = frame - 1;
	distM = t - timeOfFrame[image][frame-1];

	frameP = frame;
	distP = timeOfFrame[image][frame] - t;
}

void VolumeData::display3D(int timeframe, string str)
{
	cout << "VolumeData::display3D()" << timeframe << endl;

	nifti_1_header hdr;
	nifti1_extender pad={0,0,0,0};
	FILE *fp;
	int ret;

	/********** fill in the minimal default header fields */
	bzero((void *)&hdr, sizeof(hdr));
	hdr.sizeof_hdr = MIN_HEADER_SIZE;
	hdr.dim[0] = 3;
	hdr.dim[1] = size[0];
	hdr.dim[2] = size[1];
	hdr.dim[3] = size[2];
	hdr.datatype = NIFTI_TYPE_FLOAT32;
	hdr.bitpix = 32;
	hdr.pixdim[1] = 1;
	hdr.pixdim[2] = 1;
	hdr.pixdim[3] = 1;
	hdr.vox_offset = (float) NII_HEADER_SIZE;
	hdr.scl_slope = 0;
	hdr.xyzt_units = NIFTI_UNITS_MM;
	strncpy(hdr.magic, "n+1\0", 4);


	/********** write first 348 bytes of header   */
	ostringstream name;
	name << "VolumeData";
	if ( str.compare("") ) name << "_" << str;
	name << ".nii";
	fp = fopen(name.str().c_str(),"w");
	if (fp == NULL) {
		fprintf(stderr, "\nError opening header file phi.nii for write\n");
		exit(1);
	}
	ret = fwrite(&hdr, MIN_HEADER_SIZE, 1, fp);
	if (ret != 1) {
		fprintf(stderr, "\nError writing header file phi.nii\n");
		exit(1);
	}


	/********** write extender pad and image data   */
	ret = fwrite(&pad, 4, 1, fp);
	if (ret != 1) {
		fprintf(stderr, "\nError writing header file extension pad phi.nii\n");
		exit(1);
	}

	float* buffer = new float[nbPixVol];
	int i;

#pragma omp parallel for default(shared) private(i) schedule(dynamic) num_threads(NB_THREADS)
	for (i=0; i<nbPixVol; i++)
	{
		if (volume[timeframe][i] == VAL_INVALID) buffer[i] = 0;
		else buffer[i] = volume[timeframe][i];
	}

	ret = fwrite(buffer, (size_t)(hdr.bitpix/8), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], fp);
	if (ret != hdr.dim[1]*hdr.dim[2]*hdr.dim[3]) {
		fprintf(stderr, "\nError writing data to phi.nii\n");
		exit(1);
	}

	fclose(fp);
}

const vector< pair<int, int> >& VolumeData::get_selectedImages()
{
	return selectedImages;
}

const pair<int, int>& VolumeData::getImageSize(int serie)
{
	return seriesImageSizes[serie];
}

float VolumeData::getPixelValue(int serie, int slice, int frame, int pixel)
{
	return seriesDataValues[serie][slice][frame][pixel];
}

float VolumeData::getPixelValue(int image, int frame, int pixel)
{
	return seriesDataValues[selectedImages[image].first][selectedImages[image].second][frame][pixel];
}

float VolumeData::getTimeOfFrame(int image, int timeframe)
{
	return timeOfFrame[image][timeframe];
}

float VolumeData::getVoxelValue(int voxel, int timeframe)
{
	return volume[timeframe][voxel];
}

const int* VolumeData::getVolumeSize()
{
	return size;
}

bool VolumeData::is_motion_periodic()
{
	return periodic_motion;
}

int VolumeData::getNbFrames()
{
	return nbFrames;
}

int VolumeData::getNbImages()
{
	return nbImages;
}

void VolumeData::printImage(int numSerie, int slice, int timeframe, string win_name, int range)
{
	int Rows = seriesImageSizes[numSerie].first, Columns = seriesImageSizes[numSerie].second;

	if (!range)
	{
		for(int r=0; r<Rows; r++)
			for(int c=0; c<Columns; c++)
			{
				if (seriesDataValues[numSerie][slice][timeframe][c + Columns * r] > range) range = seriesDataValues[numSerie][slice][timeframe][c + Columns * r];
			}
	}

	cv::Mat imgDisp(Rows, Columns, CV_64FC1);

	for(int r=0; r<Rows; r++)
		for(int c=0; c<Columns; c++)
		{
			imgDisp.at<double>(r, c) = (double)seriesDataValues[numSerie][slice][timeframe][c + Columns * r] / (float)range;

			if (imgDisp.at<double>(r, c) > 1) imgDisp.at<double>(r, c) = 1;
		}

	cv::namedWindow(win_name, CV_WINDOW_AUTOSIZE);
	cv::imshow(win_name, imgDisp);

	cv::waitKey(10);
}

int VolumeData::getNbSlices(int serie)
{
	return seriesNbSlices[serie];
}

int VolumeData::getImageNumber(int serie, int slice)
{
	return serieFirstImageIndex[serie] + slice;
}

void VolumeData::move_spatial_sequence(int serie, Vector3 d_translation, Vector3 d_rotation)
{
	int last_index = serieFirstImageIndex[serie] + seriesNbSlices[serie];

	for (int image = serieFirstImageIndex[serie]; image<last_index; image++)
	{
		translationVect[image] += d_translation;
		rotationAngles[image] += d_rotation;

		pixelCoordsInVolume[image].clear();
		voxelCoordsInImageVolume[image].clear();
	}
}

void VolumeData::move_temporal_sequence(int image, Vector3 d_translation, Vector3 d_rotation)
{
	translationVect[image] += d_translation;
	rotationAngles[image] += d_rotation;

	pixelCoordsInVolume[image].clear();
	voxelCoordsInImageVolume[image].clear();
}

Vector3 VolumeData::get_rotation_angles(int image)
{
	return rotationAngles[image];
}

Vector3 VolumeData::get_translation_vector(int image)
{
	return translationVect[image];
}

bool VolumeData::comprise_temporal_sequences()
{
	if (nbFrames == 1) return false;
	return true;
}

//vectIni: Vector3(row, col, 0)
Vector3 VolumeData::imageToVolumeVect(int image, const Vector3& vectIni)
{
	int serie = selectedImages[image].first;

	Vector3 colVector(seriesImageOrientation[serie][0],
			seriesImageOrientation[serie][1],
			seriesImageOrientation[serie][2]);

	Vector3 rowVector(seriesImageOrientation[serie][3],
			seriesImageOrientation[serie][4],
			seriesImageOrientation[serie][5]);

	Vector3 normalVector = vectProduct(colVector, rowVector);

	Matrice3x3 mat(rowVector[0], colVector[0], normalVector[0],
			rowVector[1], colVector[1], normalVector[1],
			rowVector[2], colVector[2], normalVector[2]);

	Vector3 vectPatient = mat * vectIni;
	if ( isnan( vectPatient[0] ) || isnan( vectPatient[1] ) || isnan( vectPatient[2] ) )
		throw string("VolumeData::imageToVolumeVect() Erreur : nan #0 !");

	Vector3 colVector2(seriesImageOrientation[serie_ref][0],
			seriesImageOrientation[serie_ref][1],
			seriesImageOrientation[serie_ref][2]);

	Vector3 rowVector2(seriesImageOrientation[serie_ref][3],
			seriesImageOrientation[serie_ref][4],
			seriesImageOrientation[serie_ref][5]);

	Vector3 normalVector2 = vectProduct(colVector2, rowVector2);

	Matrice3x3 mat2(rowVector2[0], colVector2[0], normalVector2[0],
			rowVector2[1], colVector2[1], normalVector2[1],
			rowVector2[2], colVector2[2], normalVector2[2]);

	Matrice3x3 matInv = mat2.inv();

	Vector3 res = matInv * vectPatient;
	if ( isnan( res[0] ) || isnan( res[1] ) || isnan( res[2] ) ) throw string("VolumeData::imageToVolumeVect() Erreur : nan #1 !");

	if ( rotationAngles[image][2] != 0 )
	{
		Matrice3x3 matz = Matrice3x3( cos( -rotationAngles[image][2] ), -sin( -rotationAngles[image][2] ), 0,
				sin( -rotationAngles[image][2] ), cos( -rotationAngles[image][2] ), 0,
				0, 0, 1 );
		res = matz * res;
	}

	if ( rotationAngles[image][1] )
	{
		Matrice3x3 maty = Matrice3x3( cos( -rotationAngles[image][1] ), 0, sin( -rotationAngles[image][1] ),
				0, 1, 0,
				-sin( -rotationAngles[image][1] ), 0, cos( -rotationAngles[image][1] ) );
		res = maty * res;
	}

	if ( rotationAngles[image][0] )
	{
		Matrice3x3 matx = Matrice3x3( 1, 0, 0,
				0, cos( -rotationAngles[image][0] ), -sin( -rotationAngles[image][0] ),
				0, sin( -rotationAngles[image][0] ), cos( -rotationAngles[image][0] ) );
		res = matx * res;
	}

	return res;
}


/// VolumeData_dicom ///

VolumeData_dicom::VolumeData_dicom(XMLreader& _xml_reader):
	VolumeData(_xml_reader)
{}

void VolumeData_dicom::listSeries()
{
	cout << "listing series..." << endl;

	/// get the filenames

	DIR *rep = opendir(repertory.c_str());
	list<string> filenames;

	if (rep != NULL)
	{
		listFilenames(rep, filenames, repertory + "/");

		closedir (rep);
		rep = NULL;
	}
	else
	{
		cout << repertory << endl;
		throw(string("Error: cannot open the directory"));
	}

	/// read the images
	list<string>::const_iterator lit(filenames.begin()), lend(filenames.end());

	DcmFileFormat fileformat;
	OFCondition status;

	string SerieDescription;

	for(; lit != lend; ++lit)
	{
		//cout << *lit << endl;

		status = fileformat.loadFile(lit->c_str());

		bool stop = false;

		if (status.good())
		{
			/**read serie attributes**/
			unsigned short int Rows, Columns;
			double ImageOrientationPatient[6], PixelSpacing[2];
			OFString ImageType, SerieDescription, DerivationDescription;

			if (!fileformat.getDataset()->findAndGetOFStringArray(DCM_ImageType, ImageType).good())
			{
				cout << *lit << endl;
				throw(string("Error: cannot access ImageType!"));
			}

			//cout << ImageType << endl;

			if (!ImageType.compare("DERIVED\\SECONDARY\\OTHER\\CSA BLACK IMAGE")) continue;
			if (!ImageType.compare("ORIGINAL\\SECONDARY\\SCREEN SAVE")) continue;

			if (!fileformat.getDataset()->findAndGetUint16(DCM_Rows, Rows).good())
			{
				cout << *lit << endl;
				throw(string("Error: cannot access Rows!"));
			}

			//cout << Rows << endl;

			if (!fileformat.getDataset()->findAndGetUint16(DCM_Columns, Columns).good())
			{
				cout << *lit << endl;
				throw(string("Error: cannot access Columns!"));
			}

			//cout << Columns << endl;

			if (!fileformat.getDataset()->findAndGetOFStringArray(DCM_SeriesDescription, SerieDescription).good())
				//throw(string("Error: cannot access SeriesDescription!"));
				cout << "Warning: cannot access SeriesDescription!" << endl;

			//cout << SerieDescription << endl;

			if (!ImageType.compare(0, 7, "DERIVED"))
				if (!fileformat.getDataset()->findAndGetOFStringArray(DCM_DerivationDescription, DerivationDescription).good())
				{
					//throw(string("Error: cannot access DerivationDescription!"));
					cout << "Warning: cannot access DerivationDescription!" << endl;
				}

			//cout << DerivationDescription << endl;

			for (int i=0; i<6; i++)
			{
				/*OFString tmpString;

				  if (!fileformat.getDataset()->findAndGetOFString(DCM_ImageOrientationPatient, tmpString, i).good())
				  throw(string("Error: cannot access ImageOrientationPatient!"));

				  cout << tmpString << endl;*/

				if (!fileformat.getDataset()->findAndGetFloat64(DCM_ImageOrientationPatient, ImageOrientationPatient[i], i).good())
				{
					cout << "Warning: cannot access ImageOrientationPatient!" << endl;
					stop = true;
					break;
				}

				//cout << ImageOrientationPatient[i] << endl;
			}

			if (stop) continue;

			for (int i=0; i<2; i++)
			{
				if (!fileformat.getDataset()->findAndGetFloat64(DCM_PixelSpacing, PixelSpacing[i], i).good())
				{
					cout << *lit << endl;
					throw(string("Error: cannot access PixelSpacing!"));
					//PixelSpacing[i] = 1;
				}

				//cout << PixelSpacing[i] << endl;
			}

			/**known serie?**/
			bool knownSerie = false;
			if (!seriesDescription.empty())
			{
				int listeSeriesSize = seriesDescription.size();
				for (int n=0; n<listeSeriesSize; n++)
				{
					bool sameOrientation = true;
					for (int i=0; i<6; i++)
					{
						//if (ImageOrientationPatient[i] != seriesImageOrientationPatient[n][i])
						if (fabs(ImageOrientationPatient[i] - seriesImageOrientation[n][i]) >= 1e-6)
						{
							sameOrientation = false;
							break;
						}
					}
					//cout << "sameOrientation = " << sameOrientation << endl;

					bool sameSpacing = true;
					//if (PixelSpacing[0] != seriesPixelSpacing[n].first || PixelSpacing[1] != seriesPixelSpacing[n].second)
					if ( fabs(PixelSpacing[1] - seriesPixelSpacing[n].first) >= 1e-6 ||
							fabs(PixelSpacing[0] - seriesPixelSpacing[n].second) >= 1e-6 )
						sameSpacing = false;

					if (	!SerieDescription.compare(OFString(seriesDescription[n].c_str())) &&
							Rows == seriesImageSizes[n].first &&
							Columns == seriesImageSizes[n].second &&
							!ImageType.compare(OFString(seriesImageType[n].c_str())) &&
							sameOrientation &&
							sameSpacing )
					{						
						seriesImagesFilenames[n].push_back(*lit);

						knownSerie = true;
						break;
					}
				}
			}

			/**unknown serie? -> add it to the list**/
			if (!knownSerie)
			{
				seriesDescription.push_back(SerieDescription.c_str());
				seriesImageSizes.push_back( make_pair(Rows, Columns) );
				seriesImageType.push_back(ImageType.c_str());

				vector<double> tmp;
				for (int i=0; i<6; i++) tmp.push_back(ImageOrientationPatient[i]);
				seriesImageOrientation.push_back( tmp );

				seriesPixelSpacing.push_back( make_pair(PixelSpacing[1], PixelSpacing[0]) );

				seriesImagesFilenames.push_back( list<string>(1, *lit) );

				if (!ImageType.compare(0, 8, "ORIGINAL")) seriesDerivationDescription.push_back("NA");
				else seriesDerivationDescription.push_back(DerivationDescription.c_str());
			}

			//cout << endl;
		}
		else
		{
			cout << *lit << endl;

			ostringstream message;
			message << "Error: cannot read DICOM file (" << status.text() << ")";
			throw(message.str());
		}
	} // all the listed files

	/// nb timeframes ///

	status = fileformat.loadFile(seriesImagesFilenames.front().front().c_str());
	Sint32 value;

	if (status.good())
	{
		if (!fileformat.getDataset()->findAndGetSint32(DCM_CardiacNumberOfImages, value).good())
		{
			cout << "Error: cannot access CardiacNumberOfImages!" << endl;
			cout << "Enter the number of timeframes:" << endl;
			cin >> nbFrames;
		}
		else if (value == 0)
		{
			cout << "Error: CardiacNumberOfImages = 0 !" << endl;
			cout << "Enter the number of timeframes:" << endl;
			cin >> nbFrames;
		}
		else
			nbFrames = value;

		cout << "number of timeframes: " << nbFrames << endl;
	}
	else
	{
		cout << seriesImagesFilenames.front().front() << endl;

		ostringstream message;
		message << "Error: cannot read DICOM file (" << status.text() << ")";
		throw(message.str());
	}

	seriesNbSlices.resize( seriesDescription.size(), 0 ); //will be filled by VolumeData_dicom::readSerie()
}

void VolumeData_dicom::printListedSeries()
{
	cout << "listed series:" << endl;

	cout << left;
	cout << setw(15) << "Serie number" << setw(25) << "Serie description" << setw(40) << "Image type and Derivation description" << endl;
	cout << setw(15) << "Nb of images" << setw(25) << "Images size (rows x cols)" << "Pixel spacing (rows x cols)" << endl;
	cout << "Images orientation" << endl;
	cout << endl;

	int listeSeriesSize = seriesDescription.size();
	for (int n=0; n<listeSeriesSize; n++)
	{
		cout << setw(15) << n << setw(25) << seriesDescription[n] << setw(40) << seriesImageType[n] << seriesDerivationDescription[n] << endl;
		cout << setw(15) << seriesImagesFilenames[n].size() << setw(25) << seriesPixelSpacing[n].first << " x " << seriesPixelSpacing[n].second << endl;
		for (int i=0; i<6; i++) cout << seriesImageOrientation[n][i] << " ";
		cout << endl << endl;
	}
}

void VolumeData_dicom::printSelectedSeries()
{
	cout << "selected series:" << endl;

	cout << left;
	cout << setw(15) << "Serie number" << setw(25) << "Serie description" << setw(40) << "Image type and Derivation description" << endl;
	cout << setw(15) << "Nb of images" << setw(25) << "Images size (rows x cols)" << "Pixel spacing (rows x cols)" << endl;
	cout << "Images orientation" << endl;
	cout << endl;

	set<int>::const_iterator lit(selectedSeries.begin()), lend(selectedSeries.end());

	for (; lit != lend; ++lit)
	{
		cout << setw(15) << *lit << setw(25) << seriesDescription[*lit] << setw(40) << seriesImageType[*lit] << seriesDerivationDescription[*lit] << endl;
		cout << setw(15) << seriesImagesFilenames[*lit].size() << setw(25) << seriesPixelSpacing[*lit].first << " x " << seriesPixelSpacing[*lit].second << endl;
		for (int i=0; i<6; i++) cout << seriesImageOrientation[*lit][i] << " ";
		cout << endl << endl;
	}
}

void VolumeData_dicom::readSerie(int serieNumber)
{
	cout << "reading serie " << serieNumber << "..." << endl;

	DcmFileFormat fileformat;
	OFCondition status;

	Vector3 rowVector(seriesImageOrientation[serieNumber][0],
			seriesImageOrientation[serieNumber][1],
			seriesImageOrientation[serieNumber][2]);

	Vector3 colVector(seriesImageOrientation[serieNumber][3],
			seriesImageOrientation[serieNumber][4],
			seriesImageOrientation[serieNumber][5]);

	Vector3 normalVector = vectProduct(rowVector, colVector);

	/** list the possible slice positions and sort them **/
	/** some slice positions may have been acquired several times -> treat them as different slices an choose later which one to use **/
	Vector3 ImagePositionPatient;
	double distSlice;
	pair<double, string> dist_time_slice;
	set< pair<double, string> > dists; // distance of the slice to the origine of the coordinates system along the normal axis: used to sort the slices
	set< pair<double, string> >::const_iterator sitdists, senddists;
	OFString ofstr;

	map<double, Vector3> dists_posPatient;

	list<string>::const_iterator lit(seriesImagesFilenames[serieNumber].begin()), lend(seriesImagesFilenames[serieNumber].end());
	for (; lit!=lend; ++lit)
	{
		status = fileformat.loadFile( (*lit).c_str() );

		if (status.good())
		{
			//which slice is it?
			for (int i=0; i<3; i++)
			{
				if (!fileformat.getDataset()->findAndGetFloat64(DCM_ImagePositionPatient, ImagePositionPatient[i], i).good())
					throw(string("Error: cannot access ImagePositionPatient!"));
			}

			distSlice = normalVector * ImagePositionPatient;

			fileformat.getDataset()->findAndGetOFString(DCM_SeriesTime, ofstr);
			dist_time_slice = make_pair(distSlice, string(ofstr.c_str()));

			dists.insert(dist_time_slice);

			dists_posPatient[distSlice] = ImagePositionPatient;
		}
		else
		{
			ostringstream message;
			message << "Error: cannot read DICOM file (" << status.text() << ")";
			throw(message.str());
		}
	}

	int nbSlicesTmp = dists.size();

	/**sitdists = dists.begin();
	  senddists = dists.end();
	  for (; sitdists != senddists; ++sitdists)
	  cout << setprecision(10) << sitdists->first << endl;

	  map<double, Vector3>::const_iterator sitTest(dists_posPatient.begin()), sendTest(dists_posPatient.end());
	  for (; sitTest !=sendTest; ++sitTest)
	  cout << setprecision(8) << sitTest->first << " " << sitTest->second << endl;

	  cout << ((--(dists.end()))->first - dists.begin()->first) / (dists.size()-1) << endl;

	  cv::waitKey()**/

	/** list the possible trigger times and sort them **/
	int posSlice;
	double TriggerTime;

	senddists = dists.end();

	seriesTriggerTimes[serieNumber].resize(nbSlicesTmp);

	for (lit = seriesImagesFilenames[serieNumber].begin(); lit!=lend; ++lit)
	{
		status = fileformat.loadFile( (*lit).c_str() );

		if (status.good())
		{
			//which slice is it?
			for (int i=0; i<3; i++)
			{
				if (!fileformat.getDataset()->findAndGetFloat64(DCM_ImagePositionPatient, ImagePositionPatient[i], i).good())
					throw(string("Error: cannot access ImagePositionPatient!"));
			}

			distSlice = normalVector * ImagePositionPatient;

			fileformat.getDataset()->findAndGetOFString(DCM_SeriesTime, ofstr);
			dist_time_slice = make_pair(distSlice, string(ofstr.c_str()));

			for (posSlice=0, sitdists=dists.begin(); sitdists!=senddists; ++sitdists, posSlice++)
				if (*sitdists == dist_time_slice) break;

			//which timeframe is it?
			if (!fileformat.getDataset()->findAndGetFloat64(DCM_TriggerTime, TriggerTime).good())
			{
				if (nbFrames != 1)
					throw(string("Error: cannot access TriggerTime!"));
				else
				{
					cout << "Warning: cannot access TriggerTime!" << endl;
					TriggerTime = 0;
				}
			}

			seriesTriggerTimes[serieNumber][posSlice].insert(TriggerTime);
		}
		else
		{
			ostringstream message;
			message << "Error: cannot read DICOM file (" << status.text() << ")";
			throw(message.str());
		}
	}


	/** read the images **/
	int posTime, Rows = seriesImageSizes[serieNumber].first, Columns = seriesImageSizes[serieNumber].second;
	int nbFrameTmp;
	unsigned long count;
	const unsigned short int* img;
	set<double>::const_iterator sittimes, sendtimes = seriesTriggerTimes[serieNumber][posSlice].end();




	seriesDataValues[serieNumber].resize(nbSlicesTmp);
	seriesImagePosition[serieNumber].resize(nbSlicesTmp);
	seriesImagesFilenamesSorted[serieNumber].resize(nbSlicesTmp);
	for (int i=0; i<nbSlicesTmp; i++)
	{
		nbFrameTmp = seriesTriggerTimes[serieNumber][posSlice].size();

		seriesDataValues[serieNumber][i].resize(nbFrameTmp);
		seriesImagesFilenamesSorted[serieNumber][i].resize(nbFrameTmp, string(""));
		for (int j=0; j<nbFrameTmp; j++)
			seriesDataValues[serieNumber][i][j].resize(Rows * Columns);
	}

	for (lit = seriesImagesFilenames[serieNumber].begin(); lit != lend; ++lit)
	{
		status = fileformat.loadFile( (*lit).c_str() );
		//cout << (*lit).c_str() << endl;

		if (status.good())
		{
			//which slice is it?
			for (int i=0; i<3; i++)
			{
				if (!fileformat.getDataset()->findAndGetFloat64(DCM_ImagePositionPatient, ImagePositionPatient[i], i).good())
					throw(string("Error: cannot access ImagePositionPatient!"));
			}

			distSlice = normalVector * ImagePositionPatient;

			fileformat.getDataset()->findAndGetOFString(DCM_SeriesTime, ofstr);
			dist_time_slice = make_pair(distSlice, string(ofstr.c_str()));

			for (posSlice=0, sitdists=dists.begin(); sitdists!=senddists; ++sitdists, posSlice++)
				if (*sitdists == dist_time_slice) break;

			//which timeframe is it?
			if (!fileformat.getDataset()->findAndGetFloat64(DCM_TriggerTime, TriggerTime).good())
			{
				if (nbFrames != 1)
					throw(string("Error: cannot access TriggerTime!"));
				else
				{
					//cout << "Warning: cannot access TriggerTime!" << endl;
					TriggerTime = 0;
				}
			}

			for (posTime=0, sittimes=seriesTriggerTimes[serieNumber][posSlice].begin(); sittimes!=sendtimes; ++sittimes, posTime++)
				if (*sittimes == TriggerTime) break;

			if (posTime == 0) seriesImagePosition[serieNumber][posSlice] = Vector3(ImagePositionPatient);

			if (seriesImagesFilenamesSorted[serieNumber][posSlice][posTime].compare(string("")))
				throw string("Error: several images at the same space-time position!");

			//read and save the image and its properties
			if (!fileformat.getDataset()->findAndGetUint16Array(DCM_PixelData, img, &count).good())
				throw(string("Error: cannot access PixelData!"));

			if (count != Rows * Columns)
				throw(string("Error: wrong pixel format!"));

			for(int r=0; r<Rows; r++)
				for(int c=0; c<Columns; c++)
				{
					float val = (float)(img[c + Columns * r]);

					if (val > 32767) val -= 65536;

					seriesDataValues[serieNumber][posSlice][posTime][c + Columns * r] = val;
				}

			seriesImagesFilenamesSorted[serieNumber][posSlice][posTime] = *lit;




		}
		else
		{
			ostringstream message;
			message << "Error: cannot read DICOM file (" << status.text() << ")";
			throw(message.str());
		}
	}







	/// resolve conflicts: several acquisitions for one slice ///

	map<double, list<int> > dists_posSlice; //for each distance, the list of corresponding indexes in dists and the other arrays
	map<double, list<int> >::const_iterator mit, mend;
	list<int>::const_iterator litPosSlice, lendPosSlice;

	list<int> rejectedSlices;

	for (posSlice=0, sitdists=dists.begin(); sitdists != senddists; ++sitdists, posSlice++)
		dists_posSlice[sitdists->first].push_back(posSlice);

	listDistances[serieNumber].clear();

	mit = dists_posSlice.begin();
	mend = dists_posSlice.end();

	for (; mit != mend; ++mit)
	{
		if (mit->second.size() != 1)
		{
			cout << endl;
			cout << "Warning: several acqusitions for slice position " << mit->first << ": " << mit->second.size() << " different acquisition times" << endl;

			vector<int> rejectedSlicesTmp(mit->second.size());

			litPosSlice = mit->second.begin();
			lendPosSlice = mit->second.end();

			//display the first frame of all acquisitions to help choosing
			vector<myargs> userData(mit->second.size());
			for (int i=0; litPosSlice != lendPosSlice; ++litPosSlice, i++)
			{
				ostringstream titre;
				titre << "Choice #" << i;

				userData[i].ptr = this;
				userData[i].serie = serieNumber;
				userData[i].slice = *litPosSlice;
				userData[i].frame = 0;
				userData[i].win_name = titre.str();

				float maxVal = 0;
				int Rows = seriesImageSizes[serieNumber].first, Columns = seriesImageSizes[serieNumber].second;
				for(int r=0; r<Rows; r++)
					for(int c=0; c<Columns; c++)
					{
						if (seriesDataValues[serieNumber][*litPosSlice][0][c + Columns * r] > maxVal)
							maxVal = seriesDataValues[serieNumber][*litPosSlice][0][c + Columns * r];
					}

				int valueIni = (int)(maxVal / 2);
				printImage(serieNumber, *litPosSlice, 0, titre.str(), valueIni);
				cv::createTrackbar("range", titre.str(), &valueIni, (int)maxVal, updateRange, &userData[i]);

				rejectedSlicesTmp[i] = *litPosSlice;
			}

			int reponse;
			bool stop = false;

			while (!stop)
			{
				cout << endl;
				cout << "which slice should be used?" << endl;
				cout << "-1: display settings" << endl;

				cin >> reponse;

				if(reponse == -1)
				{
					cout << "Set the position and grey value range of the images, then press a key (while images selected) to continue" << endl;
					cv::waitKey();
					cout << "display setting done" << endl;
				}
				else if (reponse >= 0 && reponse < mit->second.size())
				{
					for (int i=0; i<rejectedSlicesTmp.size(); i++)
					{
						if (i == reponse) continue;
						else rejectedSlices.push_back( rejectedSlicesTmp[i] );
					}

					stop = true;
				}
				else
				{
					cout << "Please, choose between 0 and " << mit->second.size() - 1 << " or enter -1 for display settings" << endl;
				}
			}

			//close windows
			for (int i=0; i < mit->second.size(); i++)
			{
				ostringstream titre;
				titre << "Choice #" << i;
				cvDestroyWindow( titre.str().c_str() );
				cv::waitKey(50);
			}
		} //if more than one acquisition

		listDistances[serieNumber].insert(mit->first);

	} //each slice

	// erase unused images
	vector< set< double> >::iterator vitTT;
	vector< vector< vector<float> > >::iterator vitDV;
	vector<Vector3>::iterator vitPP;
	vector< vector<string> >::iterator vitFS;
	vector< vector<int> >::iterator vitPV;

	list<int>::const_iterator litRS(rejectedSlices.end()), lendRS(rejectedSlices.begin());

	bool stop3;
	if (litRS == lendRS) stop3 = true;
	else stop3 = false;

	while (!stop3)
	{
		--litRS;

		vitTT = seriesTriggerTimes[serieNumber].begin();
		for (int i=0; i<*litRS; i++) ++vitTT;
		seriesTriggerTimes[serieNumber].erase(vitTT);

		vitDV = seriesDataValues[serieNumber].begin();
		for (int i=0; i<*litRS; i++) ++vitDV;
		seriesDataValues[serieNumber].erase(vitDV);

		vitPP = seriesImagePosition[serieNumber].begin();
		for (int i=0; i<*litRS; i++) ++vitPP;
		seriesImagePosition[serieNumber].erase(vitPP);

		vitFS = seriesImagesFilenamesSorted[serieNumber].begin();
		for (int i=0; i<*litRS; i++) ++vitFS;
		seriesImagesFilenamesSorted[serieNumber].erase(vitFS);

		if (litRS == lendRS) stop3 = true;
	}


	map<int, map<int, string> > seriesNumberMap;
	map<int, vector< vector< string > > >::iterator msit = seriesImagesFilenamesSorted.begin();
	while(msit != seriesImagesFilenamesSorted.end())
	{
		int series = msit->first;
		vector<vector<string> > slices = msit->second;

		for(unsigned int i = 0; i < slices.size(); i++)
		{
			// get an image 
			string fname = slices[i].front();
			status = fileformat.loadFile( fname.c_str() );

			// get the series number 
			const char *finalSeriesNumber;
			if (!fileformat.getDataset()->findAndGetString(DCM_SeriesNumber, finalSeriesNumber).good())
				throw(string("Error: cannot access Series number"));

			seriesNumberMap[series][i] = finalSeriesNumber;



		}

		++msit;
	}


	// print out the look up 
	ofstream outFile;
	outFile.open("series_lookup.txt");
	map<int, map<int, string> >::iterator sit = seriesNumberMap.begin();
	while(sit != seriesNumberMap.end())
	{
		map<int, string>::iterator sit2 = sit->second.begin();
		int seriesId = sit->first;
		while(sit2 != sit->second.end())
		{
			int sliceId = sit2->first;
			string serNum = sit2->second;

			outFile << "image (" << seriesId << "," << sliceId << ")" << " = " << serNum << endl;
			sit2++;
		}

		sit++;
	}

	outFile.close();



	/// nb slices

	seriesNbSlices[serieNumber] = listDistances[serieNumber].size();

	cout << "serie read" << endl;
}

void VolumeData_dicom::add_slices_to_serie(int serieNumber, set< pair<double, pair<int, int> > >& dists)
{
	VolumeData::add_slices_to_serie(serieNumber, dists);

	set< pair<double, pair<int, int> > >::const_iterator sitdists(dists.begin()), senddists(dists.end());
	int serie, slice;

	for (int s=0; sitdists!=senddists; ++sitdists, s++)
	{
		serie = sitdists->second.first;
		slice = sitdists->second.second;

		if (serie == serieNumber) continue;

		for (int f=0; f<nbFrames; f++) seriesImagesFilenames[serieNumber].push_back( seriesImagesFilenamesSorted[serie][slice][f] );

		seriesTriggerTimes[serieNumber].insert(seriesTriggerTimes[serieNumber].begin()+s, seriesTriggerTimes[serie][slice]); //serie, slice, timeframe, triggerTime
		seriesImagesFilenamesSorted[serieNumber].insert(seriesImagesFilenamesSorted[serieNumber].begin()+s, seriesImagesFilenamesSorted[serie][slice]); //serie, slice, timeframe, name
	}
}


/// VolumeData_dicom_cardiac ///

VolumeData_dicom_cardiac::VolumeData_dicom_cardiac(XMLreader& _xml_reader):
	VolumeData_dicom(_xml_reader)
{
	periodic_motion = true;
}

void VolumeData_dicom_cardiac::selectSeries()
{
	//cout << "VolumeData_dicom_cardiac::selectSeries()" << endl;

	selectedSeries.clear();

	int listeSeriesSize = seriesDescription.size();
	for (int n=0; n<listeSeriesSize; n++)
	{
		if (!seriesDescription[n].compare(0, 7, "Cine_SA") || !seriesDescription[n].compare(0, 12, "Cine_ipat SA") ||
				!seriesDescription[n].compare(0, 18, "Cine_ipat_retro_SA") || !seriesDescription[n].compare(0, 7, "cine_SA") )
		{
			if (!selectedSeries.empty()) cout << "Warning: several SA series detected" << endl;

			selectedSeries.insert(n);
			seriesSA.push_back(n);
		}
		else
			if (!seriesDescription[n].compare(0, 5, "Cine_") || !seriesDescription[n].compare(0, 5, "cine_"))
				selectedSeries.insert(n);
	}

	VolumeData::selectSeries();
}

void VolumeData_dicom_cardiac::choose_serie_ref()
{
	if ( seriesSA.empty() ) return VolumeData::choose_serie_ref();

	if (seriesSA.size() == 1)
	{
		serie_ref = seriesSA.front();
		return;
	}

	//several series named "SA" have been found, so we let the user decide which one to use
	list<int>::iterator lit, lend, litComp;

	bool compatible = false;

	while (!compatible)
	{
		cout << endl;
		cout << "Several SA series found!" << endl;

		cout << left;
		cout << setw(15) << "Serie number" << setw(25) << "Serie description" << setw(40) << "Image type" << endl;
		cout << setw(15) << "Nb of images" << setw(25) << "Images size (rows x cols)" << "Pixel spacing (rows x cols)" << endl;
		cout << "Images orientation" << endl << endl;

		vector<myargs> userData( seriesSA.size() );
		vector<int> ini_range( seriesSA.size() );

		lit = seriesSA.begin();
		lend = seriesSA.end();
		for (int n=0; lit != lend; ++lit, n++)
		{
			cout << "Choice #" << n << endl;
			cout << setw(15) << *lit << setw(25) << seriesDescription[*lit] << setw(40) << seriesImageType[*lit] << endl;

			ostringstream tmp;
			tmp << seriesImageSizes[*lit].first << " x " << seriesImageSizes[*lit].second;

			cout << setw(15) << seriesNbSlices[*lit]*nbFrames << setw(25) << tmp.str() << seriesPixelSpacing[*lit].first << " x " << seriesPixelSpacing[*lit].second << endl;
			for (int i=0; i<6; i++) cout << seriesImageOrientation[*lit][i] << " ";
			cout << endl << endl;

			ostringstream nameDisp;
			nameDisp << "Choice #" << n;

			int sliceDisp = seriesNbSlices[*lit] / 2;

			userData[n].ptr = this;
			userData[n].serie = *lit;
			userData[n].slice = sliceDisp;
			userData[n].frame = 0;
			userData[n].win_name = nameDisp.str();

			float maxVal = 0;
			int Rows = seriesImageSizes[*lit].first, Columns = seriesImageSizes[*lit].second;
			for(int r=0; r<Rows; r++)
				for(int c=0; c<Columns; c++)
				{
					if (seriesDataValues[*lit][sliceDisp][0][c + Columns * r] > maxVal)
						maxVal = seriesDataValues[*lit][sliceDisp][0][c + Columns * r];
				}

			ini_range[n] = (int)maxVal;
			printImage(*lit, sliceDisp, 0, nameDisp.str(), ini_range[n]);
			cv::createTrackbar("range", nameDisp.str(), &ini_range[n], (int)maxVal, updateRange, &userData[n]);
		}

		bool stop;
		int reponse;

		lit = seriesSA.begin();
		for (int n=0; lit != lend; n++)
		{
			stop=false;

			while(!stop)
			{
				cout << endl;
				cout << "What should we do with serie #" << n << "?" << endl;
				cout << "1: keep it as SA serie" << endl;
				cout << "2: make it a LA serie" << endl;
				cout << "0: discard it" << endl;
				cout << "9: show display settings" << endl;

				cin >> reponse;

				switch(reponse)
				{
					case 1:
						cout << "Serie #"<< n <<" will be a SA serie" << endl;
						++lit;
						stop = true;
						break;

					case 2:
						cout << "Serie #"<< n <<" becomes now a LA serie" << endl;

						lit = seriesSA.erase(lit);
						lend = seriesSA.end();

						stop = true;
						break;

					case 0:
						cout << "Discarding serie #" << n << endl;

						selectedSeries.erase(*lit);
						lit = seriesSA.erase(lit);
						lend = seriesSA.end();

						stop = true;
						break;

					case 9:
						cout << "Set the position and grey value range of the images, then press a key (while images selected) to continue" << endl;
						cv::waitKey();
						cout << "display setting done" << endl;
						break;

					default:
						cout << "Please, choose between #0, #1 and #2 or reset display #9" << endl;
						break;
				} //switch
			} //while !stop (user interaction)
		} //foreach serieSA

		int taille = userData.size();

		for (int n=0; n<taille; n++)
		{
			ostringstream nameDisp;
			nameDisp << "Choice #" << n;
			cvDestroyWindow( nameDisp.str().c_str() );
			cv::waitKey(50);
		}
		cv::waitKey(10);

		//check that remaining SA series are compatible with each others before fusion
		compatible = true;

		lit = seriesSA.begin();
		lend = seriesSA.end();

		pair<int, int> dim = seriesImageSizes[*lit];
		vector<double> orientation = seriesImageOrientation[*lit];
		pair<double, double> pixelSp = seriesPixelSpacing[*lit];
		lit++;

		for (; lit != lend; ++lit)
		{
			if ( seriesPixelSpacing[*lit] != pixelSp )
			{
				cout << "Warning: different pixel spacing" << endl;
				compatible = false;
				break;
			}

			if ( seriesImageSizes[*lit] != dim )
			{
				cout << "Warning: different image sizes" << endl;
				compatible = false;
				break;
			}

			for (int i=0; i<6; i++)
			{
				if ( fabs(seriesImageOrientation[*lit][i] - orientation[i]) > 1e-5 )
				{
					cout << "Warning: different orientations" << endl;
					compatible = false;
					break;
				}
			}
		}

		if (compatible == false) cout << "SA series are not compatible!" << endl;

	} //while incompatible SA

	if ( seriesSA.empty() ) return VolumeData::choose_serie_ref();

	// fusion of remaining SA series

	set<int>::iterator lit_ini, lit_add, lend_selectedSeries(selectedSeries.end());

	while (seriesSA.size() != 1)
	{
		int serie_ini = seriesSA.front();
		int serie_to_add = seriesSA.back();

		lit_ini = selectedSeries.find(serie_ini);
		lit_add = selectedSeries.find(serie_to_add);

		if (lit_ini == lend_selectedSeries || lit_add == lend_selectedSeries)
			throw string("VolumeData_dicom_cardiac::choose_serie_ref() Error: unknown serie");

		fuse_series(lit_ini, lit_add);

		seriesSA.pop_back();
	}

	serie_ref = seriesSA.front();
}

void VolumeData_dicom_cardiac::buildListTimeframes()
{
	//cout << "VolumeData_dicom_cardiac::buildListTimeframes()" << endl;

	if ( xml_reader.use_first_timeframe_only() ) nbFrames = 1;

	timeOfFrame = new float*[nbImages];

	int image;

	for (image = 0; image < nbImages; image++)
		timeOfFrame[image] = new float[nbFrames];

	if (nbFrames <= 1)
	{
#pragma omp parallel for default(shared) private(image) schedule(dynamic) num_threads(NB_THREADS)
		for (image=0; image<nbImages; image++)
			timeOfFrame[image][0] = 0;

		size[3] = 1;

		return;
	}

	/// fill array ///
	//time = frame number, we don't care about the triggerTime because the frames are synchronised by ECG, so all the slices of a timeframe should show the same moment in the cardiac cycle

	int f;

#ifdef _OPENMP
#if GCC_VERSION >= 40400 && _OPENMP >= 200805

#pragma omp parallel for default(shared) private(image, f) schedule(dynamic) collapse(2) num_threads(NB_THREADS)
	for (image=0; image<nbImages; image++)
		for (f=0; f<nbFrames; f++)
			timeOfFrame[image][f] = f;

#else

#pragma omp parallel for default(shared) private(image, f) schedule(dynamic) num_threads(NB_THREADS)
	for (image=0; image<nbImages; image++)
		for (f=0; f<nbFrames; f++)
			timeOfFrame[image][f] = f;

#endif
#else
	for (image=0; image<nbImages; image++)
		for (f=0; f<nbFrames; f++)
			timeOfFrame[image][f] = f;
#endif

	size[3] = nbFrames;
}

bool VolumeData_dicom_cardiac::comprise_temporal_sequences()
{
	return true;
}

/// VolumeData_nifti ///

VolumeData_nifti::VolumeData_nifti(XMLreader& _xml_reader):
	VolumeData(_xml_reader)
{}

void VolumeData_nifti::listSeries()
{
	/** get the filenames **/
	DIR *rep = opendir (repertory.c_str());
	list<string> filenames;

	if (rep != NULL)
	{
		listFilenames(rep, filenames, repertory + "/");

		closedir (rep);
		rep = NULL;
	}
	else
	{
		cout << repertory << endl;
		throw(string("Error: cannot open the directory"));
	}

	/** read the images **/
	list<string>::const_iterator lit(filenames.begin()), lend(filenames.end());

	nifti_1_header hdr;
	FILE *fp;
	int ret;

	int nbDims, Rows, Columns, nbSlices;
	float spacing_r, spacing_c, spacing_z;

	for(;lit!=lend;++lit)
	{
		cout << *lit << endl;

		if ( lit->compare(lit->size()-4, 4, ".nii") ) continue;

		/// serie description
		string file = lit->substr(lit->find_last_of('/')+1);
		string description = file.substr(0, file.find('.'));
		seriesDescription.push_back(description);

		cout << description << endl;

		/// slice info

		///// read first 348 bytes of header   //

		fp = fopen(lit->c_str(),"r");
		if (fp == NULL)
		{
			fprintf(stderr, "\nError opening header file volumeData.nii for read\n");
			exit(1);
		}

		ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
		if (ret != 1)
		{
			fprintf(stderr, "\nError writing header file volumeData.nii\n");
			exit(1);
		}

		nbDims = hdr.dim[0];

		Rows = hdr.dim[1];
		Columns = hdr.dim[2];
		nbSlices = hdr.dim[3];

		if (nbDims == 4) nbFrames = hdr.dim[4];
		else nbFrames = 1;

		seriesImageSizes.push_back( make_pair(Rows, Columns) );
		seriesNbSlices.push_back(nbSlices);

		vector<double> tmp(6);
		tmp[0] = 1;
		tmp[1] = 0;
		tmp[2] = 0;
		tmp[3] = 0;
		tmp[4] = 1;
		tmp[5] = 0;
		seriesImageOrientation.push_back( tmp );

		spacing_r = hdr.pixdim[1];
		spacing_c = hdr.pixdim[2];
		spacing_z = hdr.pixdim[3];
		seriesPixelSpacing.push_back( make_pair(spacing_r, spacing_c) );

		list<string> filenames;
		filenames.push_back(*lit);
		seriesImagesFilenames.push_back(filenames);
	} // all the listed files
}

void VolumeData_nifti::readSerie(int serieNumber)
{
	cout << "VolumeData_nifti::readSerie() " << serieNumber << endl;

	nifti_1_header hdr;
	nifti1_extender pad={0,0,0,0};
	FILE *fp;
	int ret;

	///// read first 348 bytes of header   //

	string filename = seriesImagesFilenames[serieNumber].front();
	cout << filename << endl;

	fp = fopen(filename.c_str(),"r");
	if (fp == NULL)
	{
		fprintf(stderr, "\nError opening header file volumeData.nii for read\n");
		exit(1);
	}

	ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
	if (ret != 1)
	{
		fprintf(stderr, "\nError writing header file volumeData.nii\n");
		exit(1);
	}

	int nbSlices = hdr.dim[3];

	float spacing_s, spacing_t;

	spacing_s = hdr.pixdim[3];

	if (hdr.dim[0] == 4) spacing_t = hdr.pixdim[4];
	else spacing_t = 0;

	seriesImagePosition[serieNumber].resize(nbSlices);

	int ind;
	int Rows = seriesImageSizes[serieNumber].first;
	int Columns = seriesImageSizes[serieNumber].second;
	int pitch = Rows*Columns;

	for (int s=0; s<nbSlices; s++)
	{
		seriesImagePosition[serieNumber][s] = Vector3(0, 0, s*spacing_s);
		listDistances[serieNumber].insert(s*spacing_s);
	}

	///// read extender pad and image data   //
	ret = fread(&pad, 4, 1, fp);
	if (ret != 1)
	{
		fprintf(stderr, "\nError reading header file extension pad volumeData.nii\n");
		exit(1);
	}

	short* buffer = new short[hdr.dim[1]*hdr.dim[2]*hdr.dim[3]];

	ret = fread(&buffer[0], (size_t)(hdr.bitpix/8), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], fp);
	if (ret != hdr.dim[1] * hdr.dim[2] * hdr.dim[3])
	{
		fprintf(stderr, "\nError reading data from volumeData.nii\n");
		exit(1);
	}

	//values
	seriesDataValues[serieNumber].resize(nbSlices);
	seriesImagesFilenamesSorted[serieNumber].resize(nbSlices);

	for (int s=0; s<nbSlices; s++)
	{
		seriesDataValues[serieNumber][s].resize(nbFrames);
		seriesImagesFilenamesSorted[serieNumber][s].resize(nbFrames);

		for (int t=0; t<nbFrames; t++)
		{
			seriesDataValues[serieNumber][s][t].resize(Rows*Columns);

			ostringstream name;
			name << filename.c_str() << "_slice_" << s;
			seriesImagesFilenamesSorted[serieNumber][s][t] = name.str();

			for (int r=0; r<Rows; r++)
				for (int c=0; c<Columns; c++)
				{
					ind = c + Columns * r;

					seriesDataValues[serieNumber][s][t][ind] = (float)buffer[ind + (hdr.dim[3]/2) * pitch];
				}
		}
	}

	delete[] buffer;

	fclose(fp);
}

void VolumeData_nifti::add_slices_to_serie(int serieNumber, set< pair<double, pair<int, int> > >& dists)
{
	VolumeData::add_slices_to_serie(serieNumber, dists);

	set< pair<double, pair<int, int> > >::const_iterator sitdists(dists.begin()), senddists(dists.end());
	int serie, slice;

	for (int s=0; sitdists!=senddists; ++sitdists, s++)
	{
		serie = sitdists->second.first;
		slice = sitdists->second.second;

		if (serie == serieNumber) continue;

		seriesImagesFilenames[serieNumber].push_back( seriesImagesFilenames[serie].front() );
	}
}


/// VolumeData_nrrd ///

VolumeData_nrrd::VolumeData_nrrd(XMLreader& _xml_reader):
	VolumeData(_xml_reader)
{}

void VolumeData_nrrd::listSeries()
{
	/** get the filenames **/
	DIR *rep = opendir (repertory.c_str());
	list<string> filenames;

	if (rep != NULL)
	{
		listFilenames(rep, filenames, repertory + "/");

		closedir (rep);
		rep = NULL;
	}
	else
	{
		cout << repertory << endl;
		throw(string("Error: cannot open the directory"));
	}

	/** read the images **/
	list<string>::const_iterator lit(filenames.begin()), lend(filenames.end());

	int numDim, Rows, Columns, nbSlices;
	float spacing_r, spacing_c, spacing_z;

	for(;lit!=lend;++lit)
	{
		cout << *lit << endl;

		if ( lit->compare(lit->size()-5, 5, ".nrrd") ) continue;

		/// serie description
		string file = lit->substr(lit->find_last_of('/')+1);
		string description = file.substr(0, file.find('.'));
		seriesDescription.push_back(description);

		cout << description << endl;

		char* err;

		/** create a nrrd; at this point this is just an empty container **/
		Nrrd *nin = nrrdNew();

		/** read in the nrrd from file **/
		if (nrrdLoad(nin, lit->c_str(), NULL))
		{
			err = biffGetDone(NRRD);
			fprintf(stderr, "VolumeData_nrrd::listSeries(): trouble reading \"%s\":\n%s", lit->c_str(), err);
			free(err);
			return;
		}

		numDim = nin->dim;
		cout << numDim << " dimensions" << endl;
		if (numDim < 3 || numDim > 4) throw string("VolumeData_Medical(): Erreur: mauvais nombre de dimensions");

		size_t info[4];
		nrrdAxisInfoGet_nva(nin, nrrdAxisInfoSize, &info);
		Columns = info[0];
		Rows = info[1];
		nbSlices = info[2];

		if (numDim == 3) nbFrames = 1;
		else if (numDim == 4) nbFrames = info[3];

		seriesImageSizes.push_back( make_pair(Rows, Columns) );
		seriesNbSlices.push_back(nbSlices);

		/*if (numDim == 3)
		  {
		  nrrdAxisInfoGet_va(nin, nrrdAxisInfoSpacing, &spacing_c, &spacing_r, &spacing_z);
		  }
		  else if (numDim == 4)
		  {
		  nrrdAxisInfoGet_va(nin, nrrdAxisInfoSpacing, &spacing_c, &spacing_r, &spacing_z, &spacing_t);
		  }*/

		double info_spacing[4];
		nrrdAxisInfoGet_nva(nin, nrrdAxisInfoSpacing, &info_spacing);

		spacing_c = info_spacing[0];
		spacing_r = info_spacing[1];
		spacing_z = info_spacing[2];

		cout << "spacings: " << spacing_r << " " << spacing_c << " " << spacing_z << endl;

		seriesPixelSpacing.push_back( make_pair(spacing_r, spacing_c) );

		seriesImageOrientation.resize(1);
		seriesImageOrientation[0].resize(6);
		seriesImageOrientation[0][0] = 0;
		seriesImageOrientation[0][1] = 1;
		seriesImageOrientation[0][2] = 0;
		seriesImageOrientation[0][3] = 1;
		seriesImageOrientation[0][4] = 0;
		seriesImageOrientation[0][5] = 0;

		list<string> filenames;
		filenames.push_back(*lit);
		seriesImagesFilenames.push_back(filenames);
	}
}

void VolumeData_nrrd::readSerie(int serieNumber)
{
	cout << "VolumeData_nifti::readSerie() " << serieNumber << endl;

	string filename = seriesImagesFilenames[serieNumber].front();
	cout << filename << endl;

	char* err;

	/* create a nrrd; at this point this is just an empty container */
	Nrrd *nin = nrrdNew();

	/* read in the nrrd from file */
	if (nrrdLoad(nin, filename.c_str(), NULL))
	{
		err = biffGetDone(NRRD);
		fprintf(stderr, "VolumeData_nifti::readSerie(): trouble reading \"%s\":\n%s", filename.c_str(), err);
		free(err);
		return;
	}

	int ind;
	int Rows = seriesImageSizes[serieNumber].first;
	int Columns = seriesImageSizes[serieNumber].second;

	double info_spacing[4];
	nrrdAxisInfoGet_nva(nin, nrrdAxisInfoSpacing, &info_spacing);

	float spacing_z = info_spacing[2];

	seriesImagePosition[0].resize(seriesNbSlices[0]);

	Vector3 colVector(seriesImageOrientation[0][0],
			seriesImageOrientation[0][1],
			seriesImageOrientation[0][2]);

	Vector3 rowVector(seriesImageOrientation[0][3],
			seriesImageOrientation[0][4],
			seriesImageOrientation[0][5]);

	Vector3 normalVector = vectProduct(colVector, rowVector);

	for (int s=0; s<seriesNbSlices[0]; s++)
	{
		float coord_z = s * spacing_z;

		seriesImagePosition[0][s] = Vector3(0, 0, -coord_z);
		//seriesImagePosition[0][s] = Vector3(-coord_z, 0, 0);

		listDistances[0].insert(coord_z);
	}

	/// values

	double (*lup)(const void *, size_t I);
	size_t I;
	double val;

	lup = nrrdDLookup[nin->type];

	seriesDataValues[0].resize(seriesNbSlices[0]);
	seriesImagesFilenamesSorted[0].resize(seriesNbSlices[0]);

	for (int s=0; s<seriesNbSlices[0]; s++)
	{
		seriesDataValues[0][s].resize(nbFrames);
		seriesImagesFilenamesSorted[0][s].resize(nbFrames);

		for (int t=0; t<nbFrames; t++)
		{
			seriesDataValues[0][s][t].resize(Rows*Columns);

			ostringstream name;
			name << filename.c_str() << "_slice_" << s;
			seriesImagesFilenamesSorted[0][s][t] = name.str();

			bool NANfound = false;

			for (int r=0; r<Rows; r++)
				for (int c=0; c<Columns; c++)
				{
					I = c + Columns * (r + Rows * (s + seriesNbSlices[0] * t));
					ind = c + Columns * r;

					val = lup(nin->data, I);

					if (isinf(val) || isnan(val))
					{
						if (NANfound == false) cout << "VolumeData_nifti::readSerie() nrrd Error: slice " << s << " frame " << t << " : " << val << " -> value changed into VAL_INVALID" << endl;

						val = VAL_INVALID;
						NANfound = true;
					}

					seriesDataValues[0][s][t][ind] = (float)val;
				}
		}
	}

	nrrdNuke(nin);
}

void VolumeData_nrrd::add_slices_to_serie(int serieNumber, set< pair<double, pair<int, int> > >& dists)
{
	VolumeData::add_slices_to_serie(serieNumber, dists);

	set< pair<double, pair<int, int> > >::const_iterator sitdists(dists.begin()), senddists(dists.end());
	int serie, slice;

	for (int s=0; sitdists!=senddists; ++sitdists, s++)
	{
		serie = sitdists->second.first;
		slice = sitdists->second.second;

		if (serie == serieNumber) continue;

		seriesImagesFilenames[serieNumber].push_back( seriesImagesFilenames[serie].front() );
	}
}


/// utilities ///

void updateRange(int value, void* arg)
{
	myargs *arguments = (myargs *)arg;
	arguments->ptr->printImage(arguments->serie, arguments->slice, arguments->frame, arguments->win_name, value);
}

void mouse_callback( int event, int x, int y, int flags, void* param )
{
	mousse_args_box* arguments = (mousse_args_box*)param;

	if (arguments->actif == false) return;

	switch( event )
	{
		case CV_EVENT_MOUSEMOVE:
			if( arguments->drawing_box )
			{
				arguments->box.width = x-arguments->box.x;
				arguments->box.height = y-arguments->box.y;
			}
			break;

		case CV_EVENT_LBUTTONDOWN:
			arguments->drawing_box = true;
			arguments->box = cvRect( x, y, 0, 0 );
			break;

		case CV_EVENT_LBUTTONUP:
			arguments->drawing_box = false;
			if( arguments->box.width < 0 )
			{
				arguments->box.x += arguments->box.width;
				arguments->box.width *= -1;
			}
			if( arguments->box.height < 0 )
			{
				arguments->box.y += arguments->box.height;
				arguments->box.height *= -1;
			}

			cv::Mat tmp(*(arguments->img));
			arguments->img->copyTo(tmp);

			cv::rectangle(tmp, cvPoint(arguments->box.x, arguments->box.y), cvPoint(arguments->box.x+arguments->box.width,arguments->box.y+arguments->box.height), cvScalar(255, 0, 0) );

			cv::imshow("ROI selection", tmp);

			bool stop = false;
			char choice;

			while (!stop)
			{
				choice = cv::waitKey();

				switch (choice)
				{
					case 'v':
						arguments->validated = true;
						arguments->actif = false;
						stop = true;
						break;

					case 'a':
						arguments->cancelled = true;
						arguments->actif = false;
						stop = true;
						break;

					default:
						cout << "Choisissez 'v' ou 'a'" << endl;
						break;
				}
			}

			break;
	}
}

