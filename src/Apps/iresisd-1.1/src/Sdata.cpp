#include "Sdata.h"

using namespace std;

Sdata::Sdata(XMLreader* _xml_reader, VolumeData* _volumeData, vector<Phi*>& _listPhi, const int* _size_fft, const vector<int>& _borders):
listPhi(_listPhi), borders(_borders)
{
	xml_reader = _xml_reader;
	volumeData = _volumeData;
	size_fft = _size_fft;
	
	volume_size = volumeData->getVolumeSize();
	nbPixVol = volume_size[0] * volume_size[1] * volume_size[2];
	
	const set<int>& series = volumeData->get_selectedSeries();
	set<int>::const_iterator sit_series, send_series(series.end());
	
	/// create the SpeedFunctions
	
	int nbSF = xml_reader->get_nb_SFdata();
	
	if (nbSF == 0) // interactive user choice of SFs
	{
		query_speed_functions();
		
		nbSF = sf_type.size();
		sf_pointers.resize(nbSF);
	}
	else
	{
		sf_pointers.resize(nbSF);
		sf_type.resize(nbSF);
		
		for (int sf=0; sf < nbSF; sf++)
		{
			xml_reader->get_SFdata_params(sf, sf_type[sf], sf_params[sf]);
			cout << sf_type[sf] << " : " << sf_params[sf].size() << " params" << endl;
			
			float weight;
			if (nbSF == 1) weight = 1;
			else weight = xml_reader->get_SFdata_weight(sf);
			
			if (weight != -1) //all series have the same speed functions
			{
				for (sit_series=series.begin(); sit_series != send_series; ++sit_series)
				{
					sf_of_serie[*sit_series].push_back(make_pair(sf, weight));
				}
			}
		}
	}
	
	for (int sf=0; sf < nbSF; sf++) create_SF(sf);
	
	// associate SFs to series
	
	vector< pair<int, float> >::iterator vit, vend;
	
	for (sit_series=series.begin(); sit_series != send_series; ++sit_series)
	{
		if ( sf_of_serie[*sit_series].empty() ) query_speed_functions_of_serie(*sit_series);
	
		vit = sf_of_serie[*sit_series].begin();
		vend = sf_of_serie[*sit_series].end();
		
		for (; vit != vend; ++vit) sf_pointers[vit->first]->add_serie(*sit_series);
	}
	
	/// allocate the arrays for storing S^n and S
	
	nbImages = volumeData->getNbImages();
	nbFrames = volumeData->getNbFrames();

	S_image = new map<int, float>*[nbImages];
	
	for (int i=0; i<nbImages; i++)
		S_image[i] = new map<int, float>[nbFrames];
	
	nbPixFFT = size_fft[0] * size_fft[1] * size_fft[2] * size_fft[3];
	
	fft_array = (double*)fftw_malloc(sizeof(double) * nbPixFFT);

	const set<int>& listSeries = volumeData->get_selectedSeries();
	set<int>::const_iterator sit(listSeries.begin()), send(listSeries.end());
	int maxNbPixImg = 0;
	nbPixelsTot = 0;
	for (; sit != send; ++sit)
	{
		const pair<int, int>& imgSize = volumeData->getImageSize(*sit);
		int nbPixImg = imgSize.first * imgSize.second;
		if (nbPixImg > maxNbPixImg) maxNbPixImg = nbPixImg;
		nbPixelsTot += nbPixImg;
	}
	S_image_pixel = new float[maxNbPixImg];
	
	list_pixels = new set<int>*[nbImages];
	
	for (int i=0; i<nbImages; i++)
		list_pixels[i] = new set<int>[nbFrames];
}

Sdata::~Sdata()
{
	for (int i=0; i<nbImages; i++) delete[] S_image[i];
	delete[] S_image;

	delete[] S_image_pixel;
	
	for (int i=0; i<nbImages; i++) delete[] list_pixels[i];
	delete[] list_pixels;
	
	fftw_free(fft_array);
	
	vector<SpeedFunction_Data*>::iterator vit(sf_pointers.begin()), vend(sf_pointers.end());
	for (; vit != vend; ++vit) delete *vit;
}

void Sdata::query_speed_functions()
{
	bool all_functions_defined = false;
	
	while (!all_functions_defined)
	{
		cout << "Please choose speed function to create amongst the following:" << endl;
		cout << "0: Piecewise constant model" << endl << "1: GVF" << endl << "-1 to finish" << endl;
	
		bool stop = false;
		int sf;
		vector<float> params;
	
		while (!stop)
		{
			cin >> sf;
			
			switch (sf)
			{
				case -1:
					if ( !sf_type.empty() ) stop = true;
					else cout << "Please choose at least one speed function" << endl;
					break;
				case 0:
					sf_type.push_back("CV");
					break;
				case 1:
					sf_type.push_back("GVF");
					
					params.clear();
					params.push_back(query_sf_param("baloon_force"));
					sf_params[sf] = params;
					
					break;
				default:
					break;
			} //switch
		} //while
	} //while
}

float Sdata::query_sf_param(string param_name)
{
	cout << "Please enter value for param " << param_name << endl;

	float value;
	cin >> value;
	
	return value;
}

void Sdata::query_speed_functions_of_serie(int serie)
{
	cout << "Please choose speed functions and their weights for serie " << serie <<
		": [speed function number] [weight]" << endl << endl;
		
	int nbSlices = volumeData->getNbSlices(serie);
	volumeData->printImage(serie, nbSlices/2, 0, "Central slice");
	
	cout << "Speed functions availables:" << endl;
	
	int nbSF = sf_type.size();
	vector<float>::const_iterator vit, vend;
	
	for (int sf=0; sf<nbSF; sf++)
	{
		cout << sf << ": " << sf_type[sf];
		if ( !sf_params[sf].empty() )
		{
			cout << " with params ";
			vit = sf_params[sf].begin();
			vend = sf_params[sf].end();
			for (; vit != vend; ++vit) cout << *vit << " ";
		}
	}
	
	bool stop = false;
	int sf;
	float weight;
	
	while (!stop)
	{
		cin >> sf >> weight;
		
		switch (sf)
		{
			case -1:
				if ( !sf_of_serie[serie].empty() ) stop = true;
				else cout << "Please choose at least one speed function" << endl;
				break;
				
			default:
				if (sf >= nbSF) cout << "Please choose amongst available functions" << endl;
				else sf_of_serie[serie].push_back( make_pair(sf, weight) );
				break;
		} //switch
	} //while
}

void Sdata::create_SF(int sf)
{
	if ( !sf_type[sf].compare("CV") ) sf_pointers[sf] = new SpeedFunction_CV(volumeData, listPhi);
	else if ( !sf_type[sf].compare("GVF") )
	{
		if (sf_params[sf].size() != 1) throw string("Sdata::create_SF() Error: no baloon force provided for GVF");
		sf_pointers[sf] = new SpeedFunction_GVF(volumeData, listPhi, sf_params[sf][0]);
	}
	else throw string("Sdata() Error: unknown SpeedFunction type");
}

double* Sdata::get_S()
{
	return fft_array;
}

const map<int, float>& Sdata::get_S_image(int image, int timeframe)
{
	return S_image[image][timeframe];
}

void Sdata::initialise_iteration()
{
	vector<SpeedFunction_Data*>::const_iterator vit(sf_pointers.begin()), vend(sf_pointers.end());
	for (; vit != vend; ++vit) (*vit)->prepare_iteration();
}

void Sdata::compute_S(int numPhi, bool global_variant, bool select_connected_regions)
{
	cout << "Sdata::compute_S()" << endl;
	
	/// initialise the fft array
	
	int i;
	
	#pragma omp parallel for default(shared) private(i) schedule(dynamic) num_threads(NB_THREADS)
	for (i=0; i<nbPixFFT; i++)
		fft_array[i] = 0;
	
	/// compute the speeds S^n produced by the images
	
	int contourSize = listPhi[numPhi]->getContourSize();
	
	if (!global_variant && contourSize < nbPixelsTot) get_pixels_on_contour(numPhi);
	
	int display_f = (int)ceil(nbFrames/2);
	
	for (int image=0; image < nbImages; image++)
	{
		for (int f=0; f<nbFrames; f++)
		{
			S_image[image][f].clear();
			
			if (global_variant)
				compute_S_in_image_global_variant(image, f, numPhi, select_connected_regions);
			else
			{
				if (contourSize < nbPixelsTot)
				{
					if ( list_pixels[image][f].empty() ) continue;
					compute_S_on_contour_in_image(image, f, numPhi);
				}
				else compute_S_in_image_local_variant(image, f, numPhi);
			}
			
			if ( xml_reader->allow_display() == 2 && f == display_f ) display_S_image(image);
			
			compute_S_in_volume(image, f, numPhi, global_variant);
		}
	}
}

void Sdata::compute_S_in_image_local_variant(int image, int frame, int numPhi)
{
	waitingList.clear();

	int serie = volumeData->get_selectedImages()[image].first;
	const pair<int, int>& imgSize = volumeData->getImageSize(serie);
	int nbPixImg = imgSize.first * imgSize.second;
	int indImg;
	vector< pair<int, float> >::const_iterator vitSF, vendSF(sf_of_serie[serie].end());
	bool usable_image  = false;
	
	for (vitSF=sf_of_serie[serie].begin(); vitSF != vendSF; ++vitSF)
	{
		if ( sf_pointers[vitSF->first]->is_image_usable(image, frame) )
		{
			usable_image = true;
			break;
		}
	}
	
	if (!usable_image)
	{
		cout << "Sdata::compute_S_in_image_local_variant(): image " << image << " not usable (S undefined)" << endl;
		
		#pragma omp parallel for default(shared) private(indImg) schedule(dynamic) num_threads(NB_THREADS)
		for (indImg=0; indImg<nbPixImg; indImg++)
		{
			S_image_pixel[indImg] = 0;
		}
		
		return;
	}
	
	float time = volumeData->getTimeOfFrame(image, frame);
	Point3<float> point;
	int xM, yM, zM, tM, xP, yP, zP, tP, x, y, z, t;
	bool need_S;
	float speed, sum_speeds, count_speeds;
	S_max = 0;
	
	#pragma omp parallel for default(shared) private(indImg, point, xM, yM, zM, tM, xP, yP, zP, tP, x, y, z, t, need_S, speed, sum_speeds, count_speeds, vitSF) schedule(dynamic) num_threads(NB_THREADS)
	for (indImg=0; indImg<nbPixImg; indImg++)
	{
		S_image_pixel[indImg] = 0;
		
		point = volumeData->getImagePointCoordsInVolume(image, indImg);

		xM = (int)floor(point.x());
		yM = (int)floor(point.y());
		zM = (int)floor(point.z());
		tM = (int)floor(time);

		xP = (int)ceil(point.x());
		yP = (int)ceil(point.y());
		zP = (int)ceil(point.z());
		tP = (int)ceil(time);

		need_S = false;

		for (x=xM; x<=xP; x++)
		for (y=yM; y<=yP; y++)
		for (z=zM; z<=zP; z++)
		for (t=tM; t<=tP; t++)
		{
			if ( !volumeData->isPointInVolume( Point3<int>(x, y, z), t ) ) continue;
			
			if ( listPhi[numPhi]->isContourPoint( computeIndex(volume_size, x, y, z), t ) )
			{
				#pragma omp critical
				waitingList.insert( Point4<int>(x, y, z, t) );
				need_S = true;
			}
		}

		if (!need_S) continue;

		/// compute S using all SpeedFunctions
		
		sum_speeds = 0;
		count_speeds = 0;

		for (vitSF=sf_of_serie[serie].begin(); vitSF != vendSF; ++vitSF)
		{
			if ( !sf_pointers[vitSF->first]->is_image_usable(image, frame) ) continue;
			
			speed = sf_pointers[vitSF->first]->get_speed(image, frame, indImg, numPhi);
			
			if (speed == VAL_INVALID) continue;

			sum_speeds += speed * vitSF->second;
			count_speeds += vitSF->second;
		}
		
		if (!count_speeds) continue;
		
		#pragma omp critical
		if ( !isfinite(sum_speeds) )
		{
			cout << "Sdata::compute_S_in_image_local_variant() Error: S = " << sum_speeds << endl;
			int wait;
			cin >> wait;
		}
		
		S_image_pixel[indImg] = sum_speeds / count_speeds;
		
		#pragma omp critical
		S_max = max(fabs( S_image_pixel[indImg] ), S_max);
	}
}

void Sdata::compute_S_on_contour_in_image(int image, int frame, int numPhi)
{
	waitingList.clear();

	int serie = volumeData->get_selectedImages()[image].first;
	const pair<int, int>& imgSize = volumeData->getImageSize(serie);
	int nbPixImg = imgSize.first * imgSize.second;
	int indImg;
	
	#pragma omp parallel for default(shared) private(indImg) schedule(dynamic) num_threads(NB_THREADS)
	for (indImg=0; indImg<nbPixImg; indImg++)
	{
		S_image_pixel[indImg] = 0;
	}
	
	bool usable_image  = false;
	vector< pair<int, float> >::const_iterator vitSF, vendSF(sf_of_serie[serie].end());
	
	for (vitSF=sf_of_serie[serie].begin(); vitSF != vendSF; ++vitSF)
	{
		if ( sf_pointers[vitSF->first]->is_image_usable(image, frame) )
		{
			usable_image = true;
			break;
		}
	}
	
	if (!usable_image)
	{
		cout << "Sdata::compute_S_on_contour_in_image(): image " << image << " not usable (S undefined)" << endl;
		return;
	}
	
	float time = volumeData->getTimeOfFrame(image, frame);
	vector<int> vector_pixels(list_pixels[image][frame].begin(), list_pixels[image][frame].end()); //convert into a vector for older versions of OpenMP
	int nbPix = vector_pixels.size();
	int i;
	Point3<float> point;
	int xM, yM, zM, tM, xP, yP, zP, tP, x, y, z, t;
	bool need_S;
	float speed, sum_speeds, count_speeds;
	
	S_max = 0;
	
	#pragma omp parallel for default(shared) private(i, indImg, point, xM, yM, zM, tM, xP, yP, zP, tP, x, y, z, t, need_S, speed, sum_speeds, count_speeds, vitSF) schedule(dynamic) num_threads(NB_THREADS)
	for (i=0; i<nbPix; i++)
	{
		indImg = vector_pixels[i];
		
		S_image_pixel[indImg] = 0;
		
		point = volumeData->getImagePointCoordsInVolume(image, indImg);

		xM = (int)floor(point.x());
		yM = (int)floor(point.y());
		zM = (int)floor(point.z());
		tM = (int)floor(time);

		xP = (int)ceil(point.x());
		yP = (int)ceil(point.y());
		zP = (int)ceil(point.z());
		tP = (int)ceil(time);

		need_S = false;

		for (x=xM; x<=xP; x++)
		for (y=yM; y<=yP; y++)
		for (z=zM; z<=zP; z++)
		for (t=tM; t<=tP; t++)
		{
			if ( !volumeData->isPointInVolume( Point3<int>(x, y, z), t ) ) continue;
			
			if ( listPhi[numPhi]->isContourPoint( computeIndex(volume_size, x, y, z), t ) )
			{
				#pragma omp critical
				waitingList.insert( Point4<int>(x, y, z, t) );
				need_S = true;
			}
		}

		if (!need_S) continue;

		/// compute S using all SpeedFunctions
		
		sum_speeds = 0;
		count_speeds = 0;

		for (vitSF=sf_of_serie[serie].begin(); vitSF != vendSF; ++vitSF)
		{
			if ( !sf_pointers[vitSF->first]->is_image_usable(image, frame) ) continue;
			
			speed = sf_pointers[vitSF->first]->get_speed(image, frame, indImg, numPhi);
			
			if (speed == VAL_INVALID) continue;

			sum_speeds += speed * vitSF->second;
			count_speeds += vitSF->second;
		}
		
		if (!count_speeds) continue;
		
		#pragma omp critical
		if ( !isfinite(sum_speeds) )
		{
			cout << "Sdata::compute_S_in_image_local_variant() Error: S = " << sum_speeds << endl;
			int wait;
			cin >> wait;
		}
		
		S_image_pixel[indImg] = sum_speeds / count_speeds;
		
		#pragma omp critical
		S_max = max(fabs( S_image_pixel[indImg] ), S_max);
	}
}

void Sdata::get_pixels_on_contour(int numPhi)
{
	int image, f;
	
	#pragma omp parallel for default(shared) private(image, f) schedule(dynamic) num_threads(NB_THREADS)
	for (image=0; image < nbImages; image++)
	for (f=0; f<nbFrames; f++)
		list_pixels[image][f].clear();
	
	const map<int, float>* contour = listPhi[numPhi]->getContour();
	map<int, float>::const_iterator mit, mend;
	set<int>::const_iterator sit, send;
	Point3<float> imageCoord;
	int xM, xP, yM, yP;
	int frameM, frameP;
	float distM, distP;
	
	for (int t=0; t<volume_size[3]; t++)
	{
		mit = contour[t].begin();
		mend = contour[t].end();
		
		for (; mit != mend; ++mit)
		{
			const set<int>& list_images = volumeData->get_images_at_voxel(mit->first);
			
			sit = list_images.begin();
			send = list_images.end();
			
			for (; sit != send; ++sit)
			{
				imageCoord = volumeData->getVoxelCoordsInImageVolume(mit->first, *sit);
				
				xM = (int)floor(imageCoord.x());
				xP = (int)ceil(imageCoord.x());
				
				yM = (int)floor(imageCoord.y());
				yP = (int)ceil(imageCoord.y());
				
				const pair<int, int>& imgSize = volumeData->getImageSize( volumeData->get_selectedImages()[*sit].first );

				if (xM < 0) xM = 0;
				if (xP < 0) xP = 0;
				if (xM >= imgSize.first) xM = imgSize.first - 1;
				if (xP >= imgSize.first) xP = imgSize.first - 1;
				if (yM < 0) yM = 0;
				if (yP < 0) yP = 0;
				if (yM >= imgSize.second) yM = imgSize.second - 1;
				if (yP >= imgSize.second) yP = imgSize.second - 1;
				
				volumeData->getClosestFrames(t, *sit, frameM, distM, frameP, distP);
				
				if (distM < 1)
				{
					list_pixels[*sit][frameM].insert( computeIndex(imgSize, xM, yM) );
					list_pixels[*sit][frameM].insert( computeIndex(imgSize, xM, yP) );
					list_pixels[*sit][frameM].insert( computeIndex(imgSize, xP, yM) );
					list_pixels[*sit][frameM].insert( computeIndex(imgSize, xP, yP) );
				}
				
				if (distP < 1)
				{
					list_pixels[*sit][frameP].insert( computeIndex(imgSize, xM, yM) );
					list_pixels[*sit][frameP].insert( computeIndex(imgSize, xM, yP) );
					list_pixels[*sit][frameP].insert( computeIndex(imgSize, xP, yM) );
					list_pixels[*sit][frameP].insert( computeIndex(imgSize, xP, yP) );
				}
			}
		}
	}
}

void Sdata::compute_S_in_image_global_variant(int image, int frame, int numPhi, bool do_selection_of_connected_regions)
{
	waitingList.clear();

	int serie = volumeData->get_selectedImages()[image].first;
	const pair<int, int>& imgSize = volumeData->getImageSize(serie);
	int nbPixImg = imgSize.first * imgSize.second;
	float time = volumeData->getTimeOfFrame(image, frame);
	
	/// first, compute S on the contour or in Omega_diff

	int indImg;
	Point3<float> point;
	int xM, yM, zM, tM, xP, yP, zP, tP, x, y, z, t;
	vector< pair<int, float> >::const_iterator vitSF, vendSF(sf_of_serie[serie].end());
	float speed, sum_speeds, count_speeds;
	bool isInVolume;
	int indVol;
	float phiVal;
	
	bool usable_image  = false;
	for (vitSF=sf_of_serie[serie].begin(); vitSF != vendSF; ++vitSF)
	{
		if ( sf_pointers[vitSF->first]->is_image_usable(image, frame) )
		{
			usable_image = true;
			break;
		}
	}
	
	if (!usable_image)
	{
		cout << "Sdata::compute_S_in_image_global_variant(): image " << image << " not usable (S undefined)" << endl;
		
		#pragma omp parallel for default(shared) private(indImg) schedule(dynamic) num_threads(NB_THREADS)
		for (indImg=0; indImg<nbPixImg; indImg++)
		{
			S_image_pixel[indImg] = 0;
		}
		
		return;
	}
	
	map<int, set< Point4<int> > > positiveSpeeds; //list of pixels with (possibly several) associated voxels
	
	S_max = 0;
	
	#pragma omp parallel for default(shared) private(indImg, point, isInVolume, indVol, xM, yM, zM, tM, xP, yP, zP, tP, x, y, z, t, vitSF, speed, sum_speeds, count_speeds, phiVal) schedule(dynamic) num_threads(NB_THREADS)
	for (indImg=0; indImg<nbPixImg; indImg++)
	{
		S_image_pixel[indImg] = 0;

		point = volumeData->getImagePointCoordsInVolume(image, indImg);

		xM = (int)floor(point.x());
		yM = (int)floor(point.y());
		zM = (int)floor(point.z());
		tM = (int)floor(time);

		xP = (int)ceil(point.x());
		yP = (int)ceil(point.y());
		zP = (int)ceil(point.z());
		tP = (int)ceil(time);

		isInVolume = false;

		for (x=xM; x<=xP; x++)
		for (y=yM; y<=yP; y++)
		for (z=zM; z<=zP; z++)
		{
			if ( !volumeData->isPointInVolume( Point3<int>(x, y, z) ) ) continue;

			isInVolume = true;
		}
		
		if (!isInVolume) continue;
		
		/// compute S using all SpeedFunctions

		sum_speeds = 0;
		count_speeds = 0;

		for (vitSF=sf_of_serie[serie].begin(); vitSF != vendSF; ++vitSF)
		{
			if ( !sf_pointers[vitSF->first]->is_image_usable(image, frame) ) continue;
			
			speed = sf_pointers[vitSF->first]->get_speed(image, frame, indImg, numPhi);
			
			if (speed == VAL_INVALID) continue;

			sum_speeds += speed * vitSF->second;
			count_speeds += vitSF->second;
		}
		
		if (!count_speeds) continue;
		
		speed = sum_speeds / count_speeds;
		
		for (x=xM; x<=xP; x++)
		for (y=yM; y<=yP; y++)
		for (z=zM; z<=zP; z++)
		{
			for (t=tM; t<=tP; t++)
			{
				indVol = computeIndex(volume_size, x, y, z);
				
				/// keep all speeds on the contour
				
				if ( listPhi[numPhi]->isContourPoint(indVol, t) )
				{
					#pragma omp critical
					{
						waitingList.insert( Point4<int>(x, y, z, t) );
						S_max = max(fabs(speed), S_max);
					}
					continue;
				}
				
				/// discard other speeds that are not in Omega_diff
				
				phiVal = listPhi[numPhi]->get_phi_value(indVol + t * nbPixVol);
				if (phiVal * speed >= 0) continue; //this is not Omega_diff
				
				if (do_selection_of_connected_regions)
				{
					// keep all remaining negative speeds...
					if (speed <= 0)
					{
						#pragma omp critical
						{
							waitingList.insert( Point4<int>(x, y, z, t) );
							S_max = max(fabs(speed), S_max);
						}
					}
					else // ... add positive speeds to the list of positive speeds waiting to be selected (for connected regions selection)
					{
						#pragma omp critical
						positiveSpeeds[indImg].insert( Point4<int>(x, y, z, t) );
					}
				}
				else
				{
					#pragma omp critical
					{
						waitingList.insert( Point4<int>(x, y, z, t) );
						S_max = max(fabs(speed), S_max);
					}
				}
			}
		}

		S_image_pixel[indImg] = speed;
	}
	
	/// select regions connected to the contour
	
	map<int, set< Point4<int> > >::const_iterator mit(positiveSpeeds.begin()), mend(positiveSpeeds.end());
	set<int> checked;
	set< Point4<int> >::const_iterator sit, send;
	bool contour_found;
	
	for (; mit != mend; ++mit)
	{
		if ( checked.count(mit->first) ) continue;
		
		sit = mit->second.begin();
		send = mit->second.end();
		contour_found = false;
		
		for (; sit != send; ++sit)
		{
			indVol = computeIndex(volume_size, sit->x(), sit->y(), sit->z());
			phiVal = listPhi[numPhi]->get_phi_value(indVol + sit->t() * nbPixVol);
			if (phiVal > 0) contour_found = true;
			if (contour_found) break;
		}
		
		if (contour_found) // select all positive speeds connected to this point
		{
			select_connected_region(positiveSpeeds, checked, mit, imgSize.second, nbPixImg);
		}
	}
}

void Sdata::select_connected_region(map<int, set< Point4<int> > >& positiveSpeeds, set<int>& checked, map<int, set< Point4<int> > >::const_iterator mit, int Nc, int nbPixImg)
{
	#pragma omp critical
	checked.insert(mit->first);
	
	/// add the corresponding 4D points to waitingList
	
	set< Point4<int> >::const_iterator sit(mit->second.begin()), send(mit->second.end());
	
	for (; sit != send; ++sit)
	{
		#pragma omp critical
		waitingList.insert(*sit);
	}
	
	float speed = S_image_pixel[mit->first];
	
	#pragma omp critical
	S_max = max(fabs(speed), S_max);
	
	/// check if the neighbouring pixels have positive speeds
	
	list<int> neighbours;
	
	neighbours.push_back(mit->first - 1);
	neighbours.push_back(mit->first - 1 - Nc);
	neighbours.push_back(mit->first - 1 + Nc);
	neighbours.push_back(mit->first - Nc);
	neighbours.push_back(mit->first + Nc);
	neighbours.push_back(mit->first + 1);
	neighbours.push_back(mit->first + 1 - Nc);
	neighbours.push_back(mit->first + 1 + Nc);
	
	map<int, set< Point4<int> > >::const_iterator mit2;
	list<int>::const_iterator lit(neighbours.begin()), lend(neighbours.end());
	
	for (; lit != lend; ++lit)
	{
		if (*lit < 0 || *lit >= nbPixImg) continue;
		
		if ( checked.count(*lit) ) continue;
		
		mit2 = positiveSpeeds.find(*lit);
		
		if (mit2 != positiveSpeeds.end())
		{
			select_connected_region(positiveSpeeds, checked, mit2, Nc, nbPixImg);
		}
	}
}

void Sdata::compute_S_in_volume(int image, int frame, int numPhi, bool global_variant)
{
	vector< Point4<int> > waitingList_vect(waitingList.begin(), waitingList.end());
	waitingList.clear();

	int serie = volumeData->get_selectedImages()[image].first;
	const pair<int, int>& imgSize = volumeData->getImageSize(serie);
	float time_of_frame = volumeData->getTimeOfFrame(image, frame);
	
	int indList, list_size = waitingList_vect.size();
	int indV;
	Point3<float> imageCoord;
	int x, y, z, t;
	int xM, yM, xP, yP;
	float ax, ay, az, at;
	float valMM, valMP, valPM, valPP, valCM, valCP, res;
	int indfft;
	
	int framePitch = volume_size[0] * volume_size[1] * volume_size[2];
	
	if (S_max == 0) S_max = 1;
	
	#pragma omp parallel for default(shared) private(indList, x, y, z, t, indV, imageCoord, xM, yM, xP, yP, ax, ay, az, at, res, valMM, valMP, valPM, valPP, valCM, valCP, indfft) schedule(dynamic) num_threads(NB_THREADS)
	for (indList=0; indList < list_size; indList++)
	{
		x = waitingList_vect[indList].x();
		y = waitingList_vect[indList].y();
		z = waitingList_vect[indList].z();
		t = waitingList_vect[indList].t();
		
		indV = computeIndex(volume_size, x, y, z);
		
		imageCoord = volumeData->getVoxelCoordsInImageVolume(indV, image);
		
		/// bi-linear interpolation of the pixels' S at imageCoord
		
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

		ax = imageCoord.x() - xM;
		ay = imageCoord.y() - yM;
		
		if (ax == 0 && ay == 0) // the projection of the point on the image is a pixel, so no interpolation needed
		{
			res = S_image_pixel[ computeIndex(imgSize, xM, yM) ];
		}
		else
		{
			valMM = S_image_pixel[ computeIndex(imgSize, xM, yM) ];
			valMP = S_image_pixel[ computeIndex(imgSize, xM, yP) ];
			valPM = S_image_pixel[ computeIndex(imgSize, xP, yM) ];
			valPP = S_image_pixel[ computeIndex(imgSize, xP, yP) ];
			
			valCM = ax * valPM + (1. - ax) * valMM;
			valCP = ax * valPP + (1. - ax) * valMP;
			
			res = ay * valCP + (1. - ay) * valCM;
			
			if ( !isfinite(res) )
			{
				cout << "Sdata::compute_S_in_volume() Error: S_1 = " << res << endl;
				int wait;
				cin >> wait;
			}
		}

		/// linear interpolation in the z direction, using the value at z=0 (computed above) and at z=+-1 (=0)

		az = fabs(imageCoord.z());
		res *= (1. - az);
		
		/// linear interpolation in the t direction

		at = fabs(time_of_frame - t);
		res *= (1. - at);
		
		if ( !isfinite(res) )
		{
			cout << "Sdata::compute_S_in_volume() Error: S_2 = " << res << endl;
			int wait;
			cin >> wait;
		}

		/// normalisation
		// just to make sure that all images have the same weight in the computation of S

		res /= S_max;
		
		if ( !isfinite(res) )
		{
			cout << "Sdata::compute_S_in_volume() Error: S = " << res << endl;
			int wait;
			cin >> wait;
		}
		
		/// save interpolation result
		
		if (global_variant)
		{
			if ( listPhi[numPhi]->isContourPoint(indV, t) )
			{
				indfft = (x + borders[0]) + size_fft[0] * ((y + borders[1]) + size_fft[1] * ((z + borders[2]) + size_fft[2] * (t + borders[3])));
		
				#pragma omp critical
				{
					S_image[image][frame][indV + t * framePitch] = res;
					fft_array[indfft] += res;
				}
			}
			else //only keep the speed for the computation of the registration
			{
				#pragma omp critical
				S_image[image][frame][indV + t * framePitch] = res;
			}
		}
		else
		{
			indfft = (x + borders[0]) + size_fft[0] * ((y + borders[1]) + size_fft[1] * ((z + borders[2]) + size_fft[2] * (t + borders[3])));
			
			#pragma omp critical
			{
				S_image[image][frame][indV + t * framePitch] = res;
				fft_array[indfft] += res;
			}
		}
	}
}

void Sdata::compute_S_everywhere(int numPhi, bool global_variant, bool select_connected_regions)
{
	cout << "Sdata::compute_S_everywhere()" << endl;
	
	/// initialise the fft array
	
	int i;
	
	#pragma omp parallel for default(shared) private(i) schedule(dynamic) num_threads(NB_THREADS)
	for (i=0; i<nbPixFFT; i++)
		fft_array[i] = 0;
	
	/// compute the speeds S^n produced by the images
	
	for (int image=0; image < nbImages; image++)
	{
		cout << image << endl;
		
		for (int f=0; f<nbFrames; f++)
		{
			if (!global_variant || select_connected_regions)
			{
				S_image[image][f].clear();
				compute_S_in_image_global_variant(image, f, numPhi, false);
			}
			
			compute_S_in_volume(image, f, numPhi, false);
		}
	}
	
	/// only keep the positive values at voxels where images are present
	
	
}

void Sdata::display_S_image(int image)
{
	const pair<int, int>& serieSlice = volumeData->get_selectedImages()[image];
	pair<int, int> imgSize = volumeData->VolumeData::getImageSize(serieSlice.first);
	int Rows = imgSize.first;
	int Columns = imgSize.second;
	
	cv::Mat imgDisp(Rows, Columns, CV_64FC1);
	
	float val;
	float valMax = -VAL_INVALID;
	int indImg;
	
	for (int i=0; i<Rows; i++)
	for (int j=0; j<Columns; j++)
	{
		indImg = computeIndex(imgSize, i, j);
		
		val = S_image_pixel[indImg];
		imgDisp.at<double>(i, j) = val;
		
		if (fabs(val) > valMax) valMax = fabs(val);
	}
	
	imgDisp /= 2.*valMax;
	imgDisp += 0.5;
	
	ostringstream name;
	name << "S_image_" << image;
	
	cv::namedWindow(name.str(), CV_WINDOW_AUTOSIZE);
	cv::imshow(name.str(), imgDisp);
	cv::waitKey(10);
}
