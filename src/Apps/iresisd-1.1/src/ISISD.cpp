#include "ISISD.h"

using namespace std;


ISISD::ISISD(Sdata* _S_data, Sgeom* _S_geom, vector<Phi*>& _listPhi, const vector<int>& _borders, VolumeData* _volumeData, int _display): listPhi(_listPhi), borders(_borders)
{
	S_data = _S_data;
	S_geom = _S_geom;
	
	volumeData = _volumeData;
	
	display = _display;
	
	speeds_geom = S_geom->get_S();
	
	size = listPhi[0]->get_volume_size();
	nbPixVol = size[0] * size[1] * size[2];
	nbPix = nbPixVol * size[3];
	
	int nbPhi = listPhi.size();
	oscillationsNumber = 10;
	contourLengths = new list<float>[nbPhi];
	maxTranslation = new float[nbPhi];
	contourStoppedGrowing = new bool[nbPhi];
	dynamic_maxTranslation = true;
	
	for (int i=0; i<nbPhi; i++)
	{
		maxTranslation[i] = 1;
		contourStoppedGrowing[i] = false;
	}
}

ISISD::~ISISD()
{
	delete[] contourLengths;
	delete[] maxTranslation;
	delete[] contourStoppedGrowing;
}

bool ISISD::compute_dt_data(int numPhi)
{
	int x, y, z, indfft;
	float dphival, normGrad, speed;
	float speedMax = 0;
	
	const map<int, float>* contour = listPhi[numPhi]->getContour();
	map<int, float>::const_iterator mit, mend;
	
	for (int t=0; t<size[3]; t++)
	{
		mit = contour[t].begin();
		mend = contour[t].end();
		
/*#ifdef _OPENMP
	#if GCC_VERSION >= 40400 && _OPENMP >= 200805
		#pragma omp parallel for default(shared) private(mit, x, y, z, indfft, dphival, normGrad, speed) schedule(dynamic) num_threads(NB_THREADS)
	#endif
#endif*/
		for (; mit != mend; ++mit)
		{
			computeCoords(mit->first, size, x, y, z);
			indfft = (x + borders[0]) + size_fft[0] * ((y + borders[1]) + size_fft[1] * ((z + borders[2]) + size_fft[2] * (t + borders[3])));
			
			dphival = speeds_data[indfft] / (float)nbPixFFT;
			
			normGrad = listPhi[numPhi]->computeNormGradientInNormalDirection(x, y, z, t);
			if (normGrad == 0) normGrad = 1;

			speed = dphival / normGrad;

			if ( !isfinite(speed) )
			{
				cout << "ISISD::compute_dt_data() Error: speed = " << speed << endl << "dphival = " << dphival << " normGrad = " << normGrad << endl;
				int reponse;
				cin >> reponse;
			}
			
			#pragma omp critical
			speedMax = max(speedMax, fabs(speed));
		}
	}
	
	if (speedMax == 0) return false;
	
	dt_data = maxTranslation[numPhi] / speedMax;
	
	cout << "timestep data : " << dt_data << " (maxTranslation = " << maxTranslation[numPhi] << ")" << endl;
	
	return true;
}

void ISISD::compute_dt_geom(int numPhi)
{
	if (speeds_geom == NULL) return;
	
	float Smax = 0;
	for (int i=0; i<nbPix; i++)
	{
		#pragma omp critical
		Smax = max(Smax, fabs(speeds_geom[i]));
	}
	
	dt_geom = 0.9 / Smax;
	
	if (dt_geom > 2) dt_geom = 2;
	cout << "Warning: dt_geom > 2 (" << dt_geom << ") -> setting to 2" << endl;
}

bool ISISD::auto_reinitialisation(int numPhi, bool global_variant, bool select_connected_regions)
{
	cout << "Warning: contour stopped moving" << endl;
	cout << "Possible reason: no more contour. Solution proposed:" << endl;
	cout << "1) reinitialise the contour automatically on the object (if object can be found)" << endl;
	cout << "2) stop the program" << endl;
	
	int choice;
	bool stop = false;
	
	while (!stop)
	{
		cin >> choice;
		switch (choice)
		{
			case 1:
				stop = true;
				break;
				
			case 2:
				cout << "Stopping the program..." << endl;
				return true;
				
			default:
				cout << "Please choose between options 1 and 2" << endl;
				break;
		}
	}
	
	/// very crude re-initialisation of phi, by replacing its values by S at voxels where images are present
	
	S_data->compute_S_everywhere(numPhi, global_variant, select_connected_regions);
	
	
	int x, y, z, t, ind_phi, indfft;
	
#ifdef _OPENMP
	#if GCC_VERSION >= 40400 && _OPENMP >= 200805
	#pragma omp parallel for default(shared) private(x, y, z, t, indfft, ind_phi) schedule(dynamic) collapse(4) num_threads(NB_THREADS)
	for (x=0; x<size[0]; x++)
	for (y=0; y<size[1]; y++)
	for (z=0; z<size[2]; z++)
	for (t=0; t<size[3]; t++)
	#else
	#pragma omp parallel for default(shared) private(x, y, z, t, indfft, ind_phi) schedule(dynamic) num_threads(NB_THREADS)
	for (ind_phi=0; ind_phi<nbPix; ind_phi++)
	#endif
#else
	for (x=0; x<size[0]; x++)
	for (y=0; y<size[1]; y++)
	for (z=0; z<size[2]; z++)
	for (t=0; t<size[3]; t++)
#endif
	{
#ifdef _OPENMP
	#if GCC_VERSION >= 40400 && _OPENMP >= 200805
		ind_phi = computeIndex(size, x, y, z, t);
	#else
		computeCoords(ind_phi, size, x, y, z, t);
	#endif
#else
		ind_phi = computeIndex(size, x, y, z, t);
#endif
		const set<int>& list_images = volumeData->get_images_at_voxel(ind_phi);
		if (list_images.size() <= 1) continue;
		
		indfft = (x + borders[0]) + size_fft[0] * ((y + borders[1]) + size_fft[1] * ((z + borders[2]) + size_fft[2] * (t + borders[3])));
		if (speeds_data[indfft] < 0) continue;
		
		listPhi[numPhi]->set_phi_value(ind_phi, speeds_data[indfft]);
	}
	
	
	listPhi[numPhi]->computeContour(max_dirac_width);
	
	int current_contourLength = listPhi[numPhi]->getContourSize();
	cout << "contour size: " << current_contourLength << endl;
	
	contourLengths[numPhi].push_back(current_contourLength);
	if (contourLengths[numPhi].size() > oscillationsNumber+1) contourLengths[numPhi].pop_front();
	
	checkContourStability(numPhi);
	
	return false;
}

bool ISISD::updatePhi(int numPhi)
{
	compute_dphi_data(numPhi);
	bool test_success = compute_dt_data(numPhi);
	
	if (!test_success) return false;
	
	if (speeds_geom != NULL) compute_dt_geom(numPhi);
	
	/// update phi by adding dphi_data and dphi_geom
	
	int x, y, z, t, ind_phi, indfft;
	float newVal;
	
	dt_data /= (float)nbPixFFT; //speeds_data needs to be divided by nbPixFFT, so we do it here
	
#ifdef _OPENMP
	#if GCC_VERSION >= 40400 && _OPENMP >= 200805
	#pragma omp parallel for default(shared) private(x, y, z, t, indfft, ind_phi, newVal) schedule(dynamic) collapse(4) num_threads(NB_THREADS)
	for (x=0; x<size[0]; x++) //size_fft = size + 2*spatial_frame + padding
	for (y=0; y<size[1]; y++)
	for (z=0; z<size[2]; z++)
	for (t=0; t<size[3]; t++)
	#else
	#pragma omp parallel for default(shared) private(x, y, z, t, indfft, ind_phi, newVal) schedule(dynamic) num_threads(NB_THREADS)
	for (ind_phi=0; ind_phi<nbPix; ind_phi++)
	#endif
#else
	for (x=0; x<size[0]; x++)
	for (y=0; y<size[1]; y++)
	for (z=0; z<size[2]; z++)
	for (t=0; t<size[3]; t++)
#endif
	{
#ifdef _OPENMP
	#if GCC_VERSION >= 40400 && _OPENMP >= 200805
		ind_phi = computeIndex(size, x, y, z, t);
	#else
		computeCoords(ind_phi, size, x, y, z, t);
	#endif
#else
		ind_phi = computeIndex(size, x, y, z, t);
#endif
		indfft = (x + borders[0]) + size_fft[0] * ((y + borders[1]) + size_fft[1] * ((z + borders[2]) + size_fft[2] * (t + borders[3])));
		
		newVal = listPhi[numPhi]->get_phi_value(ind_phi) + speeds_data[indfft] * dt_data;
		if (speeds_geom != NULL) newVal += speeds_geom[ind_phi] * dt_geom;
		
		listPhi[numPhi]->set_phi_value(ind_phi, newVal);
	}
	
	
	listPhi[numPhi]->computeContour(max_dirac_width);
	
	int current_contourLength = listPhi[numPhi]->getContourSize();
	cout << "contour size: " << current_contourLength << " pixels" << endl;
	
	contourLengths[numPhi].push_back(current_contourLength);
	if (contourLengths[numPhi].size() > oscillationsNumber+1) contourLengths[numPhi].pop_front();
	
	checkContourStability(numPhi);
	
	return true;
}

void ISISD::checkContourStability(int numPhi)
{
	if (contourLengths[numPhi].size() < oscillationsNumber+1) //not enough measurements
		return;
	
	/// detect oscillations of the contour and adjust the maxTranslation step accordingly
	
	list<float>::const_iterator lit(contourLengths[numPhi].begin()), lend(contourLengths[numPhi].end()), litPrec;
	litPrec = lit;
	++lit;
	vector<float> diff(oscillationsNumber);
	
	for (int i=0; lit != lend; ++lit, ++litPrec, i++)
	{
		diff[i] = *lit - *litPrec;
	}
	
	bool sameSign = true; //diffs are always of the same sign
	bool diffSign = true; //consecutive diffs are always of different signs
	
	for (int i=1; i<oscillationsNumber; i++)
	{
		if (diff[i] * diff[i-1] > 0)
		{
			diffSign = false;
		}
		else if (diff[i] * diff[i-1] < 0)
		{
			sameSign = false;
		}
		else if (diff[i] == 0 && diff[i-1] == 0)
		{
			diffSign = false;
			sameSign = false;
			break;
		}
	}
	
	if (diffSign) //diffs are always of different signs
	{
		// if variation is very small, don't do anything
		/*for (lit=contourLengths[numPhi].begin(); lit != lend; ++lit)
		{
			if (fabs(contourLengths[numPhi].back() - *lit) > contourLengths[numPhi].back() * 0.05)
			{
				maxTranslation -= maxTranslation * 0.05;
				cout << "oscillations detected -> maxTranslation = " << maxTranslation << endl;
				contourLengths[numPhi].clear();
				
				break;
			}
		}*/
		
		if (dynamic_maxTranslation)
		{
			maxTranslation[numPhi] -= maxTranslation[numPhi] / 10.;
			cout << "oscillations detected -> maxTranslation = " << maxTranslation[numPhi] << endl;
		}
		contourLengths[numPhi].clear();
		
		contourStoppedGrowing[numPhi] = true;
	}
	
	if (sameSign) //consecutive diffs are always of the same sign
	{
		if (dynamic_maxTranslation)
		{
			maxTranslation[numPhi] += 0.1; //= 1
			cout << "very stable contour -> maxTranslation = " << maxTranslation[numPhi] << endl;
		}
		contourLengths[numPhi].clear();
		
		contourStoppedGrowing[numPhi] = false;
	}
	
	if (maxTranslation[numPhi] < 0.1) maxTranslation[numPhi] = 0.1;
	if (maxTranslation[numPhi] > 1) maxTranslation[numPhi] = 1;
}

void ISISD::set_dynamic_maxTranslation(bool newValue)
{
	dynamic_maxTranslation = newValue;
	
	int nbPhi = listPhi.size();
	for (int p=0; p<nbPhi; p++) maxTranslation[p] = 1;
}

bool ISISD::doesContourGrow(int numPhi)
{
	return !contourStoppedGrowing[numPhi];
}

void ISISD::display_S()
{
	cv::Mat imgDisp(size[0], size[1], CV_64FC1);
	
	int x, y;
	int z = size[2]/2;
	int t = size[3]/2;
        double max = 0, min = VAL_INVALID;
        int indfft;
         
        for (x=0; x<size[0]; x++)
        for (y=0; y<size[1]; y++)
	{
		indfft = (x + borders[0]) + size_fft[0] * ((y + borders[1]) + size_fft[1] * ((z + borders[2]) + size_fft[2] * (t + borders[3])));
		imgDisp.at<double>(x, y) = speeds_data[indfft];
		
		if (fabs(imgDisp.at<double>(x, y)) > max) max = fabs(imgDisp.at<double>(x, y));
		if (fabs(imgDisp.at<double>(x, y)) < min) min = fabs(imgDisp.at<double>(x, y));
        }
        imgDisp /= 2.*max;
        imgDisp += 0.5;
         
        ostringstream name;
        name << "S data " << t;
         
        cv::namedWindow(name.str(), CV_WINDOW_AUTOSIZE);
        cv::imshow(name.str(), imgDisp);
         
        if (size[2]!= 1)
        {
            cv::Mat imgDisp2(size[0], size[2], CV_64FC1);
             
            y = size[1]/2;
            max = 0;
            min = VAL_INVALID;
             
            for (x=0; x<size[0]; x++)
            for (z=0; z<size[2]; z++)
            {
                indfft = (x + borders[0]) + size_fft[0] * ((y + borders[1]) + size_fft[1] * ((z + borders[2]) + size_fft[2] * (t + borders[3])));
		imgDisp2.at<double>(x, z) = speeds_data[indfft];
                                 
                if (fabs(imgDisp2.at<double>(x, z)) > max) max = fabs(imgDisp2.at<double>(x, z));
                if (fabs(imgDisp2.at<double>(x, z)) < min) min = fabs(imgDisp2.at<double>(x, z));
            }
            imgDisp2 /= 2.*max;
            imgDisp2 += 0.5;
             
            ostringstream name2;
            name2 << "S data vertical " << t;
             
            cv::namedWindow(name2.str(), CV_WINDOW_AUTOSIZE);
            cv::imshow(name2.str(), imgDisp2);
        }
     
    	cv::waitKey(10);
}


/// ISISD_simpleScheme

ISISD_simpleScheme::ISISD_simpleScheme(Sdata* _S_data, Sgeom* _S_geom, vector<Phi*>& _listPhi, RBF* _psi, float _dirac_width, VolumeData* _volumeData, int _display):
 ISISD(_S_data, _S_geom, _listPhi, _psi->get_borders(), _volumeData, _display)
{
	dirac_width = _dirac_width;
	psi = _psi;
	
	max_dirac_width = dirac_width;
	speeds_data = S_data->get_S();
	
	size_fft = psi->get_size_fft();
	nbPixFFT = size_fft[0] * size_fft[1] * size_fft[2] * size_fft[3];
}

void ISISD_simpleScheme::compute_dphi_data(int numPhi)
{
	if (display == 2) display_S();
	
	/// multiply by dirac
	
	const map<int, float>* contour = listPhi[numPhi]->getContour();
	map<int, float>::const_iterator mit, mend;
	int x, y, z, indfft;
	
	for (int f=0; f<size[3]; f++)
	{
		mit = contour[f].begin();
		mend = contour[f].end();
		
		for (; mit != mend; ++mit)
		{
			computeCoords(mit->first, size, x, y, z);
			indfft = (x + borders[0]) + size_fft[0] * ((y + borders[1]) + size_fft[1] * ((z + borders[2]) + size_fft[2] * (f + borders[3])));
			speeds_data[indfft] *= (1. + cos(M_PI * mit->second / dirac_width)) / (2. * dirac_width);
		}
	}
	
	///compute d_alpha
	
	psi->do_fft(speeds_data, 1);
	psi->multiply_with_RBF(speeds_data, speeds_data);
	
	///compute d_phi
	
	psi->multiply_with_RBF(speeds_data, speeds_data);
	psi->do_fft(speeds_data, -1);
}

void ISISD_simpleScheme::change_psi(RBF* _psi, float _dirac_width)
{
	psi = _psi;
	size_fft = psi->get_size_fft();
	
	dirac_width = _dirac_width;
	max_dirac_width = dirac_width;
}


/// ISISD_mixedScheme

ISISD_mixedScheme::ISISD_mixedScheme(Sdata* _S_data, Sgeom* _S_geom, vector<Phi*>& _listPhi, vector<RBF*>& _listPsi, vector<float>& _dirac_widths, VolumeData* _volumeData, int _display):
 ISISD(_S_data, _S_geom, _listPhi, _listPsi[0]->get_borders(), _volumeData, _display), dirac_widths(_dirac_widths), listPsi(_listPsi)
{
	size_fft = listPsi[0]->get_size_fft();
	nbPixFFT = size_fft[0] * size_fft[1] * size_fft[2] * size_fft[3];
	
	max_dirac_width = *max_element(dirac_widths.begin(), dirac_widths.end());
	
	unique_width = true;
	
	vector<float>::const_iterator vit(dirac_widths.begin()), vend(dirac_widths.end());
	float current_width = *vit;
	++vit;
	
	for (; vit != vend; ++vit)
	{
		if (*vit != current_width)
		{
			unique_width = false;
			break;
		}
	}
	
	dalpha = (double*) fftw_malloc(sizeof(double) * nbPixFFT);
	dphi = new double[nbPixFFT];
	
	/// map of the domains of the different psi
	domains_psi = new unsigned int[nbPix];

	int ind;

	#pragma omp parallel for default(shared) private(ind) schedule(dynamic) num_threads(NB_THREADS)
	for (ind=0; ind<nbPix; ind++) domains_psi[ind] = VAL_INVALID;

	const vector< pair<int, int> >& images = volumeData->get_selectedImages();
	int nbImages = volumeData->getNbImages();
	int x, y, z, t, indVol, indImg;
	float time;
	Point3<float> point;

	for (int image; image < nbImages; image++)
	{
		const pair<int, int>& imgSize = volumeData->getImageSize(images[image].first);
		int nbPixImg = imgSize.first * imgSize.second;

		for (int f=0; f<size[3]; f++)
		{
			time = volumeData->getTimeOfFrame(image, f);
			t = (int)round(time);
			
			if ( t < 0 || t >= size[3] ) continue;

			#pragma omp parallel for default(shared) private(x, y, z, indVol, indImg, point) schedule(dynamic) num_threads(NB_THREADS)
			for (indImg=0; indImg < nbPixImg; indImg++)
			{
				#pragma omp critical
				point = volumeData->getImagePointCoordsInVolume(image, indImg);

				x = (int)round(point.x());
				y = (int)round(point.y());
				z = (int)round(point.z());
				
				if ( x < 0 || x >= size[0] ||
					y < 0 || y >= size[1] ||
					z < 0 || z >= size[2] )
					continue;

				indVol = computeIndex(size, x, y, z, t);
				
				#pragma omp critical
				domains_psi[indVol] = 0;
			}
		}
	}

	int xM, yM, zM, tM, xP, yP, zP, tP, indN, i;
	unsigned int newValue;
	bool test_value;
	
	nbPsi = listPsi.size();

#ifdef _OPENMP
	#if GCC_VERSION >= 40400 && _OPENMP >= 200805
	#pragma omp parallel for default(shared) private(x, y, z, t, ind, xM, yM, zM, tM, i, indN, newValue, test_value) schedule(dynamic) num_threads(NB_THREADS) collapse(4)
	for (x=0; x<size[0]; x++)
	for (y=0; y<size[1]; y++)
	for (z=0; z<size[2]; z++)
	for (t=0; t<size[3]; t++)
	#else
	#pragma omp parallel for default(shared) private(x, y, z, t, ind, xM, yM, zM, tM, i, indN, newValue, test_value) schedule(dynamic) num_threads(NB_THREADS)
	for (ind=0; ind<nbPix; ind++)
	#endif
#else
	for (x=0; x<size[0]; x++)
	for (y=0; y<size[1]; y++)
	for (z=0; z<size[2]; z++)
	for (t=0; t<size[3]; t++)
#endif
	{
#ifdef _OPENMP
	#if GCC_VERSION >= 40400 && _OPENMP >= 200805
		ind = computeIndex(size, x, y, z, t);
	#else
		computeCoords(ind, size, x, y, z, t);
	#endif
#else
		ind = computeIndex(size, x, y, z, t);
#endif
		#pragma omp critical
		test_value = !domains_psi[ind];
		
		if (test_value) continue;

		xM = x-1 >= 0 ? x-1 : 0;
		yM = y-1 >= 0 ? y-1 : 0;
		zM = z-1 >= 0 ? z-1 : 0;
		tM = t-1 >= 0 ? t-1 : 0;

		int xN2[4] = {xM, x, x, x};
		int yN2[4] = {y, yM, y, y};
		int zN2[4] = {z, z, zM, z};
		int tN2[4] = {t, t, t, tM};

		for (i=0; i<4; i++)
		{
			indN = computeIndex(size, xN2[i], yN2[i], zN2[i], tN2[i]);

			newValue = domains_psi[indN] + 1;

			if (newValue >= nbPsi) newValue = nbPsi - 1;

			#pragma omp critical
			test_value = domains_psi[ind] > newValue;

			if (test_value)
			{
				#pragma omp critical
				domains_psi[ind] = newValue;

			}
		}
	}

#ifdef _OPENMP
	#if GCC_VERSION >= 40400 && _OPENMP >= 200805
	#pragma omp parallel for default(shared) private(x, y, z, t, ind, xP, yP, zP, tP, i, indN, newValue, test_value) schedule(dynamic) num_threads(NB_THREADS) collapse(4)
	for (x=size[0]-1; x>=0; x--)
	for (y=size[1]-1; y>=0; y--)
	for (z=size[2]-1; z>=0; z--)
	for (t=size[3]-1; t>=0; t--)
	#else
	#pragma omp parallel for default(shared) private(x, y, z, t, ind, xP, yP, zP, tP, i, indN, newValue, test_value) schedule(dynamic) num_threads(NB_THREADS)
	for (ind=nbPix-1; ind>=0; ind--)
	#endif
#else
	for (x=size[0]-1; x>=0; x--)
	for (y=size[1]-1; y>=0; y--)
	for (z=size[2]-1; z>=0; z--)
	for (t=size[3]-1; t>=0; t--)
#endif
	{
#ifdef _OPENMP
	#if GCC_VERSION >= 40400 && _OPENMP >= 200805
		ind = computeIndex(size, x, y, z, t);
	#else
		computeCoords(ind, size, x, y, z, t);
	#endif
#else
		ind = computeIndex(size, x, y, z, t);
#endif
		#pragma omp critical
		test_value = !domains_psi[ind];
		
		if (test_value) continue;

		xP = x+1 < size[0] ? x+1 : size[0]-1;
		yP = y+1 < size[1] ? y+1 : size[1]-1;
		zP = z+1 < size[2] ? z+1 : size[2]-1;
		tP = t+1 < size[3] ? t+1 : size[3]-1;

		int xN2[4] = {xP, x, x, x};
		int yN2[4] = {y, yP, y, y};
		int zN2[4] = {z, z, zP, z};
		int tN2[4] = {t, t, t, tP};

		for (i=0; i<4; i++)
		{
			indN = computeIndex(size, xN2[i], yN2[i], zN2[i], tN2[i]);

			newValue = domains_psi[indN] + 1;

			if (newValue >= nbPsi) newValue = nbPsi - 1;

			#pragma omp critical
			test_value = domains_psi[ind] > newValue;
			
			if (test_value)
			{
				#pragma omp critical
				domains_psi[ind] = newValue;
			}
		}
	}
	
	if (display == 2) display_psi_domains();
}

ISISD_mixedScheme::~ISISD_mixedScheme()
{
	fftw_free(dalpha);
	delete[] dphi;
	delete[] domains_psi;
}

void ISISD_mixedScheme::compute_dphi_data(int numPhi)
{
	speeds_data = S_data->get_S();
	
	memory_S.clear();
	
	/// 1) multiply speeds_data by dirac, and keep speeds_data in memory if needed (ie if multiple dirac widths)
	
	float current_gamma = listPsi[0]->get_gamma();
	float current_epsilon = dirac_widths[current_gamma];
	
	const map<int, float>* contour = listPhi[numPhi]->getContour();
	map<int, float>::const_iterator mit, mend;
	
	int x, y, z, t, indfft;
	
	for (t=0; t<size[3]; t++)
	{
		mit = contour[t].begin();
		mend = contour[t].end();
		
		for (; mit != mend; ++mit)
		{
			computeCoords(mit->first, size, x, y, z);
			
			indfft = (x + borders[0]) + size_fft[0] * ((y + borders[1]) + size_fft[1] * ((z + borders[2]) + size_fft[2] * (t + borders[3])));
			
			if (!unique_width) memory_S[indfft] = make_pair(speeds_data[indfft], mit->second);
			
			if (mit->second > current_epsilon) speeds_data[indfft] = 0;
			else speeds_data[indfft] *= (1. + cos(M_PI * mit->second / current_epsilon)) / (2. * current_epsilon);
		}
	}
	
	if (display == 2) display_S();
	
	//go to frequency domain
	listPsi[0]->do_fft(speeds_data, 1);
	
	/// 2) compute dphi for each gamma
	
	int ind;
	map<int, pair<float, float> >::const_iterator mit_mem, mend_mem(memory_S.end());
	
	int xP, yP, zP, tP;
	int lim_x = size[0] + 2 * borders[0];
	int dx = (int)ceil( (float)lim_x/2. );
	int dy = (int)ceil( (float)size_fft[1]/2. );
	int dz = (int)ceil( (float)size_fft[2]/2. );
	int dt = (int)ceil( (float)size_fft[3]/2. );
	
	#pragma omp parallel for default(shared) private(ind) schedule(dynamic) num_threads(NB_THREADS)
	for (ind=0; ind<nbPixFFT; ind++)
		dphi[ind] = 0;
	
	for (int num_psi=0; num_psi < nbPsi; num_psi++)
	{
		if ( current_epsilon != dirac_widths[num_psi] )
		{
			current_epsilon = dirac_widths[num_psi];
			
			#pragma omp parallel for default(shared) private(ind) schedule(dynamic) num_threads(NB_THREADS)
			for (ind=0; ind<nbPixFFT; ind++)
				speeds_data[ind] = 0;
			
			for (mit_mem=memory_S.begin(); mit_mem != mend_mem; ++mit_mem)
			{
				if (mit_mem->second.second <= current_epsilon)
					speeds_data[mit_mem->first] = mit_mem->second.first * (1. + cos(M_PI * mit_mem->second.second / current_epsilon)) / (2. * current_epsilon);
			}
			
			display_S();
			
			//go to frequency domain
			listPsi[num_psi]->do_fft(speeds_data, 1);
		}
		
		// 1rst step: compute d_alpha_n/dt then go back to spatial domain
		
		listPsi[num_psi]->multiply_with_RBF(speeds_data, dalpha);
		listPsi[num_psi]->do_fft(dalpha, -1);
		
		// 2nd step: apply mask to d_alpha_n/dt
		
#ifdef _OPENMP
	#if GCC_VERSION >= 40400 && _OPENMP >= 200805
		#pragma omp parallel for default(shared) private(x, y, z, t, xP, yP, zP, tP, indfft, ind) schedule(dynamic) collapse(4) num_threads(NB_THREADS)
		for (x=0; x<size_fft[0]; x++) //size_fft = size + 2*spatial_frame + padding
		for (y=0; y<size_fft[1]; y++)
		for (z=0; z<size_fft[2]; z++)
		for (t=0; t<size_fft[3]; t++)
	#else
		#pragma omp parallel for default(shared) private(x, y, z, t, xP, yP, zP, tP, indfft, ind) schedule(dynamic) num_threads(NB_THREADS)
		for (indfft=0; indfft<nbPixFFT; indfft++)
	#endif
#else
		for (x=0; x<size_fft[0]; x++) //size_fft = size + 2*spatial_frame + padding
		for (y=0; y<size_fft[1]; y++)
		for (z=0; z<size_fft[2]; z++)
		for (t=0; t<size_fft[3]; t++)
#endif
		{
#ifdef _OPENMP
	#if GCC_VERSION >= 40400 && _OPENMP >= 200805
			indfft = computeIndex(size_fft, x, y, z, t);
	#else
			computeCoords(indfft, size, x, y, z, t);
	#endif
#else
			indfft = computeIndex(size_fft, x, y, z, t);
#endif
			
			if (x >= lim_x) //we are in the padding area
			{
				dalpha[indfft] = 0;
				continue;
			}
			
			xP = x - dx;
			yP = y - dy;
			zP = z - dz;
			tP = t - dt;
			
			if (xP < 0) xP += lim_x;
			if (yP < 0) yP += size_fft[1];
			if (zP < 0) zP += size_fft[2];
			if (tP < 0) tP += size_fft[3];
			
			xP -= borders[0];
			yP -= borders[1];
			zP -= borders[2];
			tP -= borders[3];
			
			if ( xP < 0 || xP >= size[0] ||
				yP < 0 || yP >= size[1] ||
				zP < 0 || zP >= size[2] ||
				tP < 0 || tP >= size[3] )
			{
				dalpha[indfft] = 0;
				continue;
			}
			else
			{
				ind = computeIndex(size, xP, yP, zP, tP);
				if (domains_psi[ind] != num_psi)
				{
					dalpha[indfft] = 0;
					continue;
				}
				else dalpha[indfft] /= (float)nbPixFFT;
			}

			if (isnan(dalpha[indfft]))
			{
				cout << "ISISD_mixedScheme::compute_dphi_data() " << dalpha[indfft] << " error #-1" << endl;
				int reponse;
				cin >> reponse;
			}
		}
		
		// 3rd step: convolve with psi_n
		
		listPsi[num_psi]->do_fft(dalpha, 1);
		listPsi[num_psi]->multiply_with_RBF(dalpha, dalpha);
		listPsi[num_psi]->do_fft(dalpha, -1);
		
		// 4th step: add to dphi/dt
		
		// normalisation of dphi/dt|n so that all levels n have an equal weight in the final sum
		double maxVal = 0;
		
		#pragma omp parallel for default(shared) private(ind) schedule(dynamic) num_threads(NB_THREADS)
		for (ind=0; ind<nbPixFFT; ind++)
		{
			#pragma omp critical
			if (fabs(dalpha[ind]) > maxVal) maxVal = fabs(dalpha[ind]);
			
			if (isnan(dalpha[ind]))
			{
				cout << "ISISD_mixedScheme::compute_dphi_data() " << dalpha[ind] << " error #-2" << endl;
				int reponse;
				cin >> reponse;
			}
		}

		if (maxVal == 0) maxVal = 1;
		
		#pragma omp parallel for default(shared) private(ind) schedule(dynamic) num_threads(NB_THREADS)
		for (ind=0; ind<nbPixFFT; ind++)
		{
			dphi[ind] += dalpha[ind] / maxVal;
		}

		if (display == 2) display_dphi_tmp();
	} // each RBF
	
	// just so that dphi is in the same array for all ISISD classes
	speeds_data = dphi;
}

void ISISD_mixedScheme::display_psi_domains()
{
	cv::Mat imgDisp(size[0], size[1], CV_64FC1);

	int x, y;
	int z = size[2]/2;
	int t = size[3]/2;
	int ind;

	for (x=0; x<size[0]; x++)
	for (y=0; y<size[1]; y++)
	{
		ind = computeIndex(size, x, y, z, t);
		imgDisp.at<double>(x, y) = domains_psi[ind];
	}
	imgDisp /= (float)(nbPsi - 1);

	cv::namedWindow("psi domains - central xy plane", CV_WINDOW_AUTOSIZE);
	cv::imshow("psi domains - central xy plane", imgDisp);

	if (size[2]!= 1)
	{
		cv::Mat imgDisp2(size[0], size[2], CV_64FC1);

		y = size[1]/2;

		for (x=0; x<size[0]; x++)
		for (z=0; z<size[2]; z++)
		{
			ind = computeIndex(size, x, y, z, t);

			imgDisp2.at<double>(x, z) = domains_psi[ind];
		}
		imgDisp2 /= (float)(nbPsi - 1);

		cv::namedWindow("psi domains - central xz plane", CV_WINDOW_AUTOSIZE);
		cv::imshow("psi domains - central xz plane", imgDisp2);
	}

	cv::waitKey(10);
}

void ISISD_mixedScheme::display_dphi_tmp()
{
	cv::Mat imgDisp(size[0], size[1], CV_64FC1);
	
	int x, y;
	int z = size[2]/2;
	int t = size[3]/2;
        double max = 0, min = VAL_INVALID;
        int indfft;
         
        for (x=0; x<size[0]; x++)
        for (y=0; y<size[1]; y++)
	{
		indfft = (x + borders[0]) + size_fft[0] * ((y + borders[1]) + size_fft[1] * ((z + borders[2]) + size_fft[2] * (t + borders[3])));
		imgDisp.at<double>(x, y) = dphi[indfft];
		
		if (fabs(imgDisp.at<double>(x, y)) > max) max = fabs(imgDisp.at<double>(x, y));
		if (fabs(imgDisp.at<double>(x, y)) < min) min = fabs(imgDisp.at<double>(x, y));
        }
        imgDisp /= 2.*max;
        imgDisp += 0.5;
         
        ostringstream name;
        name << "dphi data tmp " << t;
         
        cv::namedWindow(name.str(), CV_WINDOW_AUTOSIZE);
        cv::imshow(name.str(), imgDisp);
         
        if (size[2]!= 1)
        {
            cv::Mat imgDisp2(size[0], size[2], CV_64FC1);
             
            y = size[1]/2;
            max = 0;
            min = VAL_INVALID;
             
            for (x=0; x<size[0]; x++)
            for (z=0; z<size[2]; z++)
            {
                indfft = (x + borders[0]) + size_fft[0] * ((y + borders[1]) + size_fft[1] * ((z + borders[2]) + size_fft[2] * (t + borders[3])));
		imgDisp2.at<double>(x, z) = dphi[indfft];
                                 
                if (fabs(imgDisp2.at<double>(x, z)) > max) max = fabs(imgDisp2.at<double>(x, z));
                if (fabs(imgDisp2.at<double>(x, z)) < min) min = fabs(imgDisp2.at<double>(x, z));
            }
            imgDisp2 /= 2.*max;
            imgDisp2 += 0.5;
             
            ostringstream name2;
            name2 << "dphi data tmp vertical " << t;
             
            cv::namedWindow(name2.str(), CV_WINDOW_AUTOSIZE);
            cv::imshow(name2.str(), imgDisp2);
        }
     
    	cv::waitKey(10);
}
