#include "SpeedFunction.h"

SpeedFunction::SpeedFunction(vector<Phi*>& _listPhi): listPhi(_listPhi)
{}

SpeedFunction::~SpeedFunction()
{}


/// SpeedFunction_Data ///

SpeedFunction_Data::SpeedFunction_Data(VolumeData* _volumeData, vector<Phi*>& _listPhi): SpeedFunction(_listPhi)
{
	volumeData = _volumeData;
}

void SpeedFunction_Data::add_serie(int serie)
{
	list_series.insert(serie);
	
	int nb_frames = volumeData->getNbFrames();
	int nbSlices = volumeData->getNbSlices(serie);
	
	for (int s=0; s<nbSlices; s++)
		image_can_be_used[ volumeData->getImageNumber(serie, s) ].resize(nb_frames, false);
}

bool SpeedFunction_Data::is_image_usable(int image, int timeframe)
{
	return image_can_be_used[image][timeframe];
}


/// SpeedFunction_CV ///

SpeedFunction_CV::SpeedFunction_CV(VolumeData* _volumeData, vector<Phi*>& _listPhi): SpeedFunction_Data(_volumeData, _listPhi)
{
	nb_phases = pow(2, listPhi.size());
	
	if (nb_phases > 4) throw string("SpeedFunction_CV() Error: number of phases > 4 not supported");
	
	c.resize(nb_phases);
	count.resize(nb_phases);
	
	if ( volumeData->comprise_temporal_sequences() ) slice_wise = true;
	else slice_wise = false;
}

void SpeedFunction_CV::update_c(int image, int timeframe, int sequence)
{
	const pair<int, int>& imgSize = volumeData->getImageSize( volumeData->get_selectedImages()[image].first );
	int nbPixImg = imgSize.first * imgSize.second;
	float time = volumeData->getTimeOfFrame(image, timeframe);

	int indImg;
	float phiVal0, phiVal1;
	Point3<float> point;
	float pixVal;

	#pragma omp parallel for default(shared) private(indImg, phiVal0, phiVal1, point, pixVal) schedule(dynamic) num_threads(NB_THREADS)
	for (indImg = 0; indImg < nbPixImg; indImg++)
	{
		pixVal = volumeData->getPixelValue(image, timeframe, indImg);
		if (pixVal == VAL_INVALID) continue;
		
		point = volumeData->getImagePointCoordsInVolume(image, indImg);
		if ( !volumeData->isPointInVolume(point) ) continue;
		
		phiVal0 = listPhi[0]->get_phi_value(point.x(), point.y(), point.z(), time);
		if (nb_phases > 2) phiVal1 = listPhi[1]->get_phi_value(point.x(), point.y(), point.z(), time);
		
		switch(nb_phases)
		{
			case 2:
				if (phiVal0 >= 0)
				{
					#pragma omp critical
					c[0][sequence] += pixVal;
					#pragma omp critical
					count[0][sequence] += 1;
				}
				else if (phiVal0 < 0)
				{
					#pragma omp critical
					c[1][sequence] += pixVal;
					#pragma omp critical
					count[1][sequence] += 1;
				}
				
				break;
				
			case 4:
				if (phiVal0 > 0 && phiVal1 > 0)
				{
					#pragma omp critical
					c[0][sequence] += pixVal;
					#pragma omp critical
					count[0][sequence] += 1;
				}
				else if (phiVal0 > 0 && phiVal1 < 0)
				{
					#pragma omp critical
					c[1][sequence] += pixVal;
					#pragma omp critical
					count[1][sequence] += 1;
				}
				else if (phiVal0 < 0 && phiVal1 > 0)
				{
					#pragma omp critical
					c[2][sequence] += pixVal;
					#pragma omp critical
					count[2][sequence] += 1;
				}
				else if (phiVal0 < 0 && phiVal1 < 0)
				{
					#pragma omp critical
					c[3][sequence] += pixVal;
					#pragma omp critical
					count[3][sequence] += 1;
				}
				
				break;

			default:
				throw( string("SpeedFunction_CV::update_c() Error: nb of phases not supported") );
				break;
		}
	}
}

int SpeedFunction_CV::compute_sequence_number(int serie, int slice)
{
	if (slice_wise) return volumeData->getImageNumber(serie, slice);
	
	return serie;
}

void SpeedFunction_CV::prepare_iteration()
{
	for (int zone=0; zone < nb_phases; zone++)
	{
		c[zone].clear();
		count[zone].clear();
	}
	
	int nb_frames = volumeData->getNbFrames();
	set<int>::const_iterator sit(list_series.begin()), send(list_series.end());
	
	for (; sit != send; ++sit)
	{
		int nb_slices = volumeData->getNbSlices(*sit);
		
		for (int s=0; s<nb_slices; s++)
		for (int f=0; f<nb_frames; f++)
			update_c(volumeData->getImageNumber(*sit, s), f, compute_sequence_number(*sit, s));
	}
	
	map<int, int> image_usable;
	
	for (int zone=0; zone < nb_phases; zone++)
	{
		map<int, float>::iterator mit_c(c[zone].begin());
		map<int, int>::const_iterator mit(count[zone].begin()), mend(count[zone].end());
		
		for (; mit != mend; ++mit, ++mit_c)
		{
			if ( !mit->second ) mit_c->second = VAL_INVALID;
			else
			{
				mit_c->second /= (float)mit->second;
				
				cout << "phase " << zone << " sequence " << mit_c->first << ": c = " << mit_c->second << endl;
				
				if (slice_wise) image_usable[mit->first] += 1;
				else
				{
					int nb_slices = volumeData->getNbSlices(mit->first);
	
					for (int s=0; s<nb_slices; s++)
					{
						image_usable[ volumeData->getImageNumber(mit->first, s) ] += 1;
					}
				}
			}
		}
	}
	
	map<int, vector<bool> >::iterator mit_u(image_can_be_used.begin()), mend_u(image_can_be_used.end());
	
	for (; mit_u != mend_u; ++mit_u)
	{
		if ( image_usable[mit_u->first] > 1 )
			for (int f=0; f<nb_frames; f++) mit_u->second[f] = true;
		else
			for (int f=0; f<nb_frames; f++) mit_u->second[f] = false;
	}
}

float SpeedFunction_CV::get_speed(int image, int timeframe, int pixel, int numPhi)
{
		if ( !is_image_usable(image, timeframe) ) return VAL_INVALID;
		
		int sequence;
		if (slice_wise) sequence = image;
		else sequence = volumeData->get_selectedImages()[image].first + list_series.size()*timeframe; //serie+timeframe*nbSeries
		
	float pixVal = volumeData->getPixelValue(image, timeframe, pixel);
	if (pixVal == VAL_INVALID) return VAL_INVALID;
	
	float time = volumeData->getTimeOfFrame(image, timeframe);

	float phiVal;
	
	if (nb_phases == 4)
	{
		Point3<float> point = volumeData->getImagePointCoordsInVolume(image, pixel);

		if (numPhi == 0) phiVal = listPhi[1]->get_phi_value(point.x(), point.y(), point.z(), time);
		else phiVal = listPhi[0]->get_phi_value(point.x(), point.y(), point.z(), time);
	}

	float diff0, diff1, diff2, diff3, res;

	diff0 = pixVal - c[0][sequence];
	diff0 *= diff0;

	diff1 = pixVal - c[1][sequence];
	diff1 *= diff1;

	switch(nb_phases)
	{
		case 2:
			if (c[1][sequence] == VAL_INVALID || c[0][sequence] == VAL_INVALID) res = VAL_INVALID;
			else res = diff1 - diff0;
		
			break;
		
		case 4:
			diff2 = pixVal - c[2][sequence];
			diff2 *= diff2;
		
			diff3 = pixVal - c[3][sequence];
			diff3 *= diff3;
		
			if (numPhi == 0)
			{
				if (phiVal >= 0)
				{
					if (c[2][sequence] == VAL_INVALID || c[0][sequence] == VAL_INVALID)
						res = VAL_INVALID;
					else res = diff2 - diff0;
				}
				else
				{
					if (c[3][sequence] == VAL_INVALID || c[1][sequence] == VAL_INVALID)
						res = VAL_INVALID;
					else res = diff3 - diff1;
				}
			}
			else
			{
				if (phiVal >= 0)
				{
					if (c[1][sequence] == VAL_INVALID || c[0][sequence] == VAL_INVALID)
						res = VAL_INVALID;
					else res = diff1 - diff0;
				}
				else
				{
					if (c[3][sequence] == VAL_INVALID || c[2][sequence] == VAL_INVALID)
						res = VAL_INVALID;
					else res = diff3 - diff2;
				}
			}
		
			break;
			
		default:
			throw( string("SpeedFunction_CV::get_speed() Error: number of phases not supported") );
			break;
	}
	
	return res;
}


/// SpeedFunction_GVF ///

SpeedFunction_GVF::SpeedFunction_GVF(VolumeData* _volumeData, vector<Phi*>& _listPhi, float _baloonForce): SpeedFunction_Data(_volumeData, _listPhi)
{
	cout << "constructor SpeedFunction_GVF" << endl;
	
	baloonForce = _baloonForce;
}

void SpeedFunction_GVF::add_serie(int serie)
{
	SpeedFunction_Data::add_serie(serie);
	
	/// compute g and grad_g
	
	int rKernel, cKernel;
	int infr, supr, infc, supc;
	float tmp, tmpc, tmpr;
	float normGrad2;
	
	float kernel[9] = {1./16., 2./16., 1./16.,
			2./16., 4./16., 2./16.,
			1./16., 2./16., 1./16.};
	
	int nbFrames = volumeData->getNbFrames();
	
	const pair<int, int>& imgSize = volumeData->getImageSize(serie);
	int Rows = imgSize.first;
	int Columns = imgSize.second;
	int nbPixImg = Rows * Columns;
	int nbSlices = volumeData->getNbSlices(serie);
	
	float* smoothedImage = new float[nbPixImg];
	
        for (int s=0; s<nbSlices; s++)
	{
		int image = volumeData->getImageNumber(serie, s);
		
		gValues[image].resize(nbFrames);
		gradgValues[image].resize(nbFrames);
		
		for (int f=0; f<nbFrames; f++)
		{
			gValues[image][f].resize(nbPixImg, 0);
			gradgValues[image][f].resize(nbPixImg, VECT_NULL);
		
			// smooth image
			for (int r=0; r<Rows; r++)
			for (int c=0; c<Columns; c++)
			{
				tmp = 0;
		
				for (int i=0; i<3; i++)
				for (int j=0; j<3; j++)
				{
					rKernel = r + (i - 1) >= 0 ? r + (i - 1) : 0;
					rKernel = rKernel < Rows ? rKernel : Rows - 1;
			
					cKernel = c + (j - 1) >= 0 ? c + (j - 1) : 0;
					cKernel = cKernel < Columns ? cKernel : Columns - 1;
			
					if (volumeData->getPixelValue(image, f, cKernel + Columns * rKernel) == -VAL_INVALID) continue;
						
					tmp += kernel[i + 3 * j] * volumeData->getPixelValue(image, f, cKernel + Columns * rKernel);
				}
				
				smoothedImage[c + Columns * r] = tmp;
			} //smooth image
			
			cv::Mat imgDisp_g(Rows, Columns, CV_64FC1);
			
			// compute gradient of smoothed image, then g
			for (int r=0; r<Rows; r++)
			for (int c=0; c<Columns; c++)
			{
				if (volumeData->getPixelValue(image, f, c + Columns * r) == -VAL_INVALID)
					gValues[image][f][c + Columns * r] = 1;
				else
				{
					infr = r - 1 >= 0 ? r - 1 : 0;
					supr = r + 1 < Rows ? r + 1 : Rows - 1;
				
					infc = c - 1 >= 0 ? c - 1 : 0;
					supc = c + 1 < Columns ? c + 1 : Columns - 1;
				
					tmpc = (smoothedImage[supc + Columns * r] - smoothedImage[infc + Columns * r]) / 2.;
					tmpr = (smoothedImage[c + Columns * supr] - smoothedImage[c + Columns * infr]) / 2.;
					
					normGrad2 = tmpc * tmpc + tmpr * tmpr;
					gValues[image][f][c + Columns * r] = 1. / (1. + normGrad2);
				}
				
				imgDisp_g.at<double>(r, c) = gValues[image][f][c + Columns * r];
			} //g
			
			cv::namedWindow("g", CV_WINDOW_AUTOSIZE);
			cv::imshow("g", imgDisp_g);
			cv::waitKey(10);
			
			// compute grad_g
			for (int r=0; r<Rows; r++)
			for (int c=0; c<Columns; c++)
			{
				infr = r - 1 >= 0 ? r - 1 : 0;
				supr = r + 1 < Rows ? r + 1 : Rows - 1;
			
				infc = c - 1 >= 0 ? c - 1 : 0;
				supc = c + 1 < Columns ? c + 1 : Columns - 1;
			
				tmpc = (gValues[image][f][supc + Columns * r] - gValues[image][f][infc + Columns * r]) / 2.;
				tmpr = (gValues[image][f][c + Columns * supr] - gValues[image][f][c + Columns * infr]) / 2.;
				
				gradgValues[image][f][c + Columns * r] = Vector3(tmpr, tmpc, 0);
				
				imgDisp_g.at<double>(r, c) = gradgValues[image][f][c + Columns * r].norm();
			} //grad_g
			
			cv::namedWindow("grad_g", CV_WINDOW_AUTOSIZE);
			cv::imshow("grad_g", imgDisp_g);
			cv::waitKey(10);
			
			image_can_be_used[image][f] = true;
		} //all frames
	} //all slices
		
	delete[] smoothedImage;
}

void SpeedFunction_GVF::prepare_iteration()
{}

float SpeedFunction_GVF::get_speed(int image, int timeframe, int pixel, int numPhi)
{
	Vector3 gradg = volumeData->imageToVolumeVect(image, gradgValues[image][timeframe][pixel]);
	
	Point3<float> point = volumeData->getImagePointCoordsInVolume(image, pixel);
	Vector3 n = listPhi[numPhi]->getSpatialNormal(point.x(), point.y(), point.z(), volumeData->getTimeOfFrame(image, timeframe));
	
	return gValues[image][timeframe][pixel] * baloonForce - gradg * n;
}


/// SpeedFunction_Geom ///

SpeedFunction_Geom::SpeedFunction_Geom(vector<Phi*>& _listPhi, const int* _size, bool _periodic_motion): SpeedFunction(_listPhi)
{
	size = _size;
	nbPixVol = size[0] * size[1] * size[2];
	pitch.resize(3);
	pitch[0] = 1;
	pitch[1] = size[0];
	pitch[2] = size[0] * size[1];
	periodic_motion = _periodic_motion;
}

double SpeedFunction_Geom::computeNormGradientEntropy2ndOrder(int voxel, int time, float speed, int numPhi)
{
	Vector4 dphiP = computeForwardDiff(voxel, time, numPhi);
	Vector4 dphiM = computeBackwardDiff(voxel, time, numPhi);
	
	Vector4 dphiPM = computeForwardBackwardDiff(voxel, time, numPhi);
	Vector4 dphiPP = computeForwardForwardDiff(voxel, time, numPhi);
	Vector4 dphiMM = computeBackwardBackwardDiff(voxel, time, numPhi);
	
	double A = dphiM[0] + m(dphiMM[0], dphiPM[0]) / 2.;
	double B = dphiP[0] - m(dphiPP[0], dphiPM[0]) / 2.;
	double C = dphiM[1] + m(dphiMM[1], dphiPM[1]) / 2.;
	double D = dphiP[1] - m(dphiPP[1], dphiPM[1]) / 2.;
	double E = dphiM[2] + m(dphiMM[2], dphiPM[2]) / 2.;
	double F = dphiP[2] - m(dphiPP[2], dphiPM[2]) / 2.;
	double G = dphiM[3] + m(dphiMM[3], dphiPM[3]) / 2.;
	double H = dphiP[3] - m(dphiPP[3], dphiPM[3]) / 2.;

	Vector4 dPhi2 = VECT_NULL_4;
	if (speed < 0)
	{
		dPhi2[0] = pow( max(A, 0.), 2 ) + pow( min(B, 0.), 2 );
		dPhi2[1] = pow( max(C, 0.), 2 ) + pow( min(D, 0.), 2 );
		dPhi2[2] = pow( max(E, 0.), 2 ) + pow( min(F, 0.), 2 );
		dPhi2[3] = pow( max(G, 0.), 2 ) + pow( min(H, 0.), 2 );
	}
	else if (speed > 0)
	{
		dPhi2[0] = pow( max(B, 0.), 2 ) + pow( min(A, 0.), 2 );
		dPhi2[1] = pow( max(D, 0.), 2 ) + pow( min(C, 0.), 2 );
		dPhi2[2] = pow( max(F, 0.), 2 ) + pow( min(E, 0.), 2 );
		dPhi2[3] = pow( max(H, 0.), 2 ) + pow( min(G, 0.), 2 );
	}
	
	return sqrt(dPhi2[0] + dPhi2[1] + dPhi2[2] + dPhi2[3]);
}

double SpeedFunction_Geom::m(double x, double y)
{
	if (x * y < 0) return 0;
	
	if (abs(x) <= abs(y)) return x;
	return y;
}

Vector4 SpeedFunction_Geom::computeForwardDiff(int voxel, int time, int numPhi)
{
	int voxel_xP = voxel + pitch[0] < nbPixVol ? voxel + pitch[0] : nbPixVol - 1;
	int voxel_yP = voxel + pitch[1] < nbPixVol ? voxel + pitch[1] : nbPixVol - 1;
	int voxel_zP = voxel + pitch[2] < nbPixVol ? voxel + pitch[2] : nbPixVol - 1;
	int tP;
	if (periodic_motion) tP = time + 1 < size[3] ? time + 1 : time + 1 - size[3];
	else tP = time + 1 < size[3] ? time + 1 : size[3] - 1;
	
	double phiVal = listPhi[numPhi]->get_phi_value(voxel + time * nbPixVol);
	
	double dphixP = listPhi[numPhi]->get_phi_value(voxel_xP + time * nbPixVol) - phiVal;
	double dphiyP = listPhi[numPhi]->get_phi_value(voxel_yP + time * nbPixVol) - phiVal;
	double dphizP = listPhi[numPhi]->get_phi_value(voxel_zP + time * nbPixVol) - phiVal;
	double dphitP = listPhi[numPhi]->get_phi_value(voxel + tP * nbPixVol) - phiVal;
	
	return Vector4(dphixP, dphiyP, dphizP, dphitP);
}

Vector4 SpeedFunction_Geom::computeBackwardDiff(int voxel, int time, int numPhi)
{
	int voxel_xM = voxel - pitch[0] >= 0 ? voxel - pitch[0] : 0;
	int voxel_yM = voxel - pitch[1] >= 0 ? voxel - pitch[1] : 0;
	int voxel_zM = voxel - pitch[2] >= 0 ? voxel - pitch[2] : 0;
	int tM;
	if (periodic_motion) tM = time - 1 >= 0 ? time - 1 : time - 1 + size[3];
	else tM = time - 1 >= 0 ? time -1 : 0;
	
	double phiVal = listPhi[numPhi]->get_phi_value(voxel + time * nbPixVol);
	
	double dphixM = phiVal - listPhi[numPhi]->get_phi_value(voxel_xM + time * nbPixVol);
	double dphiyM = phiVal - listPhi[numPhi]->get_phi_value(voxel_yM + time * nbPixVol);
	double dphizM = phiVal - listPhi[numPhi]->get_phi_value(voxel_zM + time * nbPixVol);
	double dphitM = phiVal - listPhi[numPhi]->get_phi_value(voxel + tM * nbPixVol);
	
	return Vector4(dphixM, dphiyM, dphizM, dphitM);
}

Vector4 SpeedFunction_Geom::computeForwardBackwardDiff(int voxel, int time, int numPhi)
{
	int voxel_xP = voxel + pitch[0] < nbPixVol ? voxel + pitch[0] : nbPixVol - 1;
	int voxel_yP = voxel + pitch[1] < nbPixVol ? voxel + pitch[1] : nbPixVol - 1;
	int voxel_zP = voxel + pitch[2] < nbPixVol ? voxel + pitch[2] : nbPixVol - 1;
	int tP;
	if (periodic_motion) tP = time + 1 < size[3] ? time + 1 : time + 1 - size[3];
	else tP = time + 1 < size[3] ? time + 1 : size[3] - 1;
	
	int voxel_xM = voxel - pitch[0] >= 0 ? voxel - pitch[0] : 0;
	int voxel_yM = voxel - pitch[1] >= 0 ? voxel - pitch[1] : 0;
	int voxel_zM = voxel - pitch[2] >= 0 ? voxel - pitch[2] : 0;
	int tM;
	if (periodic_motion) tM = time - 1 >= 0 ? time - 1 : time - 1 + size[3];
	else tM = time - 1 >= 0 ? time -1 : 0;
	
	double phiVal = listPhi[numPhi]->get_phi_value(voxel + time * nbPixVol);
	
	double dphixPM = listPhi[numPhi]->get_phi_value(voxel_xP + time * nbPixVol) - 2. * phiVal + listPhi[numPhi]->get_phi_value(voxel_xM + time * nbPixVol);
	double dphiyPM = listPhi[numPhi]->get_phi_value(voxel_yP + time * nbPixVol) - 2. * phiVal + listPhi[numPhi]->get_phi_value(voxel_yM + time * nbPixVol);
	double dphizPM = listPhi[numPhi]->get_phi_value(voxel_zP + time * nbPixVol) - 2. * phiVal + listPhi[numPhi]->get_phi_value(voxel_zM + time * nbPixVol);
	double dphitPM = listPhi[numPhi]->get_phi_value(voxel + tP * nbPixVol) - 2. * phiVal + listPhi[numPhi]->get_phi_value(voxel + tM * nbPixVol);
	
	return Vector4(dphixPM, dphiyPM, dphizPM, dphitPM);
}

Vector4 SpeedFunction_Geom::computeForwardForwardDiff(int voxel, int time, int numPhi)
{
	int voxel_xP = voxel + pitch[0] < nbPixVol ? voxel + pitch[0] : nbPixVol - 1;
	int voxel_yP = voxel + pitch[1] < nbPixVol ? voxel + pitch[1] : nbPixVol - 1;
	int voxel_zP = voxel + pitch[2] < nbPixVol ? voxel + pitch[2] : nbPixVol - 1;
	int tP;
	if (periodic_motion) tP = time + 1 < size[3] ? time + 1 : time + 1 - size[3];
	else tP = time + 1 < size[3] ? time + 1 : size[3] - 1;
	
	int voxel_xPP = voxel_xP + pitch[0] < nbPixVol ? voxel_xP + pitch[0] : nbPixVol - 1;
	int voxel_yPP = voxel_yP + pitch[1] < nbPixVol ? voxel_yP + pitch[1] : nbPixVol - 1;
	int voxel_zPP = voxel_zP + pitch[2] < nbPixVol ? voxel_zP + pitch[2] : nbPixVol - 1;
	int tPP;
	if (periodic_motion) tPP = tP + 1 < size[3] ? tP + 1 : tP + 1 - size[3];
	else tPP = tP + 1 < size[3] ? tP + 1 : size[3] - 1;
	
	double phiVal = listPhi[numPhi]->get_phi_value(voxel + time * nbPixVol);
	
	double dphixPP = listPhi[numPhi]->get_phi_value(voxel_xPP + time * nbPixVol) - 2. * listPhi[numPhi]->get_phi_value(voxel_xP + time * nbPixVol) + phiVal;
	double dphiyPP = listPhi[numPhi]->get_phi_value(voxel_yPP + time * nbPixVol) - 2. * listPhi[numPhi]->get_phi_value(voxel_yP + time * nbPixVol) + phiVal;
	double dphizPP = listPhi[numPhi]->get_phi_value(voxel_zPP + time * nbPixVol) - 2. * listPhi[numPhi]->get_phi_value(voxel_zP + time * nbPixVol) + phiVal;
	double dphitPP = listPhi[numPhi]->get_phi_value(voxel + tPP * nbPixVol) - 2. * listPhi[numPhi]->get_phi_value(voxel + tP * nbPixVol) + phiVal;
	
	return Vector4(dphixPP, dphiyPP, dphizPP, dphitPP);
}

Vector4 SpeedFunction_Geom::computeBackwardBackwardDiff(int voxel, int time, int numPhi)
{
	int voxel_xM = voxel - pitch[0] >= 0 ? voxel - pitch[0] : 0;
	int voxel_yM = voxel - pitch[1] >= 0 ? voxel - pitch[1] : 0;
	int voxel_zM = voxel - pitch[2] >= 0 ? voxel - pitch[2] : 0;
	int tM;
	if (periodic_motion) tM = time - 1 >= 0 ? time - 1 : time - 1 + size[3];
	else tM = time - 1 >= 0 ? time -1 : 0;
	
	int voxel_xMM = voxel_xM - pitch[0] >= 0 ? voxel_xM - pitch[0] : 0;
	int voxel_yMM = voxel_yM - pitch[1] >= 0 ? voxel_yM - pitch[1] : 0;
	int voxel_zMM = voxel_zM - pitch[2] >= 0 ? voxel_zM - pitch[2] : 0;
	int tMM;
	if (periodic_motion) tMM = tM - 1 >= 0 ? tM - 1 : tM - 1 + size[3];
	else tMM = tM - 1 >= 0 ? tM -1 : 0;
	
	double phiVal = listPhi[numPhi]->get_phi_value(voxel + time * nbPixVol);
	
	double dphixMM = phiVal - 2. * listPhi[numPhi]->get_phi_value(voxel_xM + time * nbPixVol) + listPhi[numPhi]->get_phi_value(voxel_xMM + time * nbPixVol);
	double dphiyMM = phiVal - 2. * listPhi[numPhi]->get_phi_value(voxel_yM + time * nbPixVol) + listPhi[numPhi]->get_phi_value(voxel_yMM + time * nbPixVol);
	double dphizMM = phiVal - 2. * listPhi[numPhi]->get_phi_value(voxel_zM + time * nbPixVol) + listPhi[numPhi]->get_phi_value(voxel_zMM + time * nbPixVol);
	double dphitMM = phiVal - 2. * listPhi[numPhi]->get_phi_value(voxel + tM * nbPixVol) + listPhi[numPhi]->get_phi_value(voxel + tMM * nbPixVol);
	
	return Vector4(dphixMM, dphiyMM, dphizMM, dphitMM);
}


/// SpeedFunction_normalisation ///

SpeedFunction_normalisation::SpeedFunction_normalisation(vector<Phi*>& _listPhi, const int* _size, bool _periodic_motion): SpeedFunction_Geom(_listPhi, _size, _periodic_motion)
{
	border_values.insert(0);
	border_values.insert(size[0]-1);
	border_values.insert(size[0] * (size[1]-1));
	border_values.insert(size[0]-1 + size[0] * (size[1]-1));
	border_values.insert(size[0] * size[1] * (size[2]-1));
	border_values.insert(size[0]-1 + size[0] * size[1] * (size[2]-1));
	border_values.insert(size[0] * (size[1]-1 + size[1] * (size[2]-1) ));
	border_values.insert(size[0]-1 + size[0] * (size[1]-1) + size[1] * (size[2]-1));
}

void SpeedFunction_normalisation::prepare_iteration()
{
	max_speed = 0;
}

float SpeedFunction_normalisation::get_speed(int voxel, int time, int numPhi)
{
	float phiVal = listPhi[numPhi]->get_phi_value(voxel + time * nbPixVol);
	float signPhi;
	
	if (phiVal > 0) signPhi = 1;
	else if (phiVal == 0) signPhi = 0;
	else signPhi = -1;
	
	double normGradPhi = computeNormGradientEntropy2ndOrder(voxel, time, -signPhi, numPhi);
	
	if ( normGradPhi == 0 && border_values.count(voxel) )
		normGradPhi = 1;
	
	if (phiVal == 0) signPhi = 0;
	else signPhi = phiVal / sqrt((double)(phiVal * phiVal) + normGradPhi);
	
	float res = signPhi * (1. - normGradPhi);
	
	#pragma omp critical
	max_speed = max( fabs(res), max_speed );
	
	return res;
}

float SpeedFunction_normalisation::get_dt()
{
	if (!max_speed) return 0;
	
	float dt = 0.9 / max_speed;
	
	if (dt > 2)
	{
		dt = 2;
		cout << "SpeedFunction_normalisation::get_dt(): dt too high, limiting to 2" << endl;
	}
	
	return dt;
}
