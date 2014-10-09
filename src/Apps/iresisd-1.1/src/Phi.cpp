#include "Phi.h"

typedef float MY_DATATYPE;

#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352

using namespace std;


Phi::Phi(const int* s)
{
	cout << "constructor Phi" << endl;

	size = s;

	nbPixVol = size[0] * size[1] * size[2];
	nbPix = nbPixVol * size[3];

	contour = new map<int, float>[ size[3] ];
	data = new float[nbPix];
	for (int i=0; i<nbPix; i++) data[i] = -VAL_INVALID;
	
	cout << "end constructor Phi" << endl;
}

Phi::~Phi()
{
	cout << "destructor Phi" << endl;
	
	delete[] contour;
	delete[] data;
	
	cout << "end destructor Phi" << endl;
}

void Phi::initialiseSpheres(const list< Point4<int> >& centreIni, const list<float>& radiusIni)
{
	cout << "Phi::initialiseSpheres()" << endl;

	int ind, x, y, z, t;
	float distcentreTmp, newValue, valTmp;
	
	/** compute the distances to the sphere **/
	
	list< Point4<int> >::const_iterator lit, lend(centreIni.end());
	list<float>::const_iterator litR;
	
#ifdef _OPENMP
	#if GCC_VERSION >= 40400 && _OPENMP >= 200805
	#pragma omp parallel for default(shared) private(x, y, z, t, ind, newValue, lit, litR, distcentreTmp, valTmp) schedule(dynamic) num_threads(NB_THREADS) collapse(4)
	for (x=0; x<size[0]; x++)
	for (y=0; y<size[1]; y++)
	for (z=0; z<size[2]; z++)
	for (t=0; t<size[3]; t++)
	#else
	#pragma omp parallel for default(shared) private(x, y, z, t, ind, newValue, lit, litR, distcentreTmp, valTmp) schedule(dynamic) num_threads(NB_THREADS)
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
		calculCoords(ind, size, x, y, z, t);
	#endif
#else
		ind = computeIndex(size, x, y, z, t);
#endif
		newValue = data[ind];

		for (lit = centreIni.begin(), litR = radiusIni.begin(); lit != lend; ++lit, ++litR) //for each sphere
		{
			//distance to centre
			distcentreTmp = sqrt( pow((float)(x - lit->x()), 2) +
					pow((float)(y - lit->y()), 2) +
					pow((float)(z - lit->z()), 2) +
					pow((float)(t - lit->t()), 2) );
			
			//comparison with radius
			valTmp = *litR - distcentreTmp;
			
			//comparison with current value
			if (valTmp > newValue)
			{
				newValue = valTmp;
			}
		}
		
		/*if (newValue > 3) newValue = 3;
		if (newValue < -3) newValue = -3;*/
			
		data[ind] = newValue;
	}
	
	//display_phi();
}

void Phi::initialise_from_nii(string nii_file)
{
	cout << "Phi::initialise_from_nii()" << endl;

	nifti_1_header hdr;
	nifti1_extender pad={0,0,0,0};
	FILE *fp;
	int ret;
	
	///// read first 348 bytes of header   //
	
	fp = fopen(nii_file.c_str(),"r");
	if (fp == NULL)
	{
		fprintf(stderr, "\nError opening header file phi.nii for read\n");
		exit(1);
	}
	
	ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
	if (ret != 1)
	{
		fprintf(stderr, "\nError writing header file phi.nii\n");
		exit(1);
	}
	
	int nbDims = hdr.dim[0];
	
	if (size[0] != hdr.dim[1] || size[1] != hdr.dim[2] || size[2] != hdr.dim[3]) throw string("Phi::initialise_from_nii() Error: wrong volume size");
	if (nbDims == 4) if (size[3] != hdr.dim[4]) throw string("Phi::initialise_from_nii() Error: wrong volume size");
	
	///// read extender pad and image data   //
	ret = fread(&pad, 4, 1, fp);
	if (ret != 1)
	{
		fprintf(stderr, "\nError reading header file extension pad phi.nii\n");
		exit(1);
	}
	
	int size_read = hdr.dim[1]*hdr.dim[2]*hdr.dim[3];
	if (nbDims == 4) size_read *= hdr.dim[4];
		
	ret = fread(&data, (size_t)(hdr.bitpix/8), size_read, fp);
	if (ret != size_read)
	{
		fprintf(stderr, "\nError reading data from phi.nii\n");
		exit(1);
	}
	
	fclose(fp);
	
	cout << "Phi::initialise_from_nii() done" << endl;
}

const int* Phi::get_volume_size()
{
	return size;
}

void Phi::set_phi_value(int ind_phi, float newVal)
{
	data[ind_phi] = newVal;
}

float Phi::get_phi_value(int x, int y, int z, int t) const
{
	if (x < 0) x = 0;
	if (x > size[0] - 1) x = size[0] - 1;
	if (y < 0) y = 0;
	if (y > size[1] - 1) y = size[1] - 1;
	if (z < 0) z = 0;
	if (z > size[2] - 1) z = size[2] - 1;
	if (t < 0) t = 0;
	if (t > size[3] - 1) t = size[3] - 1;
	
	return data[ computeIndex(size, x, y, z, t) ];
}

float Phi::get_phi_value(Point4<int> point) const
{
	return get_phi_value(point.x(), point.y(), point.z(), point.t());
}

float Phi::get_phi_value(int ind_phi) const
{
	if (ind_phi < 0) ind_phi = 0;
	if (ind_phi > nbPix - 1) ind_phi = nbPix - 1;
	
	return data[ind_phi];
}

float Phi::get_phi_value(Point4<float> point) const
{
	return get_phi_value(point.x(), point.y(), point.z(), point.t());
}

float Phi::get_phi_value(float x, float y, float z, float t) const
{
	if (x < 0) x = 0;
	if (x > size[0] - 1) x = size[0] - 1;
	if (y < 0) y = 0;
	if (y > size[1] - 1) y = size[1] - 1;
	if (z < 0) z = 0;
	if (z > size[2] - 1) z = size[2] - 1;
	
	if (t <= -1 || t >= size[3])
	{
		cout << size[3] << endl;
		cout << x << " " << y << " " << z << " " << t << endl;
		throw string("Phi::get_phi_value(float x, float y, float z, float t) Error : t out of bounds");
	}
	
	if (t < 0) t = size[3] - t; //periodic signal: time 0 = time size[3]
	
	
	int xM = (int)floor(x);
	int xP = (int)ceil(x);
	
	int yM = (int)floor(y);
	int yP = (int)ceil(y);
	
	int zM = (int)floor(z);
	int zP = (int)ceil(z);
	
	int tM = (int)floor(t);
	int tP = (int)ceil(t);
	if (tP == size[3]) tP = 0;
	
	/// tri-linear interpolation of tM
	
	float ax = x - xM;
	float ay = y - yM;
	float az = z - zM;
	float at = t - tM;
	
	if (!at && !az && !ay && !ax) return get_phi_value(xM, yM, zM, tM);
	
	float valtM, valtP;
	
	//plane zM
	
	float valMMM = get_phi_value(xM, yM, zM, tM);
	float valPMM = get_phi_value(xP, yM, zM, tM);
	float valMPM = get_phi_value(xM, yP, zM, tM);
	float valPPM = get_phi_value(xP, yP, zM, tM);
	
	float valCMM = ax * valPMM + (1. - ax) * valMMM;
	float valCPM = ax * valPPM + (1. - ax) * valMPM;
	
	float valCCM = ay * valCPM + (1. - ay) * valCMM;
	
	if (!az)
		valtM = valCCM;
	else
	{
		//plane zP
		
		float valMMP = get_phi_value(xM, yM, zP, tM);
		float valPMP = get_phi_value(xP, yM, zP, tM);
		float valMPP = get_phi_value(xM, yP, zP, tM);
		float valPPP = get_phi_value(xP, yP, zP, tM);
		
		float valCMP = ax * valPMP + (1. - ax) * valMMP;
		float valCPP = ax * valPPP + (1. - ax) * valMPP;
		
		float valCCP = ay * valCPP + (1. - ay) * valCMP;
		
		valtM = az * valCCP + (1. - az) * valCCM;
	}
	
	if (!at) return valtM;
	
	/// tri-linear interpolation of tP
	
	valMMM = get_phi_value(xM, yM, zM, tP);
	valPMM = get_phi_value(xP, yM, zM, tP);
	valMPM = get_phi_value(xM, yP, zM, tP);
	valPPM = get_phi_value(xP, yP, zM, tP);
	
	//plane zM
	
	valCMM = ax * valPMM + (1. - ax) * valMMM;
	valCPM = ax * valPPM + (1. - ax) * valMPM;
	
	valCCM = ay * valCPM + (1. - ay) * valCMM;
	
	if (!az)
		valtP = valCCM;
	else
	{
		//plane zP
		
		float valMMP = get_phi_value(xM, yM, zP, tP);
		float valPMP = get_phi_value(xP, yM, zP, tP);
		float valMPP = get_phi_value(xM, yP, zP, tP);
		float valPPP = get_phi_value(xP, yP, zP, tP);
		
		float valCMP = ax * valPMP + (1. - ax) * valMMP;
		float valCPP = ax * valPPP + (1. - ax) * valMPP;
		
		float valCCP = ay * valCPP + (1. - ay) * valCMP;
		
		valtP = az * valCCP + (1. - az) * valCCM;
	}
	
	/// linear interpolation of t
	
	return at * valtP + (1. - at) * valtM;
}

void Phi::computeContour(float epsilon)
{
	cout << "Phi::computeContour() " << epsilon << endl;
	
	int epsilon_int = (int)ceil(epsilon);
	
	contourSize = 0;
	
	int x, y, z, f, ind;
	
	for (f=0; f<size[3]; f++)
	{
		contour[f].clear();
		
#ifdef _OPENMP
	#if GCC_VERSION >= 40400 && _OPENMP >= 200805
		#pragma omp parallel for default(shared) private(x, y, z) schedule(dynamic) num_threads(NB_THREADS) collapse(3)
		for (x=0; x<size[0]; x++)
		for (y=0; y<size[1]; y++)
		for (z=0; z<size[2]; z++)
	#else
		#pragma omp parallel for default(shared) private(x, y, z, ind) schedule(dynamic) num_threads(NB_THREADS)
		for (ind=0; ind<nbPixVol; ind++)
	#endif
#else
		for (x=0; x<size[0]; x++)
		for (y=0; y<size[1]; y++)
		for (z=0; z<size[2]; z++)
#endif
		{
#ifdef _OPENMP
	#if GCC_VERSION < 40400 || _OPENMP < 200805
			calculCoords(ind, size, x, y, z);
	#endif
#endif
			bool isContour = false;
			float dist = VAL_INVALID;
			
			float val = get_phi_value(x, y, z, f);

			if (val == 0) isContour = true;
			else
			{
				int supX = x+epsilon_int < size[0] ? x+epsilon_int : size[0]-1;
				int infX = x-epsilon_int >= 0 ? x-epsilon_int : 0;
				int supY = y+epsilon_int < size[1] ? y+epsilon_int : size[1]-1;
				int infY = y-epsilon_int >= 0 ? y-epsilon_int : 0;
				int supZ = z+epsilon_int < size[2] ? z+epsilon_int : size[2]-1;
				int infZ = z-epsilon_int >= 0 ? z-epsilon_int : 0;
				int lP = f+epsilon_int < size[3] ? f+epsilon_int : size[3]-1;
				int lM = f-epsilon_int >= 0 ? f-epsilon_int : 0;
				
				int xn[8] = {infX, supX, x, x, x, x, x, x};
				int yn[8] = {y, y, infY, supY, y, y, y, y};
				int zn[8] = {z, z, z, z, infZ, supZ, z, z};
				int fn[8] = {f, f, f, f, f, f, lP, lM};
				
				for (int in=0; in<8; in++)
				{
					if (val * get_phi_value(xn[in], yn[in], zn[in], fn[in]) <= 0)
					{
						isContour = true;
						break;
					}
				}
			}
			
			if (!isContour) continue;
			
			if (val == 0) dist = 0;
			else
			{
				/// check neighbours
				
				float valN, tmp, t, xC, yC, zC, fC;
				
				int iP = x+epsilon_int < size[0] ? x+epsilon_int : size[0]-1;
				int iM = x-epsilon_int >= 0 ? x-epsilon_int : 0;
				int jP = y+epsilon_int < size[1] ? y+epsilon_int : size[1]-1;
				int jM = y-epsilon_int >= 0 ? y-epsilon_int : 0;
				int kP = z+epsilon_int < size[2] ? z+epsilon_int : size[2]-1;
				int kM = z-epsilon_int >= 0 ? z-epsilon_int : 0;
				int lP = f+epsilon_int < size[3] ? f+epsilon_int : size[3]-1;
				int lM = f-epsilon_int >= 0 ? f-epsilon_int : 0;
				
				float di, dj, dk, dl;
				
				for (int i=iM; i<=iP; i++)
				for (int j=jM; j<=jP; j++)
				for (int k=kM; k<=kP; k++)
				for (int l=lM; l<=lP; l++)
				{
					if (i == x && j == y && k == z && l == f) continue;
					
					valN = get_phi_value(i, j, k, l);
					
					if (val * valN <= 0)
					{
						t = -val / (valN - val);
						xC = x + t * (i - x);
						yC = y + t * (j - y);
						zC = z + t * (k - z);
						fC = f + t * (l - f);
						
						di = x - xC;
						dj = y - yC;
						dk = z - zC;
						dl = f - fC;
						tmp = di*di + dj*dj + dk*dk + dl*dl;
						
						if (tmp < dist) dist = tmp;
					}
				}

				if (dist == VAL_INVALID)
				{
					cout << "Phi::computeContour() error : no  distance to contour found!" << endl;
					#pragma omp critical
					cv::waitKey();
				}
				else if (dist) dist = sqrt(dist);
			}
			
			if (dist > epsilon) continue;
			
#ifdef _OPENMP
	#if GCC_VERSION >= 40400 && _OPENMP >= 200805
			#pragma omp critical
			contour[f][ computeIndex(size, x, y, z) ] = dist;
	#else
			#pragma omp critical
			contour[f][ind] = dist;
	#endif
#else
			#pragma omp critical
			contour[f][ computeIndex(size, x, y, z) ] = dist;
#endif			
		} // voxels
		
		contourSize += contour[f].size();
	} // timeframes
}

bool Phi::isContourPoint(int voxel, int time)
{
	map<int, float>::const_iterator mit = contour[time].find(voxel);

	if (mit == contour[time].end()) return false;
	return true;
}

const map<int, float>* Phi::getContour()
{
	return contour;
}

int Phi::getContourSize()
{
	return contourSize;
}

float Phi::getDistToContour(int voxel, int time)
{
	map<int, float>::const_iterator mit = contour[time].find(voxel);
	if (mit == contour[time].end()) return VAL_INVALID;
	return mit->second;
}

void Phi::savePhi(int numPhi, int t, string str)
{
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
	hdr.dim[4] = size[3];
	hdr.datatype = NIFTI_TYPE_FLOAT32;
	hdr.bitpix = 32;
	hdr.pixdim[1] = 1;
	hdr.pixdim[2] = 1;
	hdr.pixdim[3] = 1;
	hdr.pixdim[4] = 1;
	hdr.vox_offset = (float) NII_HEADER_SIZE;
	hdr.xyzt_units = NIFTI_UNITS_MM;
	strncpy(hdr.magic, "n+1\0", 4);
	
	
	/********** write first 348 bytes of header   */
	ostringstream name4;
	name4 << "phi" << numPhi << "_" << str << ".nii";
	fp = fopen(name4.str().c_str(),"w");
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
	
	ret = fwrite(&data, (size_t)(hdr.bitpix/8), hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.dim[4], fp);
	if (ret != hdr.dim[1] * hdr.dim[2] * hdr.dim[3] * hdr.dim[4])
	{
		fprintf(stderr, "\nError writing data to phi.nii\n");
		exit(1);
	}
	
	fclose(fp);
}

void Phi::saveContour(int numPhi, int t, string str)
{
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
	hdr.xyzt_units = NIFTI_UNITS_MM;
	strncpy(hdr.magic, "n+1\0", 4);
	
	
	/********** write first 348 bytes of header   */
	ostringstream name;
	name << "phiContour" << numPhi << "_t=" << t << "_" << str << ".nii";
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
	
	float* buffer = new float[size[0] * size[1] * size[2]];

	int ind;
	float pix, dist_cont;
	
	for (int x=0; x<size[0]; x++)
	for (int y=0; y<size[1]; y++)
	for (int z=0; z<size[2]; z++)
	{
		ind = computeIndex(size, x, y, z);

		pix = get_phi_value(x, y, z, t);
		dist_cont = getDistToContour(ind, 0);

		if (dist_cont == 0) pix = 0;
		else if (dist_cont <= 5)
		{
			if (pix < 0) pix = -dist_cont;
			else if (pix > 0) pix = dist_cont;
		}
		else
		{
			if (pix < 0) pix = -5;
			else if (pix > 0) pix = 5;
		}
		
		buffer[ind] = pix;
	}
	
	ret = fwrite(&buffer[0], (size_t)(hdr.bitpix/8), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], fp);
	if (ret != hdr.dim[1] * hdr.dim[2] * hdr.dim[3])
	{
		fprintf(stderr, "\nError writing data to phi.nii\n");
		exit(1);
	}
	
	delete[] buffer;

	fclose(fp);
}

Vector3 Phi::getSpatialNormal(int x, int y, int z, int f)
{
	int xP = x + 1 < size[0] ? x + 1 : size[0] - 1;
	int yP = y + 1 < size[1] ? y + 1 : size[1] - 1;
	int zP = z + 1 < size[2] ? z + 1 : size[2] - 1;
	
	int xM = x - 1 >= 0 ? x - 1 : 0;
	int yM = y - 1 >= 0 ? y - 1 : 0;
	int zM = z - 1 >= 0 ? z - 1 : 0;
	
	float val = get_phi_value(x, y, z, f);
	float valxP = get_phi_value(xP, y, z, f);
	float valyP = get_phi_value(x, yP, z, f);
	float valzP = get_phi_value(x, y, zP, f);
	float valxM = get_phi_value(xM, y, z, f);
	float valyM = get_phi_value(x, yM, z, f);
	float valzM = get_phi_value(x, y, zM, f);
	
	float dphiP[3] = {valxP - val, valyP - val, valzP - val};
	float dphiM[3] = {val - valxM, val - valyM, val - valzM};
	
	float dphixP2 = dphiP[0] * dphiP[0];
	float dphixM2 = dphiM[0] * dphiM[0];
	float dphiyP2 = dphiP[1] * dphiP[1];
	float dphiyM2 = dphiM[1] * dphiM[1];
	float dphizP2 = dphiP[2] * dphiP[2];
	float dphizM2 = dphiM[2] * dphiM[2];
	
	Vector3 tmp(0, 0, 0);
	
	float denom = dphixP2 + dphiyP2 + dphizP2;
	if (denom)
		for (int i=0; i<3; i++) tmp[i] += dphiP[i] / sqrtf(denom);
	
	denom = dphixM2 + dphiyP2 + dphizP2;
	if (denom)
	{
		denom = sqrtf(denom);
		tmp[0] += dphiM[0] / denom;
		tmp[1] += dphiP[1] / denom;
		tmp[2] += dphiP[2] / denom;
	}
	
	denom = dphixP2 + dphiyM2 + dphizP2;
	if (denom)
	{
		denom = sqrtf(denom);
		tmp[0] += dphiP[0] / denom;
		tmp[1] += dphiM[1] / denom;
		tmp[2] += dphiP[2] / denom;
	}
	
	denom = dphixM2 + dphiyM2 + dphizP2;
	if (denom)
	{
		denom = sqrtf(denom);
		tmp[0] += dphiM[0] / denom;
		tmp[1] += dphiM[1] / denom;
		tmp[2] += dphiP[2] / denom;
	}
	
	denom = dphixP2 + dphiyP2 + dphizM2;
	if (denom)
	{
		denom = sqrtf(denom);
		tmp[0] += dphiP[0] / denom;
		tmp[1] += dphiP[1] / denom;
		tmp[2] += dphiM[2] / denom;
	}
	
	denom = dphixM2 + dphiyP2 + dphizM2;
	if (denom)
	{
		denom = sqrtf(denom);
		tmp[0] += dphiM[0] / denom;
		tmp[1] += dphiP[1] / denom;
		tmp[2] += dphiM[2] / denom;
	}
	
	denom = dphixP2 + dphiyM2 + dphizM2;
	if (denom)
	{
		denom = sqrtf(denom);
		tmp[0] += dphiP[0] / denom;
		tmp[1] += dphiM[1] / denom;
		tmp[2] += dphiM[2] / denom;
	}
	
	denom = dphixM2 + dphiyM2 + dphizM2;
	if (denom)
		for (int i=0; i<3; i++) tmp[i] += dphiM[i] / sqrtf(denom);
	
	float normN = tmp.norm();
	
	if (!normN) return VECT_NULL;
	else
	{
		tmp[0] = - tmp[0] / normN; //phi>0 inside contour -> n=-grad(phi)/|grad(phi)|
		tmp[1] = - tmp[1] / normN;
		tmp[2] = - tmp[2] / normN;
	}
	
	return tmp;
}

Vector3 Phi::getSpatialNormal(float x, float y, float z, float t)
{
	if (x < 0) x = 0;
	if (x > size[0] - 1) x = size[0] - 1;
	if (y < 0) y = 0;
	if (y > size[1] - 1) y = size[1] - 1;
	if (z < 0) z = 0;
	if (z > size[2] - 1) z = size[2] - 1;
	
	if (t < 0 || t >= size[3])
	{
		cout << size[3] << endl;
		cout << x << " " << y << " " << z << " " << t << endl;
		vector<int> test;
		cout << test[10] << endl;
		throw string("Phi::getSpatialNormal(float x, float y, float z, float t) Erreur : t out of bounds");
	}
	
	int xM = (int)floor(x);
	int xP = (int)ceil(x);
	
	int yM = (int)floor(y);
	int yP = (int)ceil(y);
	
	int zM = (int)floor(z);
	int zP = (int)ceil(z);
	
	int tM = (int)floor(t);
	int tP = (int)ceil(t);
	if (tP == size[3]) tP = 0;
	
	/// tri-linear interpolation for tM
	
	float ax = x - xM;
	float ay = y - yM;
	float az = z - zM;
	float at = t - tM;
	
	if (!at && !az && !ay && !ax) return getSpatialNormal(xM, yM, zM, tM);
	
	Vector3 valtM, valtP;
	
	//plane zM
	
	Vector3 valMMM = getSpatialNormal(xM, yM, zM, tM);
	Vector3 valPMM = getSpatialNormal(xP, yM, zM, tM);
	Vector3 valMPM = getSpatialNormal(xM, yP, zM, tM);
	Vector3 valPPM = getSpatialNormal(xP, yP, zM, tM);
	
	Vector3 valCMM = valPMM * ax + valMMM * (1. - ax);
	Vector3 valCPM = valPPM * ax + valMPM * (1. - ax);
	
	Vector3 valCCM = valCPM * ay + valCMM * (1. - ay);
	
	if (!az)
		valtM = valCCM;
	else
	{
		//plane zP
		
		Vector3 valMMP = getSpatialNormal(xM, yM, zP, tM);
		Vector3 valPMP = getSpatialNormal(xP, yM, zP, tM);
		Vector3 valMPP = getSpatialNormal(xM, yP, zP, tM);
		Vector3 valPPP = getSpatialNormal(xP, yP, zP, tM);
		
		Vector3 valCMP = valPMP * ax + valMMP * (1. - ax);
		Vector3 valCPP = valPPP * ax + valMPP * (1. - ax);
		
		Vector3 valCCP = valCPP * ay + valCMP * (1. - ay);
		
		
		valtM = valCCP * az + valCCM * (1. - az);
	}
	
	if (!at) return valtM;
	
	/// tri-linear interpolation for tP
	
	valMMM = getSpatialNormal(xM, yM, zM, tP);
	valPMM = getSpatialNormal(xP, yM, zM, tP);
	valMPM = getSpatialNormal(xM, yP, zM, tP);
	valPPM = getSpatialNormal(xP, yP, zM, tP);
	
	//plane zM
	
	valCMM = valPMM * ax + valMMM * (1. - ax);
	valCPM = valPPM * ax + valMPM * (1. - ax);
	
	valCCM = valCPM * ay + valCMM * (1. - ay);
	
	if (!az)
		valtP = valCCM;
	else
	{
		//plane zP
		
		Vector3 valMMP = getSpatialNormal(xM, yM, zP, tP);
		Vector3 valPMP = getSpatialNormal(xP, yM, zP, tP);
		Vector3 valMPP = getSpatialNormal(xM, yP, zP, tP);
		Vector3 valPPP = getSpatialNormal(xP, yP, zP, tP);
		
		Vector3 valCMP = valPMP * ax + valMMP * (1. - ax);
		Vector3 valCPP = valPPP * ax + valMPP * (1. - ax);
		
		Vector3 valCCP = valCPP * ay + valCMP * (1. - ay);
		
		
		valtP = valCCP * az + valCCM * (1. - az);
	}
	
	/// linear interpolation in t direction
	
	return valtP * at + valtM * (1. - at);
}

float Phi::computeNormGradientInNormalDirection(int x, int y, int z, int t)
{
	Vector3 normal = getSpatialNormal(x, y, z, t);
	
	float phiM = get_phi_value(x + normal[0], y + normal[1], z + normal[2], t);
	float phiP = get_phi_value(x - normal[0], y - normal[1], z - normal[2], t);
	
	if (phiM * phiP > 0) return 1;
	if (phiM >= phiP) return 1;
	
	return (phiP - phiM) / 2.;
}

void Phi::display_phi()
{
	int z = size[2]/2;
	
	cv::Mat imgDisp(size[0], size[1], CV_64FC1);
	
	double max = 0, min = VAL_INVALID;
	
	for(int i=0; i<size[0]; i++)
	for(int j=0; j<size[1]; j++)
	{
		imgDisp.at<double>(i, j) = (double)get_phi_value(i, j, z, 0);
		
		if (fabs(imgDisp.at<double>(i, j)) > max) max = fabs(imgDisp.at<double>(i, j));
		if (fabs(imgDisp.at<double>(i, j)) < min) min = fabs(imgDisp.at<double>(i, j));
	}
	
	imgDisp /= 2.*max;
	imgDisp += 0.5;
	
	cv::namedWindow("phi central horizontal slice", CV_WINDOW_AUTOSIZE);
	cv::imshow("phi central horizontal slice", imgDisp);

	if (size[2] > 1)
	{
		int x = size[0]/2;
		
		cv::Mat imgDisp2(size[2], size[1], CV_64FC1);
		
		max = 0;
		min = VAL_INVALID;
		
		for (int j=0; j<size[2]; j++)
		for (int i=0; i<size[1]; i++)
		{
			imgDisp2.at<double>(j, i) = (double)get_phi_value(x, i, j, 0);
			
			if (fabs(imgDisp2.at<double>(j, i)) > max) max = fabs(imgDisp2.at<double>(j, i));
			if (fabs(imgDisp2.at<double>(j, i)) < min) min = fabs(imgDisp2.at<double>(j, i));
		}
		
		imgDisp2 /= 2.*max;
		imgDisp2 += 0.5;
		
		cv::namedWindow("phi central vertical slice", CV_WINDOW_AUTOSIZE);
		cv::imshow("phi central vertical slice", imgDisp2);
	}
	
	cv::waitKey();
}
