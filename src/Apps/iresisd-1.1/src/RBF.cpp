#include "RBF.h"

using namespace std;

RBF::RBF(float _beta, float _gamma, const int* size_data, vector<int> _borders)
{
	beta = _beta;
	gamma = _gamma;
	borders = _borders;
	
	beta_ov_2 = beta / 2.;
	gamma_square = gamma * gamma;
	
	/// volume padding, to prevent the contour from leaking on one side and re-appearing on the other, due to the fft
	
	for (int i=0; i<4; i++) size_db[i] = size_data[i] + borders[i];
	
	
	/// allocate memory and create plans
	
	//the FFT requires the arrays to be padded in the first dimension
	size_x_pad_Fourier = (int)(floor(size_db[0] / 2.) + 1);
	size_x_pad = 2 * size_x_pad_Fourier;
	
	size_fft[0] = size_x_pad;
	size_fft[1] = size_db[1];
	size_fft[2] = size_db[2];
	size_fft[3] = size_db[3];
	nbPixPad = size_x_pad * size_db[1] * size_db[2] * size_db[3];
	cout << "volume size for fft : " << size_fft[0] << " " << size_fft[1] << " " << size_fft[2] << " " << size_fft[3] << endl;
	
	FFTWarrayPsi = (double*) fftw_malloc(sizeof(double) * nbPixPad);
	
	
	if (size_db[3] == 1)
	{
		if (size_db[2] == 1)
		{
			pF = fftw_plan_dft_r2c_2d(size_db[1], size_db[0], FFTWarrayPsi, (fftw_complex*)FFTWarrayPsi, FFTW_ESTIMATE);
			pB = fftw_plan_dft_c2r_2d(size_db[1], size_db[0], (fftw_complex*)FFTWarrayPsi, FFTWarrayPsi, FFTW_ESTIMATE);
		}
		else
		{
			pF = fftw_plan_dft_r2c_3d(size_db[2], size_db[1], size_db[0], FFTWarrayPsi, (fftw_complex*)FFTWarrayPsi, FFTW_ESTIMATE);
			pB = fftw_plan_dft_c2r_3d(size_db[2], size_db[1], size_db[0], (fftw_complex*)FFTWarrayPsi, FFTWarrayPsi, FFTW_ESTIMATE);
		}
	}
	else
	{
		int size_inv[4];
		size_inv[0] = size_db[3];
		size_inv[1] = size_db[2];
		size_inv[2] = size_db[1];
		size_inv[3] = size_db[0];
		
		pF = fftw_plan_dft_r2c(4, size_inv, FFTWarrayPsi, (fftw_complex*)FFTWarrayPsi, FFTW_ESTIMATE);
		pB = fftw_plan_dft_c2r(4, size_inv, (fftw_complex*)FFTWarrayPsi, FFTWarrayPsi, FFTW_ESTIMATE);
	}
	
	
	/// fill the arrays and pre-compute the fft
	
	float xC, yC, zC, tC;
	xC = (float)(size_db[0])/2.;
	yC = (float)(size_db[1])/2.;
	
	if (size_db[3] == 1) tC = 0;
	else tC = (float)(size_db[3])/2.;
	
	if (size_db[2] == 1) zC = 0;
	else zC = (float)(size_db[2])/2.;
	
	int x, y, z, t, ind;
	
#ifdef _OPENMP
	#if GCC_VERSION >= 40400 && _OPENMP >= 200805
	#pragma omp parallel for default(shared) private(x, y, z, t, ind) schedule(dynamic) num_threads(NB_THREADS) collapse(4)
	for (x=0; x<size_fft[0]; x++)
	for (y=0; y<size_fft[1]; y++)
	for (z=0; z<size_fft[2]; z++)
	for (t=0; t<size_fft[3]; t++)
	#else
	#pragma omp parallel for default(shared) private(ind, x, y, z, t) schedule(dynamic) num_threads(NB_THREADS)
	for (ind=0; ind<nbPixPad; ind++)
	#endif
#else
	for (x=0; x<size_fft[0]; x++)
	for (y=0; y<size_fft[1]; y++)
	for (z=0; z<size_fft[2]; z++)
	for (t=0; t<size_fft[3]; t++)
#endif
	{
#ifdef _OPENMP
	#if GCC_VERSION >= 40400 && _OPENMP >= 200805
		ind = computeIndex(size_fft, x, y, z, t);
	#else
		calculCoords(ind, size_fft, x, y, z, t);
	#endif
#else
		ind = computeIndex(size_fft, x, y, z, t);
#endif
		if (x < size_db[0])
			FFTWarrayPsi[ind] = pow( pow(xC - x, 2) + pow(yC - y, 2) + pow(zC - z, 2) + pow(tC - t, 2) + gamma_square, - beta_ov_2 );
		else FFTWarrayPsi[ind] = 0;
	}
	
	fftw_execute_dft_r2c(pF, FFTWarrayPsi, (fftw_complex*)FFTWarrayPsi);
}

RBF::~RBF()
{
	fftw_destroy_plan(pF);
	fftw_destroy_plan(pB);
	
	fftw_free(FFTWarrayPsi);
}

void RBF::do_fft(double* array, int direction)
{
	switch (direction)
	{
		case 1:
			fftw_execute_dft_r2c(pF, array, (fftw_complex*)array);
			break;
			
		case -1:
			fftw_execute_dft_c2r(pB, (fftw_complex*)array, array);
			break;
			
		default:
			throw string("RBF::do_fft() Error: unknown direction");
			break;
	}
}

void RBF::multiply_with_RBF(double* fft_array_ini, double* fft_array_res)
{
	double ini_real, ini_img, psi_real, psi_img, res_real, res_img;
	
	for (int ind=0; ind<nbPixPad; ind+=2)
	{
		ini_real = fft_array_ini[ind];
		ini_img = fft_array_ini[ind+1];
		
		psi_real = FFTWarrayPsi[ind];
		psi_img = FFTWarrayPsi[ind+1];
		
		res_real = ini_real * psi_real - ini_img * psi_img;
		res_img = ini_real * psi_img + ini_img * psi_real;
		
		fft_array_res[ind] = res_real;
		fft_array_res[ind+1] = res_img;
	}
}

const int* RBF::get_size_fft()
{
	return &(size_fft[0]);
}

const vector<int>& RBF::get_borders()
{
	return borders;
}

float RBF::get_beta()
{
	return beta;
}

float RBF::get_gamma()
{
	return gamma;
}

