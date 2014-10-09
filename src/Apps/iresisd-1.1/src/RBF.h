#ifndef RBF_H
#define RBF_H

#include <fftw3.h>

#include "AuxClasses.h"

#ifndef NO_OMP
	#include <omp.h>
#endif

#define NB_THREADS omp_get_num_procs()

#define GCC_VERSION (__GNUC__ * 10000 \
                               + __GNUC_MINOR__ * 100 \
                               + __GNUC_PATCHLEVEL__)

//inverse multiquadric RBF
class RBF
{
	protected:
		float beta, gamma; // parameters of the inverse multiquadric RBF: psi(x)=(x^2 + gamma^2)^(-beta/2)
		float beta_ov_2, gamma_square;
		
		vector<int> borders; // empty borders to be added to the data volume in case of non-periodic signals
		int size_db[4]; // new volume size after adding borders (db = data+border)
		
		int size_x_pad_Fourier, size_x_pad; // padding of the first dimension, needed for the Fourier transform
		int size_fft[4]; // size of the volume data+border with fft padding in the first dimension
		int nbPixPad; // number of voxels in size_fft
		
		double* FFTWarrayPsi; // array that contains the fft of the RBF
		fftw_plan pF, pB;
		
	public:
		RBF(float _beta, float _gamma, const int* _data_volume_size, vector<int> _borders);
		~RBF();
		
		const int* get_size_fft();
		const vector<int>& get_borders();
		
		void do_fft(double* array, int direction);
		void multiply_with_RBF(double* fft_array_ini, double* fft_array_res);

		float get_beta();
		float get_gamma();
};

#endif
