#ifndef ISISD_H
#define ISISD_H

#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <fftw3.h>
#include <list>

#include "AuxClasses.h"
#include "Sdata.h"
#include "Sgeom.h"
#include "RBF.h"
#include "Phi.h"

class ISISD
{
	protected:
		Sdata* S_data;
		Sgeom* S_geom;
		
		VolumeData* volumeData;
		
		vector<Phi*>& listPhi;
		
		const int* size;
		int nbPix, nbPixVol;
		
		const int* size_fft;
		int nbPixFFT;
		const vector<int>& borders;
		
		float max_dirac_width;
		double* speeds_data;
		float* speeds_geom;
		
		list<float>* contourLengths; // for each phi
		float* maxTranslation; // for each phi
		float oscillationsNumber;
		bool* contourStoppedGrowing; // for each phi
		bool dynamic_maxTranslation; // usually true, but it may be useful to set it to false during the second step of the two-steps scheme
		
		float dt_data, dt_geom;
		
		int display;
		
		
		/// methods ///
		
		virtual void compute_dphi_data(int numPhi) = 0;
		bool compute_dt_data(int numPhi);
		void compute_dt_geom(int numPhi);
		
		void checkContourStability(int numPhi);
		
		void display_S();
		
	public:
		ISISD(Sdata* _S_data, Sgeom* _S_geom, vector<Phi*>& _listPhi, const vector<int>& _borders, VolumeData* _volumeData, int _display);
		virtual ~ISISD();
		
		bool updatePhi(int numPhi);
		bool auto_reinitialisation(int numPhi, bool global_variant=false, bool select_connected_regions=false);
		void set_dynamic_maxTranslation(bool newValue);
		bool doesContourGrow(int numPhi);
};

class ISISD_simpleScheme: public ISISD
{
	protected:
		RBF* psi;
		float dirac_width;
		
	public:
		ISISD_simpleScheme(Sdata* _S_data, Sgeom* _S_geom, vector<Phi*>& _listPhi, RBF* _psi, float _dirac_width, VolumeData* _volumeData, int _display);
		
		void compute_dphi_data(int numPhi);
		void change_psi(RBF* _psi, float _dirac_width);
};

class ISISD_mixedScheme: public ISISD
{
	protected:
		vector<RBF*>& listPsi;
		int nbPsi;
		vector<float>& dirac_widths;
		bool unique_width;
		
		unsigned int* domains_psi; // map of the domains of the different RBFs
		
		double *dalpha, *dphi;
		map<int, pair<float, float> > memory_S; //indfft, (S, dist_cont)
		
		
		/// methods ///
		
		void display_dphi_tmp();
		
	public:
		ISISD_mixedScheme(Sdata* _S_data, Sgeom* _S_geom, vector<Phi*>& _listPhi, vector<RBF*>& _listPsi, vector<float>& _dirac_widths, VolumeData* _volumeData, int _display);
		~ISISD_mixedScheme();
		
		void compute_dphi_data(int numPhi);
		
		void display_psi_domains();
};

#endif
