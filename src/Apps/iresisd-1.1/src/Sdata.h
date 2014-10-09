#ifndef SDATA_H
#define SDATA_H

#include <vector>
#include <map>
#include <set>

#include <fftw3.h>

#include "AuxClasses.h"
#include "XMLreader.h"
#include "VolumeData.h"
#include "Phi.h"
#include "SpeedFunction.h"

class Sdata
{
	protected:
		XMLreader* xml_reader;
		
		VolumeData* volumeData;
		vector<Phi*>& listPhi;

		vector<SpeedFunction_Data*> sf_pointers; //index of sf, sf
		vector<string> sf_type;
		map<int, vector<float> > sf_params; //index of sf, params of sf
		map<int, vector< pair<int, float> > > sf_of_serie; // for each serie, list of SF index with associated weight
		
		const int* volume_size;
		int nbPixVol;
		const int* size_fft;
		int nbPixFFT;
		const vector<int>& borders;
		int nbFrames;
		int nbImages;
		int nbPixelsTot;

		/// temporary arrays, used for the computation of S^n
		float* S_image_pixel; //pixel, value of S^n
		set< Point4<int> > waitingList; //list of voxels that have a value for S (to be interpolated from S_image_pixel)
		float S_max;
		set<int>** list_pixels; //images, frame, list of pixels

		/// storage of S
		map<int, float>** S_image; //image, frame, voxel, value of S^n for the level set function phi currently in use
		double* fft_array; //value of S for the level set function phi currently in use
		
		/// methods ///
		
		void query_speed_functions();
		float query_sf_param(string param_name);
		void query_speed_functions_of_serie(int serie);
		void create_SF(int sf_index); //to be updated when new SP are defined
		
		void get_pixels_on_contour(int numPhi);
		void compute_S_on_contour_in_image(int image, int frame, int numPhi);
		void compute_S_in_image_local_variant(int image, int frame, int numPhi); //compute S^n on contour: used for local variant of registration method, or no registration
		void compute_S_in_image_global_variant(int image, int frame, int numPhi, bool do_selection_of_connected_regions); //compute S^n on contour + Omega_diff: used for global variant of registration method
		void compute_S_in_volume(int image, int frame, int numPhi, bool global_variant);
		
		void select_connected_region(map<int, set< Point4<int> > >& positiveSpeeds, set<int>& checked, map<int, set< Point4<int> > >::const_iterator mit, int Nc, int nbPixImg);
		
		void display_S_image(int image);
		
	public:
		Sdata(XMLreader* _xml_reader, VolumeData* _volumeData, vector<Phi*>& _listPhi, const int* _size_fft, const vector<int>& _borders);
		~Sdata();
		
		void initialise_iteration();
		void compute_S(int numPhi, bool global_variant=false, bool do_selection_of_connected_regions=false);
		
		double* get_S();
		const map<int, float>& get_S_image(int image, int timeframe);
		
		void compute_S_everywhere(int numPhi, bool global_variant=false, bool do_selection_of_connected_regions=false);
};

#endif

