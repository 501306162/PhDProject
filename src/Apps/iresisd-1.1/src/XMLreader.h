#ifndef XMLREADER_H
#define XMLREADER_H

#include <iostream>
#include <list>
#include <vector>

#include "tinyxml2.h"
#include "AuxClasses.h"

using namespace std;
using namespace tinyxml2;

class XMLreader
{
	protected:
		XMLDocument doc;
		
	public:
		XMLreader(string xml_file);
		
		/// methods for VolumeData ///
		
		string get_dataType();
		string get_repertory();
		bool is_motion_periodic();
		bool accurate_volume_filling();
		bool size_volume_known();
		Vector3 get_volume_size();
		Vector3 get_ROI_offset();
		bool serie_ref_defines_volume_size();
		bool ROI_selection();
		bool use_first_timeframe_only();
		Vector3 get_registration_translation(string serieDescription, int slice);
		Vector3 get_registration_rotation(string serieDescription, int slice);
		
		/// methods for RBF ///
		
		int get_nb_RBFs();
		float get_beta();
		float get_gamma(int rbf_index);
		vector<int> get_volume_borders();
		
		/// methods for Phi ///
		
		int get_nb_phi();
		string get_phi_ini_type();
		string get_phi_ini_nifti_file(int phi_index);
		void get_initial_circles(vector< list< Point4<int> > >& centre_ini, vector< list<float> >& radius_ini);
		
		/// methods for Sdata ///
		
		int get_nb_SFdata();
		void get_SFdata_params(int SF_index, string& type_sf, vector<float>& params_sf);
		float get_SFdata_weight(int SF_index);
				
		/// methods for Sgeom ///
		
		int get_nb_SFgeom();
		void get_SFgeom_params(int SF_index, string& type_sf, vector<float>& params_sf);
		float get_SFgeom_weight(int SF_index);
		
		/// methods for ISISD ///
		
		string get_schemeType();
		vector<float> get_dirac_widths();
		
		/// methods for IReSD ///
		
		bool registration_required();
		
		/// methods for Interface ///
		
		bool use_global_registration();
		bool do_selection_of_connected_regions();
		int do_slice_wise_registration();
		bool shifts_in_SA_planes_at_beginning();
		int wait_before_rotation();
		
		int allow_display();
		int get_nb_iterations(int step);
		int get_max_stable_iterations(int step);
		
};

#endif
