#ifndef IRESD_H
#define IRESD_H

#include "AuxClasses.h"
#include "Sdata.h"
#include "VolumeData.h"
#include "Phi.h"

class IReSD
{
	protected:
		VolumeData* volumeData;
		Sdata* S_data;
		vector<Phi*> listPhi;
		int nbPhi;
		
		const int* volume_size;
		int framePitch;
		
		bool global_variant, selection_of_connected_regions;
		vector< map<int, bool> > global_variant_detail; //contour, sequence (image or serie), global
		int wait_before_rotation; //number of iterations before rotation is performed. if -1: never do rotation
		
		map<int, Vector3> translations, rotations; //sequence (image or serie), amount of registration
		vector< map<int, Vector3> > translations_per_contour, rotations_per_contour; //contour, sequence (image or serie)
		Vector3 offset_t, offset_r;
		
		/// methods ///
		
		void compute_alignment_of_slice(int image, int numPhi, Vector3& shift_tmp, Vector3& rotation_tmp, int& nb_points_alignment, int& nb_S_phi_diff_signs, int& nb_S_phi_positive, int& nb_Sn_Sm_diff_signs, int& nb_intersections);
		void normalise_amount_of_registration(Vector3& shift_tmp, Vector3& rotation_tmp, float misalignment);
		void check_local_variant_condition(int image, int numPhi, int& nb_S_phi_diff_signs, int& nb_S_phi_positive);
		void compute_offsets();
		
		virtual bool same_sequence(int image_n, int image_m) = 0;
		virtual bool get_global_variant_usage(int image, int numPhi) = 0;
		virtual void unset_global_variant_usage(int image, int numPhi) = 0;
	
	public:
		IReSD(VolumeData* _volumeData, Sdata* _S_data, vector<Phi*>& _listPhi, bool _global_variant, bool _selection_of_connected_regions, int _wait_before_rotation);
		
		bool get_selection_of_connected_regions_policy();
		bool get_global_variant_usage();
		
		void initialise_iteration();
		virtual void compute_alignment_from_contour(int numPhi) = 0;
		void compute_total_alignment();
		virtual void align_data();
		
		void save_registration(string filename);
};

class IReSD_slice_wise: public IReSD
{
	protected:
		bool only_shifts_in_SA_planes;
		vector< Point3<bool> > inhibited_dir; //image
		
		bool same_sequence(int image_n, int image_m);
		bool get_global_variant_usage(int image, int numPhi);
		void unset_global_variant_usage(int image, int numPhi);
		
	public:
		IReSD_slice_wise(VolumeData* _volumeData, Sdata* _S_data, vector<Phi*>& _listPhi, bool _global_variant, bool _selection_of_connected_regions, int _wait_before_rotation, bool starts_with_only_shifts_in_SA_planes);
		
		void compute_alignment_from_contour(int numPhi);
		void align_data();
};

class IReSD_stack_wise: public IReSD
{
	protected:
		bool same_sequence(int image_n, int image_m);
		bool get_global_variant_usage(int image, int numPhi);
		void unset_global_variant_usage(int image, int numPhi);
		
	public:
		IReSD_stack_wise(VolumeData* _volumeData, Sdata* _S_data, vector<Phi*>& _listPhi, bool _global_variant, bool _selection_of_connected_regions, int _wait_before_rotation);
		
		void compute_alignment_from_contour(int numPhi);
		void align_data();
};

#endif

