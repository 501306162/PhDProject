#include "IReSD.h"
//#include <boost/concept_check.hpp>

using namespace std;

IReSD::IReSD(VolumeData* _volumeData, Sdata* _S_data, vector<Phi*>& _listPhi, bool _global_variant, bool _selection_of_connected_regions, int _wait_before_rotation): listPhi(_listPhi)
{
	volumeData = _volumeData;
	S_data = _S_data;
	nbPhi = listPhi.size();
	
	global_variant = _global_variant;
	global_variant_detail.resize(nbPhi);
	selection_of_connected_regions = _selection_of_connected_regions;
	wait_before_rotation = _wait_before_rotation;
	
	volume_size = volumeData->getVolumeSize();
	framePitch = volume_size[0] * volume_size[1] * volume_size[2];
	
	translations_per_contour.resize(nbPhi);
	rotations_per_contour.resize(nbPhi);
}

bool IReSD::get_selection_of_connected_regions_policy()
{
	return selection_of_connected_regions;
}

bool IReSD::get_global_variant_usage()
{
	return global_variant;
}

void IReSD::initialise_iteration()
{
	for (int p=0; p<nbPhi; p++)
	{
		translations_per_contour.clear();
		rotations_per_contour.clear();
	}
}

void IReSD::compute_alignment_of_slice(int image, int numPhi, Vector3& shift_tmp, Vector3& rotation_tmp, int& nb_points_alignment, int& nb_S_phi_diff_signs, int& nb_S_phi_positive, int& nb_Sn_Sm_diff_signs, int& nb_intersections)
{
	int nbFrames = volumeData->getNbFrames();
	Vector3 N, dgdalpha, dgdbeta, dgdgamma;
	int x, y, z, t;
	float S, Sm, phiVal, dist_cont;
	set<int>::const_iterator sit, send;
	bool global = get_global_variant_usage(image, numPhi);
	
	Vector3 rotationAngles = volumeData->get_rotation_angles(image);
	Matrice3x3 matx = Matrice3x3(	1, 0, 0,
					0, cos( rotationAngles[0] ), -sin( rotationAngles[0] ),
					0, sin( rotationAngles[0] ), cos( rotationAngles[0] ) );
	Matrice3x3 maty = Matrice3x3(	cos( rotationAngles[1] ), 0, sin( rotationAngles[1] ),
					0, 1, 0,
					-sin( rotationAngles[1] ), 0, cos( rotationAngles[1] ) );
	Matrice3x3 matz = Matrice3x3(	cos( rotationAngles[2] ), -sin( rotationAngles[2] ), 0,
					sin( rotationAngles[2] ), cos( rotationAngles[2] ), 0,
					0, 0, 1 );
	Matrice3x3 matxP = Matrice3x3(	0, 0, 0,
					0, -sin( rotationAngles[0] ), -cos( rotationAngles[0] ),
					0, cos( rotationAngles[0] ), -sin( rotationAngles[0] ) );
	Matrice3x3 matyP = Matrice3x3(	-sin( rotationAngles[1] ), 0, cos( rotationAngles[1] ),
					0, 0, 0,
					-cos( rotationAngles[1] ), 0, -sin( rotationAngles[1] ) );
	Matrice3x3 matzP = Matrice3x3(	-sin( rotationAngles[2] ), -cos( rotationAngles[2] ), 0,
					cos( rotationAngles[2] ), -sin( rotationAngles[2] ), 0,
					0, 0, 0 );
	int cx = volume_size[0]/2, cy = volume_size[1]/2, cz = volume_size[2]/2;
	Vector3 D = volumeData->get_translation_vector(image);
	
	for (int f=0; f<nbFrames; f++)
	{
		const map<int, float>& S_image = S_data->get_S_image(image, f);
		
		map<int, float>::const_iterator mit(S_image.begin()), mend(S_image.end()), mit2;
		
		for (; mit != mend; ++mit)
		{
			computeCoords(mit->first, volume_size, x, y, z, t);
			N = listPhi[numPhi]->getSpatialNormal(x, y, z, t);
			
			S = mit->second;
			
			if (S == VAL_INVALID || isnan(S))
				throw string("IReSD::compute_alignment_of_slice() Error: invalid S");
			
			if (global)
			{
				phiVal = listPhi[numPhi]->get_phi_value(mit->first);
				
				if (S * phiVal < 0) nb_S_phi_diff_signs++;
				if (S > 0 && phiVal > 0) nb_S_phi_positive++;
			}
			
			/// check position ///
			
			dist_cont = listPhi[numPhi]->getDistToContour(mit->first - t * framePitch, t);
			
			if (global)
			{
				if (dist_cont > 1 && S * phiVal > 0) continue;
			}
			else
			{
				if (dist_cont > 1) continue;
			}
			
			/// translation ///
			
			if (volume_size[2] > 1)
			{
				N = matx * N;
				N = maty * N;
			}
			N = matz * N;
			
			#pragma omp critical
			shift_tmp += N * S;
			
			#pragma omp critical
			nb_points_alignment++;
			
			/// rotation ///
			
			if (!wait_before_rotation)
			{
				if (volume_size[2] > 1)
				{
					dgdalpha = Vector3(x - cx, y - cy, z - cz);
					dgdalpha = matxP * dgdalpha;
					dgdalpha = maty * dgdalpha;
					dgdalpha = matz * dgdalpha;
					dgdalpha += Vector3(cx, cy, cz);
					dgdalpha += D;
					
					#pragma omp critical
					rotation_tmp[0] += S * ( dgdalpha * N );
					
					dgdbeta = Vector3(x - cx, y - cy, z - cz);
					dgdbeta = matx * dgdbeta;
					dgdbeta = matyP * dgdbeta;
					dgdbeta = matz * dgdbeta;
					dgdbeta += Vector3(cx, cy, cz);
					dgdbeta += D;
					
					#pragma omp critical
					rotation_tmp[1] += S * ( dgdbeta * N );
				}
				
				dgdgamma = Vector3(x - cx, y - cy, z - cz);
				if (volume_size[2] > 1)
				{
					dgdgamma = matx * dgdgamma;
					dgdgamma = maty * dgdgamma;
				}
				dgdgamma = matzP * dgdgamma;
				dgdgamma += Vector3(cx, cy, cz);
				dgdgamma += D;
				
				#pragma omp critical
				rotation_tmp[2] += S * ( dgdgamma * N );
			}
			
			/// check intersections ///
			
			const set<int>& all_images = volumeData->get_images_at_voxel(mit->first - t * framePitch);
			sit = all_images.begin();
			send = all_images.end();
			
			for (; sit != send; ++sit)
			{
				if ( same_sequence(image, *sit) ) continue;
				
				mit2 = S_data->get_S_image(*sit, f).find(mit->first);
				if (mit2 == S_data->get_S_image(*sit, f).end()) continue;
				
				Sm = mit2->second;
				
				if (S * Sm < 0) nb_Sn_Sm_diff_signs++;
				nb_intersections++;
			}
		}
	}
}

void IReSD::normalise_amount_of_registration(Vector3& shift_tmp, Vector3& rotation_tmp, float misalignment)
{
	if (!misalignment)
	{
		shift_tmp = VECT_NULL;
		rotation_tmp = VECT_NULL;
		return;
	}
	
	/// translation ///
	
	float limH = 0.6, limL = 0.3;
	float amount_shift;
	
	if (shift_tmp != VECT_NULL)
	{
		float norm_shift = shift_tmp.norm();
		
		if (misalignment >= limH) amount_shift = 1;
		else if (misalignment < limL)
		{
			float amount_shift_tmp = misalignment / limL * 0.1;
			
			amount_shift = min(norm_shift, amount_shift_tmp);
		}
		else amount_shift = (misalignment - limL) / (limH - limL) * (1. - 0.1) + 0.1;
		
		shift_tmp *= amount_shift / norm_shift;
	}
	
	/// rotation ///
	
	if (rotation_tmp == VECT_NULL) return;
	
	limL = 0.2;
	
	float norm_rotation = rotation_tmp.norm();
	float amount_rotation;
	float valM, valP;
	
	if (volume_size[2] == 1)
	{
		valM = 0.00174532925;
		valP = 0.0174532925;
	}
	else
	{
		valM = sqrt(3 * pow(0.00174532925, 2));
		valP = sqrt(3 * pow(0.0174532925, 2));
	}
	
	if (misalignment >= limH) amount_rotation = valP;
	else if (misalignment < limL)
	{
		float amount_rotation_tmp = misalignment / limL * valM;
		
		amount_rotation = min(norm_rotation, amount_rotation_tmp);
	}
	else amount_rotation = (misalignment - limL) / (limH - limL) * (valP - valM) + valM;
	
	rotation_tmp *= amount_rotation / norm_rotation;
}

void IReSD::check_local_variant_condition(int image, int numPhi, int& nb_S_phi_diff_signs, int& nb_S_phi_positive)
{
	if (nb_S_phi_positive > nb_S_phi_diff_signs)
	{
		unset_global_variant_usage(image, numPhi);
		
		map<int, bool>::const_iterator mit(global_variant_detail[numPhi].begin()), mend(global_variant_detail[numPhi].end());
		bool at_least_one = false;
		
		for (; mit != mend; ++mit)
		{
			if (mit->second == true)
			{
				at_least_one = true;
				break;
			}
		}
		
		if (!at_least_one)
		{
			global_variant = false;
			selection_of_connected_regions = false;
		}
	}
}

void IReSD::compute_total_alignment()
{
	map<int, Vector3>::iterator mit_t(translations.begin()), mend_t(translations.end()), mit_r(rotations.begin()), mit;
	
	for (; mit_t != mend_t; ++mit_t, ++mit_r)
	{
		mit_t->second = VECT_NULL;
		mit_r->second = VECT_NULL;
		
		int count_t = 0, count_r = 0;
		
		for (int p=0; p<nbPhi; p++)
		{
			mit = translations_per_contour[p].find(mit_t->first);
			if (mit != translations_per_contour[p].end())
			{
				mit_t->second += mit->second;
				count_t++;
			}
			
			mit = rotations_per_contour[p].find(mit_t->first);
			if (mit != rotations_per_contour[p].end())
			{
				mit_r->second += mit->second;
				count_r++;
			}
		}
		
		if (count_t) mit_t->second /= (float)count_t;
		if (count_r) mit_r->second /= (float)count_r;
	}
	
	// restrict to minimal useful transformation
	
	compute_offsets();
}

void IReSD::compute_offsets()
{
	offset_t = VECT_INVALID;
	offset_r = VECT_INVALID;
	float norm_t = offset_t.norm2(), norm_r = offset_r.norm2(), norm_tmp;
	
	map<int, Vector3>::iterator mit_t(translations.begin()), mend_t(translations.end()), mit_r(rotations.begin());
	
	for (; mit_t != mend_t; ++mit_t, ++mit_r)
	{
		norm_tmp = mit_t->second.norm2();
		if (norm_tmp < norm_t)
		{
			offset_t = mit_t->second;
			norm_t = norm_tmp;
		}
		
		norm_tmp = mit_r->second.norm2();
		if (norm_tmp < norm_r)
		{
			offset_r = mit_r->second;
			norm_r = norm_tmp;
		}
	}
}

void IReSD::save_registration(string filename)
{
	ofstream file(filename.c_str(), ios::out);
	
	if (file)
	{
		const vector< pair<int, int> >& images = volumeData->get_selectedImages();
		int nbImages = images.size();
		
		for (int i=0; i<nbImages; i++)
		{
			file << "image (" << images[i].first << ", " << images[i].second << "): translation = " << volumeData->get_translation_vector(i) << " rotation = " << volumeData->get_rotation_angles(i) << endl;
		}
		file.close();
	}
	else throw string("IReSD::save_registration() Error cannot open file for writing");
}

void IReSD::align_data()
{
	if (wait_before_rotation > 0) wait_before_rotation--;
	
	volumeData->fillVolume();
}


/// IReSD_slice_wise ///

IReSD_slice_wise::IReSD_slice_wise(VolumeData* _volumeData, Sdata* _S_data, vector<Phi*>& _listPhi, bool _global_variant, bool _selection_of_connected_regions, int _wait_before_rotation, bool starts_with_only_shifts_in_SA_planes):
IReSD(_volumeData, _S_data, _listPhi, _global_variant, _selection_of_connected_regions, _wait_before_rotation)
{
	only_shifts_in_SA_planes = starts_with_only_shifts_in_SA_planes;
	
	int nbImages = volumeData->getNbImages();
	
	for (int i=0; i<nbImages; i++)
	{
		translations[i] = VECT_NULL;
		rotations[i] = VECT_NULL;
	}
	
	for (int p=0; p<nbPhi; p++)
	{
		for (int i=0; i<nbImages; i++)
			global_variant_detail[p][i] = global_variant;
	}
	
	inhibited_dir.resize(nbImages);
	for (int i=0; i<nbImages; i++)
	{
		if (only_shifts_in_SA_planes)
		{
			//normal vector to the image
			Vector3 normal_imageRef(0, 0, 1);
			Vector3 normalVector = volumeData->imageToVolumeVect(i, normal_imageRef);
			
			//select directions
			float a_ = fabs(normalVector[0]);
			float b_ = fabs(normalVector[1]);
			float c_ = fabs(normalVector[2]);
			
			float max_val = max(a_, max(b_, c_) );
			float seuil = 90./100. * max_val;
			
			inhibited_dir[i] = Point3<bool>(a_ > seuil, b_ > seuil, c_ > seuil);
		}
		else
			inhibited_dir[i] = Point3<bool>(false, false, false);
	}
}

void IReSD_slice_wise::compute_alignment_from_contour(int numPhi)
{
	int nb_Sn_Sm_diff_signs, nb_intersections, nb_S_phi_diff_signs, nb_S_phi_positive, nb_points_alignment;
	float misalignment;
	
	Point3<bool> no_dir(false, false, false);
	
	int wait_save = wait_before_rotation;
	
	int nbImages = volumeData->getNbImages();
	
	for (int image = 0; image < nbImages; image++)
	{
		translations_per_contour[numPhi][image] = VECT_NULL;
		rotations_per_contour[numPhi][image] = VECT_NULL;
		nb_S_phi_diff_signs = 0;
		nb_S_phi_positive = 0;
		nb_points_alignment = 0;
		nb_Sn_Sm_diff_signs = 0;
		nb_intersections = 0;
		
		//only rotate when performing full translation and rotation
		if (!wait_before_rotation && inhibited_dir[image] != no_dir) wait_before_rotation = 1;
		
		compute_alignment_of_slice(image, numPhi, translations_per_contour[numPhi][image], rotations_per_contour[numPhi][image], nb_points_alignment, nb_S_phi_diff_signs, nb_S_phi_positive, nb_Sn_Sm_diff_signs, nb_intersections);
		
		if (nb_points_alignment)
		{
			translations_per_contour[numPhi][image] /= (float)nb_points_alignment;
			rotations_per_contour[numPhi][image] /= (float)nb_points_alignment;
			
			if (inhibited_dir[image] != no_dir)
			{
				if ( inhibited_dir[image].x() ) translations_per_contour[numPhi][image][0] = 0;
				if ( inhibited_dir[image].y() ) translations_per_contour[numPhi][image][1] = 0;
				if ( inhibited_dir[image].z() ) translations_per_contour[numPhi][image][2] = 0;
			}
			
			if (nb_intersections) misalignment = nb_Sn_Sm_diff_signs / (float)nb_intersections;
			else misalignment = 0;
			cout << "image " << image << ": " << nb_Sn_Sm_diff_signs << " " << nb_intersections << " estimated misalignment = " << misalignment << endl;
			
			normalise_amount_of_registration(translations_per_contour[numPhi][image], rotations_per_contour[numPhi][image], misalignment);
			
			if (global_variant_detail[numPhi][image])
				check_local_variant_condition(image, numPhi, nb_S_phi_diff_signs, nb_S_phi_positive);
		}
		
		wait_before_rotation = wait_save;
	}
}

bool IReSD_slice_wise::get_global_variant_usage(int image, int numPhi)
{
	return global_variant_detail[numPhi][image];
}

void IReSD_slice_wise::unset_global_variant_usage(int image, int numPhi)
{
	global_variant_detail[numPhi][image] = false;
}

bool IReSD_slice_wise::same_sequence(int image_n, int image_m)
{
	if (image_m == image_n) return true;
	
	return false;
}

void IReSD_slice_wise::align_data()
{
	map<int, Vector3>::iterator mit_t(translations.begin()), mend_t(translations.end()), mit_r(rotations.begin()), mit;
	
	Vector3 translation, rotation;
	
	for (; mit_t != mend_t; ++mit_t, ++mit_r)
	{
		translation = mit_t->second - offset_t;
		rotation = mit_r->second - offset_r;
		
		if (translation == VECT_NULL && rotation == VECT_NULL) continue;
		
		cout << "temporal sequence " << mit_t->first << ": translation of " << translation << " and rotation of " << rotation << endl;
		volumeData->move_temporal_sequence(mit_t->first, translation, rotation);
		cout << "new alignment vectors : " << volumeData->get_translation_vector(mit_t->first) << " " << volumeData->get_rotation_angles(mit_t->first) << endl;
	}
	
	IReSD::align_data();
}


/// IReSD_stack_wise ///

IReSD_stack_wise::IReSD_stack_wise(VolumeData* _volumeData, Sdata* _S_data, vector<Phi*>& _listPhi, bool _global_variant, bool _selection_of_connected_regions, int _wait_before_rotation):
IReSD(_volumeData, _S_data, _listPhi, _global_variant, _selection_of_connected_regions, _wait_before_rotation)
{
	const set<int>& series = volumeData->get_selectedSeries();
	set<int>::const_iterator sit, send(series.end());
	
	for (sit = series.begin(); sit != send; ++sit)
	{
		translations[*sit] = VECT_NULL;
		rotations[*sit] = VECT_NULL;
	}
	
	for (int p=0; p<nbPhi; p++)
	{
		for (sit = series.begin(); sit != send; ++sit)
			global_variant_detail[p][*sit] = global_variant;
	}
}

void IReSD_stack_wise::compute_alignment_from_contour(int numPhi)
{
	int nb_Sn_Sm_diff_signs, nb_intersections, nb_S_phi_diff_signs, nb_S_phi_positive, nb_points_alignment;
	int image;
	float misalignment;
	
	const set<int>& series = volumeData->get_selectedSeries();
	set<int>::const_iterator sit(series.begin()), send(series.end());
	
	for (; sit != send; ++sit)
	{
		translations_per_contour[numPhi][*sit] = VECT_NULL;
		rotations_per_contour[numPhi][*sit] = VECT_NULL;
		nb_S_phi_diff_signs = 0;
		nb_S_phi_positive = 0;
		nb_points_alignment = 0;
		nb_Sn_Sm_diff_signs = 0;
		nb_intersections = 0;
		
		int nbSlices = volumeData->getNbSlices(*sit);
		
		for (int slice = 0; slice < nbSlices; slice++)
		{
			image = volumeData->getImageNumber(*sit, slice);
			
			compute_alignment_of_slice(image, numPhi, translations_per_contour[numPhi][*sit], rotations_per_contour[numPhi][*sit], nb_points_alignment, nb_S_phi_diff_signs, nb_S_phi_positive, nb_Sn_Sm_diff_signs, nb_intersections);
		}
		
		if (nb_points_alignment)
		{
			translations_per_contour[numPhi][*sit] /= (float)nb_points_alignment;
			rotations_per_contour[numPhi][*sit] /= (float)nb_points_alignment;
			
			cout << nb_Sn_Sm_diff_signs << " " << nb_intersections << " ";
			if (nb_intersections) misalignment = nb_Sn_Sm_diff_signs / (float)nb_intersections;
			else misalignment = 0;
			cout << "serie " << *sit << ": estimated misalignment = " << misalignment << endl;
			
			normalise_amount_of_registration(translations_per_contour[numPhi][*sit], rotations_per_contour[numPhi][*sit], misalignment);
			
			if (global_variant_detail[numPhi][*sit])
				check_local_variant_condition(*sit, numPhi, nb_S_phi_diff_signs, nb_S_phi_positive);
		}
	}
	
	//cv::waitKey();
}

bool IReSD_stack_wise::get_global_variant_usage(int image, int numPhi)
{
	return global_variant_detail[numPhi][ volumeData->get_selectedImages()[image].first ];
}

void IReSD_stack_wise::unset_global_variant_usage(int image, int numPhi)
{
	global_variant_detail[numPhi][ volumeData->get_selectedImages()[image].first ] = false;
}

bool IReSD_stack_wise::same_sequence(int image_n, int image_m)
{
	if (image_m == image_n) return true;
	
	if (volumeData->get_selectedImages()[image_n].first == volumeData->get_selectedImages()[image_m].first)
		return true;
	
	return false;
}

void IReSD_stack_wise::align_data()
{
	map<int, Vector3>::iterator mit_t(translations.begin()), mend_t(translations.end()), mit_r(rotations.begin()), mit;
	
	for (; mit_t != mend_t; ++mit_t, ++mit_r)
	{
		cout << "spatial sequence " << mit_t->first << ": translation of " << mit_t->second - offset_t << " and rotation of " << mit_r->second - offset_r << endl;
		volumeData->move_spatial_sequence(mit_t->first, mit_t->second - offset_t, mit_r->second - offset_r);
	}
	
	IReSD::align_data();
}
