#ifndef PHI_H
#define PHI_H

#include <iostream>
#include <fstream>
#include <list>
#include <set>
#include <vector>
#include <map>

#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <nifti1.h>

#include <opencv/cv.h>
#include <opencv/highgui.h>

#ifndef NO_OMP
	#include <omp.h>
#endif

#include "AuxClasses.h"

#define NB_THREADS omp_get_num_procs()

#define GCC_VERSION (__GNUC__ * 10000 \
                               + __GNUC_MINOR__ * 100 \
                               + __GNUC_PATCHLEVEL__)

using namespace std;

class Phi
{
	protected:
		const int* size;
		int nbPixVol, nbPix;
		
		float* data;

		int contourSize;
		map<int, float>* contour; //time, (indice in 3D volume, dist to contour)
		
		void display_phi();
		
	public:
		Phi(const int* s);
		~Phi();
		
		void initialiseSpheres(const list< Point4<int> >& centreIni, const list<float>& radiusIni);
		void initialise_from_nii(string nii_file);
		
		const int* get_volume_size();
		
		void set_phi_value(int ind_phi, float newVal);

		float get_phi_value(Point4<float> point) const;
		float get_phi_value(float x, float y, float z, float t) const;
		float get_phi_value(Point4<int> point) const;
		float get_phi_value(int x, int y, int z, int t) const;
		float get_phi_value(int ind_phi) const;
		
		void computeContour(float epsilon);
		
		const map<int, float>* getContour();
		int getContourSize();
		float getDistToContour(int voxel, int time);
		bool isContourPoint(int voxel_3d, int time);
		
		void savePhi(int numPhi, int t, string str="");
		void saveContour(int numPhi, int t, string str="");
		
		Vector3 getSpatialNormal(int x, int y, int z, int f);
		Vector3 getSpatialNormal(float x, float y, float z, float t);
		float computeNormGradientInNormalDirection(int x, int y, int z, int t);
};

#endif
