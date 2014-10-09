#ifndef SGEOM_H
#define SGEOM_H

#include <vector>
#include <map>
#include <set>

#include <fftw3.h>

#include "XMLreader.h"
#include "Phi.h"
#include "SpeedFunction.h"

class Sgeom
{
	protected:
		XMLreader* xml_reader;
		vector<Phi*>& listPhi;
		
		vector<SpeedFunction_Geom*> sf_pointers;
		vector<string> sf_type;
		vector<float> sf_weight;
		map<int, vector<float> > sf_params; //index of sf, params of sf
		
		const int* volume_size;
		int nbPixVol;
		bool periodic_motion;
		
		float* S;
		
		void create_SF(int sf);
		
	public:
		Sgeom(XMLreader* _xml_reader, vector<Phi*>& listPhi, const int* vol_size, bool _periodic_motion);
		~Sgeom();
		
		void initialise_iteration();
		void compute_S(int numPhi);
		
		float* get_S();
};

#endif
