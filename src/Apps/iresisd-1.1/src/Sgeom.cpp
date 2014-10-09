#include "Sgeom.h"

using namespace std;

Sgeom::Sgeom(XMLreader* _xml_reader, vector<Phi*>& _listPhi, const int* vol_size, bool _periodic_motion): listPhi(_listPhi)
{
	xml_reader = _xml_reader;
	
	volume_size = vol_size;
	nbPixVol = volume_size[0] * volume_size[1] * volume_size[2];
	
	periodic_motion = _periodic_motion;
	
	/// create the SpeedFunctions
	
	int nbSF = xml_reader->get_nb_SFgeom();
	sf_pointers.resize(nbSF);
	sf_type.resize(nbSF);
	sf_weight.resize(nbSF);
	
	for (int sf=0; sf<nbSF; sf++)
	{
		xml_reader->get_SFgeom_params(sf, sf_type[sf], sf_params[sf]);
		sf_weight[sf] = xml_reader->get_SFgeom_weight(sf);
		
		create_SF(sf);
	}
	
	/// create the array to store S
	
	if ( !sf_pointers.empty() )
		S = new float[nbPixVol * volume_size[3]];
	else S = NULL;
}

Sgeom::~Sgeom()
{
	if (S != NULL) delete[] S;
	
	vector<SpeedFunction_Geom*>::iterator vit(sf_pointers.begin()), vend(sf_pointers.end());
	for (; vit != vend; ++vit) delete *vit;
}

void Sgeom::create_SF(int sf)
{
	if ( !sf_type[sf].compare("normalisation") ) sf_pointers[sf] = new SpeedFunction_normalisation(listPhi, volume_size, periodic_motion);
	else throw string("Sgeom() Error: unknown SpeedFunction type");
}

float* Sgeom::get_S()
{
	return S;
}

void Sgeom::initialise_iteration()
{
	vector<SpeedFunction_Geom*>::const_iterator vit(sf_pointers.begin()), vend(sf_pointers.end());
	for (; vit != vend; ++vit) (*vit)->prepare_iteration();
}

void Sgeom::compute_S(int numPhi)
{
	if (S == NULL) return;
	
	/// compute the speeds dphi/dt produced by the geometric SpeedFunctions
	
	int indV, sf;
	int nbSF = sf_pointers.size();
	float speed, sum_speeds, count_speeds;
	
	for (int t=0; t<volume_size[3]; t++)
	#pragma omp parallel for default(shared) private(indV, sf, speed, sum_speeds, count_speeds) schedule(dynamic) num_threads(NB_THREADS)
	for (indV=0; indV<nbPixVol; indV++)
	{
		sum_speeds = 0;
		count_speeds = 0;
		
		for (sf=0; sf<nbSF; sf++)
		{
			speed = sf_pointers[sf]->get_speed(indV, t, numPhi);
	
			if (speed == VAL_INVALID) continue;

			sum_speeds += speed * sf_weight[sf];
			count_speeds += sf_weight[sf];
		}
		
		if (!count_speeds)
		{
			#pragma omp critical
			S[indV + t * nbPixVol] = 0;
		}
		else
		{
			#pragma omp critical
			S[indV + t * nbPixVol] = sum_speeds / count_speeds;
		}
	}
}
