#ifndef SPEEDFUNCTION_H
#define SPEEDFUNCTION_H

#include <vector>

#include "Phi.h"
#include "VolumeData.h"

class SpeedFunction
{
	protected:
		vector<Phi*>& listPhi;

	public:
		SpeedFunction(vector<Phi*>& _listPhi);
		virtual ~SpeedFunction();
		
		virtual void prepare_iteration() = 0;
};


/// SpeedFunction_Data ///

class SpeedFunction_Data: public SpeedFunction
{
	protected:
		VolumeData* volumeData;
		set<int> list_series; //series for which this SF is required
		
		map<int, vector<bool> > image_can_be_used; // image, timeframe
		
	public:
		SpeedFunction_Data(VolumeData* _volumeData, vector<Phi*>& _listPhi);
		virtual void add_serie(int serie);
		
		bool is_image_usable(int image, int timeframe);
		virtual float get_speed(int image, int timeframe, int pixel, int numPhi) = 0;
};

class SpeedFunction_CV: public SpeedFunction_Data
{
	protected:
		int nb_phases;
		vector< map<int, float> > c; // phase, sequence (image or serie), mean intensity
		vector< map<int, int> > count;
		bool slice_wise;
		
		int compute_sequence_number(int serie, int slice);

		void update_c(int image, int timeframe, int sequence);
		
	public:
		SpeedFunction_CV(VolumeData* _volumeData, vector<Phi*>& _listPhi);
		
		void prepare_iteration();
		float get_speed(int image, int timeframe, int pixel, int numPhi);
};

class SpeedFunction_GVF: public SpeedFunction_Data
{
	protected:
		float baloonForce;
		
		map<int, vector< vector<float> > > gValues; //image, timeframe, indice c + Columns * r
		map<int, vector< vector<Vector3> > > gradgValues; //image, timeframe, indice c + Columns * r
		
		void initialise_g();
		float compute_g(float normGrad);
		
	public:
		SpeedFunction_GVF(VolumeData* _volumeData, vector<Phi*>& _listPhi, float _baloonForce);
		void add_serie(int serie);
		
		void prepare_iteration();
		float get_speed(int image, int timeframe, int pixel, int numPhi);
};


/// SpeedFunction_Geom ///

class SpeedFunction_Geom: public SpeedFunction
{
	protected:
		const int* size;
		int nbPixVol;
		bool periodic_motion;
		vector<int> pitch;
		
		/// methods ///
		
		Vector4 computeForwardDiff(int voxel, int time, int numPhi);
		Vector4 computeBackwardDiff(int voxel, int time, int numPhi);
		Vector4 computeForwardBackwardDiff(int voxel, int time, int numPhi);
		Vector4 computeForwardForwardDiff(int voxel, int time, int numPhi);
		Vector4 computeBackwardBackwardDiff(int voxel, int time, int numPhi);
		double computeNormGradientEntropy2ndOrder(int voxel, int time, float speed, int numPhi);
		double m(double x, double y);
	
	public:
		SpeedFunction_Geom(vector<Phi*>& _listPhi, const int* _size, bool _periodic_motion);
		
		virtual float get_speed(int voxel, int time, int numPhi) = 0;
		virtual float get_dt() = 0;
};

class SpeedFunction_normalisation: public SpeedFunction_Geom
{
	protected:
		set<int> border_values;
		float max_speed;
		
	public:
		SpeedFunction_normalisation(vector<Phi*>& _listPhi, const int* _size, bool _periodic_motion);
		
		void prepare_iteration();
		float get_speed(int voxel, int time, int numPhi);
		float get_dt();
};

#endif
