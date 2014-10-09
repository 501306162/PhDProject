#ifndef INTERFACE_H
#define INTERFACE_H

#include <iostream>
#include <algorithm>

#include <opencv/cv.h>
#include <opencv/highgui.h>

#include <map>

#include "XMLreader.h"
#include "VolumeData.h"
#include "Phi.h"
#include "RBF.h"
#include "Sdata.h"
#include "Sgeom.h"
#include "ISISD.h"
#include "IReSD.h"

class Interface
{
	protected:
		/// xml parameters
		XMLreader* xml_reader;
		
		/// data
		VolumeData* volumeData;
		
		/// phi
		vector<Phi*> listPhi;
		int nbPhi;
		
		/// RBF
		vector<RBF*> listPsi;
		int currentStep; //used if double scheme
		
		/// speeds
		Sdata* S_data;
		Sgeom* S_geom;
		
		/// ISISD
		ISISD* isisd;
		string scheme_type;
		int nb_steps;
		vector<float> dirac_widths;
		
		///IReSISD
		IReSD* iresd;
		bool doRegistration;
		
		/// stopping criterium: choose between one of these (from xml)
		vector<int> nb_iterations; // step (if any), total number of iterations
		vector<int> nb_stable_iterations_before_stop;  // step (if any), number of stable iterations


		/// methods ///
		
		void loadData();
		void create_psi();
		void create_phi();
		void create_isisd();
		
		void queryInitialisationPhi(vector< list< Point4<int> > >& centre_ini, vector< list<float> >& radius_ini);
		
	public:
		Interface(XMLreader* xml_reader);
		~Interface();

		void loop();
		
		void displayContour(bool showImage, bool writeImage, int time, string str);
		void displayContour(bool showImage, bool writeImage, int iamge, int frame, string str);
		void displayArea(bool showImage, bool writeImage, int time, string str);
		void displayArea(bool showImage, bool writeImage, int image, int frame, string str);
};

typedef struct _mousse_args_circle
{
	bool actif;
	bool validated;
	bool cancelled;
	bool drawing;
	cv::Point centre, perimetre;
	cv::Mat* img;
}mousse_args_circle;

void mouse_callback_initialisation( int event, int x, int y, int flags, void* param );

#endif
