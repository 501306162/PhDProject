#ifndef VOLUMEDATA_H
#define VOLUMEDATA_H

#define HAVE_CONFIG_H

#include <vector>
#include <list>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <opencv/cv.h>
#include <opencv/highgui.h>

# undef DT_UNKNOWN
#include <dirent.h>
#include <math.h>
#include <cmath>

#include <dcmtk/config/osconfig.h>
#include <dcmtk/dcmdata/dctk.h>
#include <nifti1.h>
#include <teem/nrrd.h>

#ifndef NO_OMP
	#include <omp.h>
#endif

#include "AuxClasses.h"
#include "XMLreader.h"

#define NB_THREADS omp_get_num_procs()

#define GCC_VERSION (__GNUC__ * 10000 \
                               + __GNUC_MINOR__ * 100 \
                               + __GNUC_PATCHLEVEL__)

#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352

using namespace std;


class VolumeData
{
	protected:
		XMLreader& xml_reader;
	
		string repertory;
		
		int nbFrames;
		vector<string> seriesDescription;
		vector<vector<string> > seriesNumber;

		
		
		vector< list<string> > seriesImagesFilenames;
		vector< pair<int, int> > seriesImageSizes; //serie, <rows, columns>
		vector< vector<double> > seriesImageOrientation; //OrientationPatient
		vector< pair<double, double> > seriesPixelSpacing; //serie, <row spacing, col spacing>
		vector<int> seriesNbSlices; // serie, nbSlices
		map<int, vector< Vector3 > > seriesImagePosition; //serie, slice, xyz - PositionPatient
		map<int, vector< vector< vector<float> > > > seriesDataValues; //serie, slice, timeframe, col + Columns * row
		
		map<int, set<double> > listDistances; //serie, list of distances of the slices to the origine of the coordinates system along the normal axis
		map<int, vector< vector< string> > > seriesImagesFilenamesSorted; //serie, slice, timeframe, name
		
		set<int> selectedSeries;
		int serie_ref;
		vector< pair<int, int> > selectedImages;
		int nbImages;
		map<int, int> serieFirstImageIndex;
		
		bool periodic_motion;
		
		Vector3 offset_serieRef;
		Matrice4x4 mat_serieRef, mat_serieRef_inv;
		map< pair<int, int>, Matrice4x4> mat_img_inv; //image
		map< pair<int, int>, Point3<float>* > seriesDataVolumeCoords; //image, pixel: c + Nc * r. Contains the volume coordinates of the pixels without any registration applied
		
		float** timeOfFrame; //image, frame, time of frame in 4D volume
		
		float xMin, xMax, yMin, yMax, zMin, zMax; //limits of the selected volume
		Vector3 offset_ROI;
		int size[4];
		int nbPix, nbPixVol;
		float** volume; //timeframe, index in volume.   VAL_INVALID = no data
		set<int>* imagePresence; //index in volume, list of images
		
		Vector3* translationVect; //image, translation vector
		Vector3* rotationAngles; //image, rotation angles around the (Ox, Oy, Oz) axes (in this order)
		map<int, Point3<float> >* pixelCoordsInVolume; //image, pixel for which the aligned coord has already been computed (i.e. after registration has been applied)
		map<int, Point3<float> >* voxelCoordsInImageVolume; //image, voxel for which the coord in image volume has already been computed
		
		
		/// methods ///
		
		void listFilenames(DIR *rep, list<string>& filenames, string path);
		virtual void listSeries() = 0;
		virtual void printListedSeries();
		void queryChoiceSeries();
		virtual void selectSeries();
		virtual void printSelectedSeries();
		virtual void readSerie(int serieNum) = 0;
		virtual void choose_serie_ref();
		void check_serie(set<int>::iterator& it_serie);
		void handle_similar_series(set<int>::iterator& it_serie, set<int>::iterator& it_serie_add);
		void fuse_series(set<int>::iterator& it_serie, set<int>::iterator& it_serie_add);
		virtual void add_slices_to_serie(int serieNumber, set< pair<double, pair<int, int> > >& dists);
		void computeVolumeCoords();
		virtual void buildListTimeframes();
		void buildVolume();
		void ROIselection();
		void loadRegistration();
		void replace_slice(int serieNumber, int slice, int add_serie, int add_slice, double add_distance);
		void fillVolume_closestPoint(); //quick version
		void fillVolume_interpolation(); //slow but more accurate version
		
		
	public:
		VolumeData(XMLreader& _xml_reader);
		~VolumeData();
		void loadData();
		
		virtual bool comprise_temporal_sequences();

		const set<int>& get_selectedSeries();
		const vector< pair<int, int> >& get_selectedImages();
		int getNbFrames();
		int getNbImages();
		int getNbSlices(int serie);
		std::string getSeriesNumber(int series, int slice)
		{
			return seriesNumber[series][slice];
		}
		int getImageNumber(int serie, int slice);
		
		const pair<int, int>& getImageSize(int serie);
		float getTimeOfFrame(int image, int timeframe);
		float getPixelValue(int serie, int slice, int frame, int pixel);
		float getPixelValue(int image, int frame, int pixel);
		float getVoxelValue(int voxel, int timeframe);
		
		const int* getVolumeSize();
		bool is_motion_periodic();
		
		const set<int>& get_images_at_voxel(int voxel);
		void getClosestFrames(int t, int image, int& frameM, float& distM, int& frameP, float& distP);
		
		Point3<float> getImagePointCoordsInVolume(int image, int pixel);
		Point3<float> getVoxelCoordsInImageVolume(int voxel, int image);
		
		bool isPointInVolume(Point3<int> point);
		bool isPointInVolume(Point3<int> point, int t);
		bool isPointInVolume(Point3<float> point);
		bool isImagePointInVolume(int image, int indiceImg);
		
		Vector3 imageToVolumeVect(int image, const Vector3& vectIni);
		
		void move_spatial_sequence(int serie, Vector3 d_translation=NULL, Vector3 d_rotation=NULL);
		void move_temporal_sequence(int image, Vector3 d_translation=NULL, Vector3 d_rotation=NULL);
		Vector3 get_translation_vector(int image);
		Vector3 get_rotation_angles(int image);
		void fillVolume();
		
		void printImage(int numSerie, int slice, int timeframe, string win_name, int range=0);
		void display3D(int timeframe, string str="");
};


typedef struct _myargs
{
	VolumeData* ptr;
	int serie;
	int slice;
	int frame;
	string win_name;
}myargs;

typedef struct _mousse_args_box
{
	bool actif;
	bool validated;
	bool cancelled;
	bool drawing_box;
	CvRect box;
	cv::Mat* img;
}mousse_args_box;

void updateRange(int value, void* arg);
void mouse_callback( int event, int x, int y, int flags, void* param );


class VolumeData_dicom: public VolumeData
{
	protected:
		vector<string> seriesImageType;
		vector<string> seriesDerivationDescription;
		map<int, vector< set<double> > > seriesTriggerTimes; //serie, slice, triggerTimes
		
		void listSeries();
		void printListedSeries();
		void printSelectedSeries();
		void readSerie(int serieNumber);
		void add_slices_to_serie(int serieNumber, set< pair<double, pair<int, int> > >& dists);
		
		void replace_slice(int serieNumber, int slice, int add_serie, int add_slice, double add_distance);
	public:
		VolumeData_dicom(XMLreader& _xml_reader);
};

class VolumeData_dicom_cardiac: public VolumeData_dicom
{
	protected:
		list<int> seriesSA;
		
		void selectSeries();
		void choose_serie_ref();
		void buildListTimeframes();
	
	public:
		VolumeData_dicom_cardiac(XMLreader& _xml_reader);
		
		bool comprise_temporal_sequences();
};

class VolumeData_nifti: public VolumeData
{
	protected:
		void listSeries();
		void readSerie(int serieNumber);
		void add_slices_to_serie(int serieNumber, set< pair<double, pair<int, int> > >& dists);
		void replace_slice(int serieNumber, int slice, int add_serie, int add_slice, double add_distance);
	public:
		VolumeData_nifti(XMLreader& _xml_reader);
};

class VolumeData_nrrd: public VolumeData
{
	protected:
		void listSeries();
		void readSerie(int serieNumber);
		void add_slices_to_serie(int serieNumber, set< pair<double, pair<int, int> > >& dists);
		void replace_slice(int serieNumber, int slice, int add_serie, int add_slice, double add_distance);
	public:
		VolumeData_nrrd(XMLreader& _xml_reader);
};

#endif

