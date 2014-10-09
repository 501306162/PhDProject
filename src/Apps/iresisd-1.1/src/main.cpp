/***************************************************************************
 *   Copyright (C) 2014 by A Paiement   *
 *   csatmp@bristol.ac.uk   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef NO_OMP
	#include <omp.h>
#endif
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <fftw3.h>

#include "XMLreader.h"
#include "Interface.h"

using namespace std;


int main(int argc, char *argv[])
{
try{
	/// load params
	
	string xml_file;

	if (argc < 2)
	{
		cout << "No address/name provided for the xml file -> defaulting to params.xml" << endl;
		xml_file = string("params.xml");
	}
	else
	{
		xml_file = string(argv[1]);
		cout << "loading " << xml_file << endl;
	}
	
	XMLreader xml_reader(xml_file);
	
	
	/// initialise fftw and openmp
	
	if ( !fftw_init_threads() ) throw string("Failed to initialize fftw-threads");
	
#ifndef NO_OMP
	int nb_proc = omp_get_num_procs();
	cout << nb_proc << " processing units found" << endl;
	cout << omp_get_max_threads() << " threads max possible" << endl;
	
	fftw_plan_with_nthreads(nb_proc);
	omp_set_num_threads(nb_proc);

	omp_set_dynamic(1);
#endif	
		
	/// create interface	
	
	Interface interface(&xml_reader);
	
	/// do the computations ///
	
	time_t start_time, end_time;
	double time_length;
	start_time = time(NULL);
	
	interface.loop();
	
	end_time = time(NULL);
	time_length = difftime(end_time, start_time);
	cout << "Time spent : " << time_length << " seconds, ie " << time_length/60. << " minutes or " << time_length/3600. << " hours" << endl;
	
	cv::destroyAllWindows();
	cv::waitKey(10);
	
	fftw_cleanup_threads();
	fftw_cleanup();
}
catch(string message)
{
	cout << message << endl;
	//cv::waitKey();
}
	
	return EXIT_SUCCESS;
}

