#include "Interface.h"

using namespace std;

Interface::Interface(XMLreader* _xml_reader)
{
	xml_reader = _xml_reader;
	volumeData = NULL;
	S_data = NULL;
	S_geom = NULL;
	isisd = NULL;
	iresd = NULL;
	
	currentStep = 1;
	
	scheme_type = xml_reader->get_schemeType();
	dirac_widths = xml_reader->get_dirac_widths();
	doRegistration = xml_reader->registration_required();
	
	loadData();
	create_psi();
	create_phi();
	S_data = new Sdata(xml_reader, volumeData, listPhi, listPsi[0]->get_size_fft(), listPsi[0]->get_borders() );
	S_geom = new Sgeom(xml_reader, listPhi, volumeData->getVolumeSize(), volumeData->is_motion_periodic());
	create_isisd();
	if (doRegistration)
	{
		bool global_registration = xml_reader->use_global_registration();
		bool selection_of_connected_regions = false;
		if (global_registration) selection_of_connected_regions = xml_reader->do_selection_of_connected_regions();
		int wait_before_rotation = xml_reader->wait_before_rotation();
		
		int slice_wise = xml_reader->do_slice_wise_registration();
		if (slice_wise == -1) // no instruction, defaulting based on the type of data
		{
			if ( volumeData->comprise_temporal_sequences() ) slice_wise = 1;
			else slice_wise = 0;
		}
		
		if (slice_wise) iresd = new IReSD_slice_wise(volumeData, S_data, listPhi, global_registration, selection_of_connected_regions, wait_before_rotation, xml_reader->shifts_in_SA_planes_at_beginning());
		else iresd = new IReSD_stack_wise(volumeData, S_data, listPhi, global_registration, selection_of_connected_regions, wait_before_rotation);
	}
	
	// stopping criteria
	
	for (int step = 0; step < nb_steps; step++)
	{
		nb_iterations.push_back( xml_reader->get_nb_iterations(step) );
		if (nb_iterations[step] == VAL_INVALID || nb_iterations[step] == 0)
			nb_stable_iterations_before_stop.push_back( xml_reader->get_max_stable_iterations(step) );
		else nb_stable_iterations_before_stop.push_back(VAL_INVALID);
	}
}

Interface::~Interface()
{
	if (volumeData != NULL) delete volumeData;
	
	vector<Phi*>::iterator vit(listPhi.begin()), vend(listPhi.end());
	for (; vit != vend; ++vit) delete *vit;
	
	vector<RBF*>::iterator vit2(listPsi.begin()), vend2(listPsi.end());
	for (; vit2 != vend2; ++vit2) delete *vit2;
	
	if (S_data != NULL) delete S_data;
	if (S_geom != NULL) delete S_geom;
	if (isisd != NULL) delete isisd;
	if (iresd != NULL) delete iresd;
}

void Interface::loadData()
{
	cout << "Interface::loadData()" << endl;
	
	string dataType = xml_reader->get_dataType();
	
	if ( !dataType.compare("dicom") )
		volumeData = new VolumeData_dicom(*xml_reader);
	else if ( !dataType.compare("dicom_cardiac") )
		volumeData = new VolumeData_dicom_cardiac(*xml_reader);
	else if ( !dataType.compare("nifti") )
		volumeData = new VolumeData_nifti(*xml_reader);
	else if ( !dataType.compare("nrrd") )
		volumeData = new VolumeData_nrrd(*xml_reader);
	
	volumeData->loadData();
	
	cout << "volumeData is created" << endl;
}

void Interface::create_psi()
{
	if ( !scheme_type.compare("simple") )
	{
		listPsi.resize(1);
		
		listPsi[0] = new RBF(xml_reader->get_beta(), xml_reader->get_gamma(0), volumeData->getVolumeSize(), xml_reader->get_volume_borders());
	}
	else if ( !scheme_type.compare("double") )
	{
		switch (currentStep)
		{
			case 1:
				listPsi.resize(1);
				
				// we start with only the first RBF, we will need the second one only later in the second step
				listPsi[0] = new RBF(xml_reader->get_beta(), xml_reader->get_gamma(0), volumeData->getVolumeSize(), xml_reader->get_volume_borders());
				
				break;
			
			case 2:
				delete listPsi[0];
				
				listPsi[0] = new RBF(xml_reader->get_beta(), xml_reader->get_gamma(1), volumeData->getVolumeSize(), xml_reader->get_volume_borders());
				
				break;
				
			default:
				throw string("Interface::create_psi() Error: unknown step");
		}
	}
	else if ( !scheme_type.compare("mixed") )
	{
		int nb_RBF_required = xml_reader->get_nb_RBFs();
		
		listPsi.resize(nb_RBF_required);
		
		// we need all the RBFs at the same time
		for (int i=0; i<nb_RBF_required; i++)
			listPsi[i] = new RBF(xml_reader->get_beta(), xml_reader->get_gamma(i), volumeData->getVolumeSize(), xml_reader->get_volume_borders());
	}
	else throw string("Interface::create_psi() Error: unknown scheme");
}

void Interface::create_phi()
{
	/// number of level sets ///
	
	nbPhi = xml_reader->get_nb_phi();

	if (nbPhi == 0)
	{
		while (nbPhi <= 0 || nbPhi > 2)
		{
			cout << "How many level sets should be initialised? (max 2)" << endl;
			cin >> nbPhi;
		}
	}
	
	/// create phi ///

	listPhi.resize(nbPhi);

	const int* size = volumeData->getVolumeSize();
	
	for (int i=0; i<nbPhi; i++)
	{
		listPhi[i] = new Phi(size);
	}
	
	
	/// initialise phi ///
	
	string phi_ini_type = xml_reader->get_phi_ini_type();
	
	float max_dirac_width;
	if ( !scheme_type.compare("mixed") )
		max_dirac_width = *max_element(dirac_widths.begin(), dirac_widths.end());
	else max_dirac_width = dirac_widths[currentStep-1];
	
	if ( !phi_ini_type.compare("nifti") )
	{
		for (int i=0; i<nbPhi; i++)
		{
			listPhi[i]->initialise_from_nii( xml_reader->get_phi_ini_nifti_file(i) );
			listPhi[i]->computeContour(max_dirac_width);
		}
	}
	else if ( !phi_ini_type.compare("circle") )
	{
		vector< list< Point4<int> > > centre_ini(nbPhi);
		vector< list<float> > radius_ini(nbPhi);
		
		xml_reader->get_initial_circles(centre_ini, radius_ini);
		
		bool initialised = true;
		
		for (int i=0; i<nbPhi; i++)
		{
			if ( centre_ini[i].empty() )
			{
				initialised = false;
				break;
			}
		}
		
		if (!initialised) queryInitialisationPhi(centre_ini, radius_ini);
		
		for (int i=0; i<nbPhi; i++)
		{
			listPhi[i]->initialiseSpheres(centre_ini[i], radius_ini[i]);
			listPhi[i]->computeContour(max_dirac_width);
		}
	}
	else throw string("Interface::create_phi() Error: unknown initialisation type");
}

void Interface::queryInitialisationPhi(vector< list< Point4<int> > >& centre_ini, vector< list<float> >& radius_ini)
{
	const vector< pair<int, int> >& selectedImages = volumeData->get_selectedImages();
	vector< pair<int, int> >::const_iterator vit, vend(selectedImages.end());
	int nbFrames = volumeData->getNbFrames();
	
	int centre_x, centre_y, perimetre_x, perimetre_y;
	float radius;
	
	cv::namedWindow("Initial circle/sphere", CV_WINDOW_AUTOSIZE);

	mousse_args_circle args;
	args.actif = false;
	cvSetMouseCallback("Initial circle/sphere", mouse_callback_initialisation, (void*)&args);
	
	int nbZones = (int)(pow(2., nbPhi));
	
	//vector< list< Point4<int> > > vect_centerIni(nbPhi);
	//vector< list<float> > vect_radiusIni(nbPhi);
	
	for (int p=0; p<nbZones-1; p++)
	{
		cout << "Please choose the initialisation of zone " << p << endl;
		
		list< Point4<int> > centre_tmp;
		list<float> radius_tmp;
		
		int image = 0;
		
		for (vit=selectedImages.begin(); vit != vend ; ++vit, image++)
		{
			cout << "image " << vit->first << " " << vit->second << endl;
			
			const pair<int, int>& imageSize = volumeData->getImageSize(vit->first);
			
			cv::Mat imgDisp(imageSize.first, imageSize.second, CV_64FC1);

			for (int f=0; f<nbFrames; f++)
			{
				cout << "frame " << f << endl;
				
				float max = -VAL_INVALID;
				float min = VAL_INVALID;
				
				for (int i=0; i<imageSize.first; i++)
				for (int j=0; j<imageSize.second; j++)
				{
					float val = volumeData->getPixelValue(vit->first, vit->second, f, j + imageSize.second * i);
		
					if (val != VAL_INVALID)
					{
						if (val > max) max = val;
						if (val < min) min = val;
					}
					
					imgDisp.at<double>(i, j) = (double)val;
				}
				
				if (min < -1000) min = -1000;
				//if (max < 1500) max = 1500;
				
				for (int i=0; i<imageSize.first; i++)
				for (int j=0; j<imageSize.second; j++)
				{
					if (imgDisp.at<double>(i, j) != VAL_INVALID)
						imgDisp.at<double>(i, j) = (imgDisp.at<double>(i, j) - min) / (max - min);
					else
						imgDisp.at<double>(i, j) = 0;
				}
				
				cv::imshow("Initial circle/sphere", imgDisp);
				
				args.drawing = false;
				args.img = &imgDisp;
				args.validated = false;
				args.cancelled = false;
				args.actif = true;
		
				cout << "Draw a cercle, then press 'v' to validate, or press 'c' to cancle" << endl;
				cout << "Press 'n' to go to the next image" << endl;
		
				char answer;
				bool nextSlice = false;

				while (!nextSlice)
				{
					while (!args.validated && !args.cancelled)
					{
						answer = cv::waitKey(15);
		
						if (answer == 'c')
						{
							args.cancelled = true;
							args.actif = false;
						}
						else if (answer == 'n')
						{
							args.cancelled = true;
							args.actif = false;
							nextSlice = true;
						}
					}
		
					if (args.validated)
					{
						centre_x = args.centre.y;
						centre_y = args.centre.x;
						
						perimetre_x = args.perimetre.y;
						perimetre_y = args.perimetre.x;
						
						Point3<float> point_centre = volumeData->getImagePointCoordsInVolume(image, computeIndex(imageSize, centre_x, centre_y));
						float time = volumeData->getTimeOfFrame(image, f);
						
						Point3<float> point_perimetre = volumeData->getImagePointCoordsInVolume(image, computeIndex(imageSize, perimetre_x, perimetre_y));
						radius = point_perimetre.distanceTo(point_centre);
						
						//cout << point_centre << " " << point_perimetre << " " << radius << endl;
						
						switch (p)
						{
							case 0:
								centre_ini[0].push_back( Point4<int>((int)round(point_centre.x()), (int)round(point_centre.y()), (int)round(point_centre.z()), (int)round(time)) );
								radius_ini[0].push_back(radius);
								
								cout << "new circle/sphere for phi_1 centred at " << centre_ini[0].back() << " with radius = " << radius_ini[0].back() << endl;
								
								if (nbPhi == 2)
								{
									centre_ini[1].push_back( Point4<int>((int)round(point_centre.x()), (int)round(point_centre.y()), (int)round(point_centre.z()), (int)round(time)) );
									radius_ini[1].push_back(radius);
									
									cout << "new circle/sphere for phi_2 centred at " << centre_ini[1].back() << " with radius = " << radius_ini[1].back() << endl;
								}
								
								break;
								
							case 1:
								centre_ini[0].push_back( Point4<int>((int)round(point_centre.x()), (int)round(point_centre.y()), (int)round(point_centre.z()), (int)round(time)) );
								radius_ini[0].push_back(radius);
								
								cout << "new circle/sphere for phi_1 centred at " << centre_ini[0].back() << " with radius = " << radius_ini[0].back() << endl;
								break;
								
							case 2:
								centre_ini[1].push_back( Point4<int>((int)round(point_centre.x()), (int)round(point_centre.y()), (int)round(point_centre.z()), (int)round(time)) );
								radius_ini[1].push_back(radius);
								
								cout << "new circle/sphere for phi_2 centred at " << centre_ini[1].back() << " with radius = " << radius_ini[1].back() << endl;
								break;
								
							default:
								cout << "unknown case!" << endl;
								break;
						}
		
						args.validated = false;
					}

					if (!nextSlice)
					{
						args.cancelled = false;
						args.actif = true;
					}
				}
			} //each frame
		} //each image
	} //each phi
	
	cvDestroyWindow("Initial circle/sphere");
	cv::waitKey(10);
}

void Interface::create_isisd()
{
	cout << "Interface::create_isisd()" << endl;
	
	int display = xml_reader->allow_display();
	
	if ( !scheme_type.compare("simple") )
	{
		nb_steps = 1;
		isisd = new ISISD_simpleScheme(S_data, S_geom, listPhi, listPsi[0], dirac_widths[0], volumeData, display);
	}
	else if ( !scheme_type.compare("double") )
	{
		nb_steps = 2;
		
		switch (currentStep)
		{
			case 1:
				isisd = new ISISD_simpleScheme(S_data, S_geom, listPhi, listPsi[0], dirac_widths[0], volumeData, display);
				break;
			
			case 2:
				((ISISD_simpleScheme*)isisd)->change_psi(listPsi[0], dirac_widths[1]);
				break;
				
			default:
				throw string("Interface::create_isisd() Error: unknown step");
		}
	}
	else if ( !scheme_type.compare("mixed") )
	{
		nb_steps = 1;
		
		isisd = new ISISD_mixedScheme(S_data, S_geom, listPhi, listPsi, dirac_widths, volumeData, display);
	}
	else throw string("Interface::create_isisd() Error: unknown scheme");
}

void Interface::loop()
{
	cout << "Interface::loop()" << endl;
	
	int display = xml_reader->allow_display();
	
	bool test_success, grow, stop = false;
	
	int nbFrames = volumeData->getNbFrames();
	
	for (int f=0; f<nbFrames; f++)
		displayContour(display, true, f, "start");
	
	int image;
	int nbImages = volumeData->getNbImages();
	
	
	for (int step = 0; step < nb_steps; step++)
	{
		cout << "step " << step+1 << endl;
		
		bool automatic_stop = false;
		if (nb_iterations[step] == VAL_INVALID || nb_iterations[step] == 0)
		{
			automatic_stop = true;
			cout << "loop will stop automatically when stable condition detected" << endl;
		}
		
		int nb_iterations_with_stable_contours = 0;
		int it = 0;
		
		while ( (automatic_stop && nb_iterations_with_stable_contours < nb_stable_iterations_before_stop[step]) ||
			(!automatic_stop && it != nb_iterations[step]) )
		{
			cout << "iteration " << it << endl;
			
			/// initialise the SpeedFunctions
			
			S_data->initialise_iteration();
			S_geom->initialise_iteration();
			if (doRegistration) iresd->initialise_iteration();
			
			for (int p=0; p<nbPhi; p++)
			{
				cout << "processing phi " << p << "..." << endl;
				
				/// compute S
				
				if (doRegistration)
					S_data->compute_S(p, iresd->get_global_variant_usage(), iresd->get_selection_of_connected_regions_policy());
				else S_data->compute_S(p);
				
				S_geom->compute_S(p);
				
				/// compute registration if needed
			
				if (doRegistration) iresd->compute_alignment_from_contour(p);
				
				/// update phi
				
				test_success = isisd->updatePhi(p);
				
				if ( !test_success )
				{
					if (doRegistration)
						stop = isisd->auto_reinitialisation(p, iresd->get_global_variant_usage(), iresd->get_selection_of_connected_regions_policy());
					else stop = isisd->auto_reinitialisation(p);
				}
				
				if (stop) break;
			}
			
			if (stop) break;
			
			/// apply registration from all combined contours/phi
			
			if (doRegistration)
			{
				cout << "registration..." << endl;
				
				iresd->compute_total_alignment();
				iresd->align_data();
			}
			
			/// check convergence
			
			grow = false;
			
			for (int p=0; p<nbPhi; p++)
			{
				if ( isisd->doesContourGrow(p) )
				{
					grow = true;
					break;
				}
			}
			
			if (!grow) nb_iterations_with_stable_contours++;
			else nb_iterations_with_stable_contours = 0;
			
			if (nb_iterations_with_stable_contours)
				cout << "contours have not significantly grown for " << nb_iterations_with_stable_contours << " iterations" << endl;
			
			it++;
			
			/// intermediate display ///

			if (display)
			{
				ostringstream str;
				str << "_it" << it;
				
				if (display == 2) //full display
				{
					for (image=0; image < nbImages; image++)
					{
						displayContour(display, false, image, nbFrames/2, str.str());
						if (nbPhi == 2) displayArea(display, false, image, nbFrames/2, str.str());
					}
				}
				
				displayContour(display, false, nbFrames/2, str.str());
			}
			
			cv::waitKey(10);
		} //while
		
		cout << "step " << step+1 << " complete" << endl;
		
		//cv::waitKey();
		
		if (stop) break;
		
		if ( step == 0 && !scheme_type.compare("double") )
		{
			vector<Phi*>::const_iterator vit(listPhi.begin()), vend(listPhi.end());
			const int* size = volumeData->getVolumeSize();
			
			for (int p=0; vit != vend; ++vit, p++)
			for (int t=0; t<size[3]; t++)
				(*vit)->saveContour(p, t, "end_of_first_step");
			
			currentStep++;
			create_psi();
			create_isisd();
			isisd->set_dynamic_maxTranslation(false);
		}
	} //for
	
	/// final display ///
	
	cout << "display final" << endl;
	
	for (image=0; image < nbImages; image++)
	{
		for (int f=0; f<nbFrames; f++)
		{
			displayContour(false, true, image, f, "");
			if (nbPhi == 2) displayArea(false, true, image, f, "");
		}
	}
	
	for (int f=0; f<nbFrames; f++)
		displayContour(false, true, f, "");
	
	vector<Phi*>::const_iterator vit(listPhi.begin()), vend(listPhi.end());
	const int* size = volumeData->getVolumeSize();
	
	for (int p=0; vit != vend; ++vit, p++)
	for (int t=0; t<size[3]; t++)
		(*vit)->saveContour(p, t, "final");
	
	if (doRegistration)
	{
		volumeData->display3D(0, "aligned");
		iresd->save_registration("registration.txt");
	}
	
	cout << "Interface::loop() done" << endl;
}

void Interface::displayContour(bool showImage, bool writeImage, int time, string str)
{
	const int* size = volumeData->getVolumeSize();
	
	cv::Mat imgDisp(size[0], size[1], CV_64FC1);
	
	float scaling = 1. - (nbPhi * 0.1);
	
	float pixMax = -VAL_INVALID, pixMin = VAL_INVALID;
	int planZ = size[2] / 2;
	int ind;
	float value;
	
	for (int x=0; x<size[0]; x++)
	for (int y=0; y<size[1]; y++)
	{
		ind = computeIndex(size, x, y, planZ);
		value = volumeData->getVoxelValue(ind, time);
		
		if (value != VAL_INVALID)
		{
			if ( value > pixMax ) pixMax = value;
			if ( value < pixMin ) pixMin = value;
		}
	}
	
	for (int x=0; x<size[0]; x++)
	for (int y=0; y<size[1]; y++)
	{
		ind = computeIndex(size, x, y, planZ);
		value = volumeData->getVoxelValue(ind, time);
		
		if (value != VAL_INVALID)
			imgDisp.at<double>(x, y) = (double)((value - pixMin) / (pixMax - pixMin) * scaling);
		else imgDisp.at<double>(x, y) = 0;
	}
	
	ostringstream name;
	name << "central_slice_" << time;
	
	if (showImage)
	{
		cv::namedWindow(name.str(), CV_WINDOW_AUTOSIZE);
		cv::imshow(name.str(), imgDisp);
		cv::waitKey(10);
	}
	
	if (writeImage)
	{
		name << str << ".bmp";
		imgDisp = imgDisp * 255.;
		cv::imwrite(name.str(), imgDisp);
		imgDisp = imgDisp / 255.;
	}
	
	for (int p=0; p<nbPhi; p++)
	{
		for (int x=0; x<size[0]; x++)
		for (int y=0; y<size[1]; y++)
		{
			if ( listPhi[p]->getDistToContour(computeIndex(size, x, y, planZ), time) < 1 )
				imgDisp.at<double>(x, y) = 1 - p * 0.1;
		}
	}
	
	ostringstream name2;
	name2 << "contour_central_slice_" << time;
	
	if (showImage)
	{
		cv::namedWindow(name2.str(), CV_WINDOW_AUTOSIZE);
		cv::imshow(name2.str(), imgDisp);
		cv::waitKey(10);
	}

	if (writeImage)
	{
		name2 << str << ".bmp";
		imgDisp = imgDisp * 255.;
		cv::imwrite(name2.str(), imgDisp);
	}
	
	if (size[2] > 1)
	{
		cv::Mat imgDispY(size[0], size[2], CV_64FC1);
		
		float pixMaxY = -VAL_INVALID, pixMinY = VAL_INVALID;
		int planY = size[1] / 2;
		
		for (int x=0; x<size[0]; x++)
		for (int z=0; z<size[2]; z++)
		{
			ind = computeIndex(size, x, planY, z);
			value = volumeData->getVoxelValue(ind, time);
			
			if (value != VAL_INVALID)
			{
				if ( value > pixMaxY ) pixMaxY = value;
				if ( value < pixMinY ) pixMinY = value;
			}
		}
		
		for (int x=0; x<size[0]; x++)
		for (int z=0; z<size[2]; z++)
		{
			ind = computeIndex(size, x, planY, z);
			value = volumeData->getVoxelValue(ind, time);
		
			if (value != VAL_INVALID)
				imgDispY.at<double>(x, z) = (double)((value - pixMinY) / (pixMaxY - pixMinY) * scaling);
			else imgDispY.at<double>(x, z) = 0;
		}
		
		
		ostringstream name;
		name << "central_vertical_slice_" << time;
		
		if (showImage)
		{
			cv::namedWindow(name.str(), CV_WINDOW_AUTOSIZE);
			cv::imshow(name.str(), imgDispY);
			cv::waitKey(10);
		}
		
		if (writeImage)
		{
			name << str << ".bmp";
			imgDispY = imgDispY * 255.;
			cv::imwrite(name.str(), imgDispY);
			imgDispY = imgDispY / 255.;
		}
		
		for (int p=0; p<nbPhi; p++)
		{
			for (int x=0; x<size[0]; x++)
			for (int z=0; z<size[2]; z++)
			{
				if ( listPhi[p]->getDistToContour(computeIndex(size, x, planY, z), time) < 1 )
					imgDispY.at<double>(x, z) = 1 - p * 0.1;
			}
		}
		
		ostringstream name2;
		name2 << "contour_central_vertical_slice_" << time;
		
		if (showImage)
		{
			cv::namedWindow(name2.str(), CV_WINDOW_AUTOSIZE);
			cv::imshow(name2.str(), imgDispY);
			cv::waitKey(10);
		}
		
		if (writeImage)
		{
			imgDispY = imgDispY * 255.;
			name2 << str << ".bmp";
			cv::imwrite(name2.str(), imgDispY);
		}
	}
}

void Interface::displayContour(bool showImage, bool writeImage, int image, int frame, string str)
{
	const pair<int, int>& serieSlice = volumeData->get_selectedImages()[image];
	pair<int, int> imgSize = volumeData->VolumeData::getImageSize(serieSlice.first);
	int Rows = imgSize.first;
	int Columns = imgSize.second;
	int nbPixImg = Rows * Columns;
	
	cv::Mat imgDisp(Rows, Columns, CV_64FC1);
	
	float scaling = 1. - (nbPhi * 0.1);
	
	Point3<float> point;
	Point3<int> point_int;
	
	float pixVal;
	float pixMax = -VAL_INVALID, pixMin = VAL_INVALID;
	
	for (int indImg=0; indImg<nbPixImg; indImg++)
	{
		if ( !volumeData->isImagePointInVolume(image, indImg) ) continue;
		
		pixVal = volumeData->getPixelValue(image, frame, indImg);
		
		if (pixVal == VAL_INVALID) continue;
		
		if (pixVal > pixMax) pixMax = pixVal;
		if (pixVal < pixMin) pixMin = pixVal;
	}
	
	for (int i=0; i<Rows; i++)
	for (int j=0; j<Columns; j++)
	{
		pixVal = volumeData->getPixelValue(image, frame, j + Columns * i);
		
		if (pixVal == VAL_INVALID) imgDisp.at<double>(i, j) = 0;
		else
		{
			imgDisp.at<double>(i, j) = (double)((pixVal - pixMin) / (pixMax - pixMin) * scaling);
		}
	}

	const int* size = volumeData->getVolumeSize();
	int time = (int)round(volumeData->getTimeOfFrame(image, frame));

	for (int p=0; p<nbPhi; p++)
	{
		for (int i=0; i<Rows; i++)
		for (int j=0; j<Columns; j++)
		{
			point = volumeData->getImagePointCoordsInVolume(image, j + Columns * i);
			
			point_int = Point3<int>((int)round(point.x()), (int)round(point.y()), (int)round(point.z()));
			
			if ( !volumeData->isPointInVolume(point_int) ) continue;
			
			if ( listPhi[p]->getDistToContour(computeIndex(size, point_int.x(), point_int.y(), point_int.z()), time) < 1 )
				imgDisp.at<double>(i, j) = 1 - p * 0.1;
		}
	}
	
	ostringstream name2;
	name2 << "contour_" << serieSlice.first << "_" << serieSlice.second << "_" << frame;

	if (showImage)
	{
		cv::namedWindow(name2.str(), CV_WINDOW_AUTOSIZE);
		cv::imshow(name2.str(), imgDisp);
		cv::waitKey(10);
	}

	if (writeImage)
	{
		name2 << str << ".bmp";
		imgDisp = imgDisp * 255.;
		cv::imwrite(name2.str(), imgDisp);
	}
}

void Interface::displayArea(bool showImage, bool writeImage, int time, string str)
{
	const int* size = volumeData->getVolumeSize();
	
	IplImage* imgDispColor = cvCreateImage(cvSize(size[1], size[0]), IPL_DEPTH_32F, 3);
	int step = imgDispColor->widthStep/sizeof(float);
	int channels = imgDispColor->nChannels;
	float* data = (float *)imgDispColor->imageData;
	
	float pixMax = -VAL_INVALID, pixMin = VAL_INVALID;
	int planZ = size[2] / 2;
	int ind;
	float value;
	
	for (int x=0; x<size[0]; x++)
	for (int y=0; y<size[1]; y++)
	{
		ind = computeIndex(size, x, y, planZ);
		value = volumeData->getVoxelValue(ind, time);
		
		if (value != VAL_INVALID)
		{
			if ( value > pixMax ) pixMax = value;
			if ( value < pixMin ) pixMin = value;
		}
	}
	
	for (int x=0; x<size[0]; x++)
	for (int y=0; y<size[1]; y++)
	{
		ind = computeIndex(size, x, y, planZ);
		value = volumeData->getVoxelValue(ind, time);
		
		if (value != VAL_INVALID)
		{
			data[x*step+y*channels] = (value - pixMin) / (pixMax - pixMin);
			data[x*step+y*channels+1] = (value - pixMin) / (pixMax - pixMin);
			data[x*step+y*channels+2] = (value - pixMin) / (pixMax - pixMin);
		}
		else
		{
			data[x*step+y*channels] = 0;
			data[x*step+y*channels+1] = 0;
			data[x*step+y*channels+2] = 0;
		}
	}
	
	ostringstream name;
	name << "central_slice_" << time;
	
	if (showImage)
	{
		cv::Mat mtx(imgDispColor);
		cv::namedWindow(name.str(), CV_WINDOW_AUTOSIZE);
		cv::imshow(name.str(), mtx);
		cv::waitKey(10);
	}
	
	if (writeImage)
	{
		name << str << ".bmp";
		
		for (int x=0; x<size[0]; x++)
		for (int y=0; y<size[1]; y++)
		{
			data[x*step+y*channels] *= 255.;
			data[x*step+y*channels+1] *= 255.;
			data[x*step+y*channels+2] *= 255.;
		}
		
		cv::Mat mtx(imgDispColor);
		cv::imwrite(name.str(), mtx);
		
		for (int x=0; x<size[0]; x++)
		for (int y=0; y<size[1]; y++)
		{
			data[x*step+y*channels] /= 255.;
			data[x*step+y*channels+1] /= 255.;
			data[x*step+y*channels+2] /= 255.;
		}
	}
	
	float phiVal0, phiVal1 = -VAL_INVALID;
	
	for (int x=0; x<size[0]; x++)
	for (int y=0; y<size[1]; y++)
	{
		phiVal0 = listPhi[0]->get_phi_value(x, y, planZ, time);
		if (nbPhi > 1) phiVal1 = listPhi[1]->get_phi_value(x, y, planZ, time);
		
		if (phiVal0 >= 0 && phiVal1 < 0) //red
		{
			data[x*step+y*channels+2] = 0.5 * data[x*step+y*channels+2] + 0.5;
		}
		else if (phiVal0 >= 0 && phiVal1 >= 0) //green
		{
			data[x*step+y*channels+1] = 0.5 * data[x*step+y*channels+1] + 0.5;
		}
		else if (phiVal0 < 0 && phiVal1 >= 0) //blue
		{
			data[x*step+y*channels] = 0.5 * data[x*step+y*channels] + 0.5;
		}
	}
	
	ostringstream name2;
	name2 << "areas_central_slice_" << time;
	
	if (showImage)
	{
		cv::Mat mtx(imgDispColor);
		cv::namedWindow(name2.str(), CV_WINDOW_AUTOSIZE);
		cv::imshow(name2.str(), mtx);
		cv::waitKey(10);
	}

	if (writeImage)
	{
		name << str << ".bmp";
		
		for (int x=0; x<size[0]; x++)
		for (int y=0; y<size[1]; y++)
		{
			data[x*step+y*channels] *= 255.;
			data[x*step+y*channels+1] *= 255.;
			data[x*step+y*channels+2] *= 255.;
		}
		
		cv::Mat mtx(imgDispColor);
		cv::imwrite(name2.str(), mtx);
		
		for (int x=0; x<size[0]; x++)
		for (int y=0; y<size[1]; y++)
		{
			data[x*step+y*channels] /= 255.;
			data[x*step+y*channels+1] /= 255.;
			data[x*step+y*channels+2] /= 255.;
		}
	}
	
	cvReleaseImage(&imgDispColor);
	
	if (size[2] > 1)
	{
		IplImage* imgDispColor2 = cvCreateImage(cvSize(size[2], size[0]), IPL_DEPTH_32F, 3);
		step = imgDispColor2->widthStep/sizeof(float);
		channels = imgDispColor2->nChannels;
		float* data2 = (float *)imgDispColor2->imageData;
		
		float pixMaxY = -VAL_INVALID, pixMinY = VAL_INVALID;
		int planY = size[1] / 2;
		
		for (int x=0; x<size[0]; x++)
		for (int z=0; z<size[2]; z++)
		{
			ind = computeIndex(size, x, planY, z);
			value = volumeData->getVoxelValue(ind, time);
			
			if (value != VAL_INVALID)
			{
				if ( value > pixMaxY ) pixMaxY = value;
				if ( value < pixMinY ) pixMinY = value;
			}
		}
		
		for (int x=0; x<size[0]; x++)
		for (int z=0; z<size[2]; z++)
		{
			ind = computeIndex(size, x, planY, z);
			value = volumeData->getVoxelValue(ind, time);
		
			if (value != VAL_INVALID)
			{
				data2[x*step+z*channels] = (value - pixMin) / (pixMax - pixMin);
				data2[x*step+z*channels+1] = (value - pixMin) / (pixMax - pixMin);
				data2[x*step+z*channels+2] = (value - pixMin) / (pixMax - pixMin);
			}
			else
			{
				data2[x*step+z*channels] = 0;
				data2[x*step+z*channels+1] = 0;
				data2[x*step+z*channels+2] = 0;
			}
		}
		
		
		ostringstream name;
		name << "central_vertical_slice_" << time;
		
		if (showImage)
		{
			cv::Mat mtx(imgDispColor2);
			cv::namedWindow(name.str(), CV_WINDOW_AUTOSIZE);
			cv::imshow(name.str(), mtx);
			cv::waitKey(10);
		}
		
		if (writeImage)
		{
			name << str << ".bmp";
			
			for (int x=0; x<size[0]; x++)
			for (int z=0; z<size[2]; z++)
			{
				data2[x*step+z*channels] *= 255.;
				data2[x*step+z*channels+1] *= 255.;
				data2[x*step+z*channels+2] *= 255.;
			}
			
			cv::Mat mtx(imgDispColor2);
			cv::imwrite(name.str(), mtx);
			
			for (int x=0; x<size[0]; x++)
			for (int z=0; z<size[2]; z++)
			{
				data2[x*step+z*channels] /= 255.;
				data2[x*step+z*channels+1] /= 255.;
				data2[x*step+z*channels+2] /= 255.;
			}
		}
		
		for (int x=0; x<size[0]; x++)
		for (int z=0; z<size[2]; z++)
		{
			phiVal0 = listPhi[0]->get_phi_value(x, planY, z, time);
			if (nbPhi > 1) phiVal1 = listPhi[1]->get_phi_value(x, planY, z, time);
			
			if (phiVal0 >= 0 && phiVal1 < 0) //red
			{
				data2[x*step+z*channels+2] = 0.5 * data2[x*step+z*channels+2] + 0.5;
			}
			else if (phiVal0 >= 0 && phiVal1 >= 0) //green
			{
				data2[x*step+z*channels+1] = 0.5 * data2[x*step+z*channels+1] + 0.5;
			}
			else if (phiVal0 < 0 && phiVal1 >= 0) //blue
			{
				data2[x*step+z*channels] = 0.5 * data2[x*step+z*channels] + 0.5;
			}
		}
		
		ostringstream name2;
		name2 << "areas_central_vertical_slice_" << time;
		
		if (showImage)
		{
			cv::Mat mtx(imgDispColor2);
			cv::namedWindow(name2.str(), CV_WINDOW_AUTOSIZE);
			cv::imshow(name2.str(), mtx);
			cv::waitKey(10);
		}
		
		if (writeImage)
		{
			name << str << ".bmp";
			
			for (int x=0; x<size[0]; x++)
			for (int z=0; z<size[2]; z++)
			{
				data2[x*step+z*channels] *= 255.;
				data2[x*step+z*channels+1] *= 255.;
				data2[x*step+z*channels+2] *= 255.;
			}
			
			cv::Mat mtx(imgDispColor2);
			cv::imwrite(name2.str(), mtx);
			
			for (int x=0; x<size[0]; x++)
			for (int z=0; z<size[2]; z++)
			{
				data2[x*step+z*channels] /= 255.;
				data2[x*step+z*channels+1] /= 255.;
				data2[x*step+z*channels+2] /= 255.;
			}
		}
		
		cvReleaseImage(&imgDispColor2);
	}
}

void Interface::displayArea(bool showImage, bool writeImage, int image, int frame, string str)
{
	const pair<int, int>& serieSlice = volumeData->get_selectedImages()[image];
	pair<int, int> imgSize = volumeData->VolumeData::getImageSize(serieSlice.first);
	int Rows = imgSize.first;
	int Columns = imgSize.second;
	int nbPixImg = Rows * Columns;
	
	IplImage* imgDispColor = cvCreateImage(cvSize(Columns, Rows), IPL_DEPTH_32F, 3);
	int step = imgDispColor->widthStep/sizeof(float);
	int channels = imgDispColor->nChannels;
	float* data = (float *)imgDispColor->imageData;
	
	Point3<float> point;
	Point3<int> point_int;
	
	float pixVal;
	float pixMax = -VAL_INVALID, pixMin = VAL_INVALID;
	
	for (int indImg=0; indImg<nbPixImg; indImg++)
	{
		if ( !volumeData->isImagePointInVolume(image, indImg) ) continue;
		
		pixVal = volumeData->getPixelValue(image, frame, indImg);
		
		if (pixVal == VAL_INVALID) continue;
		
		if (pixVal > pixMax) pixMax = pixVal;
		if (pixVal < pixMin) pixMin = pixVal;
	}
	
	for (int i=0; i<Rows; i++)
	for (int j=0; j<Columns; j++)
	{
		pixVal = volumeData->getPixelValue(image, frame, j + Columns * i);
		
		if (pixVal != VAL_INVALID)
		{
			data[i*step+j*channels] = (pixVal - pixMin) / (pixMax - pixMin);
			data[i*step+j*channels+1] = (pixVal - pixMin) / (pixMax - pixMin);
			data[i*step+j*channels+2] = (pixVal - pixMin) / (pixMax - pixMin);
		}
		else
		{
			data[i*step+j*channels] = 0;
			data[i*step+j*channels+1] = 0;
			data[i*step+j*channels+2] = 0;
		}
	}
	
	/*ostringstream name;
	name << "image_" << serieSlice.first << "_" << serieSlice.second << "_" << frame;
	
	if (showImage)
	{
		cv::Mat mtx(imgDispColor);
		cv::namedWindow(name.str(), CV_WINDOW_AUTOSIZE);
		cv::imshow(name.str(), mtx);
		cv::waitKey(10);
	}
	
	if (writeImage)
	{
		name << str << ".bmp";
		
		for (int i=0; i<Rows; i++)
		for (int j=0; j<Columns; j++)
		{
			data[i*step+j*channels] *= 255.;
			data[i*step+j*channels+1] *= 255.;
			data[i*step+j*channels+2] *= 255.;
		}
		
		cv::Mat mtx(imgDispColor);
		cv::imwrite(name.str(), mtx);
		
		for (int i=0; i<Rows; i++)
		for (int j=0; j<Columns; j++)
		{
			data[i*step+j*channels] /= 255.;
			data[i*step+j*channels+1] /= 255.;
			data[i*step+j*channels+2] /= 255.;
		}
	}*/

	float phiVal0, phiVal1 = -VAL_INVALID;
	float time = volumeData->getTimeOfFrame(image, frame);
	
	for (int i=0; i<Rows; i++)
	for (int j=0; j<Columns; j++)
	{
		point = volumeData->getImagePointCoordsInVolume(image, j + Columns * i);
		
		if ( !volumeData->isPointInVolume(point) ) continue;
		
		phiVal0 = listPhi[0]->get_phi_value(point.x(), point.y(), point.z(), time);
		if (nbPhi > 1) phiVal1 = listPhi[1]->get_phi_value(point.x(), point.y(), point.z(), time);
		
		if (phiVal0 >= 0 && phiVal1 < 0) //red
		{
			data[i*step+j*channels+2] = 0.5 * data[i*step+j*channels+2] + 0.5;
		}
		else if (phiVal0 >= 0 && phiVal1 >= 0) //green
		{
			data[i*step+j*channels+1] = 0.5 * data[i*step+j*channels+1] + 0.5;
		}
		else if (phiVal0 < 0 && phiVal1 >= 0) //blue
		{
			data[i*step+j*channels] = 0.5 * data[i*step+j*channels] + 0.5;
		}
	}
	
	ostringstream name2;
	name2 << "areas_" << serieSlice.first << "_" << serieSlice.second << "_" << frame;

	if (showImage)
	{
		cv::Mat mtx(imgDispColor);
		cv::namedWindow(name2.str(), CV_WINDOW_AUTOSIZE);
		cv::imshow(name2.str(), mtx);
		cv::waitKey(10);
	}

	if (writeImage)
	{
		name2 << str << ".bmp";
		
		for (int i=0; i<Rows; i++)
		for (int j=0; j<Columns; j++)
		{
			data[i*step+j*channels] *= 255.;
			data[i*step+j*channels+1] *= 255.;
			data[i*step+j*channels+2] *= 255.;
		}
		
		cv::Mat mtx(imgDispColor);
		cv::imwrite(name2.str(), mtx);
		
		for (int i=0; i<Rows; i++)
		for (int j=0; j<Columns; j++)
		{
			data[i*step+j*channels] /= 255.;
			data[i*step+j*channels+1] /= 255.;
			data[i*step+j*channels+2] /= 255.;
		}
	}
	
	cvReleaseImage(&imgDispColor);
}


/// utilities ///

void mouse_callback_initialisation( int event, int x, int y, int flags, void* param )
{
	mousse_args_circle* arguments = (mousse_args_circle*)param;

	if (arguments->actif == false) return;
	
	switch( event )
	{
		case CV_EVENT_LBUTTONDOWN:
			arguments->drawing = true;
			arguments->centre.x = x;
			arguments->centre.y = y;
			break;
			
		case CV_EVENT_LBUTTONUP:
			arguments->drawing = false;
			
			arguments->perimetre.x = x;
			arguments->perimetre.y = y;
			float radius = sqrt(pow((float)(x-arguments->centre.x), 2) + pow((float)(y-arguments->centre.y), 2));
			
			cv::Mat tmp(*(arguments->img));
			arguments->img->copyTo(tmp);
			
			cv::circle(tmp, arguments->centre, (int)round(radius), cvScalar(255, 0, 0) );
			cv::imshow("Initial circle/sphere", tmp);

			bool stop = false;
			char choice;
			
			while (!stop)
			{
				choice = cv::waitKey();
				
				switch (choice)
				{
					case 'v':
						arguments->validated = true;
						arguments->actif = false;
						stop = true;
						break;
						
					case 'c':
						arguments->cancelled = true;
						arguments->actif = false;
						stop = true;
						break;
						
					default:
						cout << "Please choose between 'v' and 'c'" << endl;
						break;
				}
			}
			
			break;
	}
}

