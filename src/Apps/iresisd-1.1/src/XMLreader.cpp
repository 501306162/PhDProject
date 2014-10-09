#include "XMLreader.h"

using namespace std;
using namespace tinyxml2;

XMLreader::XMLreader(string xml_file)
{
	doc.LoadFile(xml_file.c_str());

	if ( doc.ErrorID() )
	{
		cout << doc.ErrorID() <<endl;
		throw string("Error in XMLreader::XMLreader()");
	}
}

string XMLreader::get_dataType()
{
	XMLElement* textElement = doc.FirstChildElement("data");
	if (textElement == NULL) throw string("XMLreader::get_dataType() Error: data required");
	
	textElement = textElement->FirstChildElement("datatype");
	if (textElement == NULL) throw string("XMLreader::get_dataType() Error: datatype required");
	
	return string( textElement->GetText() );
}

string XMLreader::get_repertory()
{
	XMLElement* textElement = doc.FirstChildElement("data");
	if (textElement == NULL) throw string("XMLreader::get_repertory() Error: data required");
	
	textElement = textElement->FirstChildElement("directory");
	if (textElement == NULL) throw string("XMLreader::get_repertory() Error: directory required");
	
	return string( textElement->GetText() );
}

bool XMLreader::is_motion_periodic()
{
	bool value;
	XMLElement* element = doc.FirstChildElement("data");
	if (element == NULL) throw string("XMLreader::is_motion_periodic() Error: data required");
	
	element = element->FirstChildElement("periodic_motion");
	if (element == NULL) return false;
	
	if ( element->QueryBoolText(&value) ) return false;
	return value;
}

bool XMLreader::accurate_volume_filling()
{
	bool value;
	XMLElement* element = doc.FirstChildElement("data");
	if (element == NULL) throw string("XMLreader::accurate_volume_filling() Error: data required");
	
	element = element->FirstChildElement("accurate_volume_filling");
	if (element == NULL) return false;
	
	if ( element->QueryBoolText(&value) ) return false;
	return value;
}

bool XMLreader::size_volume_known()
{
	XMLElement* element = doc.FirstChildElement("data");
	if (element == NULL) throw string("XMLreader::size_volume_known() Error: data required");
	
	element = element->FirstChildElement("volume_size");
	if (element == NULL) return false;
	
	return true;
}

Vector3 XMLreader::get_volume_size()
{
	int sx, sy, sz;
	XMLElement* element = doc.FirstChildElement("data");
	if (element == NULL) throw string("XMLreader::get_volume_size() Error: data required");
	
	element = element->FirstChildElement("volume_size");
	if (element == NULL) throw string("XMLreader::get_volume_size() Error: volume_size element required");
	
	XMLElement* element_x = element->FirstChildElement("x");
	if (element_x == NULL) throw string("XMLreader::get_volume_size() Error: x element required");
	if ( element_x->QueryIntText(&sx) ) throw string("XMLreader::get_volume_size() Error: cannot read x");
	
	XMLElement* element_y = element->FirstChildElement("y");
	if (element_y == NULL) throw string("XMLreader::get_volume_size() Error: y element required");
	if ( element_y->QueryIntText(&sy) ) throw string("XMLreader::get_volume_size() Error: cannot read y");
	
	XMLElement* element_z = element->FirstChildElement("z");
	if (element_z == NULL) throw string("XMLreader::get_volume_size() Error: z element required");
	if ( element_z->QueryIntText(&sz) ) throw string("XMLreader::get_volume_size() Error: cannot read z");
	
	return Vector3(sx, sy, sz);
}

Vector3 XMLreader::get_ROI_offset()
{
	int ox, oy, oz;
	XMLElement* element = doc.FirstChildElement("data");
	if (element == NULL) throw string("XMLreader::get_ROI_offset() Error: data required");
	
	element = element->FirstChildElement("ROI_offset");
	if (element == NULL) throw string("XMLreader::get_ROI_offset() Error: ROI_offset element required");
	
	XMLElement* element_x = element->FirstChildElement("x");
	if (element_x == NULL) throw string("XMLreader::get_ROI_offset() Error: x element required");
	if ( element_x->QueryIntText(&ox) ) throw string("XMLreader::get_ROI_offset() Error: cannot read x");
	
	XMLElement* element_y = element->FirstChildElement("y");
	if (element_y == NULL) throw string("XMLreader::get_ROI_offset() Error: y element required");
	if ( element_y->QueryIntText(&oy) ) throw string("XMLreader::get_ROI_offset() Error: cannot read y");
	
	XMLElement* element_z = element->FirstChildElement("z");
	if (element_z == NULL) throw string("XMLreader::get_ROI_offset() Error: z element required");
	if ( element_z->QueryIntText(&oz) ) throw string("XMLreader::get_ROI_offset() Error: cannot read z");
	
	return Vector3(ox, oy, oz);
}

bool XMLreader::serie_ref_defines_volume_size()
{
	bool value;
	XMLElement* element = doc.FirstChildElement("data");
	if (element == NULL) throw string("XMLreader::serie_ref_defines_volume_size() Error: data required");
	
	element = element->FirstChildElement("serie_ref_defines_volume_size");
	if (element == NULL) return false;
	
	if ( element->QueryBoolText(&value) ) return false;
	return value;
}

bool XMLreader::ROI_selection()
{
	bool value;
	XMLElement* element = doc.FirstChildElement("data");
	if (element == NULL) throw string("XMLreader::ROI_selection() Error: data required");
	
	element = element->FirstChildElement("ROI_selection");
	if (element == NULL) return false;
	
	if ( element->QueryBoolText(&value) ) return false;
	return value;
}

bool XMLreader::use_first_timeframe_only()
{
	bool value;
	XMLElement* element = doc.FirstChildElement("data");
	if (element == NULL) throw string("XMLreader::use_first_timeframe_only() Error: data required");
	
	element = element->FirstChildElement("use_first_timeframe_only");
	if (element == NULL) return false;
	
	if ( element->QueryBoolText(&value) ) return false;
	return value;
}

Vector3 XMLreader::get_registration_translation(string serieDescription, int slice)
{
	XMLElement* element_serie = doc.FirstChildElement("data");
	if (element_serie == NULL) throw string("XMLreader::get_registration_translation() Error: data required");
	
	element_serie = element_serie->FirstChildElement("registration");
	if (element_serie == NULL) return VECT_NULL;
	
	element_serie = element_serie->FirstChildElement(serieDescription.c_str());
	if (element_serie == NULL) return VECT_NULL;
	
	ostringstream stream;
	stream << "slice_" << slice;
	XMLElement* element_slice = element_serie->FirstChildElement(stream.str().c_str());
	
	float tx, ty, tz;
	XMLElement* element;
	
	if (element_slice == NULL) element = element_serie->FirstChildElement("translation");
	else element = element_slice->FirstChildElement("translation");
	if (element == NULL) return VECT_NULL;
	
	XMLElement* element_x = element->FirstChildElement("x");
	if (element_x == NULL) throw string("XMLreader::get_registration_translation() Error: x element required for translation");
	if ( element_x->QueryFloatText(&tx) ) throw string("XMLreader::get_registration_translation() Error: cannot read x");
	
	XMLElement* element_y = element->FirstChildElement("y");
	if (element_y == NULL) throw string("XMLreader::get_registration_translation() Error: y element required for translation");
	if ( element_y->QueryFloatText(&ty) ) throw string("XMLreader::get_registration_translation() Error: cannot read y");
	
	XMLElement* element_z = element->FirstChildElement("z");
	if (element_z == NULL) throw string("XMLreader::get_registration_translation() Error: z element required for translation");
	if ( element_z->QueryFloatText(&tz) ) throw string("XMLreader::get_registration_translation() Error: cannot read z");
	
	return Vector3(tx, ty, tz);
}

Vector3 XMLreader::get_registration_rotation(string serieDescription, int slice)
{
	XMLElement* element_serie = doc.FirstChildElement("data");
	if (element_serie == NULL) throw string("XMLreader::get_registration_rotation() Error: data required");
	
	element_serie = element_serie->FirstChildElement("registration");
	if (element_serie == NULL) return VECT_NULL;
	element_serie = element_serie->FirstChildElement(serieDescription.c_str());
	if (element_serie == NULL) return VECT_NULL;
	
	ostringstream stream;
	stream << "slice_" << slice;
	XMLElement* element_slice = element_serie->FirstChildElement(stream.str().c_str());
	
	float rx, ry, rz;
	XMLElement* element;
	
	if (element_slice == NULL) element = element_serie->FirstChildElement("rotation");
	else element = element_slice->FirstChildElement("rotation");
	
	XMLElement* element_x = element->FirstChildElement("x");
	if (element_x == NULL) throw string("XMLreader::get_registration_rotation() Error: x element required for rotation");
	if ( element_x->QueryFloatText(&rx) ) throw string("XMLreader::get_registration_rotation() Error: cannot read x");
	
	XMLElement* element_y = element->FirstChildElement("y");
	if (element_y == NULL) throw string("XMLreader::get_registration_rotation() Error: y element required for rotation");
	if ( element_y->QueryFloatText(&ry) ) throw string("XMLreader::get_registration_rotation() Error: cannot read y");
	
	XMLElement* element_z = element->FirstChildElement("z");
	if (element_z == NULL) throw string("XMLreader::get_registration_rotation() Error: z element required for rotation");
	if ( element_z->QueryFloatText(&rz) ) throw string("XMLreader::get_registration_rotation() Error: cannot read z");
	
	return Vector3(rx, ry, rz);
}

int XMLreader::get_nb_RBFs()
{
	const XMLElement* element = doc.FirstChildElement("RBF");
	if (element == NULL) throw string("XMLreader::get_nb_RBFs() Error: RBF required");
	
	int count = 0;
	
	for( element = element->FirstChildElement("psi");
			 element;
			 element = element->NextSiblingElement("psi") )
	{
		count++;
	}
	
	return count;
}

float XMLreader::get_beta()
{
	float value;
	XMLElement* element = doc.FirstChildElement("RBF");
	if (element == NULL) throw string("XMLreader::get_beta() Error: RBF required");
	
	element = element->FirstChildElement("beta");
	if (element == NULL) throw string("XMLreader::get_beta() Error: beta required");
	
	if ( element->QueryFloatText(&value) ) throw string("XMLreader::get_beta() Error: cannot read beta");
	cout << "beta = " << value << endl;
	return value;
}

float XMLreader::get_gamma(int rbf_index)
{
	int nb_rbf = get_nb_RBFs();
	
	rbf_index++; //xml value starts at 1, C code value starts at 0

	float value;
	XMLElement* element = doc.FirstChildElement("RBF");
	if (element == NULL) throw string("XMLreader::get_gamma() Error: RBF required");
	
	int index = -1;
	XMLElement* element_psi, *element_index;
	
	for( element_psi = element->FirstChildElement("psi");
			 element_psi && index != rbf_index; )
	{
		element_index = element_psi->FirstChildElement("index");
		if (nb_rbf > 1 && element_index == NULL) throw string("XMLreader::get_beta() Error: index required");
		
		if (element_index == NULL) index = 1;
		else
			if ( element_index->QueryIntText(&index) ) throw string("XMLreader::get_beta() Error: cannot read index");
		
		if (index != rbf_index) element_psi = element_psi->NextSiblingElement("psi");
	}
	
	if (element_psi == NULL) throw string("XMLreader::get_gamma() Error: psi required");
	if (index != rbf_index) throw string("XMLreader::get_gamma() Error: could not find required index");
	
	element = element_psi->FirstChildElement("gamma");
	if (element == NULL) throw string("XMLreader::get_gamma() Error: gamma required");
	
	if ( element->QueryFloatText(&value) ) throw string("XMLreader::get_beta() Error: cannot read gamma");
	cout << "gamma = " << value << endl;
	return value;
}

vector<int> XMLreader::get_volume_borders()
{
	vector<int> borders(4, 0);
	
	XMLElement* element = doc.FirstChildElement("RBF");
	if (element == NULL) throw string("XMLreader::get_volume_borders() Error: RBF required");
	
	element = element->FirstChildElement("borders");
	if (element == NULL)
	{
		for (int i=0; i<4; i++) borders[i] = 20;
		return borders;
	}
	
	XMLElement* element_x = element->FirstChildElement("x");
	if (element_x == NULL) throw string("XMLreader::get_volume_borders() Error: x element required");
	if ( element_x->QueryIntText(&(borders[0])) ) throw string("XMLreader::get_volume_borders() Error: cannot read x");
	
	XMLElement* element_y = element->FirstChildElement("y");
	if (element_y == NULL) throw string("XMLreader::get_volume_borders() Error: y element required");
	if ( element_y->QueryIntText(&(borders[1])) ) throw string("XMLreader::get_volume_borders() Error: cannot read y");
	
	XMLElement* element_z = element->FirstChildElement("z");
	if (element_z == NULL) throw string("XMLreader::get_volume_borders() Error: z element required");
	if ( element_z->QueryIntText(&(borders[2])) ) throw string("XMLreader::get_volume_borders() Error: cannot read z");
	
	XMLElement* element_t = element->FirstChildElement("t");
	if (element_t == NULL) throw string("XMLreader::get_volume_borders() Error: t element required");
	if ( element_t->QueryIntText(&(borders[3])) ) throw string("XMLreader::get_volume_borders() Error: cannot read t");
	
	cout << "borders = " << borders[0] << " " << borders[1] << " " << borders[2] << " " << borders[3] << endl;
	
	return borders;
}

int XMLreader::get_nb_phi()
{
	const XMLElement* element = doc.FirstChildElement("level_set");
	if (element == NULL) throw string("XMLreader::get_nb_phi() Error: level_set required");
	
	int count = 0;
	
	for( element = element->FirstChildElement("phi");
			 element;
			 element = element->NextSiblingElement("phi") )
	{
		count++;
	}
	
	return count;
}

string XMLreader::get_phi_ini_type()
{
	XMLElement* textElement = doc.FirstChildElement("level_set");
	if (textElement == NULL) throw string("XMLreader::get_phi_ini_type() Error: level_set required");
	
	textElement = textElement->FirstChildElement("initialisation");
	if (textElement == NULL) return string("circle");
	
	return string( textElement->GetText() );
}

string XMLreader::get_phi_ini_nifti_file(int phi_index)
{
	int nb_phi = get_nb_phi();
	
	phi_index++; //xml value starts at 1, C code value starts at 0

	const XMLElement* element = doc.FirstChildElement("level_set");
	if (element == NULL) throw string("XMLreader::get_phi_ini_nifti_file() Error: level_set required");
	
	int index = -1;
	const XMLElement* element_phi, *element_index;
	
	for( element_phi = element->FirstChildElement("phi");
			 element_phi && index != phi_index;  )
	{
		element_index = element_phi->FirstChildElement("index");
		if (nb_phi > 1 && element_index == NULL) throw string("XMLreader::get_phi_ini_nifti_file() Error: index required");
		
		if (element_index == NULL) index = 1;
		else
			if ( element_index->QueryIntText(&index) ) throw string("XMLreader::get_phi_ini_nifti_file() Error: cannot read index");
		
		if (index != phi_index) element_phi = element_phi->NextSiblingElement("phi");
	}
	
	if (element_phi == NULL) throw string("XMLreader::get_phi_ini_nifti_file() Error: phi required");
	if (index != phi_index) throw string("XMLreader::get_phi_ini_nifti_file() Error: could not find required index");
	
	element = element->FirstChildElement("nifti_file");
	if (element == NULL) throw string("XMLreader::get_phi_ini_nifti_file() Error: nifti file not provided");
	
	return string( element->GetText() );
}

void XMLreader::get_initial_circles(vector< list< Point4<int> > >& centre_ini, vector< list<float> >& radius_ini)
{
	int nb_phi = get_nb_phi();
	
	int cx, cy, cz, ct;
	float r;
	list< Point4<int> > list_centre_tmp;
	list<float> list_radius_tmp;
	int index = -1;
	const XMLElement* element = doc.FirstChildElement("level_set");
	if (element == NULL) throw string("XMLreader::get_initial_circles() Error: level_set required");
	
	for( const XMLElement* element_phi = element->FirstChildElement("phi");
			 element_phi;
			 element_phi = element_phi->NextSiblingElement("phi") )
	{
		list_centre_tmp.clear();
		list_radius_tmp.clear();
		
		for ( const XMLElement* element_circle = element_phi->FirstChildElement("circle");
			 element_circle;
			 element_circle = element_circle->NextSiblingElement("circle") )
		{
			const XMLElement* element_centre = element_circle->FirstChildElement("center");
			if (element_centre == NULL) throw string("XMLreader::get_initial_circles() Error: no centre provided for the circle");
		
			const XMLElement* element_x = element_centre->FirstChildElement("x");
			if (element_x == NULL) throw string("XMLreader::get_initial_circles() Error: x element required");
			if ( element_x->QueryIntText(&cx) ) throw string("XMLreader::get_initial_circles() Error: cannot read x");
			
			const XMLElement* element_y = element_centre->FirstChildElement("y");
			if (element_y == NULL) throw string("XMLreader::get_initial_circles() Error: y element required");
			if ( element_y->QueryIntText(&cy) ) throw string("XMLreader::get_initial_circles() Error: cannot read y");
			
			const XMLElement* element_z = element_centre->FirstChildElement("z");
			if (element_z == NULL) throw string("XMLreader::get_initial_circles() Error: z element required");
			if ( element_z->QueryIntText(&cz) ) throw string("XMLreader::get_initial_circles() Error: cannot read z");
			
			const XMLElement* element_t = element_centre->FirstChildElement("t");
			if (element_t == NULL) throw string("XMLreader::get_initial_circles() Error: t element required");
			if ( element_t->QueryIntText(&ct) ) throw string("XMLreader::get_initial_circles() Error: cannot read t");
			
			const XMLElement* element_r = element_circle->FirstChildElement("radius");
			if (element_r == NULL) throw string("XMLreader::get_initial_circles() Error: no radius provided for the circle");
			if ( element_r->QueryFloatText(&r) ) throw string("XMLreader::get_initial_circles() Error: cannot read radius");
			
			list_centre_tmp.push_back( Point4<int>(cx, cy, cz, ct) );
			list_radius_tmp.push_back(r);
		}
		
		if ( list_centre_tmp.empty() ) continue;
		
		if (nb_phi == 1) index = 1;
		else
		{
			const XMLElement* element_index = element_phi->FirstChildElement("index");
			if (nb_phi > 1 && element_index == NULL) throw string("XMLreader::get_initial_circles() Error: index required");
			
			if (element_index == NULL) index = 1;
			else
				if ( element_index->QueryIntText(&index) ) throw string("XMLreader::get_initial_circles() Error: cannot read index");
			
			if (index <= 0 || index > nb_phi) throw string("XMLreader::get_initial_circles() Error: index out of bounds");
		}
		
		centre_ini[index-1] = list_centre_tmp;
		radius_ini[index-1] = list_radius_tmp;
	}
}

int XMLreader::get_nb_SFdata()
{
	int count = 0;
	const XMLElement* element = doc.FirstChildElement("speed_function");
	if (element == NULL) throw string("XMLreader::get_nb_SFdata() Error: speed_function required");
	
	for( element = element->FirstChildElement("SF_data");
			 element;
			 element = element->NextSiblingElement("SF_data") )
	{
		count++;
	}
	
	return count;
}

void XMLreader::get_SFdata_params(int SF_index, string& type_sf, vector<float>& params_sf)
{
	int nb_SF = get_nb_SFdata();
	
	SF_index++; //xml value starts at 1, C code value starts at 0
	
	float value;
	XMLElement* element = doc.FirstChildElement("speed_function");
	if (element == NULL) throw string("XMLreader::get_SFdata_params() Error: speed_function required");
	
	int index = -1;
	for( element = element->FirstChildElement("SF_data");
			 element && index != SF_index; )
	{
		XMLElement* element_index = element->FirstChildElement("index");
		if (nb_SF > 1 && element_index == NULL) throw string("XMLreader::get_SFdata_params() Error: index required");
		
		if (element_index == NULL) index = 1;
		else
			if ( element_index->QueryIntText(&index) ) throw string("XMLreader::get_SFdata_params() Error: cannot read index");
	
		if (index != SF_index) element = element->NextSiblingElement("SF_data");
	}
	
	if (element == NULL) throw string("XMLreader::get_SFdata_params() Error: SF_data required");
	if (index != SF_index) throw string("XMLreader::get_SFdata_params() Error: could not find required index");
	
	XMLElement* element_name = element->FirstChildElement("name");
	if (element_name == NULL) throw string("XMLreader::get_SFdata_params() Error: name of SF required");
	type_sf = string( element_name->GetText() );
	
	for( XMLElement* element_param = element->FirstChildElement("param");
			 element_param;
			 element_param = element_param->NextSiblingElement("param") )
	{
		if ( element_param->QueryFloatText(&value) ) throw string("XMLreader::get_SFdata_params() Error: cannot read param value");
		
		params_sf.push_back(value);
	}
}

float XMLreader::get_SFdata_weight(int SF_index)
{
	int nb_SF = get_nb_SFdata();
	
	SF_index++; //xml value starts at 1, C code value starts at 0
	
	float value;
	XMLElement* element = doc.FirstChildElement("speed_function");
	if (element == NULL) throw string("XMLreader::get_SFdata_weight() Error: speed_function required");
	
	int index = -1;
	for( element = element->FirstChildElement("SF_data");
			 element && index != SF_index; )
	{
		XMLElement* element_index = element->FirstChildElement("index");
		if (nb_SF > 1 && element_index == NULL) throw string("XMLreader::get_SFdata_weight() Error: index required");
		
		if (element_index == NULL) index = 1;
		else
			if ( element_index->QueryIntText(&index) ) throw string("XMLreader::get_SFdata_weight() Error: cannot read index");
	
		if (index != SF_index) element = element->NextSiblingElement("SF_data");
	}
	
	if (element == NULL) throw string("XMLreader::get_SFdata_weight() Error: SF_data required");
	if (index != SF_index) throw string("XMLreader::get_SFdata_weight() Error: could not find required index");
	
	element = element->FirstChildElement("weight");
	if (element == NULL) return -1;
	
	if ( element->QueryFloatText(&value) ) return -1;
	return value;
}

int XMLreader::get_nb_SFgeom()
{
	int count = 0;
	XMLElement* element = doc.FirstChildElement("speed_function");
	if (element == NULL) throw string("XMLreader::get_nb_SFgeom() Error: speed_function required");
	
	for( element = element->FirstChildElement("SF_geom");
			 element;
			 element = element->NextSiblingElement("SF_geom") )
	{
		count++;
	}
	
	return count;
}

void XMLreader::get_SFgeom_params(int SF_index, string& type_sf, vector<float>& params_sf)
{
	int nb_SF = get_nb_SFgeom();
	
	SF_index++; //xml value starts at 1, C code value starts at 0
	
	float value;
	XMLElement* element = doc.FirstChildElement("speed_function");
	if (element == NULL) throw string("XMLreader::get_SFgeom_params() Error: speed_function required");
	
	int index = -1;
	for( element = element->FirstChildElement("SF_geom");
			 element && index != SF_index; )
	{
		XMLElement* element_index = element->FirstChildElement("index");
		if (nb_SF > 1 && element_index == NULL) throw string("XMLreader::get_SFgeom_params() Error: index required");
		
		if (element_index == NULL) index = 1;
		else
			if ( element_index->QueryIntText(&index) ) throw string("XMLreader::get_SFgeom_params() Error: cannot read index");
	
		if (index != SF_index) element = element->NextSiblingElement("SF_geom");
	}
	
	if (element == NULL) throw string("XMLreader::get_SFgeom_params() Error: SF_geom required");
	if (index != SF_index) throw string("XMLreader::get_SFgeom_params() Error: could not find required index");
	
	XMLElement* element_name = element->FirstChildElement("name");
	if (element_name == NULL) throw string("XMLreader::get_SFgeom_params() Error: name of SF not provided");
	type_sf = string( element_name->GetText() );
	
	for( XMLElement* element_param = element->FirstChildElement("param");
			 element_param;
			 element_param = element_param->NextSiblingElement("param") )
	{
		if ( element_param->QueryFloatText(&value) ) throw string("XMLreader::get_SFgeom_params() Error: cannot read param value");
		
		params_sf.push_back(value);
	}
}

float XMLreader::get_SFgeom_weight(int SF_index)
{
	int nb_SF = get_nb_SFgeom();
	
	SF_index++; //xml value starts at 1, C code value starts at 0
	
	float value;
	XMLElement* element = doc.FirstChildElement("speed_function");
	if (element == NULL) throw string("XMLreader::get_SFgeom_weight() Error: speed_function required");
	
	int index = -1;
	for( element = element->FirstChildElement("SF_geom");
			 element && index != SF_index; )
	{
		XMLElement* element_index = element->FirstChildElement("index");
		if (nb_SF > 1 && element_index == NULL) throw string("XMLreader::get_SFgeom_weight() Error: index required");
		
		if (element_index == NULL) index = 1;
		else
			if ( element_index->QueryIntText(&index) ) throw string("XMLreader::get_SFgeom_weight() Error: cannot read index");
	
		if (index != SF_index) element = element->NextSiblingElement("SF_geom");
	}
	
	if (element == NULL) throw string("XMLreader::get_SFgeom_weight() Error: SF_geom required");
	if (index != SF_index) throw string("XMLreader::get_SFgeom_weight() Error: could not find required index");
	
	element = element->FirstChildElement("weight");
	if (element == NULL) return -1;
	
	if ( element->QueryFloatText(&value) ) return -1;
	return value;
}

string XMLreader::get_schemeType()
{
	XMLElement* textElement = doc.FirstChildElement("scheme");
	if (textElement == NULL) throw string("XMLreader::get_schemeType() Error: scheme required");
	
	return string( textElement->GetText() );
}

vector<float> XMLreader::get_dirac_widths()
{
	int nb_rbf = get_nb_RBFs();
	
	vector<float> tmp(nb_rbf);
	float value;
	int index;
	
	const XMLElement* element = doc.FirstChildElement("RBF");
	if (element == NULL) throw string("XMLreader::get_dirac_widths() Error: RBF required");
	
	for( element = element->FirstChildElement("psi");
			 element;
			 element = element->NextSiblingElement("psi") )
	{
		const XMLElement* element_width = element->FirstChildElement("dirac_width");
		if (element_width == NULL) value = 1;
		
		if ( element_width->QueryFloatText(&value) ) value = 1;
		
		if (nb_rbf == 1) index = 1;
		else
		{
			const XMLElement* element_index = element->FirstChildElement("index");
			if (nb_rbf > 1 && element_index == NULL) throw string("XMLreader::get_dirac_widths() Error: index required");
			
			if (element_index == NULL) index = 1;
			else
				if ( element_index->QueryIntText(&index) ) throw string("XMLreader::get_dirac_widths() Error: cannot read index");
			
			if (index <= 0 || index > nb_rbf) throw string("XMLreader::get_dirac_widths() Error: index out of bounds");
		}
		
		tmp[index-1] = value;
	}
	
	return tmp;
}

bool XMLreader::registration_required()
{
	XMLElement* element = doc.FirstChildElement("registration");
	if (element == NULL) return false;
	
	return true;
}

bool XMLreader::use_global_registration()
{
	bool value;
	XMLElement* element = doc.FirstChildElement("registration");
	if (element == NULL) throw string("XMLreader::use_global_registration() Error: registration required");
	
	element = element->FirstChildElement("global_variant");
	if (element == NULL) return false;
	
	if ( element->QueryBoolText(&value) ) return false;
	return value;
}

bool XMLreader::do_selection_of_connected_regions()
{
	bool value;
	XMLElement* element = doc.FirstChildElement("registration");
	if (element == NULL) throw string("XMLreader::do_selection_of_connected_regions() Error: registration required");
	
	element = element->FirstChildElement("selection_of_connected_regions");
	if (element == NULL) return false;
	
	if ( element->QueryBoolText(&value) ) return false;
	return value;
}

int XMLreader::do_slice_wise_registration()
{
	bool value;
	XMLElement* element = doc.FirstChildElement("registration");
	if (element == NULL) throw string("XMLreader::do_slice_wise_registration() Error: registration required");
	
	element = element->FirstChildElement("slice_wise");
	if (element == NULL) return -1;
	
	if ( element->QueryBoolText(&value) ) return -1;
	if (value) return 1;
	else return 0;
}

bool XMLreader::shifts_in_SA_planes_at_beginning()
{
	bool value;
	XMLElement* element = doc.FirstChildElement("registration");
	if (element == NULL) throw string("XMLreader::shifts_in_SA_planes_at_beginning() Error: registration required");
	
	element = element->FirstChildElement("shifts_in_SA_planes_at_beginning");
	if (element == NULL) return false;
	
	if ( element->QueryBoolText(&value) ) return false;
	return value;
}

int XMLreader::wait_before_rotation()
{
	int value;
	XMLElement* element = doc.FirstChildElement("registration");
	if (element == NULL) throw string("XMLreader::wait_before_rotation() Error: registration required");
	
	element = element->FirstChildElement("wait_before_rotation");
	if (element == NULL) return 0;
	
	if ( element->QueryIntText(&value) ) return 0;
	return value;
}

int XMLreader::allow_display()
{
	int value;
	XMLElement* element = doc.FirstChildElement("display");
	if (element == NULL) return false;
	
	if ( element->QueryIntText(&value) ) return false;
	return value;
}

int XMLreader::get_nb_iterations(int step_index)
{
	step_index++; //xml value starts at 1, C code value starts at 0
	
	int value;
	XMLElement* element = doc.FirstChildElement("stop_condition");
	if (element == NULL) throw string("XMLreader::get_nb_iterations() Error: stop_condition required");
	
	int index = -1;
	for( element = element->FirstChildElement("step");
			 element && index != step_index; )
	{
		XMLElement* element_index = element->FirstChildElement("index");
		if (step_index > 1 && element_index == NULL) throw string("XMLreader::get_nb_iterations() Error: index required");
		
		if (element_index == NULL) index = 1;
		else
			if ( element_index->QueryIntText(&index) ) throw string("XMLreader::get_nb_iterations() Error: cannot read index");
	
		if (index != step_index) element = element->NextSiblingElement("step");
	}
	
	if (element == NULL) throw string("XMLreader::get_nb_iterations() Error: step required");
	if (index != step_index) throw string("XMLreader::get_nb_iterations() Error: could not find required step index");
	
	element = element->FirstChildElement("nb_iterations");
	if (element == NULL) return 0;
	
	if ( element->QueryIntText(&value) ) throw string("XMLreader::get_nb_iterations() Error: cannot read nb_iterations");
	return value;
}

int XMLreader::get_max_stable_iterations(int step)
{
	string scheme = get_schemeType();
	
	step++; //xml value starts at 1, C code value starts at 0
	
	int value;
	XMLElement* element = doc.FirstChildElement("stop_condition");
	if (element == NULL) throw string("XMLreader::get_max_stable_iterations() Error: stop_condition required");
	
	int index = -1;
	for( element = element->FirstChildElement("step");
			 element && index != step; )
	{
		XMLElement* element_index = element->FirstChildElement("index");
		if (step > 1 && element_index == NULL) throw string("XMLreader::get_max_stable_iterations() Error: index required");
		
		if (element_index == NULL) index = 1;
		else
			if ( element_index->QueryIntText(&index) ) throw string("XMLreader::get_max_stable_iterations() Error: cannot read index");
	
		if (index != step) element = element->NextSiblingElement("step");
	}
	
	if (element == NULL) throw string("XMLreader::get_max_stable_iterations() Error: step required");
	if (index != step) throw string("XMLreader::get_max_stable_iterations() Error: could not find required step index");
	
	element = element->FirstChildElement("nb_stable_iterations");
	if (element == NULL) throw string("XMLreader::get_max_stable_iterations() Error: must provide either nb_iterations or nb_stable_iterations for each step");
	
	if ( element->QueryIntText(&value) ) throw string("XMLreader::get_max_stable_iterations() Error: must provide either nb_iterations or nb_stable_iterations for each step");
	return value;
}

