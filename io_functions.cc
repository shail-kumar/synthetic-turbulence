#include "real_field.h"
#include "json.hpp"
#include "h5si.h"
#include <hdf5.h>

using json = nlohmann::json;

using namespace blitz;
// using namespace std;

void read_parameters()
{
	ifstream paraFile("parameters.json");		// Load the parameter file 
	json container;								// Declare the json object container
	paraFile >> container;						// Parse the .json file and put the content in the container
	// cout<<container<<endl;
	N = container["shell_model"]["number_of_shells-N"];
	lambda = container["shell_model"]["shell_ratio-lambda"];
	string tmp= container["shell_model"]["viscosity-nu"];
	nu = stod(tmp);

	length_scales = container["real_space"]["number_of_length_scales"];
	gridSize = container["real_space"]["grid_size-gridSize"];

	input_path = container["read_write"]["input_dir"];
	set_path = container["read_write"]["output_dir"];
	Attempt = container["read_write"]["attempt_number"];
	starting_plane = container["read_write"]["starting_plane_number"];
}

void create_dir()
{	
	int dir;
	stringstream dir_field;
	dir_field<<set_path+"dissRate-"+to_string(Attempt);
	dir = mkdir(dir_field.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (-1==dir)
		{
			cout<<dir_field.str().c_str()<<" exists beforehand."<<endl;
		}
	else {cout<<dir_field.str().c_str()<<" created."<<endl;}
}

void write_dissRate(Array<double,3> *epsilon)
{
	h5::init();
	// -----------------------
	TinyVector<int,3> S;
	S = (*epsilon).shape();
	// cout<<S<<endl;
	string filename = set_path+"dissRate-"+to_string(Attempt)+".h5";
	h5::File f(filename, "w");

	h5::Dataset ds = f.create_dataset("eps", h5::shape(S(0),S(1),S(2)), "double");

	ds << (*epsilon).data();

	// -----------------------

	h5::finalize();
}

void write_dissRate_plane(Array<double,2> *epsilon, int plane)
{
	h5::init();
	// -----------------------
	TinyVector<int,2> S;
	S = (*epsilon).shape();
	// cout<<S<<endl;
	string filename = set_path+"dissRate-"+to_string(Attempt)+"/x-"+to_string(plane)+".h5";
	h5::File f(filename, "w");

	// string ds_name = "eps_x-"+to_string(plane);
	string ds_name = "eps";
	h5::Dataset ds = f.create_dataset(ds_name, h5::shape(S(0),S(1)), "double");

	ds << (*epsilon).data();

	// f.close();
	// -----------------------

	h5::finalize();
}