/*	This program constructs real velocity fields from the velocity represenatations obtained from shell models for incompressible flow.
	There are functions included in this program to explore the statistcs of velocity increments.

	Author: Shailendra Kumar Rathor 
	Email: skrathor@iitk.ac.in, shailkumar22@gmail.com

	This work is licensed under Creative Commons Attribution-NonCommercial 4.0 International License. If you find this work useful, please use it.
	For more information, please read the README.txt .
*/

#include "real_field.h"
using namespace blitz;
using namespace ranlib;

// ------------------- Global constructs --------------------------------------
typedef complex<double> cdouble;
typedef TinyVector<double,3> tv;
typedef TinyVector<cdouble,3> ctv;
typedef Array<ctv,1> Actv;
const cdouble  plusI = cdouble(0,1.0);			// plus I

string input_path;
string set_path;
int Attempt;
fstream field_file; ofstream real_field_file; ifstream kn_file;

int N = 24;										// # of shells
double lambda =  (1+sqrt(5))/2.0;
int P = 1;										// # of realizations of e and c.
double ld = 1.02;								// descretization ratio of real space for length scales.
int length_scales = 50;						// Number of length scales
int gridSize = 2;							// grid size for dissipation rate

int space_points = 10;						// # of space points
double xd = 1.02;								// descretization ratio of real space.
double nu;
int starting_plane = 0;
// ----------------------------------------------------------------------------

int main(int argc, char const *argv[])
{
	
	int process_id, num_procs;
	MPI_Status status;
	
	MPI_Init(NULL,NULL);
	
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
	
	read_parameters();		// Read parameters
	int procs_per_pid = gridSize/num_procs;

	if(process_id == 0)
	{
		
		create_dir();			// Create directory for dissipation rate
		eps(process_id, procs_per_pid);		// Generate dissipation rate
	
	}
	else
	{
		eps(process_id, procs_per_pid);		// Generate dissipation rate

	}
	
	
	MPI_Finalize();
	
	//~ read_parameters();		// Read parameters
	//~ create_dir();			// Create directory for dissipation rate
	//~ eps();					// Generate dissipation rate
	//~ cout<<"Task completed!"<<endl;
	return 0;
}
