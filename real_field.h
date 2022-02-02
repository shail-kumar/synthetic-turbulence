#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <string>
#include <ctime>
#include <cstdlib>
#include <blitz/array.h>
#include <random/uniform.h>
#include <time.h>
#include <sys/stat.h>
#include <mpi.h>
using namespace blitz;

typedef complex<double> cdouble;
typedef TinyVector<double,3> tv;
typedef TinyVector<cdouble,3> ctv;
typedef Array<ctv,1> Actv;

extern string input_path;
extern string set_path;
extern int Attempt;

extern const cdouble  plusI;
extern ofstream real_field_file; 
extern int N;				// # of shells
extern double nu;

extern int P; 				// # of realizations of e and c
extern double ld;			// descretization ratio of real space for length scales.
extern double xd;			// descretization ratio of real space.
extern int space_points;	// # of space points
extern double lambda;
extern int length_scales;
extern int gridSize;
extern int starting_plane;

inline double mod(tv a){return sqrt(a(0)*a(0) + a(1)*a(1) + a(2)*a(2));}
inline cdouble cmod(ctv a){return sqrt(a(0)*a(0) + a(1)*a(1) + a(2)*a(2));}

//------- Declaration of functions --------------------------------------------

void read_parameters();
void create_dir();
Array<double,1> shell_wavenumbers(Array<double,1>*, int, double);
double real_random(unsigned int);
Array<tv,1> e_tuples(int, unsigned int );
tv unit_vector(unsigned int);
Array<tv,1> c_tuples(Array<tv,1> , int , unsigned int);
// Array<ctv,1> h_tuples(Array<tv,1> , int , unsigned int);
pair<Actv, Actv> h_tuples(Array<tv,1> , int , unsigned int);
tv real_velocity(Array<cdouble,1>*,Array<double,1>*, Array<tv,1>, Array<tv,1>, tv);
void del_u_parallel(Array<cdouble,1>*, Array<double,1>*, int);
void del_u_parallel_simplified(Array<cdouble,1>*, Array<double,1>*, Array<tv,1> , Array<tv,1> , int, int);
void helical_del_u_para(Array<cdouble,1> *u, Array<cdouble,1> *v, Array<double,1> *kn);
tv helical_real_velocity(Array<cdouble,1>*,Array<cdouble,1>*, Array<double,1>*, Array<tv,1>,pair<Actv, Actv>, tv);

double norm(tv );
int number_of_descrete_values(double, double, double);
Array<double,1> descretization(double, double, double);
void distance_structure_function(Array<cdouble,1>*, Array<double,1>*, Array<tv,1>, Array<tv,1>, int, int);

void pdf_of_dsf(Array<cdouble,1>*, Array<double,1>*, Array<tv,1>, Array<tv,1>, int);

Array<double,2> structure_function_helical(Array<cdouble,1>*, Array<cdouble,1>*, Array<double,1>*, Array<tv,1>, pair<Actv, Actv>, int, int, int);
void structure_function_helical_final();

Array<double,2> derivative_matrix(Array<cdouble,1>*,Array<cdouble,1>*, Array<double,1>*, Array<tv,1> , pair<Actv, Actv> , tv );
Array<double,2> derivative_matrix(Array<cdouble,1>*, Array<double,1>*, Array<tv,1> , Array<tv,1> , tv );

double local_dissipation_rate(Array<double, 2>, double);
void eps(int, int);

void write_dissRate(Array<double,3>*);
void write_dissRate_plane(Array<double,2>*, int);
