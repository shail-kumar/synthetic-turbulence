#include "real_field.h"
using namespace blitz;
using namespace ranlib;

Array<double,1> shell_wavenumbers(Array<double,1>* k_n, int N, double lambda)
{
	double k_0 = 1.0/16;
	for (int i = 0; i < N+2; ++i)
	{
		(*k_n)(i)=k_0*pow(lambda,i-1);			//i=2 corressponds to 1st shell
		// cout<<(*k_n)(i)<<endl;
	}
	// cout<<1<<(*k_n)<<endl;
	return (*k_n);
}

double real_random(unsigned int seedint) // Generates a real random number
{
	Uniform<double> x;
	x.seed((unsigned int)time(0) + seedint);
	// cout<<x.random()<<endl;
	return x.random();
}

tv unit_vector(unsigned int seedint) // Unit vector generated randomly
{
	// unsigned int seedint = 0;
	Uniform<double> x;
	x.seed((unsigned int)time(0) + seedint);
	// seedint = 1;
	Uniform<double> y;
	y.seed((unsigned int)time(0) + seedint);
	// seedint = 2;
	Uniform<double> z;
	z.seed((unsigned int)time(0) + seedint);
	double a, b, c;
	tv e;
	a = x.random() * 2 - 1;
	b = y.random() * 2 - 1;
	c = z.random() * 2 - 1;
	e = a,b,c;
	e = e/mod(e);
	// cout<<e<<"\t"<<mod(e)<<endl;
	return e;
}

Array<tv,1> e_tuples(int N, unsigned int seed) // Array of 3-tuples of size N i.e. e_n's
{
	// unsigned int seedint = 0;
	Uniform<double> x;
	x.seed((unsigned int)time(0) + seed);
	// seedint = 1;
	Uniform<double> y;
	y.seed((unsigned int)time(0) + seed);
	// seedint = 2;
	Uniform<double> z;
	z.seed((unsigned int)time(0) + seed);

	Array<double,1> A(N);
	Array<double,1> B(N);
	Array<double,1> C(N);
	Array<tv,1> e(N); tv temp;
	int n=0;
	while(n < N){
		A(n) = x.random() * 2 - 1;
		B(n) = y.random() * 2 - 1;
		C(n) = z.random() * 2 - 1;
		temp = A(n), B(n), C(n);
		e(n) = temp/mod(temp);

		// cout<<e(n)<<"\t"<<mod(e(n))<<endl;
		n++;
	}
	// cout<<e;
	return e;
}

Array<tv,1> c_tuples(Array<tv,1> e, int N, unsigned int seed) // Array of 3-tuples of size N i.e. Cn's
{
	// unsigned int seedint = 0;
	Uniform<double> x;
	x.seed((unsigned int)time(0) + seed);

	Array<double,1> A(N);
	Array<tv,1> c(N);
	int i = 0;
	while(i<N){
		A(i) = x.random() * 2 - 1;
		double e1 = e(i)(0); double e2 = e(i)(1); double e3 = e(i)(2);
		double k1 = 1-pow(A(i),2);
		double k2 = -A(i)*e1;
		double c2, c3=0;
		double denom = 1-pow(e1,2);
		double det = pow(e3*k2, 2) - denom*(k2*k2 - k1*pow(e2,2));
		if (det > 0 && denom > 0.0001)
		{
			c3 = (e3*k2 + sqrt(det))/denom ;
			if (abs(c3) > 1)
			{
				c3 = (e3*k2 - sqrt(det))/denom ;			
			}
			c2 = (k2-c3*e3)/e2;
			c(i) = A(i), c2, c3;
			i++;
		}
	}
	return c;
}

tv real_velocity(Array<cdouble,1> *u, Array<double,1> *kn, Array<tv,1> e, Array<tv,1> c, tv x)
{
	tv U; U = 0; 	
	double phase;
	int N = (*u).numElements() - 4;
	for (int n = 0; n < N; ++n)
	{
		phase = (*kn)(n+2)*dot(e(n), x);
		U += c(n)*real((*u)(n+2)*exp(plusI*phase)) * 2;
	}
	return U;
}


void del_u_parallel(Array<cdouble,1> *u, Array<double,1> *kn, int T)
{
	double l_0 = 2*M_PI/(*kn)(N+1); 						// Define vector l = l_0 * (1.3)^(5m-1) * unit_l; m = 1,2,3,....,12 . 
	double l_max = 1.1;// 2*M_PI/(*kn)(2); 					// Characteristic length scale of the system.
	double x_max = 2*M_PI/(*kn)(2);
	int M = log(l_max/l_0)/log(ld);							// The length scale is descretised as (ld)**i .	
	int J_min = log(l_0)/log(xd);
	int J_max = log(l_max)/log(xd);
	// cout<<J_min<<endl;
	// cout<<J_max<<endl;
	for (int p = 0; p < P; ++p)								// For averaging over realizations of e and c.
	{
		int dir;
		stringstream vd;
		vd<<set_path<<"velo_disto_"<<T;
		dir = mkdir(vd.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (-1==dir)
			{
				cout<<vd.str().c_str()<<" exists beforehand."<<endl; 
				cout<<"Overwriting avoided. All processes exit..\n";
				exit(0);
			}
		else {cout<<vd.str().c_str()<<" created."<<endl;}

		double l_mod = 0;
		double x_mod = 0;
		tv unit_l; unit_l = 0;
		tv unit_x; unit_x = 0;
		tv l; l = 0;
		tv x; x = 0;
		tv du; du = 0;
		double du_parallel = 0;
		Array<tv,1> e(N); e = e_tuples(N, p);
		Array<tv,1> c(N); c = c_tuples(e,N, p);
		// -----------------------------------------------------------------------------------------------------------
		for (int m = 0; m <= M; ++m)
		{
			cout<<T<<"\t"<<m<<endl;
			stringstream file_name;
			file_name<<set_path<<"velo_disto_"<<T<<"/m_"<<m<<".d";		// vd --> velocity distribution
			real_field_file.open(file_name.str().c_str());		// Open file to write the OUTPUT.
			unit_l = unit_vector(m);
			l_mod = l_0*pow(ld,m);							// Settle this : how many shells etc..
			l = l_mod*unit_l;
			real_field_file<<"# l \n"<<l_mod<<endl;
			real_field_file<<"# du_parallel"<<endl;
			for (int i = J_min; i < J_max; ++i)						// i creates shells for the magnitude of x vector. Make the shells more fine in radius.
			{
				for (int j = 0; j < space_points; ++j)					// Tag of x vector
				{
					x_mod = x_max*real_random(j)*pow(xd, i);	// Change the 10 to 2 for making the shells finer.
					unit_x = unit_vector(j);
					x = x_mod*unit_x;
					if (mod(x+l) <= x_max)
					{
						// cout<<i<<"\t"<<j<<"\t"<<x<<endl;
						du = real_velocity(u, kn, e, c, x+l) - real_velocity(u, kn, e, c, x);
						du_parallel = dot(du, unit_l);
					}
					else if (mod(l-x) <= x_max)
					{
						// tv y = -x; cout<<i<<"\t"<<j<<"\t"<<y<<endl;
						du = real_velocity(u, kn, e, c, -x+l) - real_velocity(u, kn, e, c, -x);
						du_parallel = dot(du, unit_l);
					}
					// else{cout<<"No location of interest found !!"<<endl;}
					real_field_file<<du_parallel<<"\n";
				}
			}
			// cout<<real_velocity(u, kn, e, c, x)<<endl;
			real_field_file.close();
		}
		// ------------------------------------------------------------------------------------------------------------
	}
}

void del_u_parallel_simplified(Array<cdouble,1> *u, Array<double,1> *kn, Array<tv,1> e, Array<tv,1> c, int T, int ic)
{
	double l_0 = 1e-6; //2*M_PI/(*kn)(N+1); 						// Define vector l = l_0 * (1.3)^(5m-1) * unit_l; m = 1,2,3,....,12 . 
	double l_max = 2e2;// 2*M_PI/(*kn)(2); 					// Characteristic length scale of the system.
	double ld = 1.05;
	int M = 50;
	ld = pow(l_max/l_0, 1.0/M);							// The length scale is descretised as (ld)**i .
	// cout<<"M: "<<M<<endl;	

	/*int dir;
	stringstream vd;
	vd<<set_path<<"velo_disto_"<<T;
	dir = mkdir(vd.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (-1==dir)
		{
			cout<<vd.str().c_str()<<" exists beforehand."<<endl; 
			cout<<"Overwriting avoided. All processes exit..\n";
			exit(0);
		}
	else {cout<<vd.str().c_str()<<" created."<<endl;}*/

	double l_mod = 0;
	tv unit_l; unit_l = 0;
	tv l; l = 0;
	tv x; x = 0;
	tv du; du = 0;
	tv u0; u0 = real_velocity(u, kn, e, c, x);
	// Array<tv,1> e(N); e = e_tuples(N, T);
	// Array<tv,1> c(N); c = c_tuples(e,N, T);
	// -----------------------------------------------------------------------------------------------------------
	stringstream file_name;
	file_name<<set_path<<"realization_full-"<<ic<<"/du_parallel-"<<T<<".d";		// vd --> velocity distribution
	real_field_file.open(file_name.str().c_str());		// Open file to write the OUTPUT.
	real_field_file<<"# l_mod \t du_parallel"<<endl;
	// cout<<"# l_mod \t du_parallel"<<endl;
	int Nl = 100;
	for (int m = 0; m <= M; ++m)
	{
		l_mod = l_0*pow(ld,m);							// Settle this : how many shells etc..
		double du_parallel = 0;

		for (int il = 0; il < Nl; ++il)
		{
			unit_l = unit_vector(m);
			l = l_mod*unit_l;
			du = real_velocity(u, kn, e, c, l) - u0; 		// l is actually x+l since x = 0.
			du_parallel += abs(dot(du, unit_l));
		}
		real_field_file<<l_mod<<"\t"<<du_parallel/Nl<<"\n";
		// cout<<l_mod<<"\t"<<du_parallel/Nl<<"\n";
	}
	real_field_file.close();
	// ------------------------------------------------------------------------------------------------------------
}

int number_of_descrete_values(double min, double max, double ratio)
{
	return log(max/min) / log(ratio) + 1; // 1 is added to include min as loop index starts from 0.
}

Array<double,1> descretization(double min, double max, double ratio)
{
	// int N = number_of_descrete_values(min, max, ratio);		// Size of the set of descrete values 
	int N = log(max/min) / log(ratio) + 1;
	Array<double,1> S(N);					// S --> set of descrete values
	for (int i = 0; i < N; ++i)
	{
		S(i) = min * pow(ratio, i);
		// cout<<S(i)<<endl;
	}

	return S;
}

double norm(tv x)
{
	return pow(pow(x(0),2) + pow(x(1),2) + pow(x(2),2), 0.5);
}

void distance_structure_function(Array<cdouble,1> *u, Array<double,1> *kn, Array<tv,1> e, Array<tv,1> c, int t, int ic)
{
	stringstream file_name;
	file_name<<set_path<<"dsf-"<<ic<<"/du_vs_l-"<<t<<".d";		// vd --> velocity distribution
	real_field_file.open(file_name.str().c_str());		// Open file to write the OUTPUT.
	real_field_file<<"# du_mod \t l_mod"<<endl;
	// int N = number_of_descrete_values(1e-6, 2e2, 10.0);
	// cout<<N<<endl;
	double lambda = 1.02;
	Array<double, 1> l_fixed = descretization(1e-6, 2e2, lambda);
	Array<double, 1> du_fixed = descretization(1e-3, 1e0, lambda);
	// cout<<l<<endl;
	/*cout<<l_fixed.size()<<endl;
	cout<<du_fixed.size()<<endl;*/
	
	tv l; l = 0;
	tv x; x = 0;
	tv u0; u0 = real_velocity(u, kn, e, c, x);

	tv du; du = 0;
	double computed_du = 0;
	tv unit_l; unit_l = 0;
	int Nl = 10;
	for (unsigned int j = 0; j < du_fixed.size(); ++j)
	{
		// cout<<j<<endl;
		double computed_l = 0;
		for (int il = 0; il < Nl; ++il)
		 {

			unit_l = unit_vector(il);

			for (unsigned int i = 0; i < l_fixed.size(); ++i)
			{	
				l = l_fixed(i)*unit_l;
				du = real_velocity(u, kn, e, c, l) - u0; 		// l is actually x+l since x = 0.
				computed_du = abs(dot(du, unit_l));
				// computed_du = norm(du);
				if (computed_du > du_fixed(j))
				{
					computed_l += l_fixed(i);
					// real_field_file<< du_fixed(j) << "\t" << l_fixed(i)<<endl;
					break;
				}
			}
		}
		real_field_file<< du_fixed(j) << "\t" << computed_l/Nl <<endl;

	}

	real_field_file.close();
}

void pdf_of_dsf(Array<cdouble,1> *u, Array<double,1> *kn, Array<tv,1> e, Array<tv,1> c, int t)
{
	stringstream file_name;
	file_name<<set_path<<"dsf-1/l_for_pdf_large_norm.d";		// 
	real_field_file.open(file_name.str().c_str(), ios::app);		// Open file to write the OUTPUT.
	//~ real_field_file<<"# du_mod \t l_mod"<<t<<endl;
	// int N = number_of_descrete_values(1e-6, 2e2, 10.0);
	// cout<<N<<endl;
	double lambda = 1.02;
	Array<double, 1> l_fixed = descretization(1e-6, 2e2, lambda);
	Array<double, 1> du_fixed = descretization(0.201, 0.202, lambda);
	//~ cout<<l_fixed.size()<<endl;
	//~ cout<<du_fixed.size()<<endl;
	
	tv l; l = 0;
	tv x; x = 0;
	tv u0; u0 = real_velocity(u, kn, e, c, x);

	tv du; du = 0;
	double computed_du = 0;
	tv unit_l; unit_l = 0;
	int Nl = 1;
	for (unsigned int j = 0; j < du_fixed.size(); ++j)
	{
		//~ cout<<j<<endl;
		double computed_l = 0;
		for (int il = 0; il < Nl; ++il)
		 {

			unit_l = unit_vector(il);

			for (unsigned int i = 0; i < l_fixed.size(); ++i)
			{	
				l = l_fixed(i)*unit_l;
				du = real_velocity(u, kn, e, c, l) - u0; 		// l is actually x+l since x = 0.
				//~ computed_du = abs(dot(du, unit_l));
				computed_du = norm(du);
				if (computed_du > du_fixed(j))
				{
					computed_l += l_fixed(i);
					// real_field_file<< du_fixed(j) << "\t" << l_fixed(i)<<endl;
					break;
				}
			}
		}
		real_field_file << computed_l/Nl <<endl;

	}

	real_field_file.close();
}


pair<Actv, Actv> h_tuples(Array<tv,1> e, int N, unsigned int seed) // Random Number Tuples
{
	// unsigned int seedint = 0;
	Uniform<double> x;
	x.seed((unsigned int)time(0) + seed);
	// int M  = 2*N;
	// Array<ctv,1> h(M);
	Array<ctv,1> h_plus(N);
	Array<ctv,1> h_minus(N);
	tv e1, e2;
	// cout<<h(2)(1).real()<<endl;
	int i = 0;
	double k1 = 0;
	double k2 = 0;
	double k3 = 0;
	double norm_const = 0;
	while(i<N){
		k1 = e(i)(0); k2 = e(i)(1); k3 = e(i)(2);
		norm_const = sqrt(1 - pow(k3,2));
		e1 = k2/norm_const, -k1/norm_const, 0;
		e2 = k1*k3/norm_const, k2*k3/norm_const, -norm_const;
		// cout<<"e1: "<<e1<<endl;
		// cout<<"e2: "<<e2<<endl;
		for (int j = 0; j < 3; ++j)
		{
			// Without -std=c++11 switch
			/*h_plus(i)(j).real() = e2(j);
			h_plus(i)(j).imag() = -e1(j);
			h_minus(i)(j).real() = e2(j);
			h_minus(i)(j).imag() = e1(j);*/
			
			// With -std=c++11 switch
			h_plus(i)(j).real(e2(j));
			h_plus(i)(j).imag(-e1(j));
			h_minus(i)(j).real(e2(j));
			h_minus(i)(j).imag(e1(j));
		}
		i++;
	}
	// cout<<h_plus<<endl;
	// cout<<h_minus<<endl;
	return make_pair(h_plus, h_minus);
}

tv helical_real_velocity(Array<cdouble,1> *u,Array<cdouble,1> *v, Array<double,1> *kn, Array<tv,1> e, pair<Actv, Actv> h, tv x)
{
	ctv U_temp; 
	tv U; U = 0,0,0;
	Array<ctv,1> h_plus(N);
	Array<ctv,1> h_minus(N);
	// pair<Actv, Actv> h = h_tuples(e, N, 1);
	h_plus = h.first;
	h_minus = h.second;	
	double phase;
	int N = (*u).numElements() - 4;
	// cout<<"loop: "<<(*u)<<endl;
	for (int n = 0; n < N; ++n)
	{
		phase = (*kn)(n+2)*dot(e(n), x);
		U_temp += (h_plus(n)*(*u)(n+2) + h_minus(n)*(*v)(n+2))*exp(plusI*phase);
	}
	// cout<<U_temp<<endl;
	for (int i = 0; i < 3; ++i)
	{
		U(i)= 2 * U_temp(i).real();
	}
	// cout<<U<<endl;
	return U;
}

void helical_del_u_para(Array<cdouble,1> *u, Array<cdouble,1> *v, Array<double,1> *kn)
{
	double l_0 = 2*M_PI/(*kn)(N+1); 						// Define vector l = l_0 * (1.3)^(5m-1) * unit_l; m = 1,2,3,....,12 . 
	double l_max = 1.1;// 2*M_PI/(*kn)(2); 					// Characteristic length scale of the system.
	double x_max = 2*M_PI/(*kn)(2);
	int M = log(l_max/l_0)/log(ld) ;						// The length scale is descretised as (ld)**i .	

	for (int p = 0; p < P; ++p)								// For averaging over realizations of e and c.
	{
		int dir;
		stringstream vd;
		vd<<set_path<<"velo_disto_"<<p;
		dir = mkdir(vd.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (-1==dir)
			{
				cout<<vd.str().c_str()<<" exists beforehand."<<endl; 
				cout<<"Overwriting avoided. All processes exit..\n";
				exit(0);
			}
		else {cout<<vd.str().c_str()<<" created."<<endl;}

		double l_mod = 0;
		double x_mod = 0;
		tv unit_l; unit_l = 0;
		tv unit_x; unit_x = 0;
		tv l; l = 0;
		tv x; x = 0;
		tv du; 
		//du = 0;
		double du_parallel = 0;
		Array<tv,1> e(N); e = e_tuples(N, p);
		pair<Actv, Actv> h = h_tuples(e, N, 1);

		// Array<tv,1> c(N); e = e_tuples(N, p);
		// cout<<p<<"\t"<<e<<endl;
		// cout<<c<<endl;
		// To do: Loop over e anc c for statistical calculations.
		// -----------------------------------------------------------------------------------------------------------
		for (int m = 0; m <= M; ++m)
		{
			cout<<p<<"\t"<<m<<"\n";

			stringstream file_name;
			file_name<<set_path<<"velo_disto_"<<p<<"/m_"<<m<<".d";
			real_field_file.open(file_name.str().c_str());		// Open file to write the OUTPUT.
			unit_l = unit_vector(m);
			l_mod = l_0*pow(ld,m);							// Settle this : how many shells etc..
			l = l_mod*unit_l;
			real_field_file<<"# l \n"<<l_mod<<endl;
			real_field_file<<"# du_parallel"<<endl;
			for (int i = -5; i < 1; ++i)						// Range of i ?? i creates shells for the magnitude of x vector. Make the shells finer in radius.
			{
				for (int j = 0; j < space_points; ++j)					// Tag of x vector
				{
					x_mod = x_max*real_random(j)*pow(10, i);	// Change the 10 to 2 for making the shell spacing finer.
					unit_x = unit_vector(j);
					x = x_mod*unit_x;
					if (mod(x+l) <= x_max)
					{
						// cout<<i<<"\t"<<j<<"\t"<<x<<endl;
						du = helical_real_velocity(u,v,kn, e, h, x+l) - helical_real_velocity(u, v, kn, e, h, x);
						du_parallel = dot(du, unit_l);
					}
					else if (mod(l-x) <= x_max)
					{
						// tv y = -x; cout<<i<<"\t"<<j<<"\t"<<y<<endl;
						du = helical_real_velocity(u,v,kn, e, h, -x+l) - helical_real_velocity(u, v, kn, e, h, -x);
						du_parallel = dot(du, unit_l);
					}
					// else{cout<<"No location of interest found !!"<<endl;}
					real_field_file<<du_parallel<<"\n";
				}
			}
			// cout<<real_velocity(u, kn, e, c, x)<<endl;
			real_field_file.close();
		}
		// ------------------------------------------------------------------------------------------------------------
	}
}

Array<double,2> structure_function_helical(Array<cdouble,1> *u, Array<cdouble,1> *v, Array<double,1> *kn, Array<tv,1> e, pair<Actv, Actv> h, int T, int ic, int columns)
{
	double l_0 = 1e-6; //2*M_PI/(*kn)(N+1); 				// Define vector l = l_0 * (1.3)^(5m-1) * unit_l; m = 1,2,3,....,12 . 
	double l_max = 2e2;// 2*M_PI/(*kn)(2); 					// Characteristic length scale of the system.
	double ld = 1;
	int M = length_scales;
	ld = pow(l_max/l_0, 1.0/M);								// The length scale is descretised as (ld)**i .
	// cout<<"M: "<<M<<endl;	

	double l_mod = 0;
	tv unit_l; unit_l = 0;
	tv l; l = 0;
	tv x; x = 0;
	tv du; du = 0;
	tv u0; u0 = helical_real_velocity(u,v, kn, e, h, x);
	// Array<tv,1> e(N); e = e_tuples(N, T);
	// Array<tv,1> c(N); c = c_tuples(e,N, T);
	Array<double,2> *S;
	S = new Array<double,2>(M+1, columns);
	(*S) = 0;
	// -----------------------------------------------------------------------------------------------------------
	int Nl = 1;														// For averaging over directions of vector l
	for (int m = 0; m <= M; ++m)
	{
		l_mod = l_0*pow(ld,m);
		double du_norm = 0;

		for (int il = 0; il < Nl; ++il)
		{
			unit_l = unit_vector(m);
			l = l_mod*unit_l;
			du = helical_real_velocity(u, v, kn, e, h, l) - u0; 		// l is actually x+l since x = 0.
			// du_norm = abs(dot(du, unit_l));
			du_norm = norm(du);
			// cout<<l_mod<<"\t"<<du_norm<<endl;
			for (int q = 1; q < columns; ++q)
			{
				(*S)(m,q) += pow(du_norm,q);
				// cout<<(*S)(m,q)<<endl;
			} 
		}
		(*S)(m,Range(1,toEnd)) = (*S)(m,Range(1,toEnd))/Nl;
		(*S)(m,0) = l_mod;
		// real_field_file<<l_mod<<"\t"<<du_parallel/Nl<<"\n";			// WRITE AN ARRAY INSTEAD OF SAVING IT IN A FILE
		// cout<<l_mod<<"\t"<<du_parallel/Nl<<"\n";
	}
	// cout<<(*S)<<endl;

	// real_field_file.close();
	return (*S);
	// ------------------------------------------------------------------------------------------------------------
}


void structure_function_helical_final()
{
	int Attempt;
	cout<<"Attempt? ";
	cin>>Attempt;
	stringstream u_vs_t;
	u_vs_t<<input_path<<"velocity_field_"<<Attempt<<"/u_vs_t_1.d";
	fstream field_file;
	//============= INPUT ================
	Array<double,1> *kn;	//,FortranArray<1>());			// wave number for shell n
	kn=new Array<double,1>(N+2);
	(*kn) = shell_wavenumbers(kn, N, lambda);

	double t=0;
	Array<cdouble,1> *u;									// Field u 
	u = new Array<cdouble,1>(N+4); 
	Array<cdouble,1> *v;									// Field v 
	v= new Array<cdouble,1>(N+4); 

	//============= OUTPUT ==============
	int max_order_of_moments = 6;
	int cols = max_order_of_moments +1; // 1 is added for length scale column;
	Array<double,2> *S;
	S = new Array<double,2>(length_scales+1, cols);
	// cout<<(*S)<<endl;
	// ========== Main computation begins from here =====
	int max_ic = 1;
	int max_it = 10000;
	for (int ic = 0; ic < max_ic; ic++)
	{
		field_file.open(u_vs_t.str().c_str());
		Array<tv,1> e(N);
		e = e_tuples(N,ic);
		pair<Actv, Actv> h = h_tuples(e,N,ic);

		for (int it = 0; it < max_it; ++it)
		{
			// cout<<"==========="<<it<<"========"<<endl;
			cout<<ic<<"\t"<< it <<endl;
			field_file>>t>>(*u)>>(*v);
			(*S) += structure_function_helical(u, v, kn, e, h, it, ic, cols);
			// cout<<(*S)<<endl;


		}
		field_file.close();
	}
	(*S) = (*S)/(max_ic * max_it);

	// ********* Write the Output ***********************
	stringstream sf;
	sf<<set_path<<"sf_norm-"<<Attempt<<"_t-"<<max_it<<".d";
	ofstream sf_file;
	sf_file.open(sf.str().c_str());
	for (int m = 0; m <= length_scales; ++m)
	{
		for (int q = 0; q < cols; ++q)
		{
			sf_file<<(*S)(m,q)<<"\t";
		}
		sf_file<<endl;
	}
	// cout<<(*S)<<endl;

	delete u;
	delete v;
	delete kn;
}

Array<double,2> derivative_matrix(Array<cdouble,1> *u,Array<cdouble,1> *v, Array<double,1> *kn, Array<tv,1> e, pair<Actv, Actv> h, tv x)//, int i, int j)
{
	Array<ctv,1> h_plus(N);
	Array<ctv,1> h_minus(N);
	h_plus = h.first;
	h_minus = h.second;	

	int N = (*u).numElements() - 4;

	Array<double,2> DU(3,3);

	Array<double,1> phase(N);			// For storing phases while computing the inverse fourier transform i.e. derivative_tmp .
	for (int n = 0; n < N; ++n)
	{
		phase(n) = (*kn)(n+2)*dot(e(n), x);

	}
	// ===== Loop for the elements of the matrix DU ==== 

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			cdouble derivative_tmp  = cdouble(0,0); 
			for (int n = 0; n < N; ++n)
			{
				derivative_tmp += (h_plus(n)(i)*(*u)(n+2) + h_minus(n)(i)*(*v)(n+2)) * plusI*(*kn)(n+2)*e(n)(j) * exp(plusI*phase(n));
			}
			DU(i,j) = 2 * derivative_tmp.real();
		}
	}
	// =================================================
	// cout<<DU(0,0) + DU(1,1) + DU(2,2) <<endl; // Check if the sum of the diagonal elements is zero to ensure incompressibility. 

	return DU;
}

Array<double,2> derivative_matrix(Array<cdouble,1> *u, Array<double,1> *kn, Array<tv,1> e, Array<tv,1> c, tv x)//, int i, int j)
{
	int N = (*u).numElements() - 4;		// Number of shells

	Array<double,2> DU(3,3);			// Matrix for storing the derivatives of velocity

	Array<double,1> phase(N);			// For storing phases while computing the inverse fourier transform i.e. derivative_tmp .
	for (int n = 0; n < N; ++n)
	{
		phase(n) = (*kn)(n+2)*dot(e(n), x);

	}
	// ===== Loop for the elements of the matrix DU ==== 
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			cdouble derivative_tmp  = cdouble(0,0); 
			for (int n = 0; n < N; ++n)
			{
				derivative_tmp += c(n)(i)* (*u)(n+2) * plusI*(*kn)(n+2)*e(n)(j) * exp(plusI*phase(n));
			}
			DU(i,j) = 2 * derivative_tmp.real();
		}
	}
	// =================================================
	// cout<<DU(0,0) + DU(1,1) + DU(2,2) <<endl; // Check if the sum of the diagonal elements is zero to ensure incompressibility. 

	return DU;
}

double local_dissipation_rate(Array<double, 2> DU, double nu)
{
	double eps = 0;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			if (i==j)
			{
				// cout<<i<<","<<j<<endl;
				eps += 2 * pow(DU(i,j),2);
				// cout<<eps<<endl;
			}
			else if (i<j)
			{
				// cout<<i<<","<<j<<endl;
				eps += pow((DU(i,j) + DU(j,i)),2);
				// cout<<eps<<endl;
			}
		}
	}
	eps = nu * eps;
	return eps;
}

void eps(int pid, int procs_per_pid)
{
	stringstream u_vs_t;
	u_vs_t<<input_path<<"velocity_field_"<<Attempt<<"/u_vs_t_1.d";
	fstream field_file;
	//============= INPUT ================
	Array<double,1> *kn;	//,FortranArray<1>());			// wave number for shell n
	kn=new Array<double,1>(N+2);
	(*kn) = shell_wavenumbers(kn, N, lambda);

	double t=0;
	Array<cdouble,1> *u;									// Field u 
	u = new Array<cdouble,1>(N+4); 
	Array<cdouble,1> *v;									// Field v 
	v= new Array<cdouble,1>(N+4);
	
	field_file.open(u_vs_t.str().c_str());
	field_file>>t>>(*u)>>(*v);
	field_file.close();

	Array<tv,1> e(N);
	e = e_tuples(N,1);
	/*Array<tv,1> c(N);
	c = c_tuples(e,N,1);*/
	pair<Actv, Actv> h = h_tuples(e,N,1);
	//======================================

	Array<double,2> DU(3,3);			// Matrix for storing the derivatives of velocity
	tv x; 
	// double nu = 1e-3;
	// Array<double,3> *epsilon;
	// epsilon = new Array<double,3>(gridSize, gridSize, gridSize);
	Array<double,2> *epsilon_plane;
	epsilon_plane = new Array<double,2>(gridSize, gridSize);


	//~ for (int i = starting_plane; i < gridSize; ++i)
	for (int i = procs_per_pid * pid ; i < procs_per_pid * (pid + 1); ++i)
	{
		x(0) = 2*M_PI*i/gridSize;
		if(pid==0){	cout<<"plane: x = "<<i<<endl;}
		for (int j = 0; j < gridSize; ++j)
		{
			x(1) = 2*M_PI*j/gridSize;
			for (int k = 0; k < gridSize; ++k)
			{
				// cout<<i<<","<<j<<","<<k<<endl;
				//~ x(0) = 2*M_PI*i/gridSize; x(1) = 2*M_PI*j/gridSize; 
				x(2) = 2*M_PI*k/gridSize;
				// cout<<x<<endl;
				DU = derivative_matrix(u,v,kn,e,h,x);
				(*epsilon_plane)(j,k) = local_dissipation_rate(DU, nu);
			}
		}
		write_dissRate_plane(epsilon_plane, i);
	}
	// write_dissRate(epsilon);
}
