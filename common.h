#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>
#define PI 3.1415926535898
#define coulomb_k 0.138935444 //in elementary charges, angstroms, amu, and femtoseconds [8.987551787e9 N m^2 C^-2 * (1.60217646e-19)^2 in (amu*angstroms/femtoseconds^2) angstroms^2 C^-2]
#define water_mass 18.0153
#define abs(a)	   (((a) < 0) ? -(a) : (a))
#define round(a)   (((a) < 0) ? (int)((a) - 0.5) : (int)((a) + 0.5))
#define rf (((double)rand())/RAND_MAX)
#define nrf (1 - 2*rf)
//angstroms and radians:
#define r_min 1.5
#define r_max 15
#define theta_min 0
#define theta_max 2*PI
#define phi_min 0.0
#define phi_max 2*PI
#define rx_min 0.0
#define rx_max 2*PI
#define ry_min 0.0
#define ry_max 2*PI
#define rz_min 0.0
#define rz_max 2*PI
#define r_step ((r_max-r_min)/N)
#define theta_step ((theta_max-theta_min)/N)
#define phi_step ((phi_max-phi_min)/N)
#define rx_step ((rx_max-rx_min)/N)
#define ry_step ((ry_max-ry_min)/N)
#define rz_step ((rz_max-rz_min)/N)

#define N 20

#ifndef H2O_STRUCTS
#define H2O_STRUCTS

typedef struct {
    double x,y,z;
} vector;
typedef struct {
    vector force, torque;
} force_torque;
typedef struct {
    double m[3][3];
} matrix;
typedef struct {
	vector r;
    double q,m;
} atom;
typedef struct {
	float acceleration, tx, ty, tz;
} potential_element;
typedef struct {
	vector r, v, a, L, o;
    matrix A;
} h2o;

#endif

void print_(vector a);
void print__(matrix a);
vector add_(vector a, vector b);
vector subtract_(vector a, vector b);
matrix add__(matrix a, matrix b);
vector multiply_(vector a, double s);
matrix multiply_ms_(matrix a, double s);
vector cross_(vector u, vector v);
double len_squared_(vector a);
matrix multiply__(matrix a, matrix b);
vector rotate_(matrix T, vector v);
matrix similarity_transform__(matrix a, matrix b);
matrix tilde(vector w);
matrix reorthogonalize_(matrix A);
matrix transpose_(matrix a);
matrix inverse_(matrix t);
double determinant_(matrix t);

potential_element get_water_potential(potential_element* water_potential, int ir, int t, int p, int irx, int iry, int irz);
potential_element interpolate_water_potential(double r, double theta, double phi, double rx, double ry, double rz);

force_torque calculate_force_torque(h2o i, h2o j);
void integrate(h2o *bodies, int i, vector total_acceleration, vector total_torque, matrix initial_inverse_moment_of_inertia, double t);
matrix get_initial_inverse_moment_of_inertia();
void write_xyz_file(char* name, h2o* bodies, int n_molecules, int step);
void write_csv_file(char* name, double index, vector accel, vector torque);
double stopwatch();
