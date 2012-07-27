/*
Functions for linear algebra, rigid-body physics, and the GROK model.
*/
#include "common.h"

atom m0[] = { (atom){{0.0,0.065098,0.0},0,15.9994},
	 (atom){{-0.756950,-0.520784,0.0},0.52,1.0008},
	 (atom){{0.756950,-0.520784,0.0},0.52,1.0008},
	 (atom){{0.0,-0.084902,0.0},-1.04,0} };

void print_(vector a) {
	printf("<%.2e %.2e %.2e>\n", a.x, a.y, a.z);
}
void print__(matrix a) {
	printf("|%.2e %.2e %.2e|\n", a.m[0][0], a.m[0][1], a.m[0][2]);
	printf("|%.2e %.2e %.2e|\n", a.m[1][0], a.m[1][1], a.m[1][2]);
	printf("|%.2e %.2e %.2e|\n", a.m[2][0], a.m[2][1], a.m[2][2]);
}
vector add_(vector a, vector b) { //-O3 should inline these functions
	return (vector) {a.x+b.x, a.y+b.y, a.z+b.z};
}
vector subtract_(vector a, vector b) {
	return (vector) {a.x-b.x, a.y-b.y, a.z-b.z};
}
matrix add__(matrix a, matrix b) {
	matrix c;
	c.m[0][0] = a.m[0][0] + b.m[0][0];
	c.m[0][1] = a.m[0][1] + b.m[0][1];
	c.m[0][2] = a.m[0][2] + b.m[0][2];
	
	c.m[1][0] = a.m[1][0] + b.m[1][0];
	c.m[1][1] = a.m[1][1] + b.m[1][1];
	c.m[1][2] = a.m[1][2] + b.m[1][2];
	
	c.m[2][0] = a.m[2][0] + b.m[2][0];
	c.m[2][1] = a.m[2][1] + b.m[2][1];
	c.m[2][2] = a.m[2][2] + b.m[2][2];
	return c;
}
vector multiply_(vector a, double s) {
	return (vector) {a.x*s, a.y*s, a.z*s};
}
matrix multiply_ms_(matrix a, double s) {
	matrix c;
	c.m[0][0] = a.m[0][0]*s;
	c.m[0][1] = a.m[0][1]*s;
	c.m[0][2] = a.m[0][2]*s;
	
	c.m[1][0] = a.m[1][0]*s;
	c.m[1][1] = a.m[1][1]*s;
	c.m[1][2] = a.m[1][2]*s;
	
	c.m[2][0] = a.m[2][0]*s;
	c.m[2][1] = a.m[2][1]*s;
	c.m[2][2] = a.m[2][2]*s;
	return c;
}
vector cross_(vector u, vector v) {
	return (vector) {u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x};
}
double len_squared_(vector a) {
	return a.x*a.x + a.y*a.y + a.z*a.z;
}
matrix multiply__(matrix a, matrix b) {
	matrix c;
	c.m[0][0] = a.m[0][0]*b.m[0][0] + a.m[0][1]*b.m[1][0] + a.m[0][2]*b.m[2][0];
	c.m[0][1] = a.m[0][0]*b.m[0][1] + a.m[0][1]*b.m[1][1] + a.m[0][2]*b.m[2][1];
	c.m[0][2] = a.m[0][0]*b.m[0][2] + a.m[0][1]*b.m[1][2] + a.m[0][2]*b.m[2][2];
	
	c.m[1][0] = a.m[1][0]*b.m[0][0] + a.m[1][1]*b.m[1][0] + a.m[1][2]*b.m[2][0];
	c.m[1][1] = a.m[1][0]*b.m[0][1] + a.m[1][1]*b.m[1][1] + a.m[1][2]*b.m[2][1];
	c.m[1][2] = a.m[1][0]*b.m[0][2] + a.m[1][1]*b.m[1][2] + a.m[1][2]*b.m[2][2];
	
	c.m[2][0] = a.m[2][0]*b.m[0][0] + a.m[2][1]*b.m[1][0] + a.m[2][2]*b.m[2][0];
	c.m[2][1] = a.m[2][0]*b.m[0][1] + a.m[2][1]*b.m[1][1] + a.m[2][2]*b.m[2][1];
	c.m[2][2] = a.m[2][0]*b.m[0][2] + a.m[2][1]*b.m[1][2] + a.m[2][2]*b.m[2][2];
	return c;
}
vector rotate_(matrix T, vector v) {
	return (vector) {	v.x*T.m[0][0] + v.y*T.m[0][1] + v.z*T.m[0][2],
						v.x*T.m[1][0] + v.y*T.m[1][1] + v.z*T.m[1][2],
						v.x*T.m[2][0] + v.y*T.m[2][1] + v.z*T.m[2][2] };
}
matrix similarity_transform__(matrix a, matrix b) {
	matrix first_step = multiply__(a,b);
	double t1 = a.m[0][1], t2 = a.m[0][2], t3 = a.m[1][2];
	a.m[0][1] = a.m[1][0];
	a.m[0][2] = a.m[2][0];
	a.m[1][2] = a.m[2][1];
	a.m[1][0] = t1; a.m[2][0] = t2; a.m[2][1] = t3;
	return multiply__(first_step,a);
}
matrix tilde(vector w) {
	matrix c;
	c.m[0][0] = 0; c.m[0][1] = -w.z; c.m[0][2] = w.y;
	c.m[1][0] = w.z; c.m[1][1] = 0; c.m[1][2] = -w.x;
	c.m[2][0] = -w.y; c.m[2][1] = w.x; c.m[2][2] = 0;
	return c;
}
matrix reorthogonalize_(matrix A) { //reorthogonalize: repeatedly average the rotation matrix with its inverse transpose (itself, if it is already orthogonal). there is a more optimized way to do this, or tell if you don't have to - check paper.
	int j;
	for(j=0; j<3; j++) {
		matrix t = A;
		double t1 = t.m[0][1], t2 = t.m[0][2], t3 = t.m[1][2];
		t.m[0][1] = t.m[1][0];
		t.m[0][2] = t.m[2][0];
		t.m[1][2] = t.m[2][1];
		t.m[1][0] = t1; t.m[2][0] = t2; t.m[2][1] = t3;
		double one_over_det = 1/(t.m[0][0]*(t.m[2][2]*t.m[1][1]-t.m[2][1]*t.m[1][2]) - t.m[1][0]*(t.m[2][2]*t.m[0][1]-t.m[2][1]*t.m[0][2]) + t.m[2][0]*(t.m[1][2]*t.m[0][1]-t.m[1][1]*t.m[0][2]));
		matrix it;
		it.m[0][0] = t.m[2][2]*t.m[1][1]-t.m[2][1]*t.m[1][2];
		it.m[0][1] = t.m[2][1]*t.m[0][2]-t.m[2][2]*t.m[0][1];
		it.m[0][2] = t.m[1][2]*t.m[0][1]-t.m[1][1]*t.m[0][2];
		
		it.m[1][0] = t.m[2][0]*t.m[1][2]-t.m[2][2]*t.m[1][0];
		it.m[1][1] = t.m[2][2]*t.m[0][0]-t.m[2][0]*t.m[0][2];
		it.m[1][2] = t.m[1][0]*t.m[0][2]-t.m[1][2]*t.m[0][0];
		
		it.m[2][0] = t.m[2][1]*t.m[1][0]-t.m[2][0]*t.m[1][1];
		it.m[2][1] = t.m[2][0]*t.m[0][1]-t.m[2][1]*t.m[0][0];
		it.m[2][2] = t.m[1][1]*t.m[0][0]-t.m[1][0]*t.m[0][1];
		it = multiply_ms_(it, one_over_det);
		A = add__(A, it);
		A = multiply_ms_(A, 0.5);
	}
	return A;
}
matrix transpose_(matrix a) {
	double t1 = a.m[0][1], t2 = a.m[0][2], t3 = a.m[1][2];
	a.m[0][1] = a.m[1][0];
	a.m[0][2] = a.m[2][0];
	a.m[1][2] = a.m[2][1];
	a.m[1][0] = t1; a.m[2][0] = t2; a.m[2][1] = t3;
	return a;
}
matrix inverse_(matrix t) {
	double one_over_det = 1/(t.m[0][0]*(t.m[2][2]*t.m[1][1]-t.m[2][1]*t.m[1][2]) - t.m[1][0]*(t.m[2][2]*t.m[0][1]-t.m[2][1]*t.m[0][2]) + t.m[2][0]*(t.m[1][2]*t.m[0][1]-t.m[1][1]*t.m[0][2]));
	matrix it;
	it.m[0][0] = t.m[2][2]*t.m[1][1]-t.m[2][1]*t.m[1][2];
	it.m[0][1] = t.m[2][1]*t.m[0][2]-t.m[2][2]*t.m[0][1];
	it.m[0][2] = t.m[1][2]*t.m[0][1]-t.m[1][1]*t.m[0][2];
	
	it.m[1][0] = t.m[2][0]*t.m[1][2]-t.m[2][2]*t.m[1][0];
	it.m[1][1] = t.m[2][2]*t.m[0][0]-t.m[2][0]*t.m[0][2];
	it.m[1][2] = t.m[1][0]*t.m[0][2]-t.m[1][2]*t.m[0][0];
	
	it.m[2][0] = t.m[2][1]*t.m[1][0]-t.m[2][0]*t.m[1][1];
	it.m[2][1] = t.m[2][0]*t.m[0][1]-t.m[2][1]*t.m[0][0];
	it.m[2][2] = t.m[1][1]*t.m[0][0]-t.m[1][0]*t.m[0][1];
	it = multiply_ms_(it, one_over_det);
	return it;
}
double determinant_(matrix t) {
	return (t.m[0][0]*(t.m[2][2]*t.m[1][1]-t.m[2][1]*t.m[1][2]) - t.m[1][0]*(t.m[2][2]*t.m[0][1]-t.m[2][1]*t.m[0][2]) + t.m[2][0]*(t.m[1][2]*t.m[0][1]-t.m[1][1]*t.m[0][2]));
}
double fastsqrt(double a) { //http://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Approximations_that_depend_on_IEEE_representation, modified to use doubles
    union
    {
        long long int tmp;
        double f;
    } u;
 
    u.f = a;
 
    u.tmp -= 1LL << 52; /* Subtract 2^m. */
    u.tmp >>= 1LL; /* Divide by 2. */
    u.tmp += 1LL << (9+52); /* Add ((b + 1) / 2) * 2^m. */
 
    return u.f;
}
double fastacos(double a) {
	return acos(a);
}
double fastatan2(double a, double b) {
	//return atan2(a,b);
	return (fabs(a)>fabs(b))?(-a*b/(a*a + 0.28088*b*b) + PI/2*(a>0?1.0:-1.0)):(-PI/2*(a*b>0?1.0:-1.0) + PI/2*(a>0?1.0:-1.0) + a*b/(b*b + 0.28088*a*a));
}
potential_element get_water_potential(potential_element *water_potential, int ir, int t, int p, int irx, int iry, int irz) {
	assert(ir>=0 && t>=0 && p>=0 && irx>=0 && iry>=0 && irz>=0);
	assert(ir<N && t<N && p<N && irx<N && iry<N && irz<N);
	return *(water_potential + ir*N*N*N*N*N + t*N*N*N*N + p*N*N*N + irx*N*N + iry*N + irz);
}
potential_element interpolate_water_potential(double r, double theta, double phi, double rx, double ry, double rz) { //6D linear interpolation
	static potential_element *water_potential = NULL;
	if(water_potential==NULL) {
		FILE *input = fopen("model.dat", "rb");
		if(input==NULL) {puts("File open failed"); exit(1);}
		water_potential = malloc(sizeof(potential_element)*N*N*N*N*N*N);
		//fread(water_potential, sizeof(potential_element), N*N*N*N*N*N, input);
		int elements_read = 0;
		while(elements_read<N*N*N*N*N*N) {
			elements_read += fread(water_potential+elements_read, sizeof(potential_element), 10000, input);
		}
		if(elements_read != N*N*N*N*N*N) {printf("File read failed, read %d elements out of %d\n", elements_read, N*N*N*N*N*N); exit(1);}
		fclose(input);
		
		stopwatch(); //for timing
	}

	if(phi < 0) phi += 2*PI; //reduce to range 0 to 2*PI
	if(theta < 0) theta += 2*PI;
	if(rx < 0) rx += 2*PI;
	if(ry < 0) ry += 2*PI;
	if(rz < 0) rz += 2*PI;
	
	//arrange symmetry in rx,ry,rz here
	
	int ir = (int)((r-r_min)/r_step), t = (int)((theta-theta_min)/theta_step), p = (int)((phi-phi_min)/phi_step), irx = (int)((rx-rx_min)/rx_step), iry = (int)((ry-ry_min)/ry_step), irz = (int)((rz-rz_min)/rz_step); //low indices
	
	if(ir<0) return (potential_element) {0,0,0,0,0,0}; //too close
	else if(ir>=N) return (potential_element) {0,0,0,0,0,0}; //too far
	
	float weights[2][6] = { {0,0,0,0,0,0}, { (r-r_min-ir*r_step)/r_step, (theta-theta_min-t*theta_step)/theta_step, (phi-phi_min-p*phi_step)/phi_step, (rx-rx_min-irx*rx_step)/rx_step, (ry-ry_min-iry*ry_step)/ry_step, (rz-rz_min-irz*rz_step)/rz_step } }; //linear weights
	int i;
	for(i=0; i<6; i++) {
		weights[0][i] = 1-weights[1][i];
	}
	
	potential_element sum = (potential_element) {0,0,0,0,0,0};
	float weight_sum = 0;
	int a,b,c,d,e,f; 
	for(a=0; a<2; a++) { for(b=0; b<2; b++) { for(c=0; c<2; c++) { for(d=0; d<2; d++) { for(e=0; e<2; e++) { for(f=0; f<2; f++) {
		float w = weights[a][0]*weights[b][1]*weights[c][2]*weights[d][3]*weights[e][4]*weights[f][5];
		//if(fabs(w)<0.05) continue; //throws out low weights: since this is linear, far guesses are low-quality.
		potential_element current = get_water_potential(water_potential, (ir+a)==N?N-1:(ir+a), (t+b)==N?0:(t+b), (p+c)==N?0:(p+c), (irx+d)==N?0:(irx+d), (iry+e)==N?0:(iry+e), (irz+f)==N?0:(irz+f) ); //truncate radius and wrap angles
		sum.ax += current.ax*w;
		sum.ay += current.ay*w;
		sum.az += current.az*w;
		sum.tx += current.tx*w;
		sum.ty += current.ty*w;
		sum.tz += current.tz*w;
		weight_sum += w;
	}}}}}}
	if(weight_sum==0) { //take nearest neighbor
		int nir = round((r-r_min)/r_step), nt = round((theta-theta_min)/theta_step), np = round((phi-phi_min)/phi_step), nirx = round((rx-rx_min)/rx_step), niry = round((ry-ry_min)/ry_step), nirz = round((rz-rz_min)/rz_step);
		return get_water_potential(water_potential, nir==N?N-1:nir, nt==N?0:nt, np==N?0:np, nirx==N?0:nirx, niry==N?0:niry, nirz==N?0:nirz ); //truncate radius and wrap angles
	}
	double normalize_weights = 1/weight_sum;
	sum.ax *= normalize_weights;
	sum.ay *= normalize_weights;
	sum.az *= normalize_weights;
	sum.tx *= normalize_weights;
	sum.ty *= normalize_weights;
	sum.tz *= normalize_weights;
	
	return sum;
}
force_torque calculate_force_torque(h2o i, h2o j) {
	vector total_force = (vector) {0,0,0};
	vector total_torque = (vector) {0,0,0};
	int a, b;
	
	for(a=0; a<4; a++) { //count through atoms in TIP4P, molecule i
		vector rotated_atom_i = rotate_(i.A, m0[a].r);
		vector r_atom_i = add_(rotated_atom_i, i.r);
		for(b=0; b<4; b++) { //molecule j
			vector r_atom_j = add_(rotate_(j.A, m0[b].r), j.r);
			vector dr = subtract_(r_atom_i, r_atom_j);
			double r_squared = len_squared_(dr);
			if(r_squared > r_max*r_max) continue;
			double force_magnitude = coulomb_k*m0[a].q*m0[b].q/r_squared;
			
			if(r_squared < 10*10 && m0[a].m>15 && m0[b].m>15)
				force_magnitude += 24*6.4852e-5*(2* pow(3.16,12) * pow(r_squared,-6.5) - pow(3.16,6) * pow(r_squared,-3.5)); //Units via google [0.1550/6.0221415e23 kcal in (amu*angstroms^2) / femtoseconds^2]
			
			vector force = multiply_( dr, force_magnitude/sqrt(r_squared) );
			
			total_force = add_(total_force, force);
			total_torque = add_(total_torque, cross_(rotated_atom_i,force));
		}
	}
	
	return (force_torque) {total_force, total_torque};
}
void integrate(h2o *bodies, int i, vector total_acceleration, vector total_torque, matrix initial_inverse_moment_of_inertia, double t) {
	bodies[i].v = add_(bodies[i].v, multiply_(add_(bodies[i].a, total_acceleration), 0.5*t));
	bodies[i].a = total_acceleration;

	bodies[i].L = add_(bodies[i].L, multiply_(total_torque, t));

	matrix current_inverse_moment_of_inertia = similarity_transform__(bodies[i].A, initial_inverse_moment_of_inertia);
	bodies[i].o = rotate_(current_inverse_moment_of_inertia, bodies[i].L);
}
matrix get_initial_inverse_moment_of_inertia() {
	double moment_of_inertia[3][3] = {
		{m0[0].m*m0[0].r.y*m0[0].r.y + 2*m0[1].m*m0[1].r.y*m0[1].r.y, 0, 0},
		{0, 2*m0[1].m*m0[1].r.x*m0[1].r.x, 0},
		{0, 0, m0[0].m*(m0[0].r.x*m0[0].r.x + m0[0].r.y*m0[0].r.y) + 2*m0[1].m*(m0[1].r.x*m0[1].r.x + m0[1].r.y*m0[1].r.y)}
	};
	double inverse_moment_of_inertia[3][3] =  {
		{moment_of_inertia[1][1]*moment_of_inertia[2][2], 0, 0},
		{0, moment_of_inertia[0][0]*moment_of_inertia[2][2], 0},
		{0, 0, moment_of_inertia[0][0]*moment_of_inertia[1][1]}
	};
	double moment_of_inertia_determinant = moment_of_inertia[0][0]*moment_of_inertia[1][1]*moment_of_inertia[2][2];
	matrix initial_inverse_moment_of_inertia;	
	int i, j;
	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			initial_inverse_moment_of_inertia.m[i][j] = inverse_moment_of_inertia[i][j] / moment_of_inertia_determinant;
		}
	}
	return initial_inverse_moment_of_inertia;
}
void write_xyz_file(char* name, h2o* bodies, int n_molecules, int step) {
	static FILE *xyz = NULL;
	if(xyz==NULL) { xyz = fopen(name, "wb"); }
	fprintf(xyz, "%d\nAtoms\n", 3*n_molecules);
	int i;
	for(i=0; i<n_molecules; i++) {
		int a;
		for(a=0; a<3; a++) {
			vector r = add_(rotate_(bodies[i].A, m0[a].r), bodies[i].r);
			if(isnan(r.x)) { fclose(xyz); printf("NAN at step %d\n", step); exit(1); }
			fprintf(xyz, "%c %f %f %f\n", m0[a].m>15?'O':'H', r.x, r.y, r.z);
		}
	}
}
void write_csv_file(char* name, double index, vector v) {
	static FILE *csv = NULL;
	if(csv==NULL) { csv = fopen(name, "wb"); }
	fprintf(csv, "%e,%e,%e,%e\n", index, v.x, v.y, v.z);
}
double stopwatch() {
	static struct timeval start;
	static int started = 0;
	if(!started) {
		gettimeofday(&start, NULL);
		started = 1;
		return 0.0;
	}
	struct timeval now;
	gettimeofday(&now, NULL);
	int seconds = now.tv_sec - start.tv_sec;
    int useconds = now.tv_usec - start.tv_usec;
    if(useconds < 0) {
        --seconds;
        useconds += 1000000;
    }
    started = 0;
    return seconds + ((double)useconds)/1e6;
}
