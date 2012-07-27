/*
Rigid body dynamics!
*/
#include "common.h"

extern atom m0[4];

int main() {
	matrix initial_inverse_moment_of_inertia = get_initial_inverse_moment_of_inertia();	

	int n_molecules = 20, i, j;
	h2o bodies[n_molecules];
	srand(1);
	double wx=-8, wy=-8, wz=-8;
	for(i=0; i<n_molecules; i++) {
		bodies[i].r = (vector) {wx,wy,wz}; wx+=3.5; if(wx>8) { wx=-8; wy+=3.5;} if(wy>8) { wy=-8; wz+=3.5;}
		bodies[i].v = (vector) {0,0,0};
		bodies[i].a = (vector) {0,0,0};
		bodies[i].L = (vector) {0,0,0};
		bodies[i].o = (vector) {0,0,0};
		double rx = 2*PI*rf, ry = 2*PI*rf, rz = 2*PI*rf;
		bodies[i].A.m[0][0] = cos(rz)*cos(ry); bodies[i].A.m[0][1] = -cos(rx)*sin(rz) + cos(rz)*sin(ry)*sin(rx); bodies[i].A.m[0][2] = cos(rz)*cos(rx)*sin(ry) + sin(rz)*sin(rx);
		bodies[i].A.m[1][0] = cos(ry)*sin(rz); bodies[i].A.m[1][1] = cos(rz)*cos(rx) + sin(rz)*sin(ry)*sin(rx); bodies[i].A.m[1][2] = cos(rx)*sin(rz)*sin(ry) - cos(rz)*sin(rx);
		bodies[i].A.m[2][0] = -sin(ry); bodies[i].A.m[2][1] = cos(ry)*sin(rx); bodies[i].A.m[2][2] = cos(ry)*cos(rx);
	}
	
	vector max_r = (vector) {10,10,10}, min_r = (vector) {-10,-10,-10};
	int periodic_BC = 1;
	
	stopwatch(); //for timing

	double t = 0.1; //femtoseconds
	int step;
	for(step=0; step<10000; step++) {
		//if(step%10==0) write_xyz_file("bodies.xyz", bodies, n_molecules, step);
		
		bodies[0].r.y += 0.1;
		for(i=0; i<n_molecules; i++) { //update linear and angular positions
			bodies[i].r = add_( add_(bodies[i].r, multiply_(bodies[i].v,t)), multiply_(bodies[i].a, 0.5*t*t) );
			bodies[i].A = add__(bodies[i].A, multiply_ms_(multiply__(tilde(bodies[i].o), bodies[i].A), t));
			bodies[i].A = reorthogonalize_(bodies[i].A);
			
			if(periodic_BC) { //enforce periodic BC
				 if(bodies[i].r.x>max_r.x) bodies[i].r.x=min_r.x;
				 if(bodies[i].r.y>max_r.y) bodies[i].r.y=min_r.y;
				 if(bodies[i].r.z>max_r.z) bodies[i].r.z=min_r.z;
				 if(bodies[i].r.x<min_r.x) bodies[i].r.x=max_r.x;
				 if(bodies[i].r.y<min_r.y) bodies[i].r.y=max_r.y;
				 if(bodies[i].r.z<min_r.z) bodies[i].r.z=max_r.z;
			}
		}
		
		for(i=0; i<n_molecules; i++) { //update velocities and accelerations
			vector total_force = (vector) {0,0,0};
			vector total_torque = (vector) {0,0,0};
			for(j=0; j<n_molecules; j++) {
				if(i==j) continue;
				force_torque pair = calculate_force_torque(bodies[i], bodies[j]);
				total_force = add_(total_force, pair.force);
				total_torque = add_(total_torque, pair.torque);
			}
			//print_(total_torque);
			//printf("%e\n", len_squared_(total_torque));

			vector total_acceleration = multiply_(total_force, 1/water_mass);
			
			integrate(bodies, i, total_acceleration, total_torque, initial_inverse_moment_of_inertia, t);
			
			//if(i==0) write_csv_file("out.csv", bodies[0].r.y, total_acceleration);
		}
	}
	
	printf("%f\n", stopwatch());

	return 0;
}
