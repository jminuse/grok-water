/*
Molecular dynamics using GROK
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
		
		//printf("%f %f %f %f %f %f\n", bodies[i].r.x, bodies[i].r.y, bodies[i].r.z, rx, ry, rz);
	}

	vector max_r = (vector) {10,10,10}, min_r = (vector) {-10,-10,-10};
	int periodic_BC = 1, regulate_pressure = 0, regulate_temperature = 0;
	double desired_pressure = 1, pressure_timescale = 1000, desired_temperature = 300, temperature_timescale = 1000, friction_coefficient = 0; //atm, femtoseconds and Kelvin
	
	double t = 0.1; //femtoseconds
	int step;
	for(step=0; step<10000; step++) {
		//if(step%10==0) write_xyz_file("test.xyz", bodies, n_molecules, step);
		
		//bodies[0].r.y += 0.1;
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
			if(regulate_pressure) { //rescale atom positions
				double pressure = 1;
				pressure = desired_pressure/pressure_timescale; //dummy until real code
			}
		}
		
		for(i=0; i<n_molecules; i++) { //update velocities and accelerations
			vector total_acceleration = (vector) {0,0,0};
			vector total_torque = (vector) {0,0,0};
			matrix anti_rotate = inverse_(bodies[i].A); //it is much more accurate to use the real inverse instead of relying on inverse = transpose
			for(j=0; j<n_molecules; j++) {
				if(i==j) continue;
				vector dr = subtract_(bodies[j].r, bodies[i].r);
				if(len_squared_(dr) > (r_max-r_step)*(r_max-r_step)) continue;
				
				vector rot_dr = rotate_(anti_rotate, dr); //rotate by inverse of i.A, about i
				
				matrix dA = multiply__(anti_rotate, bodies[j].A);
				
				double r = fastsqrt(len_squared_(rot_dr));
				double theta = fastatan2(rot_dr.y,rot_dr.x);
				double phi = fastacos(rot_dr.z/r);
				
				double rx = fastatan2(dA.m[2][1], dA.m[2][2]);
				double ry = fastatan2(-dA.m[2][0], sqrt(dA.m[2][1]*dA.m[2][1] + dA.m[2][2]*dA.m[2][2]));
				double rz = fastatan2(dA.m[1][0], dA.m[0][0]);
				
				potential_element e = interpolate_water_potential(r, theta, phi, rx, ry, rz);
				vector torque = (vector){e.tx, e.ty, e.tz};
				vector acceleration = (vector){e.ax, e.ay, e.az};
				
				//Restore to global coordinates:
				torque = rotate_(bodies[i].A, torque);
				acceleration = rotate_(bodies[i].A, acceleration);
				
				total_acceleration = add_(total_acceleration, acceleration);
				total_torque = add_(total_torque, torque);
				
				if(periodic_BC) { //recalculate with j as if it were across each periodic boundary. 
					
				}
				
			}
			
			if(regulate_temperature) { //adjust acceleration, aka Nose-Hoover thermostat. Must iterate to make temperature consistent with updated velocities.
				total_acceleration = subtract_(total_acceleration, multiply_(bodies[i].v, friction_coefficient));
				double temperature = 300; //calculate this from linear and angular velocities. Count on equipartition to spread the energy out again from the linear degrees of freedom. Or, adjust bodies[i].o as well.
				friction_coefficient += t/temperature_timescale*(temperature/desired_temperature - 1);
			}
			
			integrate(bodies, i, total_acceleration, total_torque, initial_inverse_moment_of_inertia, t);
			
			//if(i==0) write_csv_file("out.csv", bodies[0].r.y, total_acceleration);
		}
	}

	printf("%f\n", stopwatch());

	return 0;
}
