/*
Parameterizes a single-point approximation of TIP4P.
*/
#include "common.h"

potential_element *water_potential;
void set_water_potential(int ir, int t, int p, int irx, int iry, int irz, double ax, double ay, double az, double tx, double ty, double tz) {
	//wrap angles to deal with getting unique values
	if(t>=N) t -= N;
	if(p>=N) p -= N;
	if(irx>=N) irx -= N;
	if(iry>=N) iry -= N;
	if(irz>=N) irz -= N;
	
	if(t<0) t += N;
	if(p<0) p += N;
	if(irx<0) irx += N;
	if(iry<0) iry += N;
	if(irz<0) irz += N;
	
	*(water_potential + ir*N*N*N*N*N + t*N*N*N*N + p*N*N*N + irx*N*N + iry*N + irz) = (potential_element){ax, ay, az, tx, ty, tz};
}

int main() {
	
	water_potential = malloc(sizeof(potential_element)*N*N*N*N*N*N);
	
	h2o stationary;
	stationary.r = (vector) {0,0,0};
	stationary.v = (vector) {0,0,0};
	stationary.a = (vector) {0,0,0};
	stationary.L = (vector) {0,0,0};
	stationary.o = (vector) {0,0,0};
	stationary.A.m[0][0] = 1; stationary.A.m[0][1] = 0; stationary.A.m[0][2] = 0;
	stationary.A.m[1][0] = 0; stationary.A.m[1][1] = 1; stationary.A.m[1][2] = 0;
	stationary.A.m[2][0] = 0; stationary.A.m[2][1] = 0; stationary.A.m[2][2] = 1;
	
	int ir, t, p, irx, iry, irz;
	for(ir=0; ir<N; ir++) {
		double r = (ir*r_step+r_min);
		for(t=0; t<N+1; t++) {
			for(p=0; p<N+1; p++) {
				double theta = t*theta_step+theta_min;
				double phi = p*phi_step+phi_min;
				vector dr = {r*sin(phi)*cos(theta), r*sin(phi)*sin(theta), r*cos(phi)};
				
				for(irx=0; irx<N+1; irx++) {
					for(iry=0; iry<N+1; iry++) {
						for(irz=0; irz<N+1; irz++) {
							double rx = irx*rx_step + rx_min, ry = iry*ry_step + ry_min, rz = irz*rz_step + rz_min;
							
							h2o mobile;
							mobile.r = dr;
							mobile.v = (vector) {0,0,0};
							mobile.a = (vector) {0,0,0};
							mobile.L = (vector) {0,0,0};
							mobile.o = (vector) {0,0,0};
							//x,y,z order rotation:
							mobile.A.m[0][0] = cos(rz)*cos(ry); mobile.A.m[0][1] = -cos(rx)*sin(rz) + cos(rz)*sin(ry)*sin(rx); mobile.A.m[0][2] = cos(rz)*cos(rx)*sin(ry) + sin(rz)*sin(rx);
							mobile.A.m[1][0] = cos(ry)*sin(rz); mobile.A.m[1][1] = cos(rz)*cos(rx) + sin(rz)*sin(ry)*sin(rx); mobile.A.m[1][2] = cos(rx)*sin(rz)*sin(ry) - cos(rz)*sin(rx);
							mobile.A.m[2][0] = -sin(ry); mobile.A.m[2][1] = cos(ry)*sin(rx); mobile.A.m[2][2] = cos(ry)*cos(rx);
							
							//get unique values:
							double unique_theta = atan2(dr.y,dr.x);
							double unique_phi = acos(dr.z/r);
							double unique_rx = atan2(mobile.A.m[2][1], mobile.A.m[2][2]);
							double unique_ry = atan2(-mobile.A.m[2][0], sqrt(mobile.A.m[2][1]*mobile.A.m[2][1] + mobile.A.m[2][2]*mobile.A.m[2][2]));
							double unique_rz = atan2(mobile.A.m[1][0], mobile.A.m[0][0]);
							//recalculate irx, iry, irz:
							int unique_t = (int)round((unique_theta-theta_min)/theta_step), unique_p = (int)round((unique_phi-phi_min)/phi_step), unique_irx = (int)round((unique_rx-rx_min)/rx_step), unique_iry = (int)round((unique_ry-ry_min)/ry_step), unique_irz = (int)round((unique_rz-rz_min)/rz_step); //low indices
							
							force_torque pair = calculate_force_torque(stationary, mobile);
							
							set_water_potential(ir, unique_t, unique_p, unique_irx, unique_iry, unique_irz, pair.force.x/water_mass, pair.force.y/water_mass, pair.force.z/water_mass, pair.torque.x, pair.torque.y, pair.torque.z);
						}
					}
				}
			}
		}
	}
	
	FILE *output = fopen("model.dat", "wb"); //check for IO errors here!
	if(output==NULL) {puts("File open failed"); exit(1);}
	int entries_written = fwrite(water_potential, sizeof(potential_element), N*N*N*N*N*N, output);
	if(entries_written!=N*N*N*N*N*N) {printf("File write failed, wrote %d elements out of %d\n", entries_written, N*N*N*N*N*N); exit(1);}
	fclose(output);
	
	return 0;
}
