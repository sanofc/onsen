#include "steam.h"

void swap(double * i, double * j){
	double * tmp = i;
	i = j;
	j = tmp;
}

void add_force(int N, double dt, double * u, double * v, double * u_prev, double * v_prev, vec * force){
	for(int i=0; i < N + 2; i++){
		for(int j=0; j < N + 2; j++){
			u[IU(i,j)] += dt * (force[IX(i-1,j)].x+force[IX(i,j)].x) * 0.5;
		}
	}
	for(int i=0; i < N + 2; i++){
		for(int j=0; j < N + 3; j++){
			v[IV(i,j)] += dt * (force[IX(i,j-1)].y+force[IX(i,j)].y) * 0.5;
		}
	}	
}

void vel_step(int N, double dt, double * u, double * v, double * u_prev, double * v_prev, vec * force){
	add_force(N,dt,u,v,u_prev,v_prev,force);
}