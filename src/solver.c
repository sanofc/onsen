#include "steam.h"

void add_force(int N, double dt, double * u, double * v, double * u_prev, double * v_prev, vec * force){
	for(int i=0; i < N + 3; i++){
		for(int j=0; j < N + 2; j++){
			u_prev[IU(i,j)] += dt * (force[IX(i-1,j)].x+force[IX(i,j)].x) * 0.5;
		}
	}
	for(int i=0; i < N + 2; i++){
		for(int j=0; j < N + 3; j++){
			v_prev[IV(i,j)] += dt * (force[IX(i,j-1)].y+force[IX(i,j)].y) * 0.5;
		}
	}	
	for(int i=0; i< N+2; i++){
		for(int j = 0; j < N+2; j++){
			force[IX(i,j)].x = force[IX(i,j)].y = 0.0;
		}
	}
}

void advect(int N, double dt, double * u, double * v, double * u_prev, double * v_prev){
	double h,x,y;
	h=1.0/(N+2);

	for(int i=1; i < N + 2; i++){
		x = i * h;
		for(int j=1; j < N + 1 ; j++){
		//int j =5;
			y = (j+0.5) * h;

			//backtrace
 			double x0 =  x - u_prev[IU(i,j)] * dt;
			double cv = (v_prev[IV(i-1,j)] + v_prev[IV(i,j)] + v_prev[IV(i-1,j+1)] + v_prev[IV(i-1,j+1)]) * 0.25;
			double y0 = y - cv * dt;

			int i0 = x0/h;
			int j0 = y0/h;



			double a = (x0 - i0 * h)/h;
			double b = (y0 - j0 * h)/h;

			if(u_prev[IU(i,j)] > 0.0){
				printf("i %d j %d i0 %d j0 %d a %f b %f  u_prev %f\n",i,j,i0,j0,a,b,u_prev[IU(i,j)]);
				fflush(stdout);
				u[IU(i,j)] = u_prev[IU(i0,j0)] * (1-a) + u_prev[IU(i0 + 1, j0)] * (a);
			}


/*
			u[IU(i,j)] = (u_prev[IU(i0,j0-1)]+u_prev[IU(i0,j0)])/2 * (1-a)  +
						 (u_prev[IU(i0+1,j0-1)]+u_prev[IU(i0+1,j0)])/2 * a  +
						 (u_prev[IU(i0,j0)]+u_prev[IU(i0,j0+1)])/2 * (1-a)  +
						 (u_prev[IU(i0+1,j0)]+u_prev[IU(i0+1,j0+1)])/2 * a ; 
						 

			if(i == 7 && j ==5){

				double t = (u_prev[IU(i0,j0-1)]+u_prev[IU(i0,j0)])/2 * (1-a) * (1-b) +
						 (u_prev[IU(i0+1,j0-1)]+u_prev[IU(i0+1,j0)])/2 * a * (1-b) +
						 (u_prev[IU(i0,j0)]+u_prev[IU(i0,j0+1)])/2 * (1-a) * b +
						 (u_prev[IU(i0+1,j0)]+u_prev[IU(i0+1,j0+1)])/2 * a * b; 
				printf("a %d %d %f %f %f %f %f %f \n",i0,j0,t,u_prev[IU(i0,j0-1)],u_prev[IU(i0,j0)],u_prev[IU(i0+1,j0-1)],u_prev[IU(i0,j0+1)],u_prev[IU(i0+1,j0+1)]);
				fflush(stdout);
			//}
			printf("a %.20f\n",a);
			//if(fabs(u_prev[IU(i,j)])>0.0){
				//printf("%f %f %f %f %f %f %f %d %d\n",u_prev[IU(i0,j0)],cv,a,b,x0,y0,h,i0,j0);
				printf("%d %d %d %d %f %f %f %f %f\n",i,j,i0,j0,a,b,u[IU(i,j)],u_prev[IU(i,j)],u_prev[IU(5,5)] );
				fflush(stdout);
			}
*/
		}
	}		
	//exit(0);
	/*
	for(int i=1; i < N + 1; i++){
		x = (i+0.5) * h;
		for(int j=1; j < N + 2 ; j++){
			y = j * h;

			//backtrace
 			double y0 =  y - v_prev[IV(i,j)] * dt;
			double cu = (u_prev[IU(i,j)] + u_prev[IU(i,j-1)] + u_prev[IU(i+1,j-1)] + u_prev[IU(i+1,j)]) * 0.25;
			double x0 = x - cu * dt;

			int i0 = (x0+0.5*h)/h;
			int j0 = y0/h;

			double a = (x0 - i0 * h -0.5*h)/h;
			double b = (y0 - j0 * h)/h;

			v[IV(i,j)] = v_prev[IV(i0,j0)] * (1-a) * (1-b) +
						 v_prev[IV(i0-1,j0)] * a * (1-b) +
						 v_prev[IV(i0,j0-1)] * (1-a) * b +
						 v_prev[IV(i0-1,j0-1)] * a * b; 

		}
	}*/
}

void swap(double **x , double **y){
	double *tmp = *x;
	*x = *y;
	*y = tmp;
}

void vel_step(int N, double dt, double * u, double * v, double * u_prev, double * v_prev, vec * force){
	add_force(N,dt,u,v,u_prev,v_prev,force);
	advect(N,dt,u,v,u_prev,v_prev);
}