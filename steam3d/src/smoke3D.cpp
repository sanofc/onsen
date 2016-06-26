/*
 *  smoke3D.cpp
 *  smoke3D
 *
 */

#include "types.h"
#include "utility.h"
#include "smoke3D.h"
#include "solver.h"
#include "advect.h"
#include "render.h"
#include "perlin.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>


#define	N			64
#define LIMIT		200
#define DT			0.2
#define SPHERE_R	0.2
#define T_AMB		10.0

FLOAT *** u[3] = { NULL, NULL, NULL };
FLOAT *** b = NULL; //boundary
FLOAT *** v = NULL; //concentration(vapor)
FLOAT *** t = NULL; //temperature
FLOAT *** s = NULL; //steam
int frame = 0;
perlin perlin_noise;

void smoke3D::init() {
	// Allocate Memory
	u[0] = alloc3D(N+1,N,N);
	u[1] = alloc3D(N,N+1,N);
	u[2] = alloc3D(N,N,N+1);
	v = alloc3D(N,N,N);
	b = alloc3D(N,N,N);
	t = alloc3D(N,N,N);
	s = alloc3D(N,N,N);

	// Mark Wall Inside A Sphere
	/*
	int w = SPHERE_R*N;
	for( int i=-w; i<=w; i++ ) {
		for( int j=-w; j<=w; j++ ) {
			for( int k=-w; k<=w; k++ ) {
				if( hypot(hypot(i,j),k) < w ) {
					b[i+N/2][j+N/2][k+N/2] = 1.0;
				}
			}
		}
	}*/

#if _OPENMP
	printf( "OpenMP Detected.\n" );
#endif
}

static void scene(){
	// Set Initial Value
	if( frame == 0 ) {

		// Temperature
		FOR_EVERY_CELL{
			t[i][j][k] = T_AMB;
		} END_FOR

	}
	
	/*
	// Set Boundary Flow Around Sphere Zero
	FOR_EVERY_CELL {
		if( b[i][j][k] ) {
			c[i][j][k] = 0.0;
			u[0][i][j][k] = u[0][i+1][j][k] = 0.0;
			u[1][i][j][k] = u[1][i][j+1][k] = 0.0;
			u[2][i][j][k] = u[2][i][j][k+1] = 0.0;
		}
	} END_FOR
	*/
	
	// Vapor 
	int w = N/2;
	for( int i=-w; i<=w; i++ ) {
		for( int j=-w; j<=w; j++ ) {
			if( hypot(i,j) < w ) {
				for( int k=0; k<1; k++ ) {
					//double noise = perlin_noise.noise((double)i*1.0/(double)(N/2),(double)j*1.0/(double)(N/2),(double)frame*1.0/(double)(N/2)) * 8.0;
					//if(noise < 0) noise = 0.0;
					double noise = 1.0;	
					v[(int)(N/2)+i][k][(int)(N/2)+j] = noise * 2.0;
					t[(int)(N/2)+i][k][(int)(N/2)+j] = T_AMB + 5.0 * noise;
				}
			}
		}
	}

	//buoyancy
	
	FOR_EVERY_CELL {

		if(j<1) continue;

		//Saturation Vapor Content
		double a = 10.0;
		double b = 30.0;
		double c = -2.3;
		double m = fmin(a * exp(-b / ((t[i][j][k]) + c)),v[i][j][k]+s[i][j][k]);

		//Phase transition
		double r = 0.7;
		double ds = r * (v[i][j][k] - m);

		s[i][j][k] += ds;
		v[i][j][k] -= ds;

		//latent heat
		t[i][j][k] += 0.5 * ds;
	
		//Ambient Temperature
		FLOAT at;

		at = T_AMB;
		/*
		if(j==1){
			at = (g_ref(t,i+1,j,k,N) + g_ref(t,i-1,j,k,N) + g_ref(t,i,j+1,k,N) + g_ref(t,i,j,k+1,N) + g_ref(t,i,j,k-1,N))/5.0;
		}else{
		    at = (g_ref(t,i+1,j,k,N) + g_ref(t,i-1,j,k,N) + g_ref(t,i,j+1,k,N) + g_ref(t,i,j-1,k,N) + g_ref(t,i,j,k+1,N) + g_ref(t,i,j,k-1,N))/6.0;
		}
		*/

	    //at = (g_ref(t,i+1,j,k,N) + g_ref(t,i-1,j,k,N) + g_ref(t,i,j,k+1,N) + g_ref(t,i,j,k-1,N))/4.0;
		//FLOAT at = (g_ref(t,i,j+1,k,N) + g_ref(t,i,j-1,k,N))/2.0;
		//FLOAT at = g_ref(t,i,j+1,k,N);

		// Give Some External Force
		FLOAT alpha = 0.002;
		FLOAT beta = 0.003;
		FLOAT buoy =  -alpha * s[i][j][k] + beta * (t[i][j][k]-at);		

		
		u[1][i][j][k] += buoy;

	} END_FOR
}

static void enforce_boundary() {
	// Set Boundary Velocity Zero
	FOR_EVERY_X_FLOW {
		if( i==0 || i==N ) u[0][i][j][k] = 0.0;
	} END_FOR
	
	FOR_EVERY_Y_FLOW {
		//if( j==0 || j==N ) u[1][i][j][k] = 0.0;
		if (j == 0) u[1][i][j][k] = 0.0;
		if (j == N) u[1][i][j][k] = 0.05;
	} END_FOR
	
	FOR_EVERY_Z_FLOW {
		if( k==0 || k==N ) u[2][i][j][k] = 0.0;
	} END_FOR

	FOR_EVERY_CELL{
		if (j == N) {
			s[i][j][k] = 0.0;
			v[i][j][k] = 0.0;
			t[i][j][k] = T_AMB;
		}
	} END_FOR

}

static void project() {
	// Cell Width
	FLOAT h = 1.0/N;
	
	// Memory Allocation
	static FLOAT *** div = alloc3D(N,N,N);
	static FLOAT *** p = alloc3D(N,N,N);
	
	// Compute Divergence
	FOR_EVERY_CELL {
		div[i][j][k] = (u[0][i+1][j][k]-u[0][i][j][k]+
						u[1][i][j+1][k]-u[1][i][j][k]+
						u[2][i][j][k+1]-u[2][i][j][k]) / h;
	} END_FOR
	
	// Solve Pressure
	solver::solve( p, div, b, N );
	
	// Subtract Pressure Gradient
	FOR_EVERY_CELL {
		if( i>0 && i<N ) u[0][i][j][k] -= (p[i][j][k]-p[i-1][j][k])/h;
		if( j>0 && j<N ) u[1][i][j][k] -= (p[i][j][k]-p[i][j-1][k])/h;
		if( k>0 && k<N ) u[2][i][j][k] -= (p[i][j][k]-p[i][j][k-1])/h;
	} END_FOR
}

static void advection() {
	// Advect
	advect::advect( u, v, s, t, N, DT );
}

static void diffuse(){

	static FLOAT ***out[2] = { alloc3D(N,N,N), alloc3D(N,N,N)};

	FLOAT Dv = DT * 0.003;
	FLOAT Dt = DT * 0.002;

	FOR_EVERY_CELL {
		out[0][i][j][k] = g_ref(v,i,j,k,N) + 
						Dv * t[i][j][k] * (g_ref(v,i+1,j,k,N) + g_ref(v,i-1,j,k,N) + 
					  	g_ref(v,i,j+1,k,N) + g_ref(v,i,j-1,k,N) +
					  	g_ref(v,i,j,k+1,N) + g_ref(v,i,j,k-1,N) - 6 * g_ref(v,i,j,k,N));

		out[1][i][j][k] = g_ref(t,i,j,k,N) + 
						Dt * out[0][i][j][k] / (double)T_AMB * (g_ref(t,i+1,j,k,N) + g_ref(t,i-1,j,k,N) + 
					  	g_ref(t,i,j+1,k,N) + g_ref(t,i,j-1,k,N) +
					  	g_ref(t,i,j,k+1,N) + g_ref(t,i,j,k-1,N) - 6 * g_ref(t,i,j,k,N));
	} END_FOR

	copy3D(v,out[0],N);
	copy3D(t,out[1],N);
}


void smoke3D::exportDensity(float ***d, int n, int i) {
	char fname[128], cmd[128];
	snprintf(fname, sizeof(fname), "output/density-%04i.vol", i);
	snprintf(cmd, sizeof(fname), "gzip -f output/density-%04i.vol&", i);
	std::cout << "    + Saving \"" << fname << "\"" << std::endl;
	std::ofstream os(fname, std::ios::out | std::ios::binary);

	int xres = n, yres = n, zres = n;
	//const vector<Density> &density = fsolver->getDensity();

	FLOAT scale = 1.0f / std::max(std::max(xres, yres), zres);

	os.write("VOL", 3);
	char version = 3;
	os.write((char *)&version, sizeof(char));
	int value = 1;
	os.write((char *)&value, sizeof(int));
	os.write((char *)&xres, sizeof(int));
	os.write((char *)&yres, sizeof(int));
	os.write((char *)&zres, sizeof(int));
	value = 1;
	os.write((char *)&value, sizeof(int));

	float minX = -xres / 2.0f*scale;
	float minY = -yres / 2.0f*scale;
	float minZ = -zres / 2.0f*scale;
	float maxX = xres / 2.0f*scale;
	float maxY = yres / 2.0f*scale;
	float maxZ = zres / 2.0f*scale;

	os.write((char *)&minX, sizeof(float));
	os.write((char *)&minY, sizeof(float));
	os.write((char *)&minZ, sizeof(float));
	os.write((char *)&maxX, sizeof(float));
	os.write((char *)&maxY, sizeof(float));
	os.write((char *)&maxZ, sizeof(float));

	for (int i = 0; i < N; i++) for (int j = N-1; j >= 0; j--) for (int k = 0; k < N ; k++) {
		{
			float value = d[i][j][k];
			if (abs(value) > 1.0) {
				std::cout << std::fixed << value << std::endl;
			}
			os.write((char *)&value, sizeof(float));
			//os.write(reinterpret_cast<const char*>(&value), sizeof(float));
		}
	}
	os.close();
	//cout << "density size:" << density.size() << endl;
	//	system(cmd);
}

void smoke3D::simulateStep() {
	scene();
	diffuse();
	project();
	advection();
	enforce_boundary();
	render::render(s,/*SPHERE_R,*/N,frame);
	//smoke3D::exportDensity(s, N, frame);
	printf( "wrote frame %d\n", frame );
	frame ++;
	if( frame > LIMIT ) exit(0);
}
