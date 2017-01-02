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

FLOAT *** u[3] = { NULL, NULL, NULL };
FLOAT *** b = NULL; //boundary
FLOAT *** v = NULL; //concentration(vapor)
FLOAT *** t = NULL; //temperature
FLOAT *** s = NULL; //steam
FLOAT *** rs = NULL; //steam for rendering
int frame = 0;
int limit = 100;
perlin perlin_noise;
vector<particle3d *> particles;


void smoke3D::init() {
	smoke3D::init(limit);
}

void smoke3D::init(int _limit) {
	// Allocate Memory
	u[0] = alloc3D(N+1,N,N);
	u[1] = alloc3D(N,N+1,N);
	u[2] = alloc3D(N,N,N+1);
	v = alloc3D(N,N,N);
	b = alloc3D(N,N,N);
	t = alloc3D(N,N,N);
	s = alloc3D(N,N,N);
	rs = alloc3D(N, N, N);
	limit = _limit;

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

	//circle
	
	for( int i=-w; i<=w; i++ ) {
		for( int j=-w; j<=w; j++ ) {
			if( hypot(i,j) < w ) {
				for( int k=0; k<1; k++ ) {
					
					double noise = 1.0;
					
					if (noise_flag) {
						noise = perlin_noise.noise((double)i*1.0 / (double)(N / 2), (double)j*1.0 / (double)(N / 2), (double)frame*1.0 / (double)(N / 2)) * 1.0;
						if (noise < 0) noise = 0.0;
					}
					v[(int)(N / 2) + i][k][(int)(N / 2) + j] = V_SUR * noise;
					t[(int)(N / 2) + i][k][(int)(N / 2) + j] = T_AMB + T_SUR *noise;
				}
			}
		}
	}
	

	//square
	/*
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			int k = 1;
			//FLOAT noise = 1;
			FLOAT noise = perlin_noise.noise((double)i*1.0 / (double)(N), (double)j*1.0 / (double)(N), (double)frame*1.0 / (double)(N)) * 1.0;
			v[i][k][j] = 5 * noise;
			t[i][k][j] = T_AMB + 3.0 * noise;
		}
	}
	*/

}

static void enforce_boundary() {
	// Set Boundary Velocity Zero
	FOR_EVERY_X_FLOW {
		if (i == 0 || i == N) {
			u[0][i][j][k] = 0.0;
			}
	} END_FOR
	
	FOR_EVERY_Y_FLOW {
		//if( j==0 || j==N ) u[1][i][j][k] = 0.0;
		if (j == 0) u[1][i][j][k] = 0.0;
		if (j == N) u[1][i][j][k] = 0.01;
	} END_FOR
	
	FOR_EVERY_Z_FLOW {
		if( k==0 || k==N ) u[2][i][j][k] = 0.0;
	} END_FOR

	FOR_EVERY_CELL{
		if (i == 0 || i == N || j == N || k == 0 || k == N) {
			s[i][j][k] = 0.0;
			v[i][j][k] = 0.0;
			t[i][j][k] = T_AMB;
		}
	} END_FOR

	auto particle = particles.begin();
	while (particle != particles.end()) {
		if ((*particle)->y() == 0 ||
			(*particle)->p[0] <= 0.0 ||(*particle)->p[0] >= M ||
			(*particle)->p[1] <= 0.0 || (*particle)->p[1] >= M ||
			(*particle)->p[2] <= 0.0 || (*particle)->p[2] >= M) {
			particle = particles.erase(particle);
		}else {
			++particle;
		}
	}

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

	FLOAT density = F_DENS;

	FOR_EVERY_CELL {
		if( i>0 && i<N ) u[0][i][j][k] -= (p[i][j][k]-p[i-1][j][k])/h * (1.0 / density);
		if( j>0 && j<N ) u[1][i][j][k] -= (p[i][j][k]-p[i][j-1][k])/h * (1.0 / density);
		if( k>0 && k<N ) u[2][i][j][k] -= (p[i][j][k]-p[i][j][k-1])/h * (1.0 / density);
	} END_FOR
}

static void grid_advection() {
	// Advect
	advect::advect( u, v, s, t, N, DT );
}

static FLOAT drag(FLOAT u_rel) {
	FLOAT b = 1.3;
	FLOAT a = 0.07;
	FLOAT drag = -a * std::powf(std::abs(u_rel), b);
	drag = std::signbit(u_rel) ? -drag : +drag;
	return drag;
}

static void particle_advection() {


	OPENMP_FOR for (auto p : particles) {
		//printf("%f %f %f\n", p->p[0],p->p[1],p->p[2]);
		FLOAT u_air = interp(u[0], N+1, N, N, N * p->p[0], N * p->p[1], N * p->p[2]);
		FLOAT v_air = interp(u[1], N, N+1 , N, N * p->p[0], N * p->p[1], N * p->p[2]);
		FLOAT w_air = interp(u[2], N, N, N + 1, N * p->p[0], N * p->p[1], N * p->p[2]);

		FLOAT u_rel = p->u[0] - u_air;
		FLOAT v_rel = p->u[1] - v_air;
		FLOAT w_rel = p->u[2] - w_air;

		FLOAT u_drag = drag(u_rel);
		FLOAT v_drag = drag(v_rel);
		FLOAT w_drag = drag(w_rel);

		//printf("%f %f %f\n", u_drag,v_drag,w_drag);	

		p->u[0] += u_drag * 1.0 / P_DENS;
		p->u[1] += v_drag * 1.0 / P_DENS -  0.001;
		p->u[2] += w_drag * 1.0 / P_DENS;

		//move particle
		
		p->p[0] += p->u[0] *DT;
		p->p[1] += p->u[1] *DT;
		p->p[2] += p->u[2] *DT;

		//move particle (particle velocity equals air velocity)
		/*
		p->p[0] += u_air * DT;
		p->p[1] += v_air * DT;
		p->p[2] += w_air * DT;
		*/
	}

}


void diffuse(){

	static FLOAT ***out[2] = { alloc3D(N,N,N), alloc3D(N,N,N)};

	FLOAT Dv = DT * 0.1;
	FLOAT Dt = DT * 0.1;

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

double KernelPoly6(const double &r, const double &h)
{
	if (r >= 0.0 && r < h) {
		double q = h*h - r*r;
		return 315.0 / (64.0*M_PI*pow(h, 9))*q*q*q;
	}
	else {
		return 0.0;
	}
}

void compute_particle_dens() {

	FOR_EVERY_CELL{
		s[i][j][k] = 0.0;
		rs[i][j][k] = 0.0;
	}END_FOR

	OPENMP_FOR for (auto particle : particles) {
		int px = particle->x();
		int py = particle->y();
		int pz = particle->z();
		s[px][py][pz] += particle->dens;

		FLOAT cp[3];
		cp[0] = cell_size * px + cell_size * 0.5;
		cp[1] = cell_size * py + cell_size * 0.5;
		cp[2] = cell_size * pz + cell_size * 0.5;

		const FLOAT r = length(cp,particle->p);
		const FLOAT h = cell_size;
		FLOAT k = KernelPoly6(r, h) *(FLOAT)P_DENS;
		rs[px][py][pz] += k * 0.00001;
	}
}

void add_particle(particle3d* p) {
	particles.push_back(p);
}



void phase_transition() {

	FOR_EVERY_CELL{

	if (j<1) continue;

	//Saturation Vapor Content
	FLOAT a =-1.0;
	FLOAT b = 10.0;
	FLOAT c = 0.3;
	FLOAT m = fmin(a * exp(-b / ((t[i][j][k]) + c)),v[i][j][k] + s[i][j][k]);
	//FLOAT m = a * exp(-b / ((t[i][j][k]) + c));



	//Phase transition
	FLOAT r = 0.9;
	FLOAT ds = r * (v[i][j][k] - m);
	FLOAT dv = 0.0;

	//find particles in grid
	vector<particle3d *> pig;
	for (auto p : particles) {
		pig.push_back(p);
	}

	if (ds > 0) {
		int p_num = ds / P_DENS;
		for (int n = 0; n < p_num; n++) {
			particle3d *p = new particle3d;
			p->p[0] = (double)i / (double)N + random() * cell_size;
			p->p[1] = (double)j / (double)N + random() * cell_size;
			p->p[2] = (double)k / (double)N + random() * cell_size;
			p->u[0] = interp(u[0], N + 1, N, N, p->p[0], p->p[1], p->p[2]);
			p->u[1] = interp(u[1], N, N + 1, N, p->p[0], p->p[1], p->p[2]);
			p->u[2] = interp(u[2], N, N, N + 1, p->p[0], p->p[1], p->p[2]);
			p->dens = P_DENS;
			add_particle(p);
		}
		dv = -ds;
	}else if (ds < 0) {
		vector<particle3d *> pig;
		for (auto particle : particles) {
			if (particle->at(i, j,k)) {
				pig.push_back(particle);
			}
		}
		for (auto particle : pig) {
			if (ds < 0) {
				if (particle->dens <= abs(ds)) {
					dv += particle->dens;
					ds += particle->dens;
					particle->dens = 0.0;
				}
				else {
					particle->dens += ds;
					dv += ds;
					ds = 0;
				}
			}
		}
	
	}

	//s[i][j][k] += ds;
	v[i][j][k] += dv;

	//latent heat
	t[i][j][k] += 0.05 * ds;

	auto particle = particles.begin();
	
	while (particle != particles.end()) {

		//printf("dens%d x%d y%d\n", (*particle)->dens, (*particle)->x(),(*particle)->y() );

		if ((*particle)->dens <= 0.0) {
			particle = particles.erase(particle);
		}
		else {
			++particle;
		}
	}

	} END_FOR

		
}

void buoyancy() {


	//Ambient Temperature
	FLOAT at;
	at = T_AMB;


	FOR_EVERY_CELL{

	if (j<1) continue;

	// Give Buoyancy 
	//FLOAT beta = 1.5;
	FLOAT buoy = BUOY * (t[i][j][k] - at);
	u[1][i][j][k] += buoy;

	} END_FOR
}


void smoke3D::simulateStep() {
	scene();

	compute_particle_dens();
	phase_transition();
	buoyancy();
	diffuse();
	project();
	grid_advection();
	particle_advection();
	enforce_boundary();
	render::render(rs,/*SPHERE_R,*/N,frame);
	//smoke3D::exportDensity(s, N, frame);
	printf("wrote frame %d %d\n", frame, particles.size() );
	frame ++;
	if( frame > limit ) exit(0);

}
