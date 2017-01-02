#include "fluid2d.h"
#include "particle.h"


double **qs,**qs_swap;					//Steam
double **qv,**qv_swap;					//Vapor
double **rs;							//Steam for Rendering
double **f[2];							//Force (0:x,1:y)
double **u, **v, **u_swap, **v_swap;	//Velocity
double **d, **p, **o;					//Divergence,Pressure,Vorticity
double **t, **t_swap;					//Temperature
vector<particle *> particles;


double ** alloc2d(int n, int m) {
	double ** ptr = (double **)malloc(sizeof(double *) * n);
	for(int i = 0; i < n; i++){
		ptr[i] = (double *)malloc(sizeof(double) * m);
	}
	return ptr;
}

void free2d(double **ptr, int n){
	for(int i=0; i < n; i++){
		free(ptr[i]);
	}
	free(ptr);
}

void fluid2d::init(){
	
	if(!qs) qs = alloc2d(X,Y);
	if (!rs) rs = alloc2d(X, Y);
	if(!qs_swap) qs_swap = alloc2d(X,Y);
	if(!qv) qv = alloc2d(X,Y);
	if(!qv_swap) qv_swap = alloc2d(X,Y);
	if(!f[0]) f[0] = alloc2d(X,Y);
	if(!f[1]) f[1] = alloc2d(X,Y);
	if(!u) u = alloc2d(X+1,Y);
	if(!v) v = alloc2d(X,Y+1);
	if(!u_swap) u_swap = alloc2d(X+1,Y);
	if(!v_swap) v_swap = alloc2d(X,Y+1);
	if(!d) d = alloc2d(X,Y);
	if(!p) p = alloc2d(X,Y);
	if(!o) o = alloc2d(X,Y);
	if(!t) t = alloc2d(X,Y);
	if(!t_swap) t_swap = alloc2d(N,N);
	
	START_FOR_X
		u[i][j] = 0.0;
		u_swap[i][j] = 0.0;
	END_FOR
	
	START_FOR_Y
		v[i][j] = 0.0;
		v_swap[i][j] = 0.0;
	END_FOR
	
	START_FOR_C
		d[i][j] = 0.0;
		p[i][j] = 0.0;
		o[i][j] = 0.0;
		qs[i][j] = 0.0;
		qs_swap[i][j] = 0.0;
		rs[i][j] = 0.0;
		qv[i][j] = 0.0;
		qv_swap[i][j] = 0.0;
		f[0][i][j] = 0.0;
		f[1][i][j] = 0.0;
		t[i][j] = A;
		t_swap[i][j] = A;
	END_FOR

	/*
	START_FOR(X, Y / 3)
		particle *p = new particle;
		p->p[0] = (double)i / (double)X;
		p->p[1] = (double)j / (double)Y;
		p->u[0] = 0.0;
		p->u[1] = 0.0;
		p->dens = 1.0;
		p->temperature = 50.0;
		p->type = WATER;
		particles.push_back(p);
	END_FOR
	*/
}

void final(){
	
	if(u)free(u);u=NULL;
	if(v)free(v);v=NULL;
	if(u_swap)free(u_swap);u_swap=NULL;
	if(v_swap)free(v_swap);v_swap=NULL;
	if(p)free(p);p=NULL;
	if(d)free(d);d=NULL;
	if(o)free(o);o=NULL;
	if(qs)free(qs);qs=NULL;
	if(qs_swap)free(qs_swap);qs_swap=NULL;
	if (rs)free(rs);rs = NULL;
	if(qv)free(qv);qv=NULL;
	if(qv_swap)free(qv_swap);qv_swap=NULL;
	if(f[0])free(f[0]);f[0]=NULL;
	if(f[1])free(f[1]);f[1]=NULL;
	if(t)free(t);t=NULL;
	if(t_swap)free(t_swap);t_swap=NULL;
	
}

int PERIOD(int src, int min, int max){
	if(src < min){
		src = max - (min-src);
	}else if(src > max){
		src = min + (src - max);
	}
	return src;
}

double g_ref(double **g, int i, int j){ 
	double ref;
	if(g==p || g==t || g==t_swap || g==qv || g==qv_swap) ref = g[CLIP(i,0,X-1)][CLIP(j,0,Y-1)];
	else if(g==u)	 ref = g[CLIP(i,0,X)][CLIP(j,0,Y-1)];
	else if(g==v)    ref = g[CLIP(i,0,X-1)][CLIP(j,0,Y)];
	return ref;
}

/*
 * 
 */
 /*
double g_ref_for_diffuse(double **g, int i, int j){ 
	double ref;
	if(g==p || g==t || g==t_swap || g==qv || g==qv_swap)  ref = g[CLIP(i,1,X-1)][CLIP(j,0,Y-1)];
	return ref;
}
*/

void add_force(int x, int y, double fx, double fy){
	f[0][x][y] += fx;
	f[1][x][y] += fy;
}

void add_steam(int x, int y, double d){

	qs[x][y] += d;
}

void add_particle(particle *p) {
	particles.push_back(p);
}

double ** get_steam(){
	return rs;
}

double ** get_vapor(){
	return qv;
}

double ** get_temperature(){
	return t;
}

double ** get_u(){
	return u;
}

double ** get_v(){
	return v;
}

vector<particle *> get_particles() {
	return particles;
}

void enforce_boundary(){
/*
	START_FOR_X
		if(i==0 || i==X) u[i][j] = 0;
	END_FOR/
	START_FOR_Y
		if(j==0 || j==Y) v[i][j] = 0;
	END_FOR
*/
	//Outflow boundary condition
	START_FOR_X
		if(i==0) u[i][j] = -0.001;
		if(i==X) u[i][j] = 0.001;
	END_FOR
	START_FOR_Y
		if (j == 0) {
			v[i][j] = 0.001;
			u[i][j] = 0.0;
		}
		if(j==Y) v[i][j] = 0.001;
	END_FOR

	START_FOR_C
		if(j==Y || i==0 || i==X) {
			qv[i][j] = 0.0;
			qs[i][j] = 0.0;
			t[i][j] = A;
		}
	END_FOR
	
	auto particle = particles.begin();
	while(particle != particles.end()) {	
		if ((*particle)->y() ==0 || (*particle)->p[0] < 0.0 || (*particle)->p[0] > M || (*particle)->p[1] < 0.0 || (*particle)->p[1] > M) {
			particle = particles.erase(particle);
		} else {
			++particle;
		}
	}
}

void vorticity_confinement(){
	
	//Compute Vorticity
	START_FOR_C
	o[i][j] = ((g_ref(v,i+1,j)-g_ref(v,i-1,j))-(g_ref(u,i,j+1)-g_ref(u,i,j-1))) * 0.5 * H;
	END_FOR
	
	START_FOR_C
	if(i==0 || i==X-1 || j==0 || j==Y-1) continue;
	
	//Compute Gradient of Vorticity Length
	double h[2]= {
		(fabs(o[i+1][j])-fabs(o[i-1][j])) * 0.5 * H,
		(fabs(o[i][j+1])-fabs(o[i][j-1])) * 0.5 * H
	};
	double h_len = hypot(h[0],h[1]);
	if(h_len > 0){
		double e = 70;
		double n[2] = { h[0]/h_len, h[1]/h_len};
		f[0][i][j] += e * H * (n[1] * o[i][j]);
		f[1][i][j] += e * H * -(n[0] * o[i][j]);
	}
	END_FOR
}


void compute_buoyancy() {
	START_FOR_C
	//calculate buoyancy

	double t_amb;
	t_amb = A;
	/*
	if(j == 1) {
	t_amb = (g_ref(t,i-1,j)+g_ref(t,i+1,j)+g_ref(t,i,j+1))/3.0;
	}else{
	t_amb = (g_ref(t,i-1,j)+g_ref(t,i+1,j)+g_ref(t,i,j-1)+g_ref(t,i,j+1))/4.0;
	}*/
	//t_amb = g_ref(t, i - 1, j) + g_ref(t, i + 1, j) / 2.0;

	double buoy = K * /*qv[i][j] **/ ((t[i][j] - t_amb) / t_amb);// -G * qs[i][j];

																 //double buoy  = G * s[i][j];
																 //double noise_x = ((double)(rand() % 100) / 100.0 - 0.5) * 0.0;
	double noise_x = 0;
	//double noise_y = ((double)(rand() % 100) / 100.0 - 0.5) * 0.0;
	double noise_y = 0;
	add_force(i, j, 0 + noise_x, buoy + noise_y);

	/*if (qs[i][j]>0.0001) {
	printf("i%d j%d ds%f s%f v%f m%f t%f buoy%f noise_x%f noise_y%f\n", i, j, ds, qs[i][j], v[i][j], m, t[i][j], buoy, noise_x, noise_y);
	}*/
	END_FOR
}


void compute_force(){
	START_FOR_C
		if(i==0)continue;
		u[i][j] += f[0][i][j] + f[0][i-1][j];
	END_FOR
	START_FOR_C
		if(j==0)continue;
		v[i][j] += f[1][i][j] + f[1][i][j-1];
	END_FOR
}

double interpolation(double **src, double x, double y, int width, int height){

	x = CLIP(x,0.0,width);
	y = CLIP(y,0.0,height);
	
	int i = MIN(x,width-2);
	int j = MIN(y,height-2);
	
	return ((i+1-x)*src[i][j]+(x-i)*src[i+1][j])*(j+1-y)+
	       ((i+1-x)*src[i][j+1]+(x-i)*src[i+1][j+1])*(y-j);
}

void semiLagragian(double **src, double **u, double **v,int width, int height, double **out){
	START_FOR(width,height)
		out[i][j]=interpolation(src, i-N*u[i][j]*DT, j-N*v[i][j]*DT, width, height);
	END_FOR
}

void compute_grid_advection(){
	
	//Compute Fluid Velocity At Each Staggerd Faces And Cell Centers
	double **uy = alloc2d(X+1,Y);
	double **vx = alloc2d(X,Y+1);
	double **cx = alloc2d(X,Y);
	double **cy = alloc2d(X,Y);
	
	START_FOR_X
		uy[i][j] = (g_ref(v,i-1,j) + g_ref(v,i-1,j+1) + g_ref(v,i,j+1) + g_ref(v,i,j)) * 0.25;
	END_FOR
	
	START_FOR_Y
		vx[i][j] = (g_ref(u,i,j-1) + g_ref(u,i,j) + g_ref(u,i+1, j) + g_ref(u,i+1,j-1)) * 0.25;
	END_FOR
	
	START_FOR_C
	    cx[i][j] = (u[i][j] + u[i+1][j]) * 0.5;
	    cy[i][j] = (v[i][j] + v[i][j+1]) * 0.5;
	END_FOR
	
	semiLagragian(u,u,uy,X+1,Y,u_swap);
	semiLagragian(v,vx,v,X,Y+1,v_swap);
	//semiLagragian(qs,cx,cy,X,Y,qs_swap);
	semiLagragian(qv,cx,cy,X,Y,qv_swap);
	semiLagragian(t,cx,cy,X,Y,t_swap);
	
	free2d(uy,X+1);
	free2d(vx,X);
	free2d(cx,X);
	free2d(cy,X);
	
	SWAP(u,u_swap);
	SWAP(v,v_swap);
	//SWAP(qs,qs_swap);
	SWAP(qv,qv_swap);
	SWAP(t,t_swap);


}

void compute_particles_advection() {
	//compute particle advection
	/*for (auto p : particles) {
		p->p[0] += p->u(u);
		p->p[1] += p->v(v);
	}*/

	double b = 1.2;
	double a = 2;

	for (auto p : particles) {
		double u_air = interpolation(u, p->p[0] * N, p->p[1] * N, X + 1, Y);
		double v_air = interpolation(v, p->p[0] * N, p->p[1] * N, X, Y + 1);
		double u_rel = p->u[0] - u_air;
		double v_rel = p->u[1] - v_air;
		/*
		double drag = - a * std::abs(u_rel) * u_rel;
		double lift = - a * std::abs(v_rel) * v_rel - G * 0.001;
		drag = isnan(drag) ? 0 : drag;
		lift = isnan(lift) ? 0 : lift;
}
		*/
		double drag = -a * std::powf(std::abs(u_rel), b);
		drag = std::signbit(u_rel) ? -drag : drag;
		double lift = -a * std::powf(std::abs(v_rel), b);
		lift = std::signbit(v_rel) ? -lift : lift - G * 0.1;
		if (abs(drag) > 0) {
			//	printf("%f %f\n", drag, lift);
		}

		p->u[0] += drag;
		p->u[1] += lift;
		//p->u[0] = u_air;
		//p->u[1] = v_air;

		p->p[0] += p->u[0] * DT;
		p->p[1] += p->u[1] * DT;
	}

}

void gausseidel(double **x, double **b){
	double h2 = H * H;
	for(int k=0; k < T; k++){
		START_FOR_C
			x[i][j] = ((g_ref(x,i+1,j)+g_ref(x,i-1,j)+g_ref(x,i,j+1)+g_ref(x,i, j-1))-h2*b[i][j])/4;
		END_FOR
	}
}

void lin_solve ( double ** x, double ** x0, double a, double c )
{
	for (int k=0 ; k<T ; k++ ) {
		START_FOR_C
			x[i][j]= (x0[i][j] + a*(g_ref(x,i-1,j)+g_ref(x,i+1,j)+g_ref(x,i,j-1)+g_ref(x,i,j+1)))/c;
		END_FOR
	}
}


void compute_diffuse(){
	/*
	double a = DT * V;
	START_FOR_X
		u_swap[i][j] = g_ref(u,i,j) + a * (g_ref(u,i-1,j) + g_ref(u,i+1,j) + g_ref(u,i,j-1) + g_ref(u,i,j+1) - 4 * g_ref(u,i,j));
	END_FOR

	START_FOR_Y
		v_swap[i][j] = g_ref(v,i,j) + a * (g_ref(v,i-1,j) + g_ref(v,i+1,j) + g_ref(v,i,j-1) + g_ref(v,i,j+1) - 4 * g_ref(v,i,j));
	END_FOR
	*/

	double Dv = DT * 0.002;
	double Dt = DT * 0.001;

	//explicit method
	START_FOR_C
		//printf("%f %f\n", g_ref(t,i,j),t_swap[i][j]);
		t_swap[i][j] = g_ref(t,i,j) + Dt * (g_ref(t,i-1,j) + g_ref(t,i+1,j) + g_ref(t,i,j-1) + g_ref(t,i,j+1) - 4 * g_ref(t,i,j));
	END_FOR

	START_FOR_C
		qv_swap[i][j] = g_ref(qv,i,j) + Dv  * t_swap[i][j] * (g_ref(qv,i-1,j) + g_ref(qv,i+1,j) + g_ref(qv,i,j-1) + g_ref(qv,i,j+1) - 4 * g_ref(qv,i,j));
	END_FOR

	/*
	// implicit method
	double a;
	a = Dt * N * N;
	lin_solve(t_swap,t,a,1.0+4.0 * a);
	a = Dv * N * N;
	lin_solve(qv_swap,qv,a,1.0+4.0*a);
	*/

	SWAP(qv,qv_swap);
	SWAP(t,t_swap);
}

void compute_divergence(){
	START_FOR_C
		d[i][j] = -1.0*((u[i+1][j]-u[i][j])+(v[i][j+1]-v[i][j]))/(double)N;
	END_FOR
}

//pressure *x
//divergence *b
void solve(double **x, double **b){
	//gausseidel(x,b);
}

void compute_pressure(){
	lin_solve(p,d,1.0/F_DENS,4);
}

void subtract_pressure(){
	START_FOR_X
		if(i>0 && i<X) u[i][j] -= (p[i][j]-p[i-1][j])/H;
	END_FOR
	START_FOR_Y
		if(j>0 && j<Y) v[i][j] -= (p[i][j]-p[i][j-1])/H;
	END_FOR
}


void initPostDisplay(){
	START_FOR_C
	f[0][i][j]=0;
	f[1][i][j]=0;
	END_FOR
}

double random() {
	return (rand() % 100) / (double)100;
}

void phase_transition() {

	
	/*
	for (auto particle : particles) {
		int px = particle->x();
		int py = particle->y();
		qs[px][py] += particle->dens;
	}
	*/

	START_FOR_C

	if (j == 0) continue;

	//Saturation Vapor Content
	double a = 1.0;
	double b = 30;
	double c = -2.0;

	//double m = MIN(a * exp(-b / ((t[i][j]) + c)), qv[i][j] + qs[i][j]);

	double m = a * exp(-b/t[i][j]+c);


	//Amount of steam generated by the phase transition
	//double r = (double)(rand()%100)/100.0; //Phase transition ratio

	double r = 1.0;
	double ds = r * (qv[i][j] - m);
	double dv = 0.0;



	/*
	if (qs[i][j] <= 0) {
		auto particle = particles.begin();
		while (particle != particles.end()) {
			if ((*particle)->at(i, j)) {
				particle = particles.erase(particle);
		}
			else {
				++particle;
			}
		}
	}else{
	*/

		
		
		if (ds > 0) {
			
			int create_particles_num = ds / P_DENS;

			for (int k = 0; k < create_particles_num; k++) {
					particle *p = new particle;
					p->p[0] = (double)i / (double)X + random() * H;
					p->p[1] = (double)j / (double)Y + random() * H;
					p->u[0] = interpolation(u, p->p[0] * N, p->p[1] * N, X + 1, Y);
					p->u[1] = interpolation(v, p->p[0] * N, p->p[1] * N, X, Y + 1);
					p->dens = P_DENS;
					add_particle(p);
			}
			dv = -ds;
		}
		else if(ds < 0){
			vector<particle *> particles_in_grid;

			for (auto particle : particles) {
				if (particle->at(i, j)) {
					particles_in_grid.push_back(particle);
				}
			}
			for (auto particle : particles_in_grid) {
				if (ds < 0) {
					if (particle->dens <= abs(ds)) {
						dv += particle->dens;
						ds += particle->dens;
						particle->dens = 0.0;
					}else {
						particle->dens += ds;
						dv += ds;
						ds = 0;
					}
				}
			}
		}
		
		//Steam density
		//qs[i][j] += dv;
		qv[i][j] += dv;
		//latent heat
		t[i][j] += 0.005*dv;


	auto particle = particles.begin();
	while (particle != particles.end()) {

		//printf("dens%d x%d y%d\n", (*particle)->dens, (*particle)->x(),(*particle)->y() );
		
		if ((*particle)->dens == 0.0) {
			particle = particles.erase(particle);
		}
		else {
			++particle;
		}
	}



	END_FOR
}


/*!
* Poly6カーネル関数値の計算(2D)
* @param[in] r 距離
* @param[in] h 有効半径
* @return 関数値
*/
double KernelPoly6_2D(const double r, const double h)
{
	if (r >= 0.0 && r < h) {
		double q = h*h - r*r;
		double p = pow(h, 8);
		double a = 4.0 /( M_PI*pow(h, 8));
		return a * pow(q,3);
	}
	else {
		return 0.0;
	}
}

/*!
* Poly6カーネル関数値の計算(3D)
* @param[in] r 距離
* @param[in] h 有効半径
* @return 関数値
*/
double KernelPoly6_3D(const double &r, const double &h)
{
	if (r >= 0.0 && r < h) {
		double q = h*h - r*r;
		return 315.0 / (64.0*M_PI*pow(h, 9))*q*q*q;
	}
	else {
		return 0.0;
	}
}

void compute_particle_density() {

	
	START_FOR_C
	qs[i][j] = 0.0;
	rs[i][j] = 0.0;
	END_FOR
	

	for (auto particle : particles) {
		int px = particle->x();
		int py = particle->y();
		qs[px][py] += particle->dens;

		double cx = H * px + H * 0.5;
		double cy = H * py + H * 0.5;
		double r = hypot(cx-particle->p[0],cy-particle->p[1]);
		double h = hypot(H, H) * 0.5;
		double k = P_DENS * KernelPoly6_2D(r, h);
		rs[px][py] += k;

	}

	
	

}


void compute_step(){

	scene();

	//interaction between grid and particles
	//heat_transfer();
	compute_particle_density();
	phase_transition();

	//grid
	//vorticity_confinement();
	compute_buoyancy();
	compute_force();
	compute_divergence();
	compute_pressure();
	subtract_pressure();
	compute_grid_advection();
	compute_diffuse();

	//particle
	compute_particles_advection();
	enforce_boundary();
}



