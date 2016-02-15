#include "fluid2d.h"

using namespace std;

double **qs,**qs_swap;					//Steam
double **qv,**qv_swap;					//Vapor
double **f[2];							//Force (0:x,1:y)
double **u, **v, **u_swap, **v_swap;	//Velocity
double **d, **p, **o;					//Divergence,Pressure,Vorticity
double **t, **t_swap;					//Temperature
list<double[2]> particles;

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

void init(){
	
	if(!qs) qs = alloc2d(X,Y);
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
		qv[i][j] = 0.0;
		qv_swap[i][j] = 0.0;
		f[0][i][j] = 0.0;
		f[1][i][j] = 0.0;
		t[i][j] = A;
		t_swap[i][j] = A;
	END_FOR
	
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
	if(g==p || g==t || g==t_swap || g==qv || g==qv_swap) ref = g[CLIP(i,1,X-1)][CLIP(j,1,Y-1)];
	else if(g==u)	 ref = g[CLIP(i,0,X)][CLIP(j,0,Y-1)];
	else if(g==v)    ref = g[CLIP(i,0,X-1)][CLIP(j,0,Y)];
	return ref;
}

double g_ref_for_diffuse(double **g, int i, int j){ 
	double ref;
	if(g==p || g==t || g==t_swap || g==qv || g==qv_swap)  ref = g[CLIP(i,1,X-1)][CLIP(j,0,Y-1)];
	return ref;
}

void add_force(int x, int y, double fx, double fy){
	f[0][x][y] += fx;
	f[1][x][y] += fy;
}

void add_steam(int x, int y, double d){
	qs[x][y] += d;
}

double ** get_steam(){
	return qs;
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


void enforce_boundary(){
/*
	START_FOR_X
		if(i==0 || i==X) u[i][j] = 0;
	END_FOR
	START_FOR_Y
		if(j==0 || j==Y) v[i][j] = 0;
	END_FOR
*/
	//Outflow boundary condition
	START_FOR_X
		if(i==0) u[i][j] = -0.01;
		if(i==X) u[i][j] = 0.01;
	END_FOR
	START_FOR_Y
		if(j==0) v[i][j] = 0.0;
		if(j==Y) v[i][j] = 0.4;
	END_FOR

	START_FOR_C
		if(j==Y || i==0 || i==X) {
			qv[i][j] = 0.0;
			qs[i][j] = 0.0;
			t[i][j] = A;
		}
	END_FOR
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
		double e = 100;
		double n[2] = { h[0]/h_len, h[1]/h_len};
		f[0][i][j] += e * H * (n[1] * o[i][j]);
		f[1][i][j] += e * H * -(n[0] * o[i][j]);
	}
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

void compute_advection(){
	
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
	semiLagragian(qs,cx,cy,X,Y,qs_swap);
	semiLagragian(qv,cx,cy,X,Y,qv_swap);
	semiLagragian(t,cx,cy,X,Y,t_swap);
	
	free2d(uy,X+1);
	free2d(vx,X);
	free2d(cx,X);
	free2d(cy,X);
	
	SWAP(u,u_swap);
	SWAP(v,v_swap);
	SWAP(qs,qs_swap);
	SWAP(qv,qv_swap);
	SWAP(t,t_swap);
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
			x[i][j]= (x0[i][j] + a*(g_ref_for_diffuse(x,i-1,j)+g_ref_for_diffuse(x,i+1,j)+g_ref_for_diffuse(x,i,j-1)+g_ref_for_diffuse(x,i,j+1)))/c;
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

	/*
	START_FOR_C
		qv_swap[i][j] = g_ref(qv,i,j) + Dv * (g_ref(qv,i-1,j) + g_ref(qv,i+1,j) + g_ref(qv,i,j-1) + g_ref(qv,i,j+1) - 4 * g_ref(qv,i,j));
	END_FOR
	START_FOR_C
		//printf("%f %f\n", g_ref(t,i,j),t_swap[i][j]);
		t_swap[i][j] = g_ref(t,i,j) + Dt * (g_ref(t,i-1,j) + g_ref(t,i+1,j) + g_ref(t,i,j-1) + g_ref(t,i,j+1) - 4 * g_ref(t,i,j));
	END_FOR
	*/


	// implicit method
	double a;
	a = Dt * N * N;
	lin_solve(t_swap,t,a,1.0+4.0 * a);
	a = Dv * N * N;
	lin_solve(qv_swap,qv,a,1.0+4.0*a);


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
	lin_solve(p,d,1,4);
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

void compute_step(){
	enforce_boundary();
	scene();
	vorticity_confinement();
	compute_force();
	compute_advection();
	compute_diffuse();
	compute_divergence();
	compute_pressure();
	subtract_pressure();
}



