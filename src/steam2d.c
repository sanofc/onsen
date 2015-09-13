#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define START_FOR(X,Y)		for(int j=0;j<Y;j++){for(int i=0;i<X;i++){
#define END_FOR				}}
#define SWAP(i,j)			{double ** tmp=i;i=j;j=tmp;}
#define MIN(i,j)			(i>j?j:i)
#define MAX(i,j)			(i<j?j:i)
#define CLIP(src,min,max)	MAX(min,MIN(max,src))

#define	DT					0.1			//Time Step
#define M					1.0			//Grid Size
#define	N					100			//Cell Numbers
#define W					512			//Window Size
#define T					200			//Iteration Numbers
#define A					0.1			//Ambient Temperature
const double H = M/N;					//Cell Size

int mx_prev,my_prev,mstat;				//Mouse
double **s,**s_swap;					//Steam
double **t,**t_swap;					//Temperature
double **f[2];							//Force (0:x,1:y)
double **u, **v, **u_swap, **v_swap;	//Velocity
double **d, **p, **o;					//Divergence,Pressure,Vorticity
bool show_g,show_v,show_f;				//Show Flag (Grid,Veclocity,Force)
int show_c;								//Show Flag (1:Pressure,2:Divergence,3:Vorticity)


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
	
	if(!s) s = alloc2d(N,N);
	if(!s_swap) s_swap = alloc2d(N,N);
	if(!f[0]) f[0] = alloc2d(N,N);
	if(!f[1]) f[1] = alloc2d(N,N);
	if(!u) u = alloc2d(N+1,N);
	if(!v) v = alloc2d(N,N+1);
	if(!u_swap) u_swap = alloc2d(N+1,N);
	if(!v_swap) v_swap = alloc2d(N,N+1);
	if(!p) p = alloc2d(N,N);
	if(!d) d = alloc2d(N,N);
	if(!o) o = alloc2d(N,N);
	if(!t) t = alloc2d(N,N);
	if(!t_swap) t_swap = alloc2d(N,N);
	
	START_FOR(N+1,N)
		u[i][j] = 0.0;
		u_swap[i][j] = 0.0;
	END_FOR
	
	START_FOR(N,N+1)
		v[i][j] = 0.0;
		v_swap[i][j] = 0.0;
	END_FOR
	
	START_FOR(N,N)
		p[i][j] = 0.0;
		d[i][j] = 0.0;
		o[i][j] = 0.0;
		s[i][j] = 0.0;
		s_swap[i][j] = 0.0;
		f[0][i][j] = 0.0;
		f[1][i][j] = 0.0;
		t[i][j] = A;
		t_swap[i][j] = A;
	END_FOR
	
	show_g = false;
	show_v = false;
	show_f = false;
	show_c = 0;
}
void scene(){
	//usleep(1000);
	for(int i = N/8; i < N-N/8; i++){
		s[i][0] += 0.02;
		t[i][0] = 40;
	}
	
	//Compute Buoyancy
	double a = 0.001;
	double b = 0.0002;
	START_FOR(N,N)
		double buoy = -a * s[i][j] + b * (t[i][j]-A);
		f[1][i][j] += buoy;
	END_FOR
}
/*
void scene(){
	usleep(1000);
}
*/


void final(){
	
	if(u)free(u);u=NULL;
	if(v)free(v);v=NULL;
	if(u_swap)free(u_swap);u_swap=NULL;
	if(v_swap)free(v_swap);v_swap=NULL;
	if(p)free(p);p=NULL;
	if(d)free(d);d=NULL;
	if(o)free(o);o=NULL;
	if(s)free(s);s=NULL;
	if(s_swap)free(s_swap);s_swap=NULL;
	if(t)free(t);t=NULL;
	if(t_swap)free(t_swap);t_swap=NULL;
	if(f[0])free(f[0]);f[0]=NULL;
	if(f[1])free(f[1]);f[1]=NULL;
	
}

double u_ref(int i, int j){
	return u[CLIP(i,0,N)][CLIP(j,0,N-1)];
}

double v_ref(int i, int j){
	return v[CLIP(i,0,N-1)][CLIP(j,0,N)];
}

double s_ref(int i, int j){
	return s[CLIP(i,0,N-1)][CLIP(j,0,N-1)];
}

double p_ref(int i, int j){
	return p[CLIP(i,0,N-1)][CLIP(j,0,N-1)];
}

void enforce_boundary(){
	START_FOR(N+1,N)
		if(i==0 || i==N) u[i][j] = 0;
	END_FOR
	START_FOR(N,N+1)
		if(j==0 || j==N) v[i][j] = 0;
	END_FOR
}

void compute_force(){
	START_FOR(N,N)
		if(i==0)continue;
		u[i][j] += f[0][i][j] + f[0][i-1][j];
	END_FOR
	START_FOR(N,N)
		if(j==0)continue;
		v[i][j] += f[1][i][j] + f[1][i][j-1];
	END_FOR

}

void init_force(){
	START_FOR(N,N)
		f[0][i][j]=0;
		f[1][i][j]=0;
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

void advection(){
	
	//Compute Fluid Velocity At Each Staggerd Faces And Cell Centers
	double **uy = alloc2d(N+1,N);
	double **vx = alloc2d(N,N+1);
	double **cx = alloc2d(N,N);
	double **cy = alloc2d(N,N);
	
	START_FOR(N+1,N)
		uy[i][j] = (v_ref(i-1,j) + v_ref(i-1,j+1) + v_ref(i,j+1) + v_ref(i,j)) * 0.25;
	END_FOR
	
	START_FOR(N,N+1)
		vx[i][j] = (u_ref(i,j-1) + u_ref(i,j) + u_ref(i+1, j) + u_ref(i+1,j-1)) * 0.25;
	END_FOR
	
	START_FOR(N,N)
	    cx[i][j] = (u[i][j] + u[i+1][j]) * 0.5;
	    cy[i][j] = (v[i][j] + v[i][j+1]) * 0.5;
	END_FOR
	
	semiLagragian(u,u,uy,N+1,N,u_swap);
	semiLagragian(v,vx,v,N,N+1,v_swap);
	semiLagragian(s,cx,cy,N,N,s_swap);
	semiLagragian(t,cx,cy,N,N,t_swap);
	
	free2d(uy,N+1);
	free2d(vx,N);
	free2d(cx,N);
	free2d(cy,N);
	
	SWAP(u,u_swap);
	SWAP(v,v_swap);
	SWAP(s,s_swap);
	SWAP(t,t_swap);
}

void compute_divergence(){
	START_FOR(N, N)
	d[i][j] = ((u[i+1][j]-u[i][j])+(v[i][j+1]-v[i][j]))/H;
	END_FOR
}

void compute_pressure(){
	double h2 = H * H;
	for(int k=0; k < T; k++){
		START_FOR(N,N)
			p[i][j] = ((p_ref(i+1,j)+p_ref(i-1,j)+p_ref(i,j+1)+p_ref(i, j-1))-h2*d[i][j])/4;
		END_FOR
	}
}

void subtract_pressure(){
	START_FOR(N+1,N)
	if(i>0 && i<N) u[i][j] -= (p[i][j]-p[i-1][j])/H;
	END_FOR
	START_FOR(N,N+1)
	if(j>0 && j<N) v[i][j] -= (p[i][j]-p[i][j-1])/H;
	END_FOR
}

void computeStep(){
	scene();
	enforce_boundary();
	vorticity_confinement();
	compute_force();
	advection();
	compute_divergence();
	compute_pressure();
	subtract_pressure();
}

void computePostDisplay(){
	init_force();
}

void idle(){
	computeStep();
	glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y){
	
	switch(key){
		case 'q':final();exit(0);break;
		case 'g':show_g=show_g?false:true;break;
		case 'v':show_v=show_v?false:true;break;
		case 'f':show_f=show_f?false:true;break;
		case 'd':show_c=(show_c==1)?0:1;break;
		case 'p':show_c=(show_c==2)?0:2;break;
		case 'o':show_c=(show_c==3)?0:3;break;
		case 'c':final();init();break;
	}
}

void mouse(int button, int state, int x, int y){
	mx_prev = x;
	my_prev = y;
	mstat = state;
}

void motion(int x, int y){
	
	if(mstat == GLUT_DOWN){
		double cx = x/(double)W;
		double cy = M-y/(double)W;
		double dx = (x-mx_prev)/(double)W;
		double dy = -(y-my_prev)/(double)W;
		int i = CLIP(cx*N, 0, N-1);
		int j = CLIP(cy*N, 0, N-1);
		f[0][i][j]+=CLIP(M*N*dx, -M/N/DT,M/N/DT);
		f[1][i][j]+=CLIP(M*N*dy, -M/N/DT,M/N/DT);
		
		int r = N/64;
		for(int u=-r; u<=r; u++){
			for(int v=-r; v<=r; v++){
				if(hypot(u,v)<=r){
					s[i+u][j+v] = 2.0;
				}
			}
		}
	}
	
	mx_prev = x;
	my_prev = y;
	
}

void drawBitmapString(const char *string){
	while(*string){
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *string++);
	}
}

void displayText(){
	glColor4f(1.0,1.0,1.0,1.0);
	glRasterPos2d(0.05,0.95);
	drawBitmapString("Press \"q\" to close window");
}

void displayConcentration(double **c){
	START_FOR(N-1,N-1)
	double pos[2]={i*H+H*0.5,j*H+H*0.5};
	double s = 2000;
	
	glBegin(GL_QUADS);
	glColor4d(c[i][j]>0,0,c[i][j]<0,fabs(c[i][j]*s));
	glVertex2d(pos[0],pos[1]);
	glColor4d(c[i+1][j]>0,0,c[i+1][j]<0,fabs(c[i+1][j]*s));
	glVertex2d(pos[0]+H,pos[1]);
	glColor4d(c[i+1][j+1]>0,0,c[i+1][j+1]<0,fabs(c[i+1][j+1]*s));
	glVertex2d(pos[0]+H,pos[1]+H);
	glColor4d(c[i][j+1]>0,0,c[i][j+1]<0,fabs(c[i][j+1]*s));
	glVertex2d(pos[0],pos[1]+H);
	glEnd();
	END_FOR
}

void displaySteam(){
	START_FOR(N-1,N-1)
	double pos[2]={i*H+H*0.5,j*H+H*0.5};
	glBegin(GL_QUADS);
	glColor4d(s[i][j],s[i][j],s[i][j],s[i][j]);
	glVertex2d(pos[0],pos[1]);
	glColor4d(s[i+1][j],s[i+1][j],s[i+1][j],s[i+1][j]);
	glVertex2d(pos[0]+H,pos[1]);
	glColor4d(s[i+1][j+1],s[i+1][j+1],s[i+1][j+1],s[i+1][j+1]);
	glVertex2d(pos[0]+H,pos[1]+H);
	glColor4d(s[i][j+1],s[i][j+1],s[i][j+1],s[i][j+1]);
	glVertex2d(pos[0],pos[1]+H);
	glEnd();
	END_FOR
}

void displayGrid(){
	glColor4d(0.4,0.4,0.4,0.5);
	for (int i=0; i<N+1; i++){
		glBegin(GL_LINES);
		glVertex2d(0.0, H*i);
		glVertex2d(1.0, H*i);
		glEnd();
	}
	for (int i=0; i<N+1; i++){
		glBegin(GL_LINES);
		glVertex2d(H*i, 0.0);
		glVertex2d(H*i, 1.0);
		glEnd();
	}
}

void displayVelocity(){
	glColor4d(1.0,1.0,0.0,1.0);
	START_FOR(N, N)
	double pos[2]={i*H+H*0.5,j*H+H*0.5};
	double vel[2]={(u[i][j]+u[i+1][j])*0.5,(v[i][j]+v[i][j+1])*0.5};
	glBegin(GL_LINES);
	glVertex2d(pos[0],pos[1]);
	glVertex2d(pos[0]+vel[0],pos[1]+vel[1]);
	glEnd();
	END_FOR
}

void displayForce(){
	glColor4d(1.0,0.0,0.0,1.0);
	START_FOR(N, N)
	double pos[2]={i*H+H*0.5,j*H+H*0.5};
	glBegin(GL_LINES);
	glVertex2d(pos[0],pos[1]);
	glVertex2d(pos[0]+f[0][i][j]*100,pos[1]+f[1][i][j]*100);
	glEnd();
	END_FOR
}

void display(){
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	displaySteam();
	if(show_c==1)displayConcentration(d);
	if(show_c==2)displayConcentration(p);
	if(show_c==3)displayConcentration(o);
	if(show_g)displayGrid();
	if(show_v)displayVelocity();
	if(show_f)displayForce();
	displayText();
	glutSwapBuffers();
	computePostDisplay();
}

void reshape(int w, int h){
	glViewport(0,0,w,h);
	glLoadIdentity();
	glOrtho(0.0,M,0.0,M,-1.0,1.0);
}

int main(int argc, char * argv[]){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(W,W);
	glutCreateWindow("steam2d");
	glutIdleFunc(idle);
	glutKeyboardFunc(keyboard);	
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	init();
	glutMainLoop();
	return 0;
}