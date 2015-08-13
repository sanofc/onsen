#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define IX(i,j) ((i)+(N+2)*(j))
#define IV(i,j) ((i)+(N+2)*(j))
#define IU(i,j) ((i)+(N+3)*(j))
#define GRID_SIZE (N+2)*(N*2)
#define VEL_SIZE (N+3)*(N*2)
#define SWAP(x0,x) {double * tmp=x0;x0=x;x=tmp;}
typedef struct _vec{
	double x;
	double y;
} vec;

extern void vel_step(int N, double dt, double * u, double * v, double * u_prev, double * v_prev, vec * force);
