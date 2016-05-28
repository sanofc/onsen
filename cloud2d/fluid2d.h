//
//  fluid2d.h
//  onsen
//
//  Created by Hiroyuki Sano on 9/17/15.
//  Copyright (c) 2015 Hiroyuki Sano. All rights reserved.
//

#ifndef onsen_fluid2d_h
#define onsen_fluid2d_h
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "scene.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define	DT					0.1			//Time Step
#define M					1.0			//Grid Size
#define	N					64			//Cell Numbers for One Grid
#define X					64			//Cell Numbers for X Direction
#define Y					64			//Cell Numbers for Y Direction
#define W					512			//Window Size
#define T					200			//Iteration Numbers
#define V					0.0001			//Viscocity
#define A					10			//Ambient Temperature


static const double H = M/N;			//Cell Size

#define START_FOR(X,Y)		for(int j=0;j<Y;j++){for(int i=0;i<X;i++){
#define START_FOR_X			START_FOR(X+1,Y)
#define START_FOR_Y			START_FOR(X,Y+1)
#define START_FOR_C			START_FOR(X,Y)
#define END_FOR				}}
#define SWAP(i,j)			{double ** tmp=i;i=j;j=tmp;}
#define MIN(i,j)			(i>j?j:i)
#define MAX(i,j)			(i<j?j:i)
#define CLIP(src,min,max)	MAX(min,MIN(max,src))

void init();
void final();
void initPostDisplay();
void compute_step();
void add_force(int x, int y, double fx, double fy);
void add_steam(int x, int y, double _s);
double ** get_steam();
double ** get_u();
double ** get_v();
double ** get_temperature();
double ** get_vapor();

#endif
