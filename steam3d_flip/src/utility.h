/*
 *  utility.h
 *  smoke3D
 *
 */

#include "types.h"

#define MAX(i,j) (i>j?i:j)
#define MIN(i,j) (i>j?j:i)

#define FOR_EVERY_X_FLOW	for( int i=0; i<N+1; i++ ) for( int j=0; j<N; j++ ) for( int k=0; k<N; k++ ) {
#define FOR_EVERY_Y_FLOW	for( int i=0; i<N; i++ ) for( int j=0; j<N+1; j++ ) for( int k=0; k<N; k++ ) {
#define FOR_EVERY_Z_FLOW	for( int i=0; i<N; i++ ) for( int j=0; j<N; j++ ) for( int k=0; k<N+1; k++ ) {
#define FOR_EVERY_CELL		for( int i=0; i<N; i++ ) for( int j=0; j<N; j++ ) for( int k=0; k<N; k++ ) {
#define END_FOR }

#ifdef _OPENMP
#include <omp.h>
#define OPENMP_FOR	__pragma("omp parallel for" )
#define OPENMP_SECTION  __pragma("omp section" )
#define OPENMP_BEGIN	__pragma("omp parallel" ) {
#define OPENMP_END		}
#define OPENMP_FOR_P	__pragma("omp for" )
#else
#define OPENMP_FOR
#define OPENMP_SECTION
#define OPENMP_BEGIN
#define OPENMP_END
#define OPENMP_FOR_P
#endif

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

FLOAT *** alloc3D( int w, int h, int d );
void free3D( FLOAT ***ptr );
void copy3D( FLOAT ***dst, FLOAT ***src, int n );
FLOAT g_ref(FLOAT ***x, int i, int j, int k, int n );
FLOAT length(FLOAT p0[3], FLOAT p1[3]);
FLOAT random();
FLOAT interp(FLOAT ***d, int width, int height, int depth, FLOAT x, FLOAT y, FLOAT z);