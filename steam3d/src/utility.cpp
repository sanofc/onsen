/*
 *  utility.cpp
 *  smoke3D
 *
 */

#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

FLOAT *** alloc3D( int w, int h, int d ) {
	FLOAT *** field = new FLOAT **[w+1];
	for( int i=0; i<w; i++ ) {
		field[i] = new FLOAT*[h+1];
		for( int j=0; j<h; j++ ) {
			field[i][j] = new FLOAT[d];
			for( int k=0; k<d; k++ ) {
				field[i][j][k] = 0.0;
			}
		}
		field[i][h] = NULL;
	}
	field[w] = NULL;	
	return field;
}

void free3D( FLOAT ***ptr ) {
	for( int i=0; ptr[i]!=NULL; i++ ) {
		for( int j=0; ptr[i][j]!=NULL; j++ ) delete [] ptr[i][j];
		delete [] ptr[i];
	}
	delete [] ptr;
}

void copy3D( FLOAT ***dst, FLOAT ***src, int N ) {
	FOR_EVERY_CELL {
		dst[i][j][k] = src[i][j][k];
	} END_FOR
}

// Clamped Fetch
FLOAT g_ref( FLOAT ***x, int i, int j, int k, int n ) {
	i = MIN(MAX(0,i),n-1);
	j = MIN(MAX(0,j),n-1);
	k = MIN(MAX(0,k),n-1);
	return x[i][j][k];
}
