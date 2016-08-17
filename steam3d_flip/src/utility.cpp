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


FLOAT length(FLOAT p0[3], FLOAT p1[3]) {
	return hypotf(hypotf(p0[0] - p1[0], p0[1] - p1[1]), p0[2] - p1[2]);
}

FLOAT random() {
	return (rand() % 100) / (double)100;
}


static FLOAT spline_cubic(const FLOAT a[4], FLOAT x) {
	int i, j;
	FLOAT alpha[4], l[4], mu[4], z[4];
	FLOAT b[4], c[4], d[4];
	for (i = 1; i < 3; i++) {
		alpha[i] = 3.0 * (a[i + 1] - a[i]) - 3.0 * (a[i] - a[i - 1]);
	}
	l[0] = 1.0;
	mu[0] = 0.0;
	z[0] = 0.0;
	for (i = 1; i < 3; i++) {
		l[i] = 4.0 - mu[i - 1];
		mu[i] = 1.0 / l[i];
		z[i] = (alpha[i] - z[i - 1]) / l[i];
	}
	l[3] = 1.0;
	z[3] = 0.0;
	c[3] = 0.0;
	for (j = 2; 0 <= j; j--) {
		c[j] = z[j] - mu[j] * c[j + 1];
		b[j] = a[j + 1] - a[j] - (c[j + 1] + 2.0 * c[j]) / 3.0;
		d[j] = (c[j + 1] - c[j]) / 3.0;
	}

	FLOAT minv = MIN(a[1], a[2]);
	FLOAT maxv = MAX(a[2], a[1]);
	return MIN(maxv, MAX(minv, (a[1] + b[1] * x + c[1] * x * x + d[1] * x * x * x)));
}

static FLOAT spline(FLOAT ***d, int width, int height, int depth, FLOAT x, FLOAT y, FLOAT z) {
	FLOAT f[16];
	FLOAT xn[4];
	FLOAT zn[4];

	x = MAX(0.0, MIN(width, x));
	y = MAX(0.0, MIN(height, y));
	z = MAX(0.0, MIN(depth, z));

	for (int k = 0; k<4; k++) {
		for (int j = 0; j<4; j++) {
			for (int i = 0; i<4; i++) {
				int h = MAX(0, MIN(width - 1, (int)x - 1 + i));
				int v = MAX(0, MIN(height - 1, (int)y - 1 + j));
				int g = MAX(0, MIN(depth - 1, (int)z - 1 + k));
				f[4 * j + i] = d[h][v][g];
			}
		}

		for (int j = 0; j<4; j++) {
			xn[j] = spline_cubic(&f[4 * j], x - (int)x);
		}
		zn[k] = spline_cubic(xn, y - (int)y);
	}

	return spline_cubic(zn, z - (int)z);
}

FLOAT interp(FLOAT ***d, int width, int height, int depth, FLOAT x, FLOAT y, FLOAT z) {
	return spline(d, width, height, depth, x, y, z);
}
