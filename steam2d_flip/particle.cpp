#include "particle.h"
#include "fluid2d.h"

// x position in grid
int particle::x() {
	return (int)(p[0] / H);
}

// y position in grid
int particle::y() {
	return (int)(p[1] / H);
}

bool particle::at(int i, int j) {
	if (x() == i && y() == j) {
		return true;
	}
	else {
		return false;
	}
}
/*
void particle::compute_position(double **grid_u, double **grid_v) {
	double px = p[0] + u(grid_u);
	double py = p[1] + v(grid_v);
	p[0] = px;
	p[1] = py;
}
*/
/*
double particle::u(double **grid_u) {
	//TODO: interpolate speed
	int px = x();
	int py = y();

	//px * H
	//px * H - p[0]

	return grid_u[px][py];
}

double particle::v(double **grid_v) {
	int px = x();
	int py = y();
	return grid_v[px][py];
}
*/