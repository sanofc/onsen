#pragma once


class particle {

public:
	double p[2];
	double u[2];
	char type;
	double temperature;
	double dens;
	int x();
	int y();
	bool at(int x, int y);
	void compute_position(double **grid_u, double **grid_v);
};
