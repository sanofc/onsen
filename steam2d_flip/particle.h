#pragma once


class particle {

public:
	double p[2];
	double u[2];
	double dens;
	int x();
	int y();
	bool at(int x, int y);
};
