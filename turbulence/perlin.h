#pragma once

#include <cmath>

class perlin{

public:
	perlin();
	double noise(double x, double y, double z);
	double noise(double x, double y);
	double noise(double x);

private:
	int p[512];
	double fade(double d);
	double grad(int hash, double x, double y,double z);
	double lerp(double t, double a, double b);

};