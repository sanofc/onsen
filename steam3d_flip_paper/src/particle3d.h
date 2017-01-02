#pragma once
#include "types.h"

class particle3d {

public:
	FLOAT p[3];
	FLOAT u[3];
	FLOAT dens;
	int x();
	int y();
	int z();
	bool at(int x, int y, int z);
};
