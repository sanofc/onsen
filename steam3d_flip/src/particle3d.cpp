#include "particle3d.h"
#include "smoke3D.h"

int particle3d::x() {
	return (int)(p[0] / cell_size);
}

int particle3d::y() {
	return (int)(p[1] / cell_size);
}

int particle3d::z() {
	return (int)(p[2] / cell_size);
}

bool particle3d::at(int i, int j, int k) {
	if (x() == i && y() == j && z() == k) {
		return true;
	}
	else {
		return false;
	}
}