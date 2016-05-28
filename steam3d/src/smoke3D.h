/*
 *  smoke3D.h
 *  smoke3D
 *
 */

#include "types.h"

namespace smoke3D {
	void init();
	void exportDensity(FLOAT ***d, int n, int i);
	void simulateStep();
}