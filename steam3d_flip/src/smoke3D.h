/*
 *  smoke3D.h
 *  smoke3D
 *
 */

#include "types.h"
#include "particle3d.h"

#define	N			32
#define DT			0.01
#define BUOY		1.5
#define SPHERE_R	0.2
#define V_SUR		1.0 + 0.0
#define T_SUR		2.0
#define T_AMB		0.1
#define M			1.0
#define P_DENS		0.1
#define F_DENS		1.0

static const bool noise_flag = false;
static const FLOAT cell_size = M / N;	//cell size


using namespace std;

namespace smoke3D {
	void init(int limit);
	void init();
	void exportDensity(FLOAT ***d, int n, int i);
	void simulateStep();
	void add_particle(particle3d *p);
}