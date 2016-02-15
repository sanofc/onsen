/*
 *  advect.h
 *  smoke3D
 *
 */

#include "types.h"

namespace advect {

	/*
	 * u : Velocity 
	 * c : Vaper Concentration
	 * s : Steam Concentration
	 * t : Temperature
	 */
	void advect( FLOAT ****u, FLOAT ***c, FLOAT ***s, FLOAT ***t, int n, FLOAT dt );
}