#include "smoke3D.h"
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char * argv[]) {
	if (argv[1] != NULL) {
		smoke3D::init(atoi(argv[1]));
	}
	else {
		smoke3D::init();
	}
	while(1) smoke3D::simulateStep();
	return 0;
}
