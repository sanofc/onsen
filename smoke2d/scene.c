//
//  scene.c
//  onsen
//
//  Created by Hiroyuki Sano on 9/20/15.
//  Copyright © 2015 Hiroyuki Sano. All rights reserved.
//

#include "scene.h"


void scene(){
	
	double **s = get_steam();
	double **t = get_temperature();
	
	for(int i = X/4; i < X-X/4; i++){
		s[i][0] = 0.5;
		t[i][0] = 12;
	}
	
	//Compute Buoyancy
	double a = 0.0001;
	double b = 0.0003;
	START_FOR_C
	double buoy = -a * s[i][j] + b * (t[i][j]-A);
	add_force(i, j, 0 , buoy);
	END_FOR
}
/*
void scene(){



}*/