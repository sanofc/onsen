#include <stdio.h>
#include "perlin.h"

using namespace std;

int main(){
	perlin p;
	for(int i=0; i<100; i++){
		//cout << p.noise(i,i,i) << endl;
		printf("%f\n",p.noise(i*0.01,0,0));
	}
}