#include "fluid2d.h"

int show_v,show_g,show_s;		//Show Flag(Velocity,Grid)
int mx_prev,my_prev,mstat;		//Mouse

void idle(){
	compute_step();
	glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y){
	switch(key){
		case 'q':final();exit(0);break;
		case 'v':show_v = !show_v;break;
		case 'g':show_g = !show_g;break;
		case 's':show_s = show_s == 2 ? 0 : show_s+1;break;
	}
}

void mouse(int button, int state, int x, int y){
	mx_prev = x;
	my_prev = y;
	mstat = state;
}

void motion(int x, int y){
	
	double wx = W * (double)X / N;
	double wy = W * (double)Y / N;
	
	if(mstat == GLUT_DOWN){
		double cx = (double)x/wx;
		double cy = M-(double)y/wy;
		double dx = (x-mx_prev)/(double)W;
		double dy = -(y-my_prev)/(double)W;
		int i = CLIP(cx*X, 0, X-1);
		int j = CLIP(cy*Y, 0, Y-1);
		double fx =CLIP(M*N*dx, -M/N/DT,M/N/DT);
		double fy =CLIP(M*N*dy, -M/N/DT,M/N/DT);
		
		add_force(i,j,fx,fy);
		
		int r = 1;
		for(int u=-r; u<=r; u++){
			for(int v=-r; v<=r; v++){
				if(hypot(u,v)<=r){
					add_steam(i+u, j+v,2.0);
				}
			}
		}
	}
	
	mx_prev = x;
	my_prev = y;
	
}

void drawBitmapString(const char *string){
	while(*string){
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *string++);
	}
}

void displayText(){
	glColor4f(1.0,1.0,1.0,1.0);
	glRasterPos2d(0.05,0.05);
	drawBitmapString("Press \"q\" to close window");
}

void displayConcentration(double **c){
	START_FOR(X-1,Y-1)
	double pos[2]={i*H+H*0.5,j*H+H*0.5};
	double s = 2000;
	
	glBegin(GL_QUADS);
	glColor4d(c[i][j]>0,0,c[i][j]<0,fabs(c[i][j]*s));
	glVertex2d(pos[0],pos[1]);
	glColor4d(c[i+1][j]>0,0,c[i+1][j]<0,fabs(c[i+1][j]*s));
	glVertex2d(pos[0]+H,pos[1]);
	glColor4d(c[i+1][j+1]>0,0,c[i+1][j+1]<0,fabs(c[i+1][j+1]*s));
	glVertex2d(pos[0]+H,pos[1]+H);
	glColor4d(c[i][j+1]>0,0,c[i][j+1]<0,fabs(c[i][j+1]*s));
	glVertex2d(pos[0],pos[1]+H);
	glEnd();
	END_FOR
}


void displayDensity(double **s, double max){
	START_FOR(X-1,Y-1)
	double pos[2]={i*H+H*0.5,j*H+H*0.5};
	glBegin(GL_QUADS);
	glColor4d(s[i][j]/max,s[i][j]/max,s[i][j]/max,s[i][j]/max);
	glVertex2d(pos[0],pos[1]);
	glColor4d(s[i+1][j]/max,s[i+1][j]/max,s[i+1][j]/max,s[i+1][j]/max);
	glVertex2d(pos[0]+H,pos[1]);
	glColor4d(s[i+1][j+1]/max,s[i+1][j+1]/max,s[i+1][j+1]/max,s[i+1][j+1]/max);
	glVertex2d(pos[0]+H,pos[1]+H);
	glColor4d(s[i][j+1]/max,s[i][j+1]/max,s[i][j+1]/max,s[i][j+1]/max);
	glVertex2d(pos[0],pos[1]+H);
	glEnd();
	END_FOR
}

void displaySteam(){
	displayDensity(get_steam(),10.0);
}

void displayVapor(){
	displayDensity(get_vapor(),1.0);
}

void displayTemperature(){
	displayDensity(get_temperature(),15.0);
}


void displayGrid(){
	glColor4d(0.4,0.4,0.4,0.5);
	for (int i=0; i<X+1; i++){
		glBegin(GL_LINES);
		glVertex2d(0.0, H*i);
		glVertex2d((double)X/N, H*i);
		glEnd();
	}
	for (int i=0; i<Y+1; i++){
		glBegin(GL_LINES);
		glVertex2d(H*i, 0.0);
		glVertex2d(H*i, (double)Y/N);
		glEnd();
	}
}
void displayVelocity(){
	
	double ** u = get_u();
	double ** v = get_v();
	
	glColor4d(1.0,1.0,0.0,1.0);
	START_FOR_C
	double pos[2]={i*H+H*0.5,j*H+H*0.5};
	double vel[2]={(u[i][j]+u[i+1][j])*0.5,(v[i][j]+v[i][j+1])*0.5};
	glBegin(GL_LINES);
	glVertex2d(pos[0],pos[1]);
	glVertex2d(pos[0]+vel[0],pos[1]+vel[1]);
	glEnd();
	END_FOR
}

void display(){
	//usleep(10000);
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	if(show_g)displayGrid();
	if(show_v)displayVelocity();
	switch(show_s){
		case 0 : displaySteam();break;
		case 1 : displayVapor();break;
		case 2 : displayTemperature();break;
	}
	displayText();
	glutSwapBuffers();
	initPostDisplay();
}

void reshape(int w, int h){
	glViewport(0,0,w,h);
	glLoadIdentity();
	glOrtho(0.0,(double)X/N,0.0,(double)Y/N,-1.0,1.0);
}

int main(int argc, char * argv[]){
	show_v = show_g = 0;
	init();
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(W*(double)X/N,W*(double)Y/N);
	glutCreateWindow("steam2d");
	glutIdleFunc(idle);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutMainLoop();
	return 0;
}