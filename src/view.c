
#include "steam.h"

static int win_x,win_y;
static int mouse_x,mouse_x_prev;
static int mouse_y,mouse_y_prev;
static int mouse_down;
static int N;

static vec * force;
static double * u, * v, * u_prev, * v_prev;
static double dt;

int allocate_data(void){

	force = (vec *)malloc(GRID_SIZE * sizeof(vec));
	u = (double *)malloc(VEL_SIZE * sizeof(double));
	v = (double *)malloc(VEL_SIZE * sizeof(double));
	u_prev = (double *)malloc(VEL_SIZE * sizeof(double));
	v_prev = (double *)malloc(VEL_SIZE * sizeof(double));

	if(!force | !u | !v){
		fprintf(stderr, "cannnot allocate data\n");
		return 0;
	}
	return 1;
}

void init_force(void){
	for(int i=0; i<GRID_SIZE; i++){
		force[i].x = force[i].y = 0.0;
	}
	force[IX(5,5)].x = 1.5;
}

void init_velocity(void){
	for(int i=0; i<VEL_SIZE; i++){
		u[i] = v[i] = u_prev[i] = v_prev[i] = 0.0;
	}
}

void init_data(void){
	init_force();
	init_velocity();
}

void free_data(void){
	if(force) free(force);
	if(u) free(u);
	if(v) free(v);
	if(u_prev) free(u_prev);
	if(v_prev) free(v_prev);
}

void draw_velocity(void){

	double h,x,y;
	h=1.0/(N+2);

	glColor3d(1.0,1.0,1.0);
	glBegin(GL_LINES);
	for(int i=0; i < N + 3; i++){
		x = i * h;
		for(int j=0; j < N + 2; j++){
			y = (j+0.5) * h;
			glVertex2d(x,y);
			glVertex2d(x + u[IU(i,j)], y);
		}
	}
	for(int i=0; i < N + 2; i++){
		x = (i+0.5) * h;
		for(int j=0; j < N + 3; j++){
			y = j * h;
			glVertex2d(x,y);
			glVertex2d(x, y + v[IV(i,j)]);
		}
	}
	//draw backtrace
	for(int i=0; i < N + 3; i++){
		x = i * h;
		for(int j=0; j < N + 2 ; j++){
			y = (j+0.5) * h;

 			double x0 =  x - u_prev[IU(i,j)] * dt;
			double cv = (v_prev[IV(i-1,j)] + v_prev[IV(i,j)] + v_prev[IV(i-1,j+1)] + v_prev[IV(i,j+1)]) * 0.25;
			double y0 = y - cv * dt;
			glColor3d(1.0,0.0,0.0);
			glVertex2d(x,y);
			glVertex2d(x0, y0);	
		}
	}
	for(int i=0; i < N + 2; i++){
		x = (i+0.5) * h;
		for(int j=0; j < N + 3; j++){
			y = j * h;
 			double y0 =  y - v_prev[IV(i,j)] * dt;
			double cu = (u_prev[IU(i,j)] + u_prev[IU(i,j-1)] + u_prev[IU(i+1,j-1)] + u_prev[IU(i+1,j)]) * 0.25;
			double x0 = x - cu * dt;
			glColor3d(1.0,0.0,0.0);
			glVertex2d(x,y);
			glVertex2d(x0, y0);	
		}
	}

	glEnd();
}

void draw_force(void){

	double h,x,y;
	h = 1.0/(N+2);

	glBegin(GL_LINES);
	glColor3d(1.0,1.0,1.0);
	for(int i = 0; i < N + 2 ; i++){
		x = (i + 0.5) * h;
		for(int j = 0; j < N + 2; j++){
			y = (j + 0.5) * h;
			glVertex2d(x,y);
			glVertex2d(x + force[IX(i,j)].x * dt, y + force[IX(i,j)].y * dt);
		}
	}
	glEnd();
}

void draw_grid(void){

	glColor3d(0.2,0.2,0.2);
	glBegin(GL_LINES);	
	for(int i = 0; i < N + 3; i++){
		glVertex2d(0.0, i * 1.0/(N+2));
		glVertex2d(1.0, i * 1.0/(N+2));
		glVertex2d(i * 1.0/(N+2), 0.0);
		glVertex2d(i * 1.0/(N+2), 1.0);
	}
	glEnd();

}


void display(void){

	glViewport(0,0,win_x,win_y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0,1.0,0.0,1.0);
	glClear(GL_COLOR_BUFFER_BIT);
	draw_grid();
	draw_force();
	draw_velocity();

	glFlush();
}

void key(unsigned char key, int x, int y){
	switch( key ){
		case 'q':
			free_data();
			exit( 0 );
			break;
	}
}

void mouse(int button, int state, int x, int y){
	mouse_x_prev = mouse_x = x;
	mouse_y_prev = mouse_y = y;
	mouse_down = state == GLUT_DOWN;
}

void motion(int x, int y){
	mouse_x = x;
	mouse_y = y;
}

void get_from_UI(void){
	int i,j;


	if(mouse_down){
		//printf("%d %d\n",mouse_x_prev,mouse_x);
		fflush(stdout);
		i = (int)((mouse_x / (float)win_x) * (N + 2));
		j = (int)(((win_y - mouse_y)/(float)win_y) * (N+2));

		if(i == 0 || i == N+1 || j ==0 || j == N+1) return;

		force[IX(i,j)].x = mouse_x - mouse_x_prev;
		force[IX(i,j)].y = mouse_y_prev - mouse_y;
	}
	mouse_x_prev = mouse_x;
	mouse_y_prev = mouse_y;
}

void idle(void){
	get_from_UI();
	//sleep(1);
	SWAP(u,u_prev);SWAP(v,v_prev);
	vel_step(N,dt,u,v,u_prev,v_prev,force);
	glutPostRedisplay();

}

int main( int argc, char ** argv )
{
	glutInit( &argc, argv);

	win_x=512;
	win_y=512;
	N=20;
	dt = 0.1;

	glutInitDisplayMode(GLUT_RGBA);
	glutInitWindowPosition(glutGet(GLUT_SCREEN_WIDTH)-win_x,0);
	glutInitWindowSize(win_x,win_y);
	glutCreateWindow(argv[0]);
	glutKeyboardFunc(key);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutIdleFunc(idle);
	glutDisplayFunc(display);
	glClearColor(0.0, 0.0, 0.0, 1.0);

	if(!allocate_data()) exit (1);
	init_data();

	glutMainLoop();
	return 0;
}