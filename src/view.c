#include <stdlib.h>
#include <stdio.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define IX(i,j) ((i)+(N+2)*(j))
#define GRID_SIZE (N+2)*(N*2)
#define VEL_SIZE (N+3)*(N*2)

typedef struct _vec{
	double x;
	double y;
} vec;

static int win_x,win_y;
static int mouse_x,mouse_x_prev;
static int mouse_y,mouse_y_prev;
static int mouse_down;
static int N;

static vec * force;
static double * u, * v;

int allocate_data(void){

	force = (vec *)malloc(GRID_SIZE * sizeof(vec));
	u = (double *)malloc(VEL_SIZE * sizeof(double));
	v = (double *)malloc(VEL_SIZE * sizeof(double));

	if(!force | !u | !v){
		fprintf(stderr, "cannnot allocate data\n");
		return 0;
	}
	return 1;
}

void init_data(void){

	for(int i=0; i<GRID_SIZE; i++){
		force[i].x = force[i].y = 0.0;
	}

	for(int i=0; i<VEL_SIZE; i++){
		u[i] = v[i] = 0.01;
	}
}

void free_data(void){
	if(force) free(force);
	if(u) free(u);
	if(v) free(v);
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
			glVertex2d(x + u[i+j*(N+3)], y + u[i+j*(N+3)]);
		}
	}
	for(int i=0; i < N + 2; i++){
		x = (i+0.5) * h;
		for(int j=0; j < N + 3; j++){
			y = j * h;
			glVertex2d(x,y);
			glVertex2d(x + u[i+j*(N+2)], y + u[i+j*(N+2)]);
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
			glVertex2d(x + force[IX(i,j)].x, y + force[IX(i,j)].y);
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


void draw_mouse_down(void){
	if(mouse_down){
		glColor3d(1.0,0.0,0.0);
		glPointSize(10);
		glBegin(GL_POINTS);
		printf("%d %d\n",mouse_x_prev,mouse_x);
		fflush(stdout);
		force->x = mouse_x - mouse_x_prev;
		force->y = mouse_y_prev - mouse_y;
		glVertex2d((double)mouse_x/win_x,(double)(win_y-mouse_y)/win_y);
		glEnd();
	}
	mouse_x_prev = mouse_x;
	mouse_y_prev = mouse_y;
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
	draw_mouse_down();
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

void idle(void){
	glutPostRedisplay();
}

int main( int argc, char ** argv )
{
	glutInit( &argc, argv);

	win_x=512;
	win_y=512;
	N=10;

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