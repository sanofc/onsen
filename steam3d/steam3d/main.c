#include <stdio.h>
#include <stdlib.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

typedef struct {
  int width;
  int height;
  char* title;
  float field_of_view_angle;
  float z_near;
  float z_far;
} glutWindow;

typedef struct {
  int down;
  int pos[2];
  float angle[2];
} glutMouse;

glutWindow win;
glutMouse m;

float g_rotation = 0;
float g_rotation_speed = 0.2f;

void display() {
  // Clear Screen and Depth Buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  
  // Define a viewing transformation
  gluLookAt(4, 2, 0, 0, 0, 0, 0, 1, 0);
  
  // Push and pop the current matrix stack.
  // This causes that translations and rotations on this matrix wont influence
  // others.
  
  glPushMatrix();
  glColor3f(1, 0, 0);
  
  // Rotate the teapot
  glRotated(90, 0.0, 1.0, 0.0);
  glRotated(m.angle[0], 0.0, 1.0, 0.0);
  glRotated(m.angle[1], 1.0, 0.0, 0.0);
  
  // Draw the teapot
  glutSolidTeapot(1);
  glPopMatrix();
  
  g_rotation += g_rotation_speed;
  glutSwapBuffers();
}

void initialize() {
  // select projection matrix
  glMatrixMode(GL_PROJECTION);
  
  // set the viewport
  glViewport(0, 0, win.width, win.height);
  
  // set matrix mode
  glMatrixMode(GL_PROJECTION);
  
  // reset projection matrix
  glLoadIdentity();
  GLfloat aspect = (GLfloat)win.width / win.height;
  
  // set up a perspective projection matrix
  gluPerspective(win.field_of_view_angle, aspect, win.z_near, win.z_far);
  
  // specify which matrix is the current matrix
  glMatrixMode(GL_MODELVIEW);
  glShadeModel(GL_SMOOTH);
  
  // specify the clear value for the depth buffer
  glClearDepth(1.0f);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_DEPTH_TEST);
  
  // specify implementation-specific hints
  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
  
  GLfloat amb_light[] = {0.1, 0.1, 0.1, 1.0};
  GLfloat diffuse[] = {0.6, 0.6, 0.6, 1};
  GLfloat specular[] = {0.7, 0.7, 0.3, 1};
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, amb_light);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);
  glShadeModel(GL_SMOOTH);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
  glEnable(GL_LIGHTING);
  glClearColor(0.0, 0.0, 0.0, 1.0);
}

void keyboard(unsigned char key, int mousePositionX, int mousePositionY) {
  switch (key) {
    case 'q':
      exit(0);
      break;
      
    default:
      break;
  }
}

void mouse(int button, int state, int x, int y) {
  if (button == GLUT_LEFT_BUTTON) {
    if (state == GLUT_DOWN) {
      m.pos[0] = x;
      m.pos[1] = y;
      m.down = 1;
    }
    if (state == GLUT_UP) {
      m.down = 0;
    }
  }
}

void motion(int x, int y) {
  int mx, my;
  
  if (m.down) {
    mx = (double)(x - m.pos[0]);
    my = (double)(y - m.pos[1]);
    m.pos[0] = x;
    m.pos[1] = y;
    m.angle[0] += mx;
    m.angle[1] += my;
    glutPostRedisplay();
  }
}

int main(int argc, char** argv) {
  // set window values
  win.width = 640;
  win.height = 480;
  win.title = "teapot";
  win.field_of_view_angle = 45;
  win.z_near = 1.0f;
  win.z_far = 500.0f;
  
  // set mouse values;
  m.down = 0;
  m.angle[0] = m.angle[1] = 0.0;
  
  // initialize and run program
  glutInit(&argc, argv);  // GLUT initialization
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);  // Display Mode
  glutInitWindowSize(win.width, win.height);                 // set window size
  glutCreateWindow(win.title);                               // create Window
  glutDisplayFunc(display);    // register Display Function
  glutIdleFunc(display);       // register Idle Function
  glutKeyboardFunc(keyboard);  // register Keyboard Handler
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  initialize();
  glutMainLoop();  // run GLUT mainloop
  return 0;
}
