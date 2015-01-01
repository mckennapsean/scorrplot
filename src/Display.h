#ifndef DISPLAY_H
#define DISPLAY_H

//include opengl stuff
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif


#include <string>

class  Display{

public:

  virtual ~Display(){};

  virtual std::string title() = 0;

  virtual void printHelp() = 0;

  virtual void init() = 0;

  virtual void reshape(int w, int h) = 0;

  virtual void display() = 0;
 
  virtual void idle(){};

  virtual void keyboard(unsigned char key, int x, int y) =0;

  virtual void special(int key, int x, int y) =0;

  virtual void mouse(int button, int state, int x, int y) =0;

  virtual void motion(int x, int y) = 0;
  
  virtual void passive(int x, int y) = 0;


};

#endif
