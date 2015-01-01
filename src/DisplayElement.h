#ifndef DISPLAYELEMENT_H
#define DISPLAYELEMENT_H

//include opengl stuff
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
//#include <FTGL/FTGL.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
//#include <FTGL/ftgl.h>
#endif


#include <string>

class  DisplayElement{

public:

  int xLeft, yTop, width, height;

  virtual ~DisplayElement(){};

  virtual void location(int x, int y, int w, int h) = 0;

  virtual void display() = 0;
 
  virtual void idle(){};
  
  virtual void reshape(int w, int h){};

  virtual void keyboard(unsigned char key, int x, int y) =0;

  virtual void special(int key, int x, int y) =0;

  virtual void mouse(int button, int state, int x, int y) =0;

  virtual void motion(int x, int y) = 0;
  
  virtual void passive(int x, int y) = 0;


  bool isInside(int x, int y){
    return x>xLeft && x < (xLeft+width) &&
           y>yTop && y < (yTop+height) ;
  };


};

#endif
