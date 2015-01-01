#ifndef GLUTSTROKEFONT_H
#define GLUTSTROKEFONT_H

#include "Font.h"

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif


class GlutStrokeFont : public Font{
private:
  double fHeight;
public:
  GlutStrokeFont():Font(){
    fHeight = 119.05;
  };


  void renderString(std::string text, int x, int y, int z=0, bool vertical=false){
    glLineWidth(1);
    glPointSize(1);
    glPushMatrix();
    glTranslatef(x,y,z);
    if(vertical){
     glRotatef(90, 0, 0, -1);
    }
    else{
     //glRotatef(180, 0, 0, -1);
    }
    double s = size/fHeight;
    glScalef(s,-s,s);
    for(int i=0; i<text.size(); i++){
      char c = text[i];
      //int w = glutStrokeWidth(GLUT_STROKE_MONO_ROMAN, c);
      glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, c);
      //glTranslatef(-w, 0, 0);
    }
    glPopMatrix();
  };

};

#endif
