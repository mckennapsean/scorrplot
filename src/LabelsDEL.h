#ifndef LABELSDEL_H
#define LABELSDEL_H

#include "Display.h"
#include "DenseVector.h"
#include "DenseMatrix.h"
#include "Rotation.h"
#include "Animator.h"

#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "Font.h"

template<typename TPrecision>
class LabelsDEL : public DisplayElement{

  private:


#define BUFSIZE 512



    Data<TPrecision> &data;
    Font &font;
    Animator &animator;
    int selected;
    int xM, yM;

  public:

    LabelsDEL(Data<TPrecision> &d, Font &f, Animator &a) 
           : data(d), font(f), animator(a){ 
       selected=-1;
    };



    void location(int xPos, int yPos, int w, int h){
      width = w;
      height = h;
      xLeft = xPos;
      yTop = yPos;
    };




    void init(){  

    };


    void display(void){
      //if(animator.isRunning()){ return; }
     
      glMatrixMode(GL_MODELVIEW); 	
      glLoadIdentity();
    
      int fs = 12;
      int x = xLeft+10;
      int y = yTop+10+fs;
      // int rIndex = data.getSelected();
      // float col[4];
      
      glPushName(-1); 

      int selectedIndex = -1;
      int index = 0;
      for(typename std::map<int, PCA<TPrecision>* >::iterator it = data.lPCA.begin();
          it!=data.lPCA.end(); ++it, ++index){
         int id = it->first;
        glLoadName(id);
        float alpha=0.35;
        if(data.labelHighlight == id || data.labelSelected(id) ){
           alpha=0.85;
        }    
     


        glColor4f(data.colLabelR(id), data.colLabelG(id), 
                  data.colLabelB(id), alpha);

         
        int w = 30;
        int h = 15;
        int x  = xLeft + index * (w+10); 
        glBegin(GL_QUADS); 
        glVertex2f(x, yTop);     
        glVertex2f(x+w, yTop);     
        glVertex2f(x+w, yTop+h);     
        glVertex2f(x, yTop+h);
        glEnd();


        glColor4f(data.colLabelR(id), data.colLabelG(id), 
                  data.colLabelB(id), 1);
        glLineWidth(2);

        glBegin(GL_LINE_LOOP); 
        glVertex2f(x, yTop);     
        glVertex2f(x+w, yTop);     
        glVertex2f(x+w, yTop+h);     
        glVertex2f(x, yTop+h);
        glEnd();



        if(id == selected){
          selectedIndex = index;
        }      
      }

     if( !isInside(x, y) ) {  return; }
    
     if(selectedIndex != -1){
        glColor4f(0.8, 0.8, 0.8, 1);
        font.setSize(15);
	font.renderString(data.labelnames[selectedIndex], xM, yM, 0); 
      } 
      

      glPopName();

    };

// resize labels to the window size (manual default values set)
void reshape(int w, int h){
  int startLocX = (int) (151. * (double) w / 1100. );
  int startLocY = (int) (650. * (double) h / 710. );
  int endLocX = (int) (630. * (double) w / 1100. );
  int endLocY = (int) (40. * (double) h / 710. );
  location(startLocX, startLocY, endLocX, endLocY);
}

    void idle(){

    };


    void keyboard(unsigned char key, int x, int y){
    };



    void special(int key, int x, int y){
    };




    void mouse(int button, int state, int x, int y){
      if( !isInside(x, y) ){
        return;
      };
      if(state == GLUT_DOWN){
        if(selected != -1){
          data.labelSelected(selected) = !data.labelSelected(selected);
          data.labelHighlight = -1; 
        }
      }
    };



    // catch mouse move events
    void motion(int x, int y){};




    void passive(int x, int y){
      xM = x;
      yM = y;
      selected = -1;
      if( !isInside(x, y) ){
        data.labelHighlight=-1;
        return;
      };



      GLint vp[4];
      glGetIntegerv(GL_VIEWPORT, vp);
      GLuint selectBuf[BUFSIZE];
      glSelectBuffer(BUFSIZE, selectBuf);
      glRenderMode(GL_SELECT);
      glInitNames();

      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glLoadIdentity();
      gluPickMatrix(x, vp[3]-y, 3, 3, vp);
      glOrtho (0, vp[2], vp[3], 0, 0, 1);

      display();
      GLint hits = glRenderMode(GL_RENDER);
      selected = -1;
      for(int i=0; i<hits; i++){
        int tmp = selectBuf[i*4 + 3];
        if(tmp != -1){
          selected = tmp;
        }
      }
      if(selected == -1){
        data.labelHighlight = -1;
      }
      else{
        data.labelHighlight = selected;
      }

      glMatrixMode(GL_PROJECTION);
      glPopMatrix();
      glutPostRedisplay();

    };




};

#endif
