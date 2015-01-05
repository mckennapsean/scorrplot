#ifndef TEXTDEL_H
#define TEXTDEL_H

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

bool allActivated = true;

template<typename TPrecision>
class TextDEL : public DisplayElement{

  private:


#define BUFSIZE 512


    int selected;
    int xM,yM;
    int pickH, pickW;
    int mod;

    Data<TPrecision> &data;
    Font &font;
    Animator &animator;


    Rotation<TPrecision> Rone;
    Rotation<TPrecision> Rtwo;



  public:

    TextDEL(Data<TPrecision> &d, Font &f, Animator &a) 
           : data(d), font(f), animator(a), Rone(d,
Rotation<TPrecision>::Primary), Rtwo(d,Rotation<TPrecision>::Secondary){ 
       pickW=5;
       pickH=5;
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
      font.setSize(fs);
      int rIndex = data.getSelected();
      float col[4];
    
      glPushName(-1); 
      int j = 0;
        for(int i=0; i<data.selectedIndicies.size(); i++){
          glLoadName(i); 
          int index = data.selectedIndicies[i]; 
	  if(index == rIndex){
	    col[3] = 1;
	  }
	  else{ 
	    col[3] = 0.75;
	  }
          data.color(index, col);
	  glColor4fv(col);

	  // create the label for all non-red (brain) elements automatically
	  if(col[0] < 0.85 || col[1] > 0.15 || col[2] > 0.15)
      font.renderString(data.rowname(index), x, y+j*(fs+2));
      
    // allow us to choose for red / brain elements
    else
    
      // only create labels if user has not pressed 'a' yet
      if(allActivated)
        font.renderString(data.rowname(index), x, y+j*(fs+2));
      else
        j--;
    j++;
        }


     glPopName();

    };

// resize text display to the window size (manual default values set)
void reshape(int w, int h){
  int startLocX = (int) (790. * (double) w / 1100. );
  int startLocY = (int) (210. * (double) h / 710. );
  int endLocX = (int) (230. * (double) w / 1100. );
  int endLocY = (int) (770. * (double) h / 710. );
  location(startLocX, startLocY, endLocX, endLocY);
}

    void idle(){};


    void keyboard(unsigned char key, int x, int y){\
      switch(key){
        case 'a':
        case 'A':
          if(allActivated)
            allActivated = false;
          else
            allActivated = true;
          break;
          
        default:
          break;
      }
    };



    void special(int key, int x, int y){
      switch(key)
      {
        case GLUT_KEY_DOWN:
          data.selectedIndex++;
          break;
        case GLUT_KEY_UP:
          data.selectedIndex--;
          break;
      }
      glutPostRedisplay();
    };




    void mouse(int button, int state, int x, int y){
      xM = x;
      yM = y;
      if(!isInside(x, y)){ return; };

      mod = glutGetModifiers();

      if (state == GLUT_DOWN ){
        if(selected >= 0){
          data.selectedIndex = selected;

	  if(button == GLUT_LEFT_BUTTON ){
	    if(mod == GLUT_ACTIVE_CTRL){
              Rtwo.setTarget(data.getSelected());
              animator.setAnimation(&Rtwo);
	    }
	    else if(mod == GLUT_ACTIVE_SHIFT){
              data.showRowname.push_back(data.getSelected());
              data.showProfile.push_back(data.getSelected());
	      glutPostRedisplay();
            }
	    else{
              Rone.setTarget(data.getSelected());
              animator.setAnimation(&Rone);
            }
	  }
        }
      }
    };



    // catch mouse move events
    void motion(int x, int y){};




    void passive(int x, int y){
      xM = x;
      yM = y;

      if( !isInside(x, y) ){ 

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
      gluPickMatrix(x, vp[3]-y, pickW, pickH, vp);
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

      glMatrixMode(GL_PROJECTION);
      glPopMatrix();
      glutPostRedisplay();
    };




};

#endif
