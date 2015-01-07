#ifndef PROJECTIONDEL_H
#define PROJECTIONDEL_H

#include "Display.h"
#include "DenseVector.h"
#include "DenseMatrix.h"
#include "Rotation.h"
#include "Animator.h"

#include <math.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>

#include "Font.h"


template<typename TPrecision>
class ProjectionDEL : public DisplayElement{

  private:


#define BUFSIZE 512


    int size; 

    int xM, yM;
    int mod, cur_button; 
    int pickH, pickW;  

    std::vector<int> select;   
    std::vector<int> selectPrevious;   

    Data<TPrecision> &data;

    Rotation<TPrecision> Rone;
    Rotation<TPrecision> Rtwo;
    Animator &animator;


    float alpha;
    float pointSize;

    Font &font;

    //displayLists
    int circleDL;
    int linesDL;


    SquaredEuclideanMetric<TPrecision> sL2;
  public:

    ProjectionDEL(Data<TPrecision> &d, Animator &a, Font &f, float &al) 
      : data(d), Rone(d, Rotation<TPrecision>::Primary),
      Rtwo(d,Rotation<TPrecision>::Secondary), animator(a), font(f){ 
        pickW = 5;
        pickH = 5;
        alpha = al;
        pointSize = 6;
      };



    void location(int xPos, int yPos, int w, int h){
      width = w;
      height = h;
      xLeft = xPos;
      yTop = yPos;
      size = std::min(w, h)/2 - 5;


      circleDL = glGenLists(2);
      glNewList(circleDL, GL_COMPILE); 
      drawCircle(0, 0, 1);
      glEndList();

      linesDL = circleDL+1;
      glNewList(linesDL, GL_COMPILE); 
      drawLines();
      glEndList();

    };



    void init(){  

    };


    void display(void){


      glMatrixMode(GL_MODELVIEW); 	
      glLoadIdentity();
      
      // draw background plot first
      glCallList(circleDL);
      //glCallList(linesDL);
      drawLines();

      //glPushName(-1); 
      glPointSize(pointSize);
      int x1;
      int x2;
      float col[4];
      col[3] = alpha;
      glBegin(GL_POINTS);
      for(int i=0; i<data.P.N(); i++){
        int lid = data.labels(i) ;
        if(data.labelSelected(lid) || lid == data.labelHighlight){
          continue;
        }
        x1 = toScreenX(data.P(0,i));
        x2 = toScreenY(data.P(1,i));


        //glLoadName(i);
        data.color(i, col);
        glColor4fv(col);
        glVertex2f(x1, x2);

      }
      glEnd();

      col[3] = 1;
      for(int i=0; i<data.P.N(); i++){
        int lid = data.labels(i) ;
        if(!data.labelSelected(lid) && lid != data.labelHighlight){
          continue;
        }
        x1 = toScreenX(data.P(0,i));
        x2 = toScreenY(data.P(1,i));


        glLoadName(i);
        data.color(i, col);
        glColor4fv(col); 
        glBegin(GL_POINTS);


        glVertex2f(x1, x2);

        glEnd();

      }


      glLoadName(-1); 
      col[3] = 1;
      font.setSize(11);
      for(int j=0; j < data.showRowname.size(); j++){
        int i = data.showRowname[j];
        x1 = toScreenX(data.P(0,i));
        x2 = toScreenY(data.P(1,i));


        //	glLoadName(i);
        data.color(i, col);
        glColor4fv(col); 
        glBegin(GL_POINTS);
        glVertex2f(x1, x2);
        glEnd();


        font.renderString(data.rownames[i], x1, x2);
      }


      /*
         glColor4f(1, 1, 1, 0.5);  
         glPointSize(9);
         glBegin(GL_POINTS);
         for(int i=0; i<data.PU.N(); i++){
         x1 = toScreenX(data.PU(0,i));
         x2 = toScreenY(data.PU(1,i));
         glVertex2f(x1, x2);
         }
         glEnd();
       */


      if(isInside(xM, yM)){
        //draw pick rectangle
        glColor4f(0.75, 0.75, 0.75, 0.75);
        glLineWidth(2);
        double pw = pickW/2.0;
        double ph = pickH/2.0;
        glBegin(GL_LINE_LOOP);
        glVertex2f(xM-pw, yM-ph);
        glVertex2f(xM+pw, yM-ph);
        glVertex2f(xM+pw, yM+ph);
        glVertex2f(xM-pw, yM+ph);
        glEnd();
      }

      glPopName();

    };

    // resize projection data to the window size (manual default values set)
    void reshape(int w, int h){
      int startLocX = (int) (151. * (double) w / 1100. );
      int startLocY = (int) (10. * (double) h / 710. );
      int endLocX = (int) (630. * (double) w / 1100. );
      int endLocY = (int) (630. * (double) h / 710. );
      location(startLocX, startLocY, endLocX, endLocY);
    }    

    void idle(){

    };


    void keyboard(unsigned char key, int x, int y){
      //std::vector<int> keep;
      if(!isInside(x, y)){ return; };
      switch(key)
      {

        case 'c':
        case 'C':
          data.showRowname.clear();  
          data.showProfile.clear();  
          break;
        case 's':
        case 'S':
          data.setSpearman(true);
          break;
        case 'p':
        case 'P':
          data.setSpearman(false);
          break;       
        case 'm':
        case 'M':
          data.colormapping =  !data.colormapping;
          break;
        case '+':
          pickW+=1;
          pickH+=1;
          break;
        case '-':
          pickW-=1;
          pickH-=1;
          if(pickW < 5){
            pickW=5;
          }
          if(pickH < 5){
            pickH=5;
          }
          break;

      }
      glutPostRedisplay();
    };



    void special(int key, int x, int y){
      switch(key)
      {
        case GLUT_KEY_LEFT:
          alpha *= 0.95;
          break;
        case GLUT_KEY_RIGHT:
          alpha *= 1.05;
          break;
        case GLUT_KEY_DOWN:
          pointSize *= 0.95;
          break;
        case GLUT_KEY_UP:
          pointSize *= 1.05;
          break;
      }
      if(alpha < 0){
        alpha=0;
      }
      else if(alpha > 1){
        alpha=1;
      }

      glutPostRedisplay();
    };




    void mouse(int button, int state, int x, int y){
      xM = x;
      yM = y;
      if(!isInside(x, y)){ return; };

      mod = glutGetModifiers();
      int index = data.getSelected();
      if (state == GLUT_DOWN && index >= 0){
        if(button == GLUT_LEFT_BUTTON ){
          if(mod == GLUT_ACTIVE_SHIFT){
            data.showRowname.push_back(index);
            data.showProfile.push_back(index);
            glutPostRedisplay();
          }
          else{
            Rone.setTarget(index);
            animator.setAnimation(&Rone);
            selectPrevious = select;
            data.selectedIndex = 0;
          }
        }else if(button == GLUT_RIGHT_BUTTON){
          selectPrevious = select;
          data.selectedIndex = 0;
          Rtwo.setTarget(index);
          animator.setAnimation(&Rtwo);
        }
      }
    };





    // catch mouse move events
    void motion(int x, int y){};




    void passive(int x, int y){
      xM = x;
      yM = y;

      if(!isInside(x, y)){ 
        if(selectPrevious.size() > 0){
          data.selectedIndicies = selectPrevious;
          selectPrevious.clear();
        }
        return; 
      };

      if(selectPrevious.size() == 0 && data.selectedIndicies.size() > 0){
        selectPrevious = data.selectedIndicies;
      }

      double xd = toDataX(x);
      double yd = toDataY(y);
      double wd = pickW/(double)size;
      double hd = pickH/(double)size;
      select.clear();
      for(int i=0; i< data.P.N(); i++){
        if(fabs(data.P(0, i) - xd) < wd && fabs(data.P(1, i) - yd ) < hd){
          select.push_back(i);
        }
      }

      if(select.size() > 0 ){
        data.clearSelection();
        data.addSelected(select);
      }

      glutPostRedisplay();


      /*
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
         select.clear();
         for(int i=0; i<hits; i++){
         GLuint selected = selectBuf[i*4 + 3];
         if(selected != -1){
         select.push_back(selected);
         }
         }
         if(select.size() > 0 ){
         data.clearSelection();
         data.addSelected(select);
         }

         glMatrixMode(GL_PROJECTION);
         glPopMatrix();
         glutPostRedisplay();
       */
    };


  private:

    double toScreenX(TPrecision x){
      return xLeft + (x+1.0) * size; 
    };
    double toScreenY(TPrecision y){
      return yTop + (y+1.0) * size; 
    };

    double toDataX(int x){
      return (x - xLeft ) / (double) size -1.0;
    };

    double toDataY(int y){
      return (y - yTop) / (double) size - 1.0;
    }; 






    void drawCircle(double x, double y, double r, int n=200){
      glLineWidth(3);
        
      glColor4f(0.2, 0.2, 0.2, 1);
      //glColor4f(0.9, 0.9, 0.9, 1);

      glBegin(GL_LINE_LOOP);
      for(int i =0; i <= n; i++){
        double angle = 2 * M_PI * i / n;
        double x = cos(angle);
        double y = sin(angle);
        glVertex2d(toScreenX(x), toScreenY(y));
      }
      glEnd();
    };



    void drawLines(int n=10){
      glLineWidth(1);
      //glColor4f(0.9, 0.9, 0.9, 0.2);
      glColor4f(0.2, 0.2, 0.2, 1);
      glBegin(GL_LINES);

      glVertex2d(toScreenX(-1) , toScreenY(0) );
      glVertex2d(toScreenX(1) , toScreenY(0) );


      for(int i=0; i < n; i++){
        double angle = acos( 1.0*i / n  );
        double x = cos(angle);
        double y = sin(angle);

        glVertex2d(toScreenX(x) , toScreenY(y) );
        glVertex2d(toScreenX(x) , toScreenY(-y) );

        glVertex2d(toScreenX(-x) , toScreenY(y) );
        glVertex2d(toScreenX(-x) , toScreenY(-y) );

      }
      glEnd();

      if(data.nPerms > 0){ 
        font.setSize(8);
        for(int i=0; i < n; i++){
          double angle = acos( 1.0*i / n );
          double x = cos(angle);
          double y = sin(angle);

          std::stringstream s1;
          std::stringstream s2;
          if(data.pValsP(i) < 0){
            s1 << "-";
            s2 << "-";
          }
          else{
            s1 << std::setiosflags(std::ios::fixed) << std::setprecision(2);
            s1 << data.pValsP(i);
            s2 << std::setiosflags(std::ios::fixed) << std::setprecision(2);
            s2 << data.pValsN(i);
          }

          font.renderString(s1.str(), toScreenX(x) + 3 , toScreenY(y) + 10 );
          font.renderString(s2.str(), toScreenX(-x) - 35 , toScreenY(y) + 10 );
        }
      }

    };



};

#endif
