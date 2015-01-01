#ifndef PCADEL_H
#define PCADEL_H

#include "Display.h"
#include "PCA.h"
#include "DenseVector.h"
#include "DenseMatrix.h"
#include "Rotation.h"
#include "Animator.h"

#include <algorithm>

template<typename TPrecision>
class PCADEL : public DisplayElement{

  private:


#define BUFSIZE 512


    int nTop;
    int xM, yM;
    int mod, cur_button; 
    int pickH, pickW;  
    float r,g,b;
    int lw;
    int selected;
    int lIndex;

    FortranLinalg::DenseVector<TPrecision> v;
    PCA<TPrecision> *pca;
   
    Data<TPrecision> &data;
    Rotation<TPrecision> Rone;
    Rotation<TPrecision> Rtwo;
    Animator &animator;


  public:

    PCADEL(int l, Data<TPrecision> &d, Animator &a, float r1, float g1, float
b1)   :  lIndex(l), data(d), Rone(d, Rotation<TPrecision>::Primary), Rtwo(d,
Rotation<TPrecision>::Secondary), animator(a){ 
      pickW = 2;
      pickH = 2;
      r=r1;
      g=g1;
      b=b1;
      lw = 4;
      setPCA();
      selected = -1;
      nTop=10;
    };


    void setPCA(){
      using namespace FortranLinalg;
      bool changed = false;
      if(lIndex == -1){
        changed = pca != data.pca;
        pca = data.pca;
      }
      else{
       changed = pca != data.lPCA[lIndex];
       pca = data.lPCA[lIndex];
      }
      if(changed){
        v.deallocate();
        v = Linalg<TPrecision>::Copy(pca->ew);
        TPrecision ma = Linalg<TPrecision>::Max(v);
        TPrecision mi = Linalg<TPrecision>::Min(v);
        Linalg<TPrecision>::Subtract(v, mi, v);
        Linalg<TPrecision>::Scale(v, 1/(ma-mi), v);
      }
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
      setPCA();
      bool active = isInside(xM, yM);

     
      glMatrixMode(GL_MODELVIEW); 	
      glLoadIdentity();
     
      
      glPushName(-1); 
      glLineWidth(lw);
      for(int i=0; i < std::min(nTop, (int)v.N()); i++){
        if(i==selected && active){
          glColor4f(r, g, b, 1);
        }
        else{
          glColor4f(r, g, b, 0.75);
        }
        glLoadName(i);
        glBegin(GL_LINES);
        glVertex2f(xLeft+(nTop-i)*(lw+1), yTop + height * (1.0 - v(v.N()-1-i)) );
        glVertex2f(xLeft+(nTop-i)*(lw+1), yTop + height );
        glEnd();
      }




      glLoadName(-1); 


if(active){
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

// due to multiple graphs, trickier to resize
// left blank for now, could be scaled later (not a huge issue)
void reshape(int w, int h){

}
  
    void idle(){

    };


    void keyboard(unsigned char key, int x, int y){
    };



    void special(int key, int x, int y){
    };




    void mouse(int button, int state, int x, int y){
      xM = x;
      yM = y;
      if(!isInside(x, y)){ return; };

      mod = glutGetModifiers();
      if (state == GLUT_DOWN && selected >= 0){
	if(button == GLUT_LEFT_BUTTON ){
	  if(mod == GLUT_ACTIVE_CTRL){
              /*//subtract current align vector and do PCA
              DenseMatrix<TPrecision> tmp = Linalg<TPrecision>::Copy(pca->data);
              for(int i=0; i<tmp.N(); i++){
                 double dot = Linalg<TPrecision>::DotColumnColumn(tmp, i, data.V, 0);
                 Linalg<TPrecision>::ColumnAddScale(tmp, i, -dot, data.V, 0);
              }
              PCA<TPrecision> pcat(tmp, tmp.M()-1, false);
              Rtwo.setTarget(pcat.ev, selected);
              pcat.cleanup();
              */
             Rtwo.setTarget(pca->ev, v.N()-1-selected);
              animator.setAnimation(&Rtwo);
	  }
	  else if(mod == GLUT_ACTIVE_SHIFT){

	  }
	  else{
              Rone.setTarget(pca->ev, v.N()-1-selected);
              animator.setAnimation(&Rone);
	  }
	}
        selected=-1;
      }
    };



    // catch mouse move events
    void motion(int x, int y){

    };




    void passive(int x, int y){
      xM = x;
      yM = y;
      
      if(!isInside(x, y)){ return; };

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
