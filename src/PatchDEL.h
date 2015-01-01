#ifndef PATCHDEL_H
#define PATCHDEL_H

#include "Display.h"
#include "PCA.h"
#include "DenseVector.h"
#include "DenseMatrix.h"
#include "Animator.h"


template<typename TPrecision>
class PatchDEL : public DisplayElement{

  private:


#define BUFSIZE 512



    int xM, yM;
    int mod, cur_button; 
    int pickH, pickW;  

   
    Data<TPrecision> &data;


  public:

    PatchDEL(Data<TPrecision> &d) : data(d){ 
      pickW = 2;
      pickH = 2;
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
      using namespace FortranLinalg;
      int sel = data.getSelected();
      if(sel >= 0){
      
       int r = sqrt( data.raw.M() );
       int index = 0;
       static int vs = width/(r+1.0);
       glBegin(GL_QUADS);
       TPrecision ma = Linalg<TPrecision>::MaxColumn(data.centered, sel); 
       TPrecision mi = Linalg<TPrecision>::MinColumn(data.centered, sel); 
       for(int i = 0; i < r; i++){
         for(int j = 0; j < r; j++){
           double col = (data.centered(index, sel)-mi)/(ma-mi); 
           glColor3f( col, col, col );
           glVertex2f( xLeft + i*vs, yTop + j*vs);
           glVertex2f( xLeft + (i+1)*vs, yTop + j*vs);
           glVertex2f( xLeft + (i+1)*vs, yTop + (j+1)*vs);
           glVertex2f( xLeft + i*vs, yTop + (j+1)*vs);
           ++index;
         }
       }  
       glEnd();

       }

    };

  
    void idle(){

    };


    void keyboard(unsigned char key, int x, int y){
    };



    void special(int key, int x, int y){
    };




    void mouse(int button, int state, int x, int y){

    };



    // catch mouse move events
    void motion(int x, int y){

    };




    void passive(int x, int y){

    };



};

#endif
