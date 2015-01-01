#ifndef SEARCHDEL_H
#define SEARCHDEL_H

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
class SearchDEL : public DisplayElement{

  private:


#define BUFSIZE 512



    Data<TPrecision> &data;
    Font &font;
    Animator &animator;
    std::string name;



  public:

    SearchDEL(Data<TPrecision> &d, Animator &a, Font &f) 
           : data(d), animator(a), font(f), name(""){ 
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
    
      int fs = 10;
      int x = xLeft;
      int y = yTop;
      font.setSize(fs);
     
      glColor4f(0.6, 0.6, 0.6, 1);
      glLineWidth(1);
      glBegin(GL_LINE_LOOP);
      glVertex2f(x, y+16);
      glVertex2f(x+width, y+16);
      glVertex2f(x+width, y+32);
      glVertex2f(x, y+32);
      glEnd();
      
      
      font.renderString("Search",x, y+14); 
      

      glColor4f(0.7, 0.7, 0.7, 1);
      font.renderString(name, x+2, y+29); 


    };


    

    void idle(){

    };


    void keyboard(unsigned char key, int x, int y){
      if(isInside(x, y)){
        if(key == 13){
          for(int i=0; i<data.rownames.size(); i++){
             if(name.compare(data.rownames[i]) == 0 ){
               data.searched.push_back(i);
               break;
             }
          }
          //name.clear();
        }
        else if(key == 127){
          name.erase(name.size()-1);
        }
        else{
          std::string tmp(1, key);
          name.append(tmp);
        }
      }
      glutPostRedisplay();
    };



    void special(int key, int x, int y){
    };




    void mouse(int button, int state, int x, int y){
    };



    // catch mouse move events
    void motion(int x, int y){};




    void passive(int x, int y){

    };




};

#endif
