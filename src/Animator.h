#ifndef ANIMATIOR_H
#define ANIMATIOR_H

#include "Animation.h"
#include "DummyAnimation.h"

class Animator{
private:
  Animation *cur;

public:

  Animator(){
    cur = new DummyAnimation();
  };

  void setAnimation(Animation *a){
    cur = a;  
    if(cur->isRunning()){
   //   glutIdleFunc(idle1);
    }
  };


  bool isRunning(){
   return cur->isRunning();
  };
  

  void step(){
    cur->step();
    if(!cur->isRunning()){
    //  glutIdleFunc(NULL);
    }
  }; 

};

#endif
