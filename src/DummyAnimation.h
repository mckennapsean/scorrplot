#ifndef DUMMYANIMATION_H
#define DUMMYANIMATION_H

class DummyAnimation : public Animation{
  public:

   bool isRunning(){
    return false;
   };


   void step(){};


};

#endif
