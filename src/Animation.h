#ifndef ANIMATION_H
#define ANIMATION_H


class Animation{

public:

  virtual ~Animation(){};

  virtual void step() = 0;

  virtual bool isRunning() = 0;
 

};

#endif
