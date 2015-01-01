// -*- C++ -*-
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//                              Color Mapper                                //
//                                                                          //
//                       Written by Kristi Potter                           //
//                           February 5, 2010                               //
//                                                                          //
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

#ifndef _COLORMAPPER_H
#define _COLORMAPPER_H

#include <vector>
#include <iostream>



template<class T>
class ColorMapper
{

 public: 

  ColorMapper<T>(T rangeMin = T(0), T rangeMax = T(1))
    : _rangeMin(rangeMin), _rangeMax(rangeMax)
  {
    _rangeMid = (_rangeMax - _rangeMin)/2.f + _rangeMin;
  }

/*
  ColorMapper<T>(const ColorMapper<T> &cm)
    : _rangeMin(cm._rangeMin), _rangeMax(cm._rangeMax)
  {}
  */


  ColorMapper<T> &operator=(const ColorMapper<T> &cm)
  {
    _rangeMid = cm._rangeMid;
    _rangeMin = cm._rangeMin;
    _rangeMax = cm._rangeMax;
    
    //default colormap
    //blue to green to red
    rmin = 0;
    rmid = 0;
    rmax = 1;
    gmin = 0;
    gmid = 1;
    gmax = 0;
    bmin = 1;
    bmid = 0;
    bmax = 0;

    return *this;
  }




  virtual ~ColorMapper() {}

      
  
  void set(T Rmin, T Rmid, T Rmax, T Gmin, T Gmid, T Gmax, T Bmin, T Bmid, T Bmax){

    rmin = Rmin;
    rmid = Rmid;
    rmax = Rmax;
    gmin = Gmin;
    gmid = Gmid;
    gmax = Gmax;
    bmin = Bmin;
    bmid = Bmid;
    bmax = Bmax;
  };

  void setRange(T rmin, T rmid, T rmax){
    _rangeMin = rmin;
    _rangeMid = rmid;
    _rangeMax = rmax;
  };

  void setRange(T rmin, T rmax){
    _rangeMin = rmin;
    _rangeMid = rmin + (rmax-rmin)/2.f;
    _rangeMax = rmax;
  };

  // Return a color for a scalar
  std::vector<T> getColor(T value)
  {
   
    T r,g,b;
    if(value < _rangeMid){
     r = affine(value, _rangeMin, _rangeMid, rmin, rmid);
     g = affine(value, _rangeMin, _rangeMid, gmin, gmid);
     b = affine(value, _rangeMin, _rangeMid, bmin, bmid);
    }
    else{     
      r = affine(value, _rangeMid, _rangeMax, rmid, rmax);
      g = affine(value, _rangeMid, _rangeMax, gmid, gmax);
      b = affine(value, _rangeMid, _rangeMax, bmid, bmax);
    }
    std::vector<T> color;
    color.push_back(r);
    color.push_back(g);
    color.push_back(b);
    return color;
  };  // Return a color for a scalar
  

  void getColor(T value, T color[3]){
    if(value < _rangeMid){
     color[0] = affine(value, _rangeMin, _rangeMid, rmin, rmid);
     color[1] = affine(value, _rangeMin, _rangeMid, gmin, gmid);
     color[2] = affine(value, _rangeMin, _rangeMid, bmin, bmid);
    }
    else{     
      color[0] = affine(value, _rangeMid, _rangeMax, rmid, rmax);
      color[1] = affine(value, _rangeMid, _rangeMax, gmid, gmax);
      color[2] = affine(value, _rangeMid, _rangeMax, bmid, bmax);
    }
  };
    
  private:

  T _rangeMin, _rangeMax, _rangeMid;
  
  inline T affine(const T x, const T a, const T A, const T b, const T B)
  {
    return T((B - b)*(x - a)/float((A-a)) + b); 
  };

  T rmin, rmid, rmax;
  T gmin, gmid, gmax;
  T bmin, bmid, bmax;

  

};

#endif
