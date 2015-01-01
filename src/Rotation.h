#ifndef ROTATION_H
#define ROTATION_H

#include "DenseVector.h"
#include "DenseMatrix.h"

#include "Animation.h"



template<typename TPrecision>
class Rotation : public Animation{

  public:
    enum RotationType{ Primary, Secondary};
  private:

    Data<TPrecision> &data;


    FortranLinalg::DenseVector<TPrecision> target;
    FortranLinalg::DenseVector<TPrecision> secant;
    FortranLinalg::DenseVector<TPrecision> v;

    int nSteps;
    int cStep;
    RotationType type;
    

    
  public:
    
    Rotation(Data<TPrecision> &d, RotationType t = Primary) : data(d), type(t){ 
      using namespace FortranLinalg;

      target = DenseVector<TPrecision>(data.raw.M());
      secant = DenseVector<TPrecision>(data.raw.M());
      v = DenseVector<TPrecision>(data.raw.M());

      nSteps = 25;
      cStep = nSteps+1;

    };



   void setTarget(int index){
     setTarget(data.centered, index);
   };

   void setTarget(FortranLinalg::DenseMatrix<TPrecision> M, int index){
      using namespace FortranLinalg;
     Linalg<TPrecision>::ExtractColumn(M, index, target);
     setTarget(target);     
   };

   void setTarget(FortranLinalg::Vector<TPrecision> &t){
      using namespace FortranLinalg;
     Linalg<TPrecision>::Copy(t, target);
     Linalg<TPrecision>::Normalize(target);
     TPrecision l = Linalg<TPrecision>::Sum(target);
     Linalg<TPrecision>::Subtract(target, l/target.N(), target);
     Linalg<TPrecision>::ExtractColumn(data.V, type, v);
     Linalg<TPrecision>::Subtract(target, v, secant);
     cStep = 0;
     if(type == Secondary){
       Linalg<TPrecision>::Copy(target, data.secondary);
     }
   };

   bool isRunning(){
    return cStep <= nSteps;
   };


   void step(){
      using namespace FortranLinalg;
     if( isRunning() ){
         double ratio = cStep/(double)nSteps;
         Linalg<TPrecision>::AddScale(target, ratio-1.0, secant, v);

         if(type==Primary){
           Linalg<TPrecision>::Normalize(v);
           data.setProjection(v, 0);
           TPrecision d = Linalg<TPrecision>::Dot(v, data.secondary);
           Linalg<TPrecision>::AddScale(data.secondary, -d, v, v);
           Linalg<TPrecision>::Normalize(v);
           data.setProjection(v, 1);
         }
         else{
           TPrecision d = Linalg<TPrecision>::DotColumn(data.V, 0, v);
           Linalg<TPrecision>::AddScale(v, -d, data.V, 0, v);
           Linalg<TPrecision>::Normalize(v);
           data.setProjection(v, 1);
         }
         data.updateProjection();
         cStep++;
         glutPostRedisplay();
      }
   };



};

#endif
