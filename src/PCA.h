#ifndef PCA_H
#define PCA_H


#include "SymmetricEigensystem.h"
#include "DenseVector.h"
#include "DenseMatrix.h"
#include "Linalg.h"

template <typename TPrecision>
class PCA{


  public:
    
    //Constructs PCA to project and unproject vectors, covariance is not known.
    PCA(FortranLinalg::DenseMatrix<TPrecision> &evecs, FortranLinalg::DenseVector<TPrecision> &evals,
        FortranLinalg::DenseVector<TPrecision> &m):data(PCA<TPrecision>::dataDummy)
    { 
      using namespace FortranLinalg;
      ev = evecs;
      ew = evals;
      mean = m;

      nProjectionDimensions = ev.N();
      
      buffer = DenseVector<TPrecision>(nProjectionDimensions);

    }

    
    PCA(FortranLinalg::DenseMatrix<TPrecision> &samples, unsigned int ndims, bool subtractMean = true) : data(samples){
      
      using namespace FortranLinalg;
      rowCov = data.M() > data.N();

      if(subtractMean){ 
        computeMean(data); 
        Linalg<TPrecision>::SubtractColumnwise(data, mean, data);
      }

      //build covariance
      computeCovariance(samples);
     
      if(ndims == 0 || ndims > covariance.N()){
        nProjectionDimensions = covariance.N();
      } 
      else{
        nProjectionDimensions = ndims;
      }
      buffer = DenseVector<TPrecision>(nProjectionDimensions);

      computePC();
    };



    FortranLinalg::DenseMatrix<TPrecision> &getCovariance(){
      return covariance;
    };

    FortranLinalg::DenseMatrix<TPrecision> &getEigenVectors(){
      return ev;
    };

    FortranLinalg::DenseVector<TPrecision> &getEigenValues(){
      return ew;
    };

    FortranLinalg::DenseVector<TPrecision> &getMean(){
      return mean;
    };
  




    //Projection methods 


    //Project column from a data matrix
    void projectColumn(FortranLinalg::DenseMatrix<TPrecision> &matrix, int index,
        FortranLinalg::DenseVector<TPrecision> &projected, bool subtractMean = false){
      using namespace FortranLinalg;
      
      if(subtractMean){
        Linalg<TPrecision>::SubtractColumn(matrix, index, mean, matrix);
      }

      Linalg<TPrecision>::MultiplyColumn(ev, matrix, index, projected, true);

    }



    FortranLinalg::DenseVector<TPrecision> projectColumn(FortranLinalg::DenseMatrix<TPrecision> &matrix, int
        index, bool subtractMean = false){
      using namespace FortranLinalg;

      DenseVector<TPrecision> projected(ev.N());
      projectColumn(matrix, index, projected, subtractMean);
      return projected;
    };


    //Project single vector
    void project(FortranLinalg::DenseVector<TPrecision> &v, FortranLinalg::DenseVector<TPrecision> &projected,
        bool subtractMean = false){
      using namespace FortranLinalg;

      if(subtractMean){
        Linalg<TPrecision>::Subtract(v, mean, v);
      }

      Linalg<TPrecision>::Multiply(ev, v, projected, true);
    };


    FortranLinalg::DenseVector<TPrecision> project(FortranLinalg::DenseVector<TPrecision> &v, bool
        subtractMean = false){ 
      using namespace FortranLinalg;

      DenseVector<TPrecision> projected(ev.N());
      project(v, projected, subtractMean);
      return projected;  

    };
   

    //Project complete data matrix
    FortranLinalg::DenseMatrix<TPrecision> project(FortranLinalg::DenseMatrix<TPrecision> &matrix, bool
        subtractMean = false){ 
      using namespace FortranLinalg;

      if(subtractMean){
        Linalg<TPrecision>::SubtractColumnwise(matrix, mean, matrix);
      }
      DenseMatrix<TPrecision> projected 
        = Linalg<TPrecision>::Multiply(ev, matrix, true, false);
      
      return projected;
    };




    //Unproject from single vector

    //unproject v -  size of v dimensions are used
    //to do the unprojection
    FortranLinalg::DenseVector<TPrecision> unproject(FortranLinalg::Vector<TPrecision> &v, bool addMean = true){
      return unproject(v, v.N(), addMean);  
    };
    //same with preallocated out of size nRows
    void unproject(FortranLinalg::Vector<TPrecision> &v, FortranLinalg::DenseVector<TPrecision> &out, 
                   bool addMean = true){
      unproject(v, v.N(), out, addMean);
    };

    //Unproject v, the first leading nDimensions principal components 
    //are used for the unprojection
    FortranLinalg::DenseVector<TPrecision> unproject(FortranLinalg::Vector<TPrecision> &v, int nDimensions, 
                                      bool addMean = true){
      using namespace FortranLinalg;

      DenseVector<TPrecision> tmp(ev.M());
      unproject(v, nDimensions, tmp, addMean);
      return tmp;    
    };

    //same with preallocated out of size nRows
    void unproject(FortranLinalg::Vector<TPrecision> &v, unsigned int nDimensions, FortranLinalg::DenseVector<TPrecision> &out, 
                   bool addMean = true){
      using namespace FortranLinalg;

      if(nDimensions > v.N()){
        nDimensions = v.N();
      }
      if(nDimensions > ev.N()){
        nDimensions = ev.N();
      }


      zeroBuffer(nDimensions);
      for(unsigned int j = 0; j < nDimensions; j++){
        buffer(buffer.N()-j-1) = v(v.N()-j-1);
      }

      Linalg<TPrecision>::Multiply(ev, buffer, out);
        
      if(addMean){
        Linalg<TPrecision>::Add(out, mean, out);
      }

    };



    //Unproject column from a matrix
    //Unproject column index from matrix data, data is assumed to be of siez
    // nEv x >index. 
    //The first leading nDimension principal components are used for projection
    // the projection is stored in out
    void unproject(FortranLinalg::Matrix<TPrecision> &matrix, int index, int
        nDimensions, FortranLinalg::DenseVector<TPrecision> &out, bool addMean = true){
      using namespace FortranLinalg;

      if(nDimensions > matrix.M()){
        nDimensions = matrix.M();
      }
      if(nDimensions > ev.N()){
        nDimensions = ev.N();
      }

        zeroBuffer(nDimensions);
        for(int j = 0; j < nDimensions; j++){
           buffer(buffer.N()-j-1) = matrix(matrix.M()-j-1, index);
        }

        Linalg<TPrecision>::Multiply(ev, buffer, out);
        
      if(addMean){
        Linalg<TPrecision>::Add(out, mean, out);
      }

    };

    
    
    FortranLinalg::DenseMatrix<TPrecision> unproject(FortranLinalg::DenseMatrix<TPrecision> &matrix, bool addMean = true){
      using namespace FortranLinalg;


      DenseMatrix<TPrecision> out; 
      if(matrix.M() < ev.N()){
        DenseMatrix<TPrecision> tmp(ev.M(), matrix.M());
        for(unsigned int i=0; i<matrix.M(); i++){
          Linalg<TPrecision>::SetColumn(tmp, tmp.N()-1-i, ev, ev.N()-1-i);
        }
        out =  Linalg<TPrecision>::Multiply(tmp, matrix);
        tmp.deallocate();

      }
      else{
        out =  Linalg<TPrecision>::Multiply(ev, matrix);
      }
      if(addMean){
        Linalg<TPrecision>::AddColumnwise(out, mean, out);
      }

      return out;
    };
    
    void cleanup(){
      covariance.deallocate();
      mean.deallocate();
      ev.deallocate();
      ew.deallocate();
      buffer.deallocate();
    };   
  



    FortranLinalg::DenseMatrix<TPrecision> data;
    FortranLinalg::DenseMatrix<TPrecision> covariance;
    FortranLinalg::DenseVector<TPrecision> mean;
    FortranLinalg::DenseMatrix<TPrecision> ev;
    FortranLinalg::DenseVector<TPrecision> ew;




  private:
    static FortranLinalg::DenseMatrix<TPrecision> dataDummy;
    FortranLinalg::DenseVector<TPrecision> buffer;
    bool rowCov;  
    int nProjectionDimensions;


    void computePC(){
      using namespace FortranLinalg;
      if(ev.N() > 0){
        return;
      }

      SymmetricEigensystem<TPrecision> eigs(covariance,
          covariance.N()-nProjectionDimensions+1, covariance.N());

      //row covariance pca
      if(rowCov){
        ev = Linalg<TPrecision>::Multiply(data, eigs.ev); 
        ew = DenseVector<TPrecision>(eigs.ew.N());

        //Convert to col eigs
        for(unsigned int i=0; i<ew.N(); i++){
          ew(i) = eigs.ew(i) / (TPrecision)(data.N()-1.0); 
        }

        //normalize
        for(unsigned int i=0; i<ev.N(); i++){
          TPrecision tmp = 0;
          for(unsigned int k=0; k< data.M(); k++){
            tmp += ev(k, i) * ev(k, i);
          }
          tmp = (TPrecision) 1.0/sqrt(tmp);
          for(unsigned int k=0; k< data.M(); k++){
            ev(k, i) *= tmp;
          }
        }
        eigs.ew.deallocate();
        eigs.ev.deallocate();
      }
      //Normal PCA
      else{
        ev = eigs.ev;
        ew = eigs.ew;
      }
      for(unsigned int i=ew.N(); i < ev.N(); i++){
        for(unsigned int k=0; k< ev.M(); k++){
          ev(k, i) =0;
        }
      }
    };


    void computeCovariance(FortranLinalg::Matrix<TPrecision> &data){
      using namespace FortranLinalg;
      if(rowCov){
        covariance = DenseMatrix<TPrecision>(data.N(), data.N());
        for(int i=0; i<data.N(); i++){
          for(int j=i; j<data.N(); j++){
            TPrecision tmp = 0;
            for(int k=0; k<data.M(); k++){
              tmp += data(k, i) * data(k, j);
            }
            covariance(i, j) = tmp;
            covariance(j, i) = tmp;
          }
        }
      }
      else{ 
        covariance = DenseMatrix<TPrecision>(data.M(), data.M());
        for(int i=0; i<data.M(); i++){
          for(int j=i; j<data.M(); j++){
            TPrecision tmp = 0;
            for(int k=0; k<data.N(); k++){
              tmp += data(i, k) * data(j, k);
            }
            covariance(i, j) = tmp;
            covariance(j, i) = tmp;
          }
        }
        TPrecision *dptr = covariance.data();
        for(int i=0; i<data.M()*data.M(); i++){
          dptr[i] /= (TPrecision)(data.N()-1);
        }
      }
    };



    void computeCovariance(FortranLinalg::DenseMatrix<TPrecision> &data){
      using namespace FortranLinalg;
      if(rowCov){
        covariance = Linalg<TPrecision>::Multiply(data, data, true, false);
      }
      else{ 
        covariance = Linalg<TPrecision>::Multiply(data, data, false, true);
        TPrecision *dptr = covariance.data();
        for(unsigned int i=0; i<data.M()*data.M(); i++){
          dptr[i] /= (TPrecision)(data.N()-1);
        }
      }
    };




    void computeMean(FortranLinalg::Matrix<TPrecision> &data){
      using namespace FortranLinalg;
        //Compute mean
        mean = DenseVector<TPrecision>(data.M());
        for(unsigned int i=0; i < data.M(); i++){
          mean(i) = 0;
          for(unsigned int j=0; j < data.N(); j++){
            mean(i) += data(i, j);
          }
        }
        for(unsigned int i=0; i < data.M(); i++){
          mean(i) /= ((TPrecision)(data.N()));
        }

    };

    void zeroBuffer(int nDimensions){
      TPrecision *ptr = buffer.data();
      for(unsigned int i=0; i<buffer.N()-nDimensions; i++){
        ptr[i] = 0;
      }
    };

};

template <typename TPrecision>
FortranLinalg::DenseMatrix<TPrecision> PCA<TPrecision>::dataDummy;

#endif
