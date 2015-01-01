#ifndef SYMMETRICEIGENSYSTEM_H
#define SYMMETRICEIGENSYSTEM_H


#include "DenseMatrix.h"
#include "SymmetricDenseMatrix.h"

#include "Vector.h"

#include "LapackDefs.h"

namespace FortranLinalg{



//Symmetric Eigensystem solver using lapack
template <typename TPrecision>
class SymmetricEigensystem{
  
  public:
    
    SymmetricEigensystem(DenseMatrix<TPrecision> &matrix, int il, int ih, bool evalsOnly = false){
      compute(evalsOnly, 'I', il, ih, 0, 0, ih-il+1, matrix);
    };

    SymmetricEigensystem(DenseMatrix<TPrecision> &matrix, bool evalsOnly = false){
      compute(evalsOnly, 'A', 0, 0, 0, 0, matrix.N(), matrix);
    };

    //Containes the ev in columnmajor order 
    //order corresponding to the eigenvectors 
    DenseMatrix<TPrecision> ev;
    //Ascending ordered eigenvaluse
    DenseVector<TPrecision> ew;

    //Return info from lapack routine
    int info;
    int nEigenvectors;
  
    void cleanup(){ 
      ev.deallocate();
      ew.deallocate();
    };



  private:
    
    void compute(bool evalsOnly, char range, int il, int ih, TPrecision vl,
        TPrecision vh, int nEvs, DenseMatrix<TPrecision> &matrix){

      char jobz = 'V';
      if(evalsOnly){
        jobz = 'N';
      }   

      char uplo = 'L';
      TPrecision tol = 0;
      nEigenvectors = 0;

      DenseVector<TPrecision> ewTmp = DenseVector<TPrecision>(matrix.N());
      ev = DenseMatrix<TPrecision>(matrix.N(), nEvs);
      int *isuppz = new int[2*nEvs];
      
      TPrecision worktmp = 0;
      int lwork = -1;
      
      int iworktmp = 0; 
      int liwork = -1;

      info = 0;
     
      int N = matrix.N();
      if(sizeof(TPrecision) == sizeof(double)){
        lapack::dsyevr_(&jobz, &range, &uplo, &N, (double *)matrix.data(),
            &N, (double *)&vl, (double *)&vh, &il, &ih, (double *)&tol,
            &nEigenvectors, (double*)ewTmp.data(), (double *)ev.data(), &N,
            isuppz, (double *)&worktmp, &lwork, &iworktmp, &liwork, &info);
      }
      else{
        lapack::ssyevr_(&jobz, &range, &uplo, &N, (float *)matrix.data(),
            &N, (float *)&vl, (float*)&vh, &il, &ih, (float *)&tol, 
            &nEigenvectors, (float *)ewTmp.data(), (float *)ev.data(), &N, 
            isuppz, (float *)&worktmp, &lwork, &iworktmp, &liwork, &info);
      }
 
      lwork = (int)worktmp;
      TPrecision *work = new TPrecision[lwork];
      liwork = (int)iworktmp;
      int *iwork = new int[liwork];

      if(sizeof(TPrecision) == sizeof(double)){
        lapack::dsyevr_(&jobz, &range, &uplo, &N, (double *)matrix.data(),
            &N, (double *)&vl, (double *)&vh, &il, &ih, (double *)&tol, 
            &nEigenvectors, (double*)ewTmp.data(), (double *)ev.data(), &N, 
            isuppz, (double *)work, &lwork, iwork, &liwork, &info);
      }
      else{
        lapack::ssyevr_(&jobz, &range, &uplo, &N, (float *)matrix.data(),
            &N, (float *)&vl, (float*)&vh, &il, &ih, (float *)&tol, 
            &nEigenvectors, (float*)ewTmp.data(), (float *)ev.data(), &N, 
            isuppz, (float *)work, &lwork, iwork, &liwork, &info);
      }

      ew = DenseVector<TPrecision>(nEigenvectors);
      for(int i=0; i<nEvs; i++){
        ew(i) = ewTmp(i);
      }
      ewTmp.deallocate();

      delete[] iwork;
      delete[] work;
      delete[] isuppz;

    };

};

}

#endif
