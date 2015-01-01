#ifndef LAPACKDEFS_H
#define LAPACKDEFS_H





//blas routines
namespace lapack{

extern "C"{

  //---- matrix matrix multiply ----//
  void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
      double *alpha, double *A, int *lda, double *B, int *ldb, double *beta,
      double *C, int *ldc);
    
  void sgemm_(char *transa, char *transb, int *m, int *n, int *k, 
      float *alpha, float *A, int *lda, float *B, int *ldb, float *beta, 
      float *C, int *ldc);
  



  //---- matrix vector multiply ----//
  void dgemv_(char *trans, int *m, int *n, double *alpha, double *A, int *lda, 
      double *x, int *incx, double *beta, double *y, int *incy);
    
  void sgemv_(char *trans, int *m, int *n, float *alpha, float *A, int *lda, 
      float *x, int *incx, float *beta, float *y, int *incy);





  //---- dot product ----/
  double ddot_(int *m, double *x, int *incx, double *y, int *incy);

  float sdot_(int *m, float *x, int *incx, float *y, int *incy);


  //---- QR decomposition ----//
  void dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);
  void sgeqrf_(int *m, int *n, float *a, int *lda, float *tau, float *work, int *lwork, int *info);

  void dorgqr_(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info);
  void sorgqr_(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *lwork, int *info);

  //---- LU decvomposition ----//
  void dgetrf_(int * m, int *n, double *a, int *lda, int *ipiv, int *info);

  void sgetrf_(int * m, int *n, float *a, int *lda, int *ipiv, int *info);
  
  //---- Cholesky decvomposition ----//
  void dpotrf_(char *u, int *n, double *a, int *lda, int *info);

  void spotrf_(char *u, int *n, float *a, int *lda, int *info);
  
  //---- inverse ---//
  void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int
      *lwork, int *info);

  void sgetri_(int *n, float *a, int *lda, int *ipiv, float *work, int
      *lwork, int *info);


  //---- inverse symmetric positive definite matrix ----/
  void dpotri_(char *uplo, int *n, double *A, int *lda, int *info);
  
  void spotri_(char *uplo, int *n, float *A, int *lda, int *info);

  //---- Solve SPD linear equations ----//
  void dposv_(char *uplo, int *n, int *nrhs, double *A, int *lda, double *B,
      int *ldb, int *info );

  void sposv_(char *uplo, int *n, int *nrhs, float *A, int *lda, float *B,
      int *ldb, int *info );

  void dposvx_( char *fact, char *uplo, int *n, int *nrhs, double *A, int *lda,
      double *af, int *ldaf, char *equed, double *s, double *b, int *ldb, double
      *x, int *ldx, double *rcond, double *ferr, double *berr, double *work,
      int *iwork, int *info );
  
  void sposvx_( char *fact, char *uplo, int *n, int *nrhs, float *A, int *lda,
      float *af, int *ldaf, char *equed, float *s, float *b, int *ldb, float
      *x, int *ldx, float *rcond, float *ferr, float *berr, float *work,
      int *iwork, int *info );

  //---- Solve linear equations ----/
  void dgesv_( int *n, int *nrhs, double *A, int *lda, int *ipiv, double *b, int
      *ldb, int *info );

  void sgesv_( int *n, int *nrhs, float *A, int *lda, int *ipiv, float *b, int
      *ldb, int *info );

  //---- Least squares  or minimum norm with QR or LQ, A is assumed to have full rank----// 
  void dgels_(char *trans, int *m, int *n, int *nhrs, double *A, int *lda, double
      *B, int *ldb, double *work, int *lwork, int *info);

  void sgels_(char *trans, int *m, int *n, int *nhrs, float *A, int *lda, float 
      *B, int *ldb, float *work, int *lwork, int *info);

  //---- Least squares  or minimum norm with SVD// 
  void dgelsd_(int *m, int *n, int *nhrs, double *A, int *lda, double
      *B, int *ldb, double *S, double *rcond, int *rank, double *work, int
      *lwork, int *iwork, int *info);

  void sgelsd_(int *m, int *n, int *nhrs, float *A, int *lda, float
      *B, int *ldb, float *S, float *rcond, int *rank, float *work, int *lwork,
      int *iwork, int *info);


  //--- Cholesky condition number ---//
  void dpocon(char *uplo, int *n, double *a, int *lda, double *anorm, double
      *rcond, double *work, int *iwork, int *info);
  
  void spocon(char *uplo, int *n, float *a, int *lda, float *anorm, float
      *rcond, float *work, int *iwork, int *info);


   //---- SVD ---//
   void dgesvd_(char *JOBU, char *JOBVT, int* M, int *N, double *A, 
               int *LDA, double *S, double *U, int* LDU, double *VT, 
               int *LDVT, double *WORK, int *LWORK, int *INFO );
   
   void sgesvd_(char *JOBU, char *JOBVT, int* M, int *N, float *A, 
               int *LDA, float *S, float *U, int* LDU, float *VT, 
               int *LDVT, float *WORK, int *LWORK, int *INFO );


   //------ lapack symmetric eigensystem routine ----//
   void dsyevr_( char *jobz, char *range, char *uplo, int *N, double *A,
                    int *LDA, double *VL, double *VU, int *IL, int *IU, 
                    double *ABSTOL, int *M, double *W, double *Z, int *LDZ, 
                    int *ISUPPZ, double *WORK, int *LWORK, int *IWORK, 
                    int *LIWORK, int *INFO );

    void ssyevr_( char *jobz, char *range, char *uplo, int *N, float *A,
                    int *LDA, float *VL, float *VU, int *IL, int *IU, 
                    float *ABSTOL, int *M, float *W, float *Z, int *LDZ, 
                    int *ISUPPZ, float *WORK, int *LWORK, int *IWORK, 
                    int *LIWORK, int *INFO );
}; 

};

#endif
