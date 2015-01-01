#ifndef SYMMETRICDENSEMATRIX_H
#define SYMMETRICDENSEMATRIX_H

#include "DenseMatrix.h"

namespace FortranLinalg{

template <typename TPrecision>
class SymmetricDenseMatrix : public DenseMatrix<TPrecision>{

  public:
    SymmetricDenseMatrix(bool rowMaj = true) : DenseMatrix<TPrecision>(rowMaj){
    };

    SymmetricDenseMatrix(int nrows, bool rowMaj = true) : DenseMatrix<TPrecision>(nrows, nrows, rowMaj){
    };

};

}


#endif
