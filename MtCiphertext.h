//
// Created by Aldo Miranda-Aguilar on 03/11/18.
//

#ifndef MATCHA_MTCIPHERTEXT_H
#define MATCHA_MTCIPHERTEXT_H

#include <stdint.h>

#include "MtMatrix.h"

class MtCiphertext : public MtMatrix
{
public:
    MtCiphertext(){}
    MtCiphertext(int theNumRows, int theNumCols, int theQ) : MtMatrix(theNumRows, theNumCols, theQ){}
    MtCiphertext(const MtMatrix& theMat) : MtMatrix(theMat){}
    virtual ~MtCiphertext(){}

    /*virtual*/ MtMatrix operator+(const MtMatrix& theMat);
    /*virtual*/ MtMatrix operator-(const MtMatrix& theMat);
    /*virtual*/ MtMatrix operator*(const MtMatrix& theMat);
    /*virtual*/ MtMatrix operator!();

    void pack(uint64_t* packMat, int* mat, int numRows, int numCols);
    void packTranspose(uint64_t* packMat, int* mat, int numRows, int numCols);

    int countOnes(uint64_t v);
};

#endif //MATCHA_MTCIPHERTEXT_H
