//
// Created by Aldo Miranda-Aguilar on 03/11/18.
//

#include <stdlib.h>
#include <math.h>
#include "MtMatrix.h"

MtMatrix::MtMatrix() { }

MtMatrix::MtMatrix(int theNumRows, int theNumCols, int theQ) {
    pNumRows = theNumRows;
    pNumCols = theNumCols;
    pQ = theQ;

    pElement = new int[pNumRows*pNumCols];
    //memset(pElement,0, pNumRows*pNumCols*sizeof(int));


    for(int i=0; i<32; i++) {
        if((pQ>>i) == 1) {
            pL = i;
            break;
        }
    }
}

MtMatrix::MtMatrix(const MtMatrix& theMat) {
    pNumRows = theMat.getNumRows();
    pNumCols = theMat.getNumCols();
    pQ = theMat.getQ();
    pL = theMat.getL();
    pElement = theMat.pElement;
}

MtMatrix::~MtMatrix() {}

// Assignment Operator
MtMatrix& MtMatrix::operator=(const MtMatrix& theMat) {
    if (&theMat == this)
        return *this;

    pNumRows = theMat.getNumRows();
    pNumCols = theMat.getNumCols();
    pQ = theMat.getQ();
    pL = theMat.getL();

    pElement = theMat.getElement();

/*
	pElement = new int[pNumRows*pNumCols];
	//memset(pElement,0, pNumRows*pNumCols*sizeof(int));

	for (int i=0; i<pNumRows; i++) {
		for (int j=0; j<pNumCols; j++) {
			(*this)(i,j) = theMat(i, j);
		}
	}
*/
    return *this;
}

// Addition of two matrices
MtMatrix MtMatrix::operator+(const MtMatrix& theMat)
{
    MtMatrix result(pNumRows, pNumCols, pQ);
    int* a = result.getElement();
    int* b = theMat.getElement();

    for (int i=0; i<pNumRows*pNumCols; i++) {
        a[i] = MOD(pElement[i]+b[i],pQ);
    }

    return result;
}

// Substraction of two matrices
MtMatrix MtMatrix::operator-(const MtMatrix& theMat)
{
    MtMatrix result(pNumRows, pNumCols, pQ);
    int* a = result.getElement();
    int* b = theMat.getElement();

    for (int i=0; i<pNumRows*pNumCols; i++) {
        a[i] = MOD(pElement[i]-b[i]+pQ,pQ);
    }

    return result;
}

MtMatrix MtMatrix::operator!()
{
    MtMatrix result(pNumRows, pNumCols, pQ);

    result = identity(pNumRows,1,pQ,pL) - *this;

    return result;
}

// Left multiplication of this matrix and another
MtMatrix MtMatrix::operator*(const MtMatrix& theMat)
{
    MtMatrix result(pNumRows, theMat.getNumCols(), pQ);
    int numCols = theMat.getNumCols();
    int* element = theMat.getElement();

    for (int i=0; i<pNumRows; i++) {
        int* a = pElement+i*pNumCols;
        for (int j=0; j<numCols; j++) {
            int* b = element+j;
            int s = 0;
            for(int k=0; k<pNumCols; k++) {
                s += a[k] * b[numCols*k];
            }
//			if(s>pQ) {
//				printf("mul : s>pQ : %d\n", s);
//			}
            result(i,j) = MOD(s,pQ);
        }
    }

    return result;
}

MtMatrix MtMatrix::identity(int theNumRows, int theInitValue, int theQ, int theL)
{
    MtMatrix result(theNumRows, theNumRows, theQ);
    for (int i=0; i<theNumRows; i++) {
        for(int j=0; j<theNumRows; j++) {
            if(i==j) result(i,i) = theInitValue;
            else	result(i,j) = 0;
        }
    }

    return result;
}



void MtMatrix::randomize(int theRange)
{
    for (int i=0; i<pNumRows; i++) {
        for (int j=0; j<pNumCols; j++) {
            (*this)(i,j) = rand() % theRange;
        }
    }
}

MtMatrix MtMatrix::flatten()
{
    MtMatrix result1 = invBitDecomp();
    return result1.bitDecomp();
}

MtMatrix MtMatrix::bitDecomp()
{
    MtMatrix result(pNumRows, pNumCols*pL, pQ);
    int* a = result.getElement();

    for (int i=0; i<pNumRows; i++) {
        for (int j=0; j<pNumCols; j++) {
            int x = pElement[i*pNumCols+j];
            int* e = a+i*pNumCols*pL+j*pL;
            for (int k=0; k<pL; k++) {
                e[k] = x & 1;
                x = x>>1;
            }
        }
    }

    return result;
}

MtMatrix MtMatrix::invBitDecomp()
{
    MtMatrix result(pNumRows, pNumCols/pL, pQ);
    int* e = pElement;
    int* a = result.getElement();

    for (int i=0; i<pNumRows; i++) {
        for (int j=0; j<pNumCols/pL; j++) {
            int x=0;
            e = pElement + i*pNumCols+j*pL;
            for(int k=0; k<pL; k++) {
                x = x + (e[k]<<k);
            }
            a[i*pNumCols/pL+j] = x;
        }
    }

    return result;
}

MtMatrix MtMatrix::getRowVector(int theRow)
{
    MtMatrix result(1, pNumCols, pQ);

    for (int j=0; j<pNumCols; j++) {
        result(0,j) = (*this)(theRow,j);
    }

    return result;
}