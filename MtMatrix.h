//
// Created by Aldo Miranda-Aguilar on 03/11/18.
//

#ifndef MATCHA_MTMATRIX_H
#define MATCHA_MTMATRIX_H

//Optimized for GSW. q = power of 2

#include <vector>

#define MOD(x,q) ((x) & ((q)-1))

class MtMatrix {
protected:
    int* pElement;
    int pNumRows;
    int pNumCols;
    int pQ; 	// modulus Q
    int pL;

public:
    MtMatrix();
    MtMatrix(int theNumRows, int theNumCols, int theQ);
    MtMatrix(const MtMatrix& theMat);
    virtual ~MtMatrix();

    // Operator overloading, for "standard" mathematical matrix operations
    MtMatrix& operator=(const MtMatrix& theMat);

    // Matrix mathematical operations
    virtual MtMatrix operator+(const MtMatrix& theMat);
    virtual MtMatrix operator-(const MtMatrix& theMat);
    virtual MtMatrix operator*(const MtMatrix& theMat);
    virtual MtMatrix operator!();

    // Access the individual elements
    int& operator()(const int& row, const int& col) { return pElement[row*pNumCols+col]; }
    const int& operator()(const int& row, const int& col) const { return pElement[row*pNumCols+col]; }

    // GSW related functions
    static MtMatrix identity(int theNumRows, int theInitValue=1,int theQ=1, int theL=1);

    void randomize(int range);

    MtMatrix flatten();
    MtMatrix bitDecomp();
    MtMatrix invBitDecomp();

    MtMatrix getRowVector(int theRow);

    // Access elements
    int* getElement() const { return pElement; }

    // Access the row and column sizes
    int getNumRows() const { return pNumRows; }
    int getNumCols() const { return pNumCols; }
    int getQ() const { return pQ; }
    int getL() const { return pL; }

    // setup parameters
    void setQ(int theQ) { pQ = theQ; }
    void setL(int theL) { pL = theL; }
};

#endif //MATCHA_MTMATRIX_H
