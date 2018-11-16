//
// Created by Aldo Miranda-Aguilar on 03/11/18.
//
#define PACK_OPTIMIZATION 1
#define POPCNT_OPTIMIZATION 1

#ifdef POPCNT_OPTIMIZATION
#include <smmintrin.h>
#endif

#include <iostream>
#include <stdint.h>
#include <stdio.h>

#include "cuda_bridge.h"

#include "MtCiphertext.h"

MtMatrix MtCiphertext::operator+(const MtMatrix& theMat)
{
    MtCiphertext result = MtMatrix::operator+(theMat);
    return result.flatten();
}

MtMatrix MtCiphertext::operator-(const MtMatrix& theMat)
{
    MtCiphertext result = MtMatrix::operator-(theMat);
    return result.flatten();
}

MtMatrix MtCiphertext::operator!()
{
    MtCiphertext result = MtMatrix::operator!();
    return result.flatten();
}

// Left multiplication of this matrix and another
// It is a rectangular matrix or NxN
MtMatrix MtCiphertext::operator*(const MtMatrix& theMat) {

    MtCiphertext result(getNumRows(), theMat.getNumCols(), getQ());
    //printf("pQ: %d", pQ);
    int num = pNumRows;
    int *element = theMat.getElement();
    int *resElement = result.getElement();
    //int *pGPUElement = new int[pNumRows*pNumCols];

    MatrixMultiplication(pElement, element, resElement, theMat.getNumCols());

/*#ifdef PACK_OPTIMIZATION
    // packed matrix
    int wSize = sizeof(uint64_t)*8;
    int packedNum = (num+wSize-1)/wSize;
    uint64_t packedA[num*packedNum];
    uint64_t packedB[num*packedNum];

//	clock_t start = clock();
    pack(packedA, pElement, num, packedNum);
    packTranspose(packedB, element, num, packedNum);
//	clock_t end = clock();

//	printf("pack time : %f\n", (float)(end-start)/CLOCKS_PER_SEC);

    for (int i=0; i<num; i++) {
        uint64_t* a = packedA+i*packedNum;
        for (int j=0; j<num; j++) {
            uint64_t* b = packedB+j*packedNum;
            int s = 0;
            for(int k=0; k<packedNum; k++) {
                uint64_t t = a[k] & b[k];
#ifdef POPCNT_OPTIMIZATION
                s += _mm_popcnt_u64(t);
#else
                s += countOnes(t);
#endif
            }
//			if(s>pQ) {
//				printf("s>pQ\n");
//			}
            result(i,j) = MOD(s,pQ);
        }
    }
#else
    int *a = pElement;
    int *b = element;

    for (int i = 0; i < num; i++) {
        printf("matmul i: %d\n", i);
        int *a = pElement + i * num;
        for (int j = 0; j < num; j++) {
            //printf("matmul j: %d\n", j);
            int *b = element + j;
            int s = 0;
            for (int k = 0; k < num; k++) {
                //printf("matmul k: %d\n", k);
                s += a[k] & b[num * k];
            }

            result(i, j) = MOD(s, pQ);
        }
    }
#endif*/
    /*int lll = pNumRows*pNumCols;
    std::cout << pNumRows << "," << pNumCols << std::endl;
    std::cout << result.getNumRows() << "," << result.getNumCols() << std::endl;
    for(int i = 0; i < lll; i++) {
        if(resElement[i] != pGPUElement[i]) {
            std::cout << "[" << i << "]" << resElement[i] << "," << pGPUElement[i] << std::endl;
            break;
        }
    }*/
    //printf("matmul end\n");
    return result.flatten();
}

void MtCiphertext::pack(uint64_t* packMat, int* mat, int num, int packedNum)
{
    int wSize = sizeof(uint64_t)*8;
    for(int i=0; i<num; i++) {
        uint64_t s = 0;
        for(int j=0; j<num; j++) {
            s = (s<<1) | mat[i*num+j];
            if((j+1)%wSize==0) {
                packMat[i*packedNum+j/wSize] = s;
                s = 0;
            }
        }
        packMat[i*packedNum+num/wSize] = s;
    }
}

void MtCiphertext::packTranspose(uint64_t* packMat, int* mat, int num, int packedNum)
{
    int wSize = sizeof(uint64_t)*8;
    for(int j=0; j<num; j++) {
        uint64_t s = 0;
        for(int i=0; i<num; i++) {
            s = (s<<1) | mat[i*num+j];
            if((i+1)%wSize==0) {
                packMat[j*packedNum + i/wSize] = s;
                s = 0;
            }
        }
        packMat[j*packedNum + num/wSize] = s;
    }
}

int MtCiphertext::countOnes(uint64_t v)
{
    int c = 0;
    for(int i=0; i<sizeof(uint64_t)*8; i++) {
        c += v&1;
        v >>= 1;
    }

    return c;
}
