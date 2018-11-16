#include <iostream>
#include <stdlib.h>
#include "MtMatrix.h"
#include "MtCiphertext.h"
#include "cuda_bridge.h"

#include <chrono>
#include <thread>

const int l = 9;
const int n = 500;	// Note that to support *, at least, N < 2^{l-2}/2
const int q = (1<<(l));
//const int m = (2*n*l);
const int m = 150;//(n+1)*(l+3);	// In BV14, m = (n+1)(l+O(1))
const int N = (n+1)*l;

// private keys
MtMatrix s(n+1,1,q);
MtMatrix v(N,1,q);

// public keys
MtMatrix A(m,n+1,q);


// the following are temporary variables. They will be moved into their functions. Currently they are moved out to the global area for debugging purpose.
//ZqMatrix B(m,n,q), e(m,1,q);
MtMatrix t(n,1,q);
MtMatrix b(m,1,q);
MtMatrix R(N,m,q);

int getZq()
{
    return rand()%q;
}

int getX()
{
    return 0;
    //return (rand()%3)-1;
}

void generateSK()
{
    s(0,0) = 1;
    int i, j;
    for(i=0; i<n;i++) {
        t(i,0) = getZq();
        s(i+1,0) = MOD(q-t(i,0),q);

    }
    for(i=0; i<n+1; i++)
        for(j=0; j<l; j++)
            v(i*l+j, 0) = MOD((s(i,0)* (1<<j)),q);
}

void generatePK()
{
    int i,j;
    // B <-- Z^{m*n}_q
    MtMatrix B(m,n,q), e(m,1,q);
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            B(i,j) = getZq();
        }
        e(i,0) = getX();
    }

    //ZqMatrix b(m,1,q);
    //ZqMatrix b = B*t+e;
    b = B*t+e;

    for(i=0;i<m; i++) {
        A(i,0) = b(i,0);
        for(j=0; j<n; j++) {
            A(i,j+1) = B(i,j);
        }
    }
}

MtMatrix generateNewRandom()
{
    int i,j;
    // R <-- Z^{N*n}_q
    MtMatrix R(N,n,q), e(N,1,q);
    for(i=0; i<N; i++) {
        for(j=0; j<n; j++) {
            R(i,j) = getZq();
        }
        e(i,0) = getX();
    }

    MtMatrix b = R*t+e;
    MtMatrix r(N,N,q);
    for(i=0;i<N; i++) {
        r(i,0) = b(i,0);
        for(j=0; j<n; j++) {
            r(i,j+1) = R(i,j);
        }
    }

    return r;
}

void genKeys()
{
    // generate sk
    printf("generate sk\n");
    generateSK();
    printf("generate pk\n");
    // generate pk
    generatePK();
}

MtMatrix enc(int mu)
{
    //printf("mu: %d\n", mu);
    //printf("N: %d\n", N);
    //printf("m: %d\n", m);
    //printf("q: %d\n", q);
    //clock_t start,end;
    //ZqMatrix R(N,m,q);

    //printf("Rand R start\n");
    R.randomize(2);
    //printf("Rand R end\n");

    //start = clock();
    //printf("Matmul R*A start\n");
    MtMatrix c = R*A;
    //printf("Matmul R*A end\n");
    //end = clock();
    //printf("c mul time : %f\n", (float)(end-start)/CLOCKS_PER_SEC);
    //printf("bit decomp start\n");
    c = c.bitDecomp();
    //printf("bit decomp end\n");

    if(mu) {
        //printf("add active bit start\n");
        c = c+MtMatrix::identity(N, mu);
        //printf("add active bit end\n");

        return c.flatten();
    } else {
        return c;
    }
}

int dec(MtMatrix& C,int index=l-1)
{
    MtMatrix C_i = C.getRowVector(index);

    MtMatrix  x_i = C_i * v;

    // we compute \lfloor x_i / v_i \rceil.
    // If v_i/2 < x_i <= 3v_i/2 then it is 1, otherwise it is 0.
    int v_i = 1<<index;
    //float t = 2.0*(x_i(0,0))/v_i;
    //if(t>1.0 && t<=3.0) return 1;
    //else return 0;

    int t = (x_i(0,0)+v_i/2) / v_i;

    return t&1;
}

MtMatrix nand(MtCiphertext& a, MtCiphertext& b) {

    //std::cout << "a dims: "<< a.getNumRows() << "," << a.getNumCols() << std::endl;
    //std::cout << "b dims: "<< b.getNumRows() << "," << b.getNumCols() << std::endl;

    return ((a * b) - MtMatrix::identity(N)).flatten();
}


int main(int argc, char* argv[]) {
    std::cout << "Matcha!" << std::endl;

    srand(1);
    genKeys();

    MtCiphertext c0 = enc(0);
    MtCiphertext c1 = enc(1);

    //MtMatrix  rnand00 = (c0*c0) + c1;
    /*MtMatrix  rnand01 = (c0*c1) + c1;
    MtMatrix  rnand10 = (c1*c0) + c1;
    MtMatrix  rnand11 = (c1*c1) + c1;*/

    MtCiphertext  rnand00 = MtCiphertext( nand(c0, c0));

    printf("dec nand00 xv1: %d\n", dec(rnand00));

    MtCiphertext  rnand00_1 = MtCiphertext(nand(c0, rnand00));
    printf("dec nand00 hop1 xv1: %d\n", dec(rnand00_1));

    MtCiphertext  rnand00_2 = MtCiphertext(nand(rnand00, rnand00_1));
    printf("dec nand00 hop2 xv0: %d\n", dec(rnand00_2));

    MtCiphertext  rnand00_3 = MtCiphertext(nand(rnand00_1, rnand00_2));
    printf("dec nand00 hop3 xv1: %d\n", dec(rnand00_3));

    MtCiphertext  rnand00_4 = MtCiphertext(nand(rnand00_2, rnand00_3));
    printf("dec nand00 hop4 xv1: %d\n", dec(rnand00_4));

    MtCiphertext  rnand00_5 = MtCiphertext(nand(rnand00_3, rnand00_4));
    printf("dec nand00 hop5 xv0: %d\n", dec(rnand00_5));

    MtCiphertext  rnand00_6 = MtCiphertext(nand(rnand00_4, rnand00_5));
    printf("dec nand00 hop6 xv1: %d\n", dec(rnand00_6));

    MtCiphertext  rnand00_7 = MtCiphertext(nand(rnand00_5, rnand00_6));
    printf("dec nand00 hop7 xv1: %d\n", dec(rnand00_7));

    MtCiphertext  rnand00_8 = MtCiphertext(nand(rnand00_6, rnand00_7));
    printf("dec nand00 hop8 xv0: %d\n", dec(rnand00_8));

    MtCiphertext  rnand00_9 = MtCiphertext(nand(rnand00_7, rnand00_8));
    printf("dec nand00 hop8 xv1: %d\n", dec(rnand00_9));

    MtCiphertext  rnand00_10 = MtCiphertext(nand(rnand00_8, rnand00_9));
    printf("dec nand00 hop8 xv1: %d\n", dec(rnand00_10));


    //printf("dec nand00 xv1: %d\n", dec(rnand00));

    //MtCiphertext  rnand00_1 = MtCiphertext(nand(c0, rnand00));
    //printf("dec nand00 hop1 xv1: %d\n", dec(rnand00_1));

    /*short a = 0, b = 0;
    MtCiphertext eA = enc(0);
    MtCiphertext eB = enc(0);

    short pNand = ((a & b) ^ 1);
    MtCiphertext eNand = MtCiphertext(nand(eA, eB));

    short pNand_m1 = 0;
    MtCiphertext eNand_m1 = enc(0);

    int count = 0;

    while(1) {
        count++;
        //if((count % 5) == 0) {
            int deceNand = dec(eNand);
            std::cout << "[" << count << "] " << "res: " << pNand << " (e)res: " << deceNand << " "
                      << (deceNand == pNand ? "PASS" : "FAIL") << std::endl;
        //}

//std::this_thread::sleep_for(std::chrono::milliseconds(1000));

        short tmp = pNand;
        MtCiphertext eTmp = eNand;

        pNand = ((pNand_m1 & pNand) ^ 1);
        eNand = MtCiphertext(nand(eNand_m1, eNand));

        pNand_m1 = tmp;
        eNand_m1 = eTmp;

    }*/

    /*printf("dec nand01 : %d\n", dec(rnand01));
    printf("dec nand10 : %d\n", dec(rnand10));
    printf("dec nand11 : %d\n", dec(rnand11));*/

  /*if (checkCmdLineFlag(argc, (const char **)argv, "help") ||
      checkCmdLineFlag(argc, (const char **)argv, "?")) {
    printf("Usage -device=n (n >= 0 for deviceID)\n");
    printf("      -wA=WidthA -hA=HeightA (Width x Height of Matrix A)\n");
    printf("      -wB=WidthB -hB=HeightB (Width x Height of Matrix B)\n");
    printf("  Note: Outer matrix dimensions of A & B matrices" \
           " must be equal.\n");

    exit(EXIT_SUCCESS);
  }*/

  /*int *h_A, *h_B, *h_C = NULL;
  int W = atoi(argv[1]);
  long len = W * W;
  h_A = new int[len];
  h_B = new int[len];
  h_C = new int[len];

  for(int i=0; i<len; i++) {
      h_A[i] = 1;
      h_B[i] = 1;
  }

  int block_size = 32;

  int dimsAx = 3;//10 * 4 * block_size;
  int dimsAy = 3;//10 * 4 * block_size;
  int dimsBx = 3;//10 * 8 * block_size;
  int dimsBy = 3;//10 * 4 * block_size;

  //printf("MatrixA(%d,%d), MatrixB(%d,%d)\n", dimsAx, dimsAy,
  //       dimsBx, dimsBy);


  std::cout << "h_A0: " << h_A[0] << std::endl;
  std::cout << "h_A1: " << h_A[1] << std::endl;
  std::cout << "h_A2: " << h_A[2] << std::endl;
  std::cout << "h_A3: " << h_A[3] << std::endl;
  std::cout << "h_A4: " << h_A[4] << std::endl;
  std::cout << "h_A5: " << h_A[5] << std::endl;
  std::cout << "h_A6: " << h_A[6] << std::endl;
  std::cout << "h_A7: " << h_A[7] << std::endl;
  std::cout << "h_A8: " << h_A[8] << std::endl;  

  std::cout << "h_B0: " << h_B[0] << std::endl;
  std::cout << "h_B1: " << h_B[1] << std::endl;
  std::cout << "h_B2: " << h_B[2] << std::endl;
  std::cout << "h_B3: " << h_B[3] << std::endl;
  std::cout << "h_B4: " << h_B[4] << std::endl;
  std::cout << "h_B5: " << h_B[5] << std::endl;
  std::cout << "h_B6: " << h_B[6] << std::endl;
  std::cout << "h_B7: " << h_B[7] << std::endl;
  std::cout << "h_B8: " << h_B[8] << std::endl;*/

  //int matrix_result = MatrixMultiply(block_size, dimsAx, dimsAy, dimsBx, dimsBx, h_A, h_B, h_C);

  //exit(matrix_result);
  //h_C[0] = 1;

  //MatrixMultiplication(h_A, h_B, h_C, W);
  //matrixMultiplication(h_A, h_B, h_C, 3);

  /*for(long i = 0; i < len; i++) {
        if(h_C[i] != (W % 512)) {
            std::cout << "[" << i << "]" << h_C[i] << std::endl;
            //break;
        }
    }*/

    //mmul();

    std::cout << "Done." << std::endl;
    return 0;
}
