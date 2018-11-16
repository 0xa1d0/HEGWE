#include <iostream>
#include <cuda_runtime.h>
#include "book.h"
#include "cuda_bridge.h"

#define BLOCK_SIZE 16

__global__ void gpu_square_matrix_mult(int *d_a, int *d_b, int *d_result, int n) 
{
    __shared__ int tile_a[BLOCK_SIZE][BLOCK_SIZE];
    __shared__ int tile_b[BLOCK_SIZE][BLOCK_SIZE];

    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    int tmp = 0;
    int idx;

    #pragma unroll
    for (int sub = 0; sub < gridDim.x; ++sub) 
    {
        idx = row * n + sub * BLOCK_SIZE + threadIdx.x;
        if(idx >= n*n)
        {
            // n may not divisible by BLOCK_SIZE
            tile_a[threadIdx.y][threadIdx.x] = 0;
        }
        else
        {
            tile_a[threadIdx.y][threadIdx.x] = d_a[idx];
        }

        idx = (sub * BLOCK_SIZE + threadIdx.y) * n + col;
        if(idx >= n*n)
        {
            tile_b[threadIdx.y][threadIdx.x] = 0;
        }  
        else
        {
            tile_b[threadIdx.y][threadIdx.x] = d_b[idx];
        }
        __syncthreads();

        for (int k = 0; k < BLOCK_SIZE; ++k) 
        {
            tmp += tile_a[threadIdx.y][k] * tile_b[k][threadIdx.x];
        }
        __syncthreads();
    }
    if(row < n && col < n)
    {
        d_result[row * n + col] = tmp % 512;
    }
}
int *Md = NULL, *Nd = NULL, *Pd = NULL;

void MatrixMultiplication(int *&M, int *&N, int *&P, int Width) {

  int size = Width * Width * sizeof(int);

  // allocate memory on the GPU
  //if(Md == NULL) {
    //std::cout << "cudaMalloc" << std::endl;
    HANDLE_ERROR( cudaMalloc((void**)&Md, size) );
  //} else {
    //std::cout << "ALLOCATED" << std::endl;
  //}
  //if(Nd == NULL)
    HANDLE_ERROR( cudaMalloc((void**)&Nd, size) );
  //if(Pd == NULL)
    HANDLE_ERROR( cudaMalloc((void**)&Pd, size) );

  // transfer M and N to device memory
  HANDLE_ERROR( cudaMemcpy(Md, M, size, cudaMemcpyHostToDevice) );
  HANDLE_ERROR( cudaMemcpy(Nd, N, size, cudaMemcpyHostToDevice) );

  unsigned int grid_rows = (Width + BLOCK_SIZE - 1) / BLOCK_SIZE;
  unsigned int grid_cols = (Width + BLOCK_SIZE - 1) / BLOCK_SIZE;
  dim3 dimGrid(grid_cols, grid_rows);
  dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);

  //std::cout << "Width: " << Width << std::endl;
  //std::cout << "Width/32: " << Width/32 << std::endl;

  //Kernel<<<dimGrid, dimBlock>>>( Md, Nd, Pd, Width);
  gpu_square_matrix_mult<<<dimGrid, dimBlock>>>(Md, Nd, Pd, Width);
  //std::cout << "Width/32: " << Width/32 << std::endl;

  // transfer P from device     
  //std::cout << "P: " << P[0] << std::endl;
  HANDLE_ERROR( cudaMemcpy(P,Pd,size,cudaMemcpyDeviceToHost) );
  cudaThreadSynchronize();

  HANDLE_ERROR( cudaFree(Md) );
  HANDLE_ERROR( cudaFree(Nd) );
  HANDLE_ERROR( cudaFree(Pd) );
}

