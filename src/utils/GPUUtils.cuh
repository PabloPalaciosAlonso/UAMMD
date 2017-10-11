/*Raul P. Pelaez 2017. Some GPU utilities

  
 */
#ifndef GPUUTILS_CUH
#define GPUUTILS_CUH


#include<iostream>
#include"utils/debugTools.cuh"


//Fill any iterator with the same value. It is much faster than cudaMemset
template<class T, class OutputIterator>
__global__ void fillWithGPU(OutputIterator array, T value, int N){
  int id = blockIdx.x*blockDim.x + threadIdx.x;
  if(id>=N) return;

  array[id] = value;	            
}

template<class T, class OutputIterator, class Iterator>
__global__ void fillWithGPU(OutputIterator array, Iterator indexIterator, T value, int N){
  int id = blockIdx.x*blockDim.x + threadIdx.x;
  if(id>=N) return;
  int i = indexIterator[id];
  array[i] = value;
}






#endif