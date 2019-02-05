/*Raul P. Pelaez 2016. BDHI Cholesky submodule.

  See BDHI_EulerMaruyama.cuh to see how it is used
  
  Any BDHI method needs to do only three things: compute M·F, sqrt(M)·dW, div(M).
  

  BDHI::Cholesky stores the full Mobility Matrix and computes the stochastic term via cholesky decomposition.

  sqrt(M)dW = B·dW -> D = B·B^T 

  It uses cuBLAS for Mv products and cuSOLVER for Cholesky decomposition


*/

#ifndef BDHI_CHOLESKY_CUH
#define BDHI_CHOLESKY_CUH

#include"System/System.h"
#include"ParticleData/ParticleData.cuh"
#include"ParticleData/ParticleGroup.cuh"
#include"global/defines.h"
#include"BDHI.cuh"
#include<cublas_v2.h>
#include<cusolverDn.h>
#include<curand.h>

namespace uammd{
  namespace BDHI{
    class Cholesky{
    public:
      using Parameters = BDHI::Parameters;
      Cholesky(shared_ptr<ParticleData> pd,
	       shared_ptr<ParticleGroup> pg,
	       shared_ptr<System> sys,
	       Parameters par);
      ~Cholesky();
      void init();
      void setup_step(              cudaStream_t st = 0);
      void computeMF(real3* MF,     cudaStream_t st = 0);    
      void computeBdW(real3* BdW,   cudaStream_t st = 0);  
      void computeDivM(real3* divM, cudaStream_t st = 0);
      void finish_step(cudaStream_t st = 0){}
    
    private:
      shared_ptr<ParticleData> pd;
      shared_ptr<ParticleGroup> pg;
      shared_ptr<System> sys;
      
      thrust::device_vector<real> mobilityMatrix; /*The full mobility matrix*/
      thrust::device_vector<real3> force3;

      bool isMup2date;
    
      /*CUBLAS*/
      cublasStatus_t status;
      cublasHandle_t handle;
      /*CUSOLVER*/
      cusolverDnHandle_t solver_handle;
      /*Cusolver temporal storage*/
      int h_work_size;
      real *d_work;
      int *d_info;

      curandGenerator_t curng;
      /*Rodne Prager Yamakawa device functions and parameters*/
      RotnePragerYamakawa rpy;
      Parameters par;
      
    };
  }
}
#include"BDHI_Cholesky.cu"
#endif
