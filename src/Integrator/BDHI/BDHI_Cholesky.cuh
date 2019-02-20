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


      real getHydrodynamicRadius(){
	return par.hydrodynamicRadius;
      }
      real getSelfMobility(){
	long double rh = par.hydrodynamicRadius;
	if(rh<0) return -1.0;
	else{
	  long double L = box.boxSize.x;
	  return  1.0l/(6.0l*M_PIl*viscosity*rh)*(1.0l
						  -2.837297l*rh/L
						  +(4.0l/3.0l)*M_PIl*pow(rh/L,3)
						  -27.4l*pow(rh/L,6.0l));
	}
      }

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
