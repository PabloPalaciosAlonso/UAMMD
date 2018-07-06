
#ifndef LBM_CUH
#define LBM_CUH

#include"uammd.cuh"
#include"Integrator/Integrator.cuh"
#include<thrust/device_vector.h>
#include<fstream>
namespace uammd{
  namespace Hydro{
    namespace LBM{    
      class D3Q19: public Integrator{
	int steps;
	thrust::device_vector<real> sourceGrid, destGrid;
	thrust::device_vector<int> cellType;
	Grid grid;
	static constexpr int numberVelocities = 19;
	real soundSpeed, relaxTime, dt, viscosity;
	std::ofstream out;
      public:
	struct Parameters{
	  Box box;
	  int3 ncells;
	  real soundSpeed;
	  real relaxTime;
	  real dt;
	  real viscosity;
	};
	D3Q19(shared_ptr<ParticleData> pd,
	      shared_ptr<System> sys,
	      Parameters par);    
	virtual void forwardTime() override;
	void write();
	void writePNG();
      };
    }
  }
}

#include"LBM.cu"
#endif
