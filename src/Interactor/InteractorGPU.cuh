/*
Raul P. Pelaez 2016. Interactor GPU kernels and callers.

Functions to compute the force acting on each particle and integrate movement

Neighbour list GPU implementation using hash short with cell index as hash.



References:
http://docs.nvidia.com/cuda/samples/5_Simulations/particles/doc/particles.pdf

More info in the .cu
*/

#ifndef INTERACTORGPU_CUH
#define INTERACTORGPU_CUH

//Stores some simulation parameters to upload as constant memory.
struct InteractorParams{
  float cellSize, invCellSize;
  int ncells;
  int ycells, xcells, zcells;  
  float L, invL, rcut, invrc2;
};

void initInteractorGPU(InteractorParams m_params, float *potDevPtr, size_t potSize,
	     uint *cellStart, uint *cellEnd, uint* particleIndex, uint ncells,
	     float4 *sortPos, uint N);


void calcCellIndex(float4 *pos, uint *cellIndex, uint *particleIndex, uint N);

void sortCellIndex(uint *&cellIndex, uint *&particleIndex, uint N);

void reorderAndFind(float4 *sortPos,
		    uint *cellIndex, uint *particleIndex, 
		    uint *cellStart, uint *cellEnd, uint ncells,
		    float4 *pos, uint N);

void computeForce(float4 *sortPos, float4 *force,
		  uint *cellStart, uint *cellEnd,
		  uint *particleIndex, 
		  uint N);


#endif








