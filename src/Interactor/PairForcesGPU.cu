/*Raul P. Pelaez 2016. Short range pair forces Interactor GPU callers and kernels.

  The Neighbour list is constructed in the GPU as follows:
  
  1-Compute a hash for each particle based on its cell. Store in particleHash, also fill particleIndex with the index of each particle (particleIndex[i] = i)
  2-Sort particleIndex based on particleHash (sort by key). This way the particles in a same cell are one after the other in particleIndex. The Morton hash also improves the memory acces patter in the GPU.
  3-Fill cellStart and cellEnd with the indices of particleIndex in which a cell starts and ends. This allows to identify where all the [indices of] particles in a cell are in particleIndex, again, one after the other.
  
  The transversal of this cell list is done by transversing, for each particle, the 27 neighbour cells of that particle's cell.
  

  Force is evaluated using table lookups (with texture memory)

  Throughout this file, all __device__ and __global__ functions have access to the params and, in the case of DPD, paramsDPD structs in constant memory. See PairForcesGPU.cuh for a list of the available information in this structures.

  You can add new information to constant params by declaring it in PairForcesGPU.cuh and initializing it either in the initGPU function  before uploading or in the PairForces constructor, before calling initGPU.

  References:
  http://docs.nvidia.com/cuda/samples/5_Simulations/particles/doc/particles.pdf

  TODO:
  100- Implement many threads per particle in force compute
  100- Make number of blocks and threads to autotune
  100- Improve the transversing of the 27 neighbour cells
  90- Implement energy and virial compute in PairForcesDPD, maybe take it to another file
  80- General functions like apply_pbc should be made global to ease development.
  10- Find a way to properly handle the alternate arrays in sortCellIntex
*/

#include<cub/cub.cuh>
#include<curand_kernel.h>
#include"PairForcesGPU.cuh"
#include"utils/helper_math.h"
#include"utils/helper_gpu.cuh"

#include<thrust/device_ptr.h>
#include<thrust/for_each.h>
#include<thrust/iterator/zip_iterator.h>
#include<thrust/sort.h>
#include<iostream>


#define BLOCKSIZE 128

using namespace thrust;
using std::cerr;
using std::endl;

namespace pair_forces_ns{
  __constant__ Params params; //Simulation parameters in constant memory, super fast access
  __constant__ ParamsDPD paramsDPD; //Simulation parameters in constant memory, super fast access
  __constant__ ullint seedGPU; //seed of the rng for DPD
  
  uint GPU_Nblocks;
  uint GPU_Nthreads;
  
  //Initialize gpu variables 
  void initGPU(Params &m_params, uint N){    
    /*Precompute some inverses to save time later*/
    m_params.invrc2 = 1.0f/(m_params.rcut*m_params.rcut);
    m_params.invrc = 1.0f/(m_params.rcut);
    m_params.invL = 1.0f/m_params.L;
    m_params.invCellSize = 1.0f/m_params.cellSize;
    m_params.getCellFactor = 0.5f*m_params.L*m_params.invCellSize;
    m_params.gridPos2CellIndex = make_int3( 1,
					    m_params.cellDim.x,
					    m_params.cellDim.x*m_params.cellDim.y);
    
    if(!m_params.texPos || !m_params.texSortPos || !m_params.texCellStart ||
       !m_params.texCellEnd || !m_params.texForce || !m_params.texEnergy){
      cerr<<"Problem setting up textures!!"<<endl;
      exit(1);

    }
    /*Upload parameters to constant memory*/
    gpuErrchk(cudaMemcpyToSymbol(params, &m_params, sizeof(Params)));

    /*Each particle is asigned a thread*/
    GPU_Nthreads = BLOCKSIZE<N?BLOCKSIZE:N;
    GPU_Nblocks  =  N/GPU_Nthreads +  ((N%GPU_Nthreads!=0)?1:0); 
  }


  
  void initDPDGPU(ParamsDPD &m_params){
    gpuErrchk(cudaMemcpyToSymbol(paramsDPD, &m_params, sizeof(ParamsDPD)));
  }

  /****************************HELPER FUNCTIONS*****************************************/
  //MIC algorithm
  // template<typename vecType>
  // inline __device__ void apply_pbc(vecType &r){
  //   float3 r3 = make_float3(r.x, r.y, r.z);
  //   float3 shift = (floorf(r3*params.invL+0.5f)*params.L); //MIC Algorithm
  //   r.x -= shift.x;
  //   r.y -= shift.y;
  //   r.z -= shift.z;    
  // }

  inline __device__ void apply_pbc(float4 &r){
    
    r -= make_float4(floorf(make_float3(r)*params.invL+0.5f)*params.L); //MIC Algorithm
  }
  inline __device__ void apply_pbc(float3 &r){
    
    r -= floorf(r*params.invL+0.5f)*params.L; //MIC Algorithm
  }
  
  
  //Get the 3D cell p is in, just pos in [0,L] divided by ncells(vector) .INT DANGER.
  template<typename vecType>
  inline __device__ int3 getCell(vecType p){
    float3 r = make_float3(p);
    apply_pbc(r); //Reduce to MIC
    // return  int( (p+0.5L)/cellSize )
    int3 cell = make_int3((r+0.5f*params.L)*params.invCellSize);
    //Anti-Traquinazo guard, you need to explicitly handle the case where a particle
    // is exactly at the box limit, AKA -L/2. This is due to the precision loss when
    // casting int from floats, which gives non-correct results very near the cell borders.
    // This is completly neglegible in all cases, except with the cell 0, that goes to the cell
    // cellDim, wich is catastrophic.
    //Doing the previous operation in double precision (by changing 0.5f to 0.5) also works, but it is a bit of a hack and the performance appears to be the same.
    if(cell.x==params.cellDim.x) cell.x = 0;
    if(cell.y==params.cellDim.y) cell.y = 0;
    if(cell.z==params.cellDim.z) cell.z = 0;
    return cell;
  }

  //Apply pbc to a cell coordinates
  inline __device__ void pbc_cell(int3 &cell){
    if(cell.x==-1) cell.x = params.cellDim.x-1;
    else if(cell.x==params.cellDim.x) cell.x = 0;

    if(cell.y==-1) cell.y = params.cellDim.y-1;
    else if(cell.y==params.cellDim.y) cell.y = 0;

    if(cell.z==-1) cell.z = params.cellDim.z-1;
    else if(cell.z==params.cellDim.z) cell.z = 0;
  }
  
  //Get linear index of a 3D cell, from 0 to ncells-1
  inline __device__ uint getCellIndex(int3 gridPos){
    return dot(gridPos, params.gridPos2CellIndex);
  }


  /*Interleave a 10 bit number in 32 bits, fill one bit and leave the other 2 as zeros.*/
  inline __device__ uint encodeMorton(uint i){
  
    uint x = i;
  
    x &= 0x3ff;
    x = (x | x << 16) & 0x30000ff;
    x = (x | x << 8) & 0x300f00f;
    x = (x | x << 4) & 0x30c30c3;
    x = (x | x << 2) & 0x9249249;
    return x;
  }
  /*Fuse three 10 bit numbers in 32 bits, producing a Z order Morton hash*/
  inline __device__ uint mortonHash(const int3 &cell){

    return encodeMorton(cell.x) | (encodeMorton(cell.y) << 1) | (encodeMorton(cell.z) << 2);
  }
  
  /****************************************************************************************/


  /*Assign a hash to each particle from its cell index*/
  __global__ void calcHashD(uint __restrict__ *particleHash, uint __restrict__ *particleIndex, 
			    const float4 __restrict__ *pos, uint N){
    uint index = blockIdx.x*blockDim.x + threadIdx.x;  
    if(index>=N) return;
    float4 p = pos[index];

    int3 cell = getCell(p);
    /*The particleIndex array will be sorted by the hashes, any order will work*/
    uint hash = mortonHash(cell);//getCellIndex(cell);
    /*Before ordering by hash the index in the array is the index itself*/
    particleIndex[index] = index;
    particleHash[index]  = hash;
  }  
  //CPU kernel caller
  void calcHash(float4 *pos, uint *particleHash, uint *particleIndex, uint N){
    calcHashD<<<GPU_Nblocks, GPU_Nthreads>>>(particleHash, particleIndex, pos, N);
    //cudaCheckErrors("Calc hash");					   
  }



  /*Sort the particleIndex list by hash*/
  // this allows to access the neighbour list of each particle in a more coalesced manner
  //Each time this is called, the pointers particleHash and particleIndex are swapped
  void sortCellHash(uint *&particleHash, uint *&particleIndex, uint N){
    //This uses the CUB API to perform a radix sort
    //CUB orders by key an array pair and copies them onto another pair
    //This function stores an internal key/value pair and switches the arrays each time
    static bool init = false;
    static void *d_temp_storage = NULL;
    static size_t temp_storage_bytes = 0; //Additional storage needed by cub
    static uint *particleHash_alt = NULL, *particleIndex_alt = NULL; //Additional key/value pair

    static cub::DoubleBuffer<uint> d_keys;
    static cub::DoubleBuffer<uint> d_values;
    /**Initialize CUB at first call**/
    if(!init){
      /*Allocate temporal value/key pair*/
      gpuErrchk(cudaMalloc(&particleHash_alt, N*sizeof(uint)));
      gpuErrchk(cudaMalloc(&particleIndex_alt, N*sizeof(uint)));
    
      /*Create this CUB like data structure*/
      d_keys = cub::DoubleBuffer<uint>(particleHash, particleHash_alt);    
      d_values = cub::DoubleBuffer<uint>(particleIndex, particleIndex_alt);
      /*On first call, this function only computes the size of the required temporal storage*/
      cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
     				      d_keys, 
     				      d_values, N);
      /*Allocate temporary storage*/
      gpuErrchk(cudaMalloc(&d_temp_storage, temp_storage_bytes));
      init = true;
    }

    /**Perform the Radix sort on the index/cell pair**/
    cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
     				    d_keys, 
     				    d_values, N); 
    /**Switch the references**/
    particleHash     = d_keys.Current();
    particleIndex = d_values.Current();

    particleHash_alt     = d_keys.Alternate();
    particleIndex_alt = d_values.Alternate();

    // thrust::stable_sort_by_key(device_ptr<uint>(particleHash),
    // 			device_ptr<uint>(particleHash+N),
    // 			device_ptr<uint>(particleIndex));
    //cudaCheckErrors("Sort hash");					   
  }

  /*This kernel fills sortPos with the positions in pos, acording to the indices in particleIndex*/
  __global__ void reorderPosD(float4 __restrict__ *sortPos,
			      const uint* __restrict__ particleIndex, uint N){
    uint i = blockIdx.x*blockDim.x + threadIdx.x;
    if(i>=N) return;

    uint sort_index = particleIndex[i]; //Coalesced

    sortPos[i] = tex1Dfetch<float4>(params.texPos, sort_index);
  }
  /*Same as above, but reordering vel aswell*/
  __global__ void reorderPosVelD(float4 *sortPos,
				 float4* sortVel,
				 float3 * vel,
				 uint* particleIndex, uint N){
    uint i = blockIdx.x*blockDim.x + threadIdx.x;
    if(i>=N) return;

    uint sort_index = particleIndex[i]; //Coalesced

    sortPos[i] = tex1Dfetch<float4>(params.texPos, sort_index);
    sortVel[i] = make_float4(vel[sort_index]);
  }
  
  /*Fill CellStart and CellEnd*/
  __global__ void fillCellListD(const float4 __restrict__ *sortPos,
				uint *cellStart, uint *cellEnd,
				uint N){
    /*A thread per particle*/
    uint i = blockIdx.x*blockDim.x + threadIdx.x;

    if(i<N){//If my particle is in range
      uint icell, icell2;
      /*Get my icell*/
      icell = getCellIndex(getCell(sortPos[i]));
      /*Get the previous part.'s icell*/
      if(i>0){ /*Shared memory target VVV*/
	icell2 = getCellIndex(getCell(sortPos[i-1]));
      }
      else
	icell2 = 0;
      //If my particle is the first or is in a different cell than the previous
      //my i is the start of a cell
      if(i ==0 || icell != icell2){
	//Then my particle is the first of my cell
	cellStart[icell] = i;
	//If my i is the start of a cell, it is also the end of the previous
	//Except if I am the first one
	if(i>0)
	  cellEnd[icell2] = i;
      }
      //If I am the last particle my cell ends 
      if(i == N-1) cellEnd[icell] = N;

      
    }

  }

  /*Reupload the parameters to constant memory*/
  void updateParams(Params m_params){
    gpuErrchk(cudaMemcpyToSymbol(params, &m_params, sizeof(Params)));
  }

  /*Create the Cell List from scratch in the GPU*/
  void makeCellList(float4 *pos, float4 *sortPos,
		    uint *&particleIndex, uint *&particleHash,
		    uint *cellStart, uint *cellEnd,
		    uint N, uint ncells){
    
    cudaMemset(cellStart, 0xffffffff, ncells*sizeof(uint));

    calcHashD<<<GPU_Nblocks, GPU_Nthreads>>>(particleHash, particleIndex, pos, N);
    
    sortCellHash(particleHash, particleIndex, N);

    reorderPosD<<<GPU_Nblocks, GPU_Nthreads>>>(sortPos, particleIndex, N);

    /*This fills cellStart and cellEnd*/
    fillCellListD<<<GPU_Nblocks, GPU_Nthreads>>>(sortPos, cellStart, cellEnd, N);
    
  }
  
  void makeCellListDPD(float4 *pos, float3* vel,  float4 *sortPos, float4 *sortVel,
		       uint *&particleIndex, uint *&particleHash,
		       uint *cellStart, uint *cellEnd,
		       uint N, uint ncells){
    
    cudaMemset(cellStart, 0xffffffff, ncells*sizeof(uint));

    calcHashD<<<GPU_Nblocks, GPU_Nthreads>>>(particleHash, particleIndex, pos, N);
    
    sortCellHash(particleHash, particleIndex, N);

    reorderPosVelD<<<GPU_Nblocks, GPU_Nthreads>>>(sortPos, sortVel, vel, particleIndex, N);

    fillCellListD<<<GPU_Nblocks, GPU_Nthreads>>>(sortPos, cellStart, cellEnd, N);
    
  }

  //TODO The naming and explanation of this function
  /*Transverses all the neighbour particles of each particle using the cell list and computes a quantity as implemented by Transversable. Each thread goes through all the neighbours of a certain particle(s)(index), transversing its 27 neighbour cells*/
  /*Computes a quantity determined by Transversable, which is a class that must implement the following methods:
    zero() -> returns the initial value of the quantity, in whatever type
    compute(float4 r1, float4 r2) -> compute the quantity depending of the pair positions/types
    set(uint index, TYPE quantity) -> sum the total quantity on particle index to global memory

    This quantity can be i.e a float4 and compute the force
                         or a float and compute the energy...
    */
  template<class Transversable>
  __global__ void transverseListD(Transversable T, 
				  const uint* __restrict__ particleIndex,
				  uint N){
    uint ii =  blockIdx.x*blockDim.x + threadIdx.x;

    //Grid-stride loop
    for(int index = ii; index<N; index += blockDim.x * gridDim.x){
      /*Compute force acting on particle particleIndex[index], index in the new order*/
      //float4 pos = tex1Dfetch<float4>(texSortPos, index);

      auto pinfoi = T.getInfo(index);//tex1Dfetch<float4>(texSortPos, index);

      /*Initial value of the quantity*/
      auto quantity = T.zero();
      
      int3 celli = getCell(pinfoi.pos);

      int x,y,z;
      int3 cellj;
      //      float4 posj;
      /**Go through all neighbour cells**/
      //For some reason unroll doesnt help here
      for(z=-1; z<=1; z++)
	for(y=-1; y<=1; y++)
	  for(x=-1; x<=1; x++){
	    cellj = celli + make_int3(x,y,z);
	    pbc_cell(cellj);

	    uint icell  = getCellIndex(cellj);
	    /*Index of the first particle in the cell's list*/ 
	    uint firstParticle = tex1Dfetch<uint>(params.texCellStart, icell);
	    if(firstParticle ==0xffFFffFF) continue;
	    /*Index of the last particle in the cell's list*/
	    uint lastParticle = lastParticle=tex1Dfetch<uint>(params.texCellEnd, icell);
	    /*Because the list is ordered, all the particle indices in the cell are coalescent!*/
	    /*If there are no particles in the cell, firstParticle=0xffffffff, the loop is not computed*/
	    /*The fetch of lastParticle eitherway reduces branch divergency and is actually faster than checking firstParticle before fetching*/	    
	    for(uint j=firstParticle; j<lastParticle; j++){
	      /*Retrieve j pos*/
	      //posj = tex1Dfetch<float4>(texSortPos, j);
	      auto pinfoj = T.getInfo(j);
	      /*Add force, i==j is handled in forceij */
	      //quantity += T.compute(pos, posj);
	      quantity += T.compute(pinfoi, pinfoj);

	    }
	    
	  }
      /*Write quantity with the original order*/
      uint pi = particleIndex[index]; //Coalesced
      T.set(pi, quantity);
    }
    
  }
    


  /***************************************FORCE*****************************/
  
  //tags: force compute force function forceij
  /*A helper class that holds the address of the force array (in device) and
    computes the force between two particles*/
  /*It also updates the global array newForce with the total force acting on particle pi*/
  /*This helper class can be passed as an argument to transverseListD, which will apply the compute
    function to every neighbouring particle pair in the system*/
  /*In order to compute any other quantity create a class like this, implementing the same functions
    but with any desired type, instead of float4 as in this case*/
  class forceTransversable{
  public:
    struct ParticleInfo{
      float4 pos;
    };

    /*I need the device pointer to force*/
    forceTransversable(float4 *newForce):newForce(newForce){
    };
    /*Compute the force between two positions*/
    inline __device__ float4 compute(const ParticleInfo &R1,const ParticleInfo &R2){

      float3 r12 = make_float3(R2.pos-R1.pos);
      apply_pbc(r12);
      /*Squared distance*/
      float r2 = dot(r12,r12);
      /*Squared distance between 0 and 1*/
      float r2c = r2*params.invrc2;
      /*Both cases handled in texForce*/
      /*Check if i==j. This way reduces warp divergence and its faster than checking i==j outside*/
      //if(r2c==0.0f) return make_float4(0.0f);  
      /*Beyond rcut..*/
      //else if(r2c>=1.0f) return make_float4(0.0f);
      /*Get the force from the texture*/
      //float fmod = tex1D(texForce, r2c);
      float fmod = tex1D<float>(params.texForce, r2c);
      // float invr2 = 1.0f/r2;
      // float invr6 = invr2*invr2*invr2;
      // float invr8 = invr6*invr2;
      // //float E =  2.0f*invr6*(invr6-1.0f);
      // float fmod = -48.0f*invr8*invr6+24.0f*invr8;
      return  make_float4(fmod*r12);
    }
    /*Update the force acting on particle pi, pi is in the normal order*/
    inline __device__ void set(uint pi, const float4 &totalForce){
      newForce[pi] += totalForce;
    }
    /*Initial value of the force, this is a trick to allow the template in transverseList
      to guess the type of my quantity, a float4 in this case. Just set it to the 0 value of 
    the type of your quantity (0.0f for a float i.e)*/
    inline __device__ float4 zero(){
      return make_float4(0.0f);
    }

    inline __device__ ParticleInfo getInfo(uint index){      
      return {tex1Dfetch<float4>(params.texSortPos, index)};
    }
  private:
    float4* newForce;
  };


  
  //CPU kernel caller
  void computePairForce(float4 *sortPos, float4 *force,
			uint *cellStart, uint *cellEnd,
			uint *particleIndex, 
			uint N){
    /*An instance of the class that holds the function that computes the force*/
    forceTransversable ft(force); //It needs the addres of the force in device memory
    /*Transverse the neighbour list for each particle, using ft to compute the force in each pair*/
    transverseListD<<<GPU_Nblocks, GPU_Nthreads>>>(ft,
						   particleIndex,
						   N);
    //cudaCheckErrors("computeForce");
  }

  /****************************ENERGY***************************************/

  /*This class is analogous to forceTransversable, see for reference*/
  //tags: energy compute energyij
  class energyTransversable{
  public:
    struct ParticleInfo{
      float4 pos;
    };

    energyTransversable(float *Energy):Energy(Energy){ };
    /*Returns the energy between two positions*/
    inline __device__ float compute(const ParticleInfo &R1,const ParticleInfo &R2){
      float3 r12 = make_float3(R2.pos-R1.pos);

      apply_pbc(r12);

      float r2 = dot(r12,r12);
      /*Squared distance between 0 and 1*/
      float r2c = r2*params.invrc2;
      /*Check if i==j. This way reduces warp divergence and its faster than checking i==j outside*/
      //if(r2c==0.0f) return 0.0f;  //Both cases handled in texForce
      /*Beyond rcut..*/
      //else if(r2c>=1.0f) return 0.0f;
      /*Get the force from the texture*/
      //float fmod = tex1D(texForce, r2c);
      //float invr2 = 1.0f/r2;
      //float invr6 = invr2*invr2*invr2;
      //TODO take from a texture*/
      //float E =  2.0f*invr6*(invr6-1.0f);
      float E = tex1D<float>(params.texEnergy, r2c);
      
      return E;
    }
    inline __device__ void set(uint pi, float energy){
      Energy[pi] = energy;
    }
    inline __device__ float zero(){
      return 0.0f;
    }
    inline __device__ ParticleInfo getInfo(uint index){
      return {tex1Dfetch<float4>(params.texSortPos, index)};
    }

  private:
    float *Energy;
  };

  
  float computePairEnergy(float4 *sortPos, float *energy,
			  uint *cellStart, uint *cellEnd,
			  uint *particleIndex, 
			  uint N){

    /*Analogous to computeForce, see for reference*/
    energyTransversable et(energy);
    transverseListD<<<GPU_Nblocks, GPU_Nthreads>>>(et,
						   particleIndex,
						   N);
    device_ptr<float> d_e(energy);
    float sum;
    sum = thrust::reduce(d_e, d_e+N, 0.0f);
    //Returns energy per particle*/
    return (sum/(float)N);
    //cudaCheckErrors("computeForce");
  }


  /****************************VIRIAL***************************************/
  /*Analogous to forceTransversable, see for reference*/
  //tags: virial compute virialij
  class virialTransversable{
  public:
    struct ParticleInfo{
      float4 pos;
    };

    virialTransversable(float *virial):Virial(virial){ };
    inline __device__ float compute(const ParticleInfo &R1,const ParticleInfo &R2){
      float3 r12 = make_float3(R2.pos-R1.pos);
      apply_pbc(r12);

      float r2 = dot(r12,r12);
      /*Squared distance between 0 and 1*/
      float r2c = r2*params.invrc2;
      //if(r2c==0.0f) return 0.0f; //No need to check i==j, tex1D(texForce, 0.0) = 0.0
      /*Beyond rcut..*/
      //if(r2c>=1.0f) return 0.0f; //Also 0 in texForce
      /*Get the force from the texture*/
      float fmod = tex1D<float>(params.texForce, r2c);
      // P = rhoKT + (1/2dV)sum_ij( Fij·rij ) //Compute only the Fij·rij, the rest is done outside
      return dot(fmod*r12,r12);
    }
    inline __device__ void set(uint pi, float virial){
      Virial[pi] = virial;
    }
    inline __device__ float zero(){
      return 0.0f;
    }

    inline __device__ ParticleInfo getInfo(uint index){
      return {tex1Dfetch<float4>(params.texSortPos, index)};
    }

  private:
    float *Virial;
  };




  //CPU kernel caller
  float computePairVirial(float4 *sortPos, float *virial,
			  uint *cellStart, uint *cellEnd,
			  uint *particleIndex, 
			  uint N){

    virialTransversable ft(virial);
    transverseListD<<<GPU_Nblocks, GPU_Nthreads>>>(ft,
						   particleIndex,
						   N);
    device_ptr<float> d_vir(virial);
    float sum;
    // P = rhoKT + (1/2dV)sum_ij( Fij·rij ) This function returns (1/2)sum_ij( Fij·rij )
    sum = thrust::reduce(d_vir, d_vir+N, 0.0f);
    return (sum/2.0f);
    //cudaCheckErrors("computeForce");
  }



  /*******************************************DPD********************************************/

  /**********************FORCE********************/


  //Random number, the seed is used to recover a certain number in the random stream
  //TODO: This is a bit awkward, probably it will be best to manually generate the number
  inline __device__ float randGPU(const ullint &seed, curandState *rng){
    curand_init(seed, 0, 0, rng);
    return curand_normal(rng);
  }

  /*Very similar to forceTransversable, but now the force depends on the velocity aswell!*/
  //tags: forceijDPD forceDPDij
  class forceDPDTransversable{
  public:
    /*A struct with all the information the force compute needs*/
    struct ParticleInfo{      
      float4 pos; //Pos has to be the first element always!
      uint pi;
    };

    /*I need the device pointer to force*/
    forceDPDTransversable(float4 *newForce):newForce(newForce){};
    
    /*Compute the force between two particles*/
    /*As in all this file, the __device__ function has access to all the parameters in params and,
      in this case, paramsDPD. See PairForcesGPU.cuh*/
    
    inline __device__ float4 compute(const ParticleInfo &R1,const ParticleInfo &R2){
      float3 r12 = make_float3(R1.pos-R2.pos);  
      apply_pbc(r12);

      float r2 = dot(r12,r12);
      /*Squared distance between 0 and 1*/
      float r2c = r2*params.invrc2;
    
      float w = 0.0f; //The intensity of the DPD thermostat 
      float rinv = 0.0f;
      if(r2c<1.0f){
	if(r2c==0.0f) return make_float4(0.0f);
	//w = r-rc -> linear
	rinv = rsqrt(r2);
	w = rinv-params.invrc;
      }
      else return make_float4(0.0f);
      


      uint i0 = R1.pi;
      uint j0 = R2.pi;

      float4 V1=tex1Dfetch<float4>(paramsDPD.texSortVel, i0);
      float4 V2=tex1Dfetch<float4>(paramsDPD.texSortVel, j0);
      
      float3 v12 = make_float3(V1-V2);      
      /*Prepare the seed for the RNG, it must be the same seed
	for pair ij and ji!*/
      if(i0>j0)
	swap(i0,j0);
      
      curandState rng;
      float randij = randGPU(i0+(ullint)params.N*j0 + seedGPU, &rng);

      //fmod = paramsDPD.A*w; //Soft force
  
      float fmod = -tex1D<float>(params.texForce, r2c); //Conservative force
      fmod -= paramsDPD.gamma*w*w*dot(r12,v12); //Damping
      fmod += paramsDPD.noiseAmp*randij*w; //Random force
      return make_float4(fmod*r12);
    }
    /*Update the force acting on particle pi, pi is in the normal order*/
    inline __device__ void set(uint pi, const float4 &totalForce){
      newForce[pi] += totalForce;
    }
    /*Initial value of the force, this is a trick to allow the template in transverseList
      to guess the type of my quantity, a float4 in this case. Just set it to the 0 value of 
    the type of your quantity (0.0f for a float i.e)*/
    inline __device__ float4 zero(){
      return make_float4(0.0f);
    }

    inline __device__ ParticleInfo getInfo(uint index){
      // ParticleInfo kk;
      // kk.pos = tex1Dfetch<float4>(params.texSortPos, index);
      // kk.vel = tex1Dfetch<float4>(paramsDPD.texSortVel, index);
      // kk.pi = index;
      
      return {tex1Dfetch<float4>(params.texSortPos, index),
	      index};
    }
  private:
    float4* newForce;
  };
  
  //CPU kernel caller
  void computePairForceDPD(float4 *force,
			   uint *particleIndex,
			   uint N, ullint seed){
    //TODO: This should not be needed somehow
    gpuErrchk(cudaMemcpyToSymbol(seedGPU, &seed, sizeof(seedGPU)));
    
    forceDPDTransversable ft(force);
    
    transverseListD<<<GPU_Nblocks, GPU_Nthreads>>>(ft, 
     						   particleIndex,
     						   N);
    //cudaCheckErrors("computeForce");
  }

}
  
