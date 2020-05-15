/* Raul P. Pelaez 2020. Neighbour Container for Cell List.

   Given a CellListData NeighbourContainer can provide for each particle a list
of its neighbours.

USAGE:

auto cl = make_shared<CellList>(pd, pg, sys);
cl->update(...);
auto nl = cl->getCellList();
auto nc = CellList_ns::NeighbourContainer(nl);

--- inside a CUDA kernel ---
//Set the container to provide the list of neighbours of particle i. i is in internal CellList indexing
ni.set(i);
//Get the group index (global index if group is All) of a particle:
int groupIndex = ni.getGroupIndexes()[i];
//There is no easy way of doing the opposite operation. That is getting the
//internal index of a particle given its global or group index. This is
//intentional and by design.
//Get the position of a particle given its internal index
const real3 pos_i = make_real3(ni.getSortedPositions()[i]);
//Iterator to the first neighbour of particle i
auto it = ni.begin();
//Note that ni.end() is not a pointer to the last neighbour, it just represents "no more neighbours" and should not be dereferenced
//Loop through neighbours
while(it){ //it will cast to false when there are no more neighbours
  auto neigh = *it++; //The iterator can only be advanced and dereferenced
  //int j = neigh.getGroupIndex();
  const real3 pj = make_real3(neigh.getPos());
  //const real3 rij = box.apply_pbc(pj-pi);
  //const real r2 = dot(rij, rij);
  //if(r2>0) f += lj(r2)*rij;
}
//Write something to a group sorted array (global sorted if group contains all particles)
force[groupIndex] += ...


--- ---


 */
#ifndef CELLLIST_NEIGHBOURCONTAINER_CUH
#define CELLLIST_NEIGHBOURCONTAINER_CUH
#include "CellListBase.cuh"
namespace uammd{

  namespace CellList_ns{
    class NeighbourContainer; //forward declaration for befriending

    class NeighbourIterator; //forward declaration for befriending

    struct Neighbour{
      __device__ Neighbour(const Neighbour &other):
	internal_i(other.internal_i),
	groupIndex(other.groupIndex),
	sortPos(other.sortPos){}

      __device__ int getInternalIndex(){return internal_i;}

      __device__ int getGroupIndex(){return groupIndex[internal_i];}

      __device__ real4 getPos(){return cub::ThreadLoad<cub::LOAD_LDG>(sortPos+internal_i);}

    private:
      int internal_i;
      const int* groupIndex;
      const real4* sortPos;
      friend class NeighbourIterator;

      Neighbour() = delete;

      __device__ Neighbour(int i, const int* gi, const real4* sp):
	internal_i(i), groupIndex(gi), sortPos(sp){}
    };

    //This forward iterator must be constructed by NeighbourContainer,
    class NeighbourIterator:
      public thrust::iterator_adaptor<
      NeighbourIterator,
      int,   Neighbour,
      thrust::any_system_tag,
      thrust::forward_device_iterator_tag,
      Neighbour, int
      >{
      friend class thrust::iterator_core_access;
      friend class NeighbourContainer;
      int particleIndex;
      int currentNeighbourIndex;
      CellListBase::CellListData nl;
      int currentCell;
      int3 celli;
      int lastParticleInCell;
      static constexpr int noMoreNeighbours = -1;
      //Take currentNeighboutIndex to the start of the next cell and return true, if no more cells remain then return false
      __device__ bool nextcell(){
	const bool is2D = nl.grid.cellDim.z<=1;
	if(currentCell >= (is2D?9:27)) return false;
	bool isCurrentCellEmpty = true;
	do{
	  int3 cellj = celli;
	  cellj.x += currentCell%3-1;
	  cellj.y += (currentCell/3)%3-1;
	  cellj.z =  is2D?0:(celli.z+currentCell/9-1);
	  cellj = nl.grid.pbc_cell(cellj);
	  const bool isPeriodicCellInNonPeriodicBox = (!nl.grid.box.isPeriodicX() and abs(cellj.x-celli.x)>1) or
	    (!nl.grid.box.isPeriodicY() and abs(cellj.y-celli.y)>1) or
	    (!nl.grid.box.isPeriodicZ() and abs(cellj.z-celli.z)>1);
	  if(!isPeriodicCellInNonPeriodicBox){
	    const int icellj = nl.grid.getCellIndex(cellj);
	    const uint cs = nl.cellStart[icellj];
	    isCurrentCellEmpty = cs<nl.VALID_CELL;
	    lastParticleInCell = isCurrentCellEmpty?noMoreNeighbours:nl.cellEnd[icellj];
	    currentNeighbourIndex = isCurrentCellEmpty?noMoreNeighbours:int(cs-nl.VALID_CELL);
	  }
	  currentCell++;
	  if(currentCell >= (is2D?9:27)) return !isCurrentCellEmpty;
	}while(isCurrentCellEmpty);
	return true;
      }

      //Go to the next neighbour
      __device__  void increment(){
	if(++currentNeighbourIndex == lastParticleInCell){
	  currentNeighbourIndex = nextcell()?currentNeighbourIndex:noMoreNeighbours;
	}
      }

      __device__ Neighbour dereference() const{
	return Neighbour(currentNeighbourIndex, nl.groupIndex, nl.sortPos);
      }

      //Can only be advanced
      __device__  void decrement() = delete;
      //Just a forward iterator
      __device__  Neighbour operator[](int i) = delete;

      __device__  bool equal(NeighbourIterator const& other) const{
	return other.particleIndex == particleIndex and other.currentNeighbourIndex==currentNeighbourIndex;
      }

      __device__ NeighbourIterator(int i, CellListBase::CellListData nl, bool begin):
	particleIndex(i),
	currentNeighbourIndex(noMoreNeighbours-1),
	nl(nl),
	currentCell(0),
	lastParticleInCell(noMoreNeighbours){
	if(begin){
	  this->celli = nl.grid.getCell(make_real3(cub::ThreadLoad<cub::LOAD_LDG>(nl.sortPos+i)));
	  increment();
	}
	else{
	  currentNeighbourIndex = noMoreNeighbours;
	}
      }

    public:
      __device__  operator bool(){ return currentNeighbourIndex != noMoreNeighbours;}

    };

    struct NeighbourContainer{
      int my_i = -1;
      CellListBase::CellListData nl;
      NeighbourContainer(CellListBase::CellListData nl): nl(nl){}

      __device__ void set(int i){this->my_i = i;}

      __device__ NeighbourIterator begin(){return NeighbourIterator(my_i, nl, true);}

      __device__ NeighbourIterator end(){  return NeighbourIterator(my_i, nl, false);}

      __host__ __device__ const real4* getSortedPositions(){return nl.sortPos;}

      __host__ __device__ const int* getGroupIndexes(){return nl.groupIndex;}

    };
  }
}
#endif
