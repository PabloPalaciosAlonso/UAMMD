/*Raul P.Pelaez 2017. Grid utility.

Given a certain box and a number of cells, subdivides the box in that number of cells and allows to:
    -Obtain the cell in which a position is located (taking into account PBC according to Box).
    -Compute the linear 1D index of a certain cell in the grid
    -Apply PBC to a cell (as in cellx=-1 goes to cellx=cellDim.x-1)

 */
#ifndef GRID_CUH
#define GRID_CUH

#include "Box.cuh"
#include "vector.cuh"
#include <algorithm>
namespace uammd{
  
  struct Grid{
    /*A magic vector that transforms cell coordinates to 1D index when dotted*/
    /*Simply: 1, ncellsx, ncellsx*ncellsy*/
    int3 gridPos2CellIndex;
    
    int3 cellDim; //ncells in each size
    real3 cellSize;
    real3 invCellSize; /*The inverse of the cell size in each direction*/
    Box box;
    Grid(): Grid(Box(), make_int3(0,0,0)){}

    Grid(Box box, real3 minCellSize):
      Grid(box, make_int3(box.boxSize/minCellSize)){}    
    Grid(Box box, real minCellSize):
      Grid(box, make_real3(minCellSize)){}
    
    Grid(Box box, int3 in_cellDim):
      box(box),
      cellDim(in_cellDim){

      if(cellDim.z == 0) cellDim.z = 1;
      cellSize = box.boxSize/make_real3(cellDim);
      invCellSize = 1.0/cellSize;
      if(box.boxSize.z == real(0.0)) invCellSize.z = 0;
	
      gridPos2CellIndex = make_int3( 1,
				     cellDim.x,
				     cellDim.x*cellDim.y);	
    }
    template<class VecType>
    inline __host__ __device__ int3 getCell(const VecType &r) const{
	// return  int( (p+0.5L)/cellSize )
	int3 cell = make_int3((box.apply_pbc(make_real3(r)) + real(0.5)*box.boxSize)*invCellSize);
	//Anti-Traquinazo guard, you need to explicitly handle the case where a particle
	// is exactly at the box limit, AKA -L/2. This is due to the precision loss when
	// casting int from floats, which gives non-correct results very near the cell borders.
	// This is completly neglegible in all cases, except with the cell 0, that goes to the cell
	// cellDim, which is catastrophic.
	//Doing the previous operation in double precision (by changing 0.5f to 0.5) also works, but it is a bit of a hack and the performance appears to be the same as this.
	//TODO: Maybe this can be skipped if the code is in double precision mode
	if(cell.x==cellDim.x) cell.x = 0;
	if(cell.y==cellDim.y) cell.y = 0;
	if(cell.z==cellDim.z) cell.z = 0;
	return cell;
    }

    inline __host__ __device__ int getCellIndex(const int3 &cell) const{
	return dot(cell, gridPos2CellIndex);
    }

    inline __host__  __device__ int3 pbc_cell(const int3 &cell) const{
	int3 cellPBC;
	cellPBC.x = pbc_cell_coord<0>(cell.x);
	cellPBC.y = pbc_cell_coord<1>(cell.y);
	cellPBC.z = pbc_cell_coord<2>(cell.z);
	return cellPBC;
    }

    template<int coordinate>
    inline __host__  __device__ int pbc_cell_coord(int cell) const{
	int ncells = 0;
	if(coordinate == 0){
	  ncells = cellDim.x;
	}
	if(coordinate == 1){
	  ncells = cellDim.y;
	}

	if(coordinate == 2){
	  ncells = cellDim.z;
	}

	if(cell <= -1) cell += ncells;
	else if(cell >= ncells) cell -= ncells;
	return cell;
    }

    inline __host__ __device__ int getNumberCells() const{ return cellDim.x*cellDim.y*cellDim.z;}
  };

  //Looks for the closest (equal or greater) number of nodes of the form 2^a*3^b*5^c*7^d*11^e
  int3 nextFFTWiseSize3D(int3 size){
    int* cdim = &size.x;

    int max_dim = std::max({size.x, size.y, size.z});
	
    int n= 14;
    int n5 = 6; //number higher than this are not reasonable...
    int n7 = 5;
    int n11 = 4;
    auto powint = [](uint64_t base, uint64_t exp){if(exp==0) return uint64_t(1); uint64_t res = base; fori(0,exp-1) res*=base; return res;};
    
    std::vector<uint64_t> tmp(n*n*n5*n7*n11, 0);
    do{
      tmp.resize(n*n*n5*n7*n11, 0);
      fori(0,n)forj(0,n)
	for(int k=0; k<n5;k++)for(int k7=0; k7<n7; k7++)for(int k11=0; k11<n11; k11++){
	      if(k11>4 or k7>5 or k>6) continue;
		
	      uint64_t id = i+n*j+n*n*k+n*n*n5*k7+n*n*n5*n7*k11;
	      tmp[id] = 0;
	      //Current fft wise size
	      uint64_t number = uint64_t(powint(2,i))*powint(3,j)*powint(5,k)*powint(7, k7);
	      //This is to prevent overflow
	      if(!(i==n-1 and j==n-1 and k==n5-1 and k7==n7-1 and k11<n11-1))
		number *= powint(11, k11);
	      else
		continue;
	      //The fastest FFTs always have at least a factor of 2
	      if(i==0) continue;
	      //I have seen empirically that factors 11 and 7 only works well with at least a factor 2 involved
	      if((k11>0 && (i==0))) continue;
	      tmp[id] = number;
	    }
      n++;
      /*Sort this array in ascending order*/
      std::sort(tmp.begin(), tmp.end());      
    }while(tmp.back()<max_dim); /*if n is not enough, include more*/

    //I have empirically seen that these sizes produce slower FFTs than they should in several platforms
    constexpr uint64_t forbiddenSizes [] = {28, 98, 150, 154, 162, 196, 242};
    /*Now look for the nearest value in tmp that is greater than each cell dimension and it is not forbidden*/
    forj(0,3){
      fori(0, tmp.size()){	    
	if(tmp[i]<uint64_t(cdim[j])) continue;
	for(int k =0;k<sizeof(forbiddenSizes)/sizeof(uint64_t); k++) if(tmp[i] == forbiddenSizes[k]) continue;
	int set = int(tmp[i]);
	if(tmp[i]>=powint(2,31)) set = -1;
	cdim[j] = set;
	break;
      }	   	  
    }
    return size;
  }
}

#endif