/*Raul P. Pelaez 2017. PairForces definition.

  PairForces Module is an interactor that computes short range forces.
  Computes the interaction between neighbour particles (pairs of particles closer tan rcut).

  For that, it uses a NeighbourList and computes the force given by Potential for each pair of particles. It sums the force for all neighbours of every particle.

  See https://github.com/RaulPPelaez/UAMMD/wiki/Pair-Forces   for more info.
*/

#ifndef PAIRFORCES_H
#define PAIRFORCES_H

#include"Interactor/Interactor.cuh"
#include"Interactor/NeighbourList/CellList.cuh"
#include"Interactor/NBody.cuh"
#include"third_party/type_names.h"

namespace uammd{
  /*This makes the class valid for any NeighbourList*/
  template<class Potential, class NeighbourList = CellList>
  class PairForces: public Interactor, public ParameterUpdatableDelegate<Potential>{
  public:
    struct Parameters{
      Box box;
      //You can provide a neighbour list from outside that will be used by PairForces
      shared_ptr<NeighbourList> nl=shared_ptr<NeighbourList>(nullptr);
    };
    PairForces(shared_ptr<ParticleData> pd,
	       shared_ptr<ParticleGroup> pg,
	       shared_ptr<System> sys,
	       Parameters par,
	       shared_ptr<Potential> pot = std::make_shared<Potential>());
    PairForces(shared_ptr<ParticleData> pd,
	       shared_ptr<System> sys,
	       Parameters par,
	       shared_ptr<Potential> pot = std::make_shared<Potential>()):
      PairForces(pd, std::make_shared<ParticleGroup>(pd, sys, "All"), sys, par, pot){
    }

    ~PairForces(){
      sys->log<System::DEBUG>("[PairForces] Destroyed.");
    }

    void updateBox(Box box) override{
      sys->log<System::DEBUG3>("[PairForces] Box updated.");
      this->box = box;
      //In case the potential wants to handle updateBox
      ParameterUpdatableDelegate<Potential>::updateBox(box);
    }
    void sumForce(cudaStream_t st) override;
    real sumEnergy() override;
    //real sumVirial() override{ return 0;}
    //sumForce and sumEnergy are defined through this function
    template<class Transverser>
    void sumTransverser(Transverser &tr, cudaStream_t st);

    void print_info(){
      sys->log<System::MESSAGE>("[PairForces] Using: %s Neighbour List.", type_name<NeighbourList>());
      //nl.print();
      sys->log<System::MESSAGE>("[PairForces] Using: %s potential.", type_name<Potential>());
    }


  private:
    shared_ptr<NeighbourList> nl;
    shared_ptr<NBody> nb;
    shared_ptr<Potential> pot;
    Box box;
    real rcut;

  };


}

#include"PairForces.cu"

#endif

