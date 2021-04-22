/*Raul P. Pelaez 2021.
A generic particle simulation code with almost everything in UAMMD.

Once compiled a lot of stuff can be chosen via a file called "data.main", which will be autogenerated with instructions if not present.
Start by compiling and running this file (which will also serve to check your CUDA environment), then read the data.main file that was generated.

If the wide variety of simulations that this code offers as-is does not meet your needs, not all hope is lost yet.
Please check customizations.cuh, were you can find several ways to extend and adapt this code.

As already mentioned, the data.main covers the kind of simulation this code can perform, but for completeness it is repeated here:
Using this code, you can fabricate via data.main a simulation using the following components:

Integration schemes (one must be chosen):
  -VerletNVT: Molecular dynamics with a thermostat.
  -VerletNVE: Molecular dynamics with constant energy (starting with a velocity distribution to match the given temperature).
  -BD: Brownian dynamics
  -BDHI: Brownian dynamics with Hydrodynamic interactions (via positively split Ewald Rotne-Prager-Yamakawa or Force Coupling Method)
  -SPH: Smooth Particle Hydrodynamics
  -DPD: Dissipative Particle Dynamics
  -FIB: Fluctuating Immersed Boundary
  -ICM: Inertial Coupling Method

Interaction modules:
  -Short range interaction: *LJ potential,
  -Bonded forces: Pairs of particles joined by harmonic bonds*.
  -Angular bonds: Groups of three particles joined by angular bonds.*
  -Torsional bonds: Groups of four particles joined by torsional bonds.*
  -External potential: An external force acting on each particle. Gravity + a wall in this example*
  -Electrostatics: Periodic electrostatics using an Ewald splitted spectral Poisson solver.

*The different potentials are taken from the file customizations.cuh accompanying this code, you can modify it to your needs.

Parameters:
Compile and run this code to generate a data.main with default parameters and explanations for each of them.
Once the data.main has been customized, run the program again.
The default data.main will run a LJ liquid MD simulation.

The code will output only the positions of the particles each print frame. With each frame separated by a line with "#".

Additional notes:
All simulations that can be performed via this file assume a periodic box. But you can put an arbitrarily large box size. Keep in mind though that some modules encode algorithms that are inherently periodic, like electrostatics and hydrodynamics. These inifinite ranged interactions make a big box not equivalent to an open system.

This code is also intended to be a skeleton code if you are trying to write some kind of UAMMD simulation code. For example, the function initialize (which starts the basic UAMMD structures and the particles) is probably going to appear everytime you write an UAMMD simulation.
If you need to create a BD integrator, you can copy paste the function createIntegratorBD.
Need to read parameters from a file? just copy paste the function readParameters.
And so on.

If you would like a more bottom up approach to UAMMD, you can surf the examples/basic folder, which will give you through UAMMD with an increasingly complex set of example codes.
Additionally, you can drop by the wiki: 
https://github.com/RaulPPelaez/UAMMD/wiki
Which has a lot of information. From basic functionality to descriptions and references for the algorithms implemented by each module.
*/

//Each UAMMD module has its own header, with uammd.cuh being the header with the basic structures
#include"uammd.cuh"
#include"Interactor/PairForces.cuh"
#include"Interactor/NeighbourList/VerletList.cuh"
#include"Interactor/Potential/Potential.cuh"
#include"Interactor/ExternalForces.cuh"
#include"Interactor/BondedForces.cuh"
#include"Interactor/AngularBondedForces.cuh"
#include"Interactor/TorsionalBondedForces.cuh"
#include"Integrator/BrownianDynamics.cuh"
#include "Integrator/BDHI/FIB.cuh"
#include "Integrator/Hydro/ICM.cuh"
#include "Integrator/BDHI/BDHI_EulerMaruyama.cuh"
#include"Integrator/BDHI/BDHI_PSE.cuh"
#include"Integrator/BDHI/BDHI_FCM.cuh"
#include"Interactor/SpectralEwaldPoisson.cuh"
#include "Integrator/Integrator.cuh"
#include "Integrator/VerletNVE.cuh"
#include "Integrator/VerletNVT.cuh"
#include "Interactor/Potential/DPD.cuh"
#include"Interactor/SPH.cuh"
#include"utils/InputFile.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <memory>
#include<random>
#include"utils/InitialConditions.cuh"
//The header with the Potentials/Parameters used in this code
#include"customizations.cuh"

using NeighbourList = VerletList;
using namespace uammd;

//Lets use this struct to pass around the basic uammd environment
struct UAMMD{
  std::shared_ptr<System> sys;
  std::shared_ptr<ParticleData> pd;
  Parameters par;
  std::shared_ptr<Integrator> integrator;
};

//Read a data.main file from "fileName" and return the contents as a Parameters struct.
Parameters readParameters(std::string fileName);

//This function reads numberParticles lines from a file and returns the contents in a std::vector
//This function will assume the format of the file is a list of objects of type "T".
//For example, in order to read 10 positions (formatted as X Y Z per particle) from file "pos.init":
//auto pos_in_file = readPropertyFromFile<real3>("pos.init", 10);
template<class T>
std::vector<T> readPropertyFromFile(std::string file, int numberParticles){
  std::ifstream in(file);
  if(not in.good()){
    System::log<System::CRITICAL>("Cannot open file %s", file);
  }
  std::istream_iterator<T> begin(in), end;
  std::vector<T> data_in_file(begin, end);
  if(not data_in_file.size()){
    System::log<System::CRITICAL>("File %s is empty", file);
  }
  if(data_in_file.size() < numberParticles){
    System::log<System::CRITICAL>("Not enough lines in file %s, expected %d, found %d", file, numberParticles, data_in_file.size());
  }
  data_in_file.resize(numberParticles);
  return data_in_file;
}

//Given a file name, reads positions from it, the format of file is a list of lines containing 3 numbers each, X Y Z for each particle.
//Only the number of particles given in the data.main will be read into positions.
void initializePositionsFromFile(std::string file, std::shared_ptr<ParticleData> pd){
  int numberParticles = pd->getNumParticles();
  auto pos_file = readPropertyFromFile<real3>(file, numberParticles);
  auto pos = pd->getPos(access::location::cpu, access::mode::write);
  std::transform(pos_file.begin(), pos_file.end(), pos.begin(), [](real3 p){return make_real4(p);});
}

//Initialize positions in an face centered cubic lattice, notice that several other lattices are available in initLattice, such as sc (simple cubic).
void initializePositionsInFCCLattice(std::shared_ptr<ParticleData> pd, real3 L){
  int numberParticles = pd->getNumParticles();
  auto pos = pd->getPos(access::location::cpu, access::mode::write);
  auto initial =  initLattice(L, numberParticles, fcc);
  std::copy(initial.begin(), initial.end(), pos.begin());
}

//Given a file name, reads the  charges from it, the format of file is a list of charges.
//The number of particles given in the data.main will be read into uammd charges.
void initializeChargesFromFile(std::string file, std::shared_ptr<ParticleData> pd){
  int numberParticles = pd->getNumParticles();
  auto charge_file = readPropertyFromFile<real>(file, numberParticles);
  auto charge = pd->getCharge(access::location::cpu, access::mode::write);
  std::copy(charge_file.begin(), charge_file.end(), charge.begin());
}

//this function will initialize particle velocities following a Boltzmann Distribution according to the given temperature
void initVelocitiesBoltzmannDistribution(std::shared_ptr<System> sys, std::shared_ptr<ParticleData> pd, real temperature){
  auto vel = pd->getVel(access::cpu, access::write);
  std::mt19937 gen(sys->rng().next());
  real mean = 0;
  real stdev = sqrt(temperature);
  std::normal_distribution<real> dis(mean, stdev);
  std::generate(vel.begin(), vel.end(), [&](){return make_real3(dis(gen),dis(gen),dis(gen));});
}

//Initialize particle properties, some properties are only initialized as a function of the parameters in the data.main
//If no initial position file is given in data.main's "readFile" positions are initialized in an FCC lattice
//Charges are only initialized if useElectrostatics is true.
//Masses are always initialzied to 1
//If a verletNVE integrator has beeen chosen, velocities are initialized following a Boltzmann distribution with the temperature in the data.main
//Execution is then handed over to furtherParticleInitialization in customizations.cuh for additional particle startup operations (if any).
std::shared_ptr<ParticleData> initializeParticles(std::shared_ptr<System> sys, Parameters par){
  auto pd = std::make_shared<ParticleData>(par.numberParticles, sys);
  if(par.readFile.empty()){
    initializePositionsInFCCLattice(pd, par.L);
  }
  else{
    initializePositionsFromFile(par.readFile, pd);
  }
  if(par.chargeReadFile.empty()){
    auto charges = pd->getCharge(access::cpu, access::write);
    std::fill(charges.begin(), charges.end(), 1);
  }
  else{
    initializeChargesFromFile(par.chargeReadFile, pd);
  }
  {
    auto mass = pd->getMass(access::cpu, access::write);
    std::fill(mass.begin(), mass.end(), 1);
  }
  if(par.integrator.compare("VerletNVE") == 0){
    initVelocitiesBoltzmannDistribution(sys, pd, par.temperature);
  }
  furtherParticleInitialization(pd, par);
  return pd;
}

//Initialize the basic UAMMD environment, including particles
//data.main will be named "data.main" unless at least an argument is passed to the program, in which case the data.main name will be the first argument
UAMMD initialize(int argc, char *argv[]){
  UAMMD sim;
  //System can optionally be handed argc/argv. In this case some command line options can be passed (see wiki for a list).
  //Things like the used GPU are selected through this.
  sim.sys = std::make_shared<System>(argc, argv);
  //System provides an uammd-wide random generator, which can be seeded like this
  std::random_device r;
  sim.sys->rng().setSeed(r());
  std::string datamain = (argc>1)?argv[1]:"data.main";
  //Process data.main
  sim.par = readParameters(datamain);
  //Initialize particle properties according to data.main
  sim.pd = initializeParticles(sim.sys, sim.par);
  return sim;
}


/*
  The family of functions below are easily copy pastable and will create and return instances of Integrators/Interactors
 */

//Brownian Dynamics
using BDMethod = BD::EulerMaruyama;
std::shared_ptr<BDMethod> createIntegratorBD(UAMMD sim){
  typename BDMethod::Parameters par;
  par.temperature = sim.par.temperature;
  par.viscosity = sim.par.viscosity;
  par.hydrodynamicRadius = sim.par.hydrodynamicRadius;
  par.dt = sim.par.dt;
  auto pg = std::make_shared<ParticleGroup>(sim.pd, sim.sys, "All");
  return std::make_shared<BDMethod>(sim.pd, pg, sim.sys, par);
}

using Verlet = VerletNVT::GronbechJensen;
std::shared_ptr<Verlet> createIntegratorVerletNVT(UAMMD sim){
  typename Verlet::Parameters par;
  par.temperature = sim.par.temperature;
  par.friction = sim.par.friction;
  par.dt = sim.par.dt;
  auto pg = std::make_shared<ParticleGroup>(sim.pd, sim.sys, "All");
  return std::make_shared<Verlet>(sim.pd, pg, sim.sys, par);
}

std::shared_ptr<Integrator> createIntegratorVerletNVE(UAMMD sim){
  typename VerletNVE::Parameters par;
  par.dt = sim.par.dt;
  par.energy = 1; //Optionally a target energy can be passed that VerletNVE will set according to velocities keep constant
  //par.initVelocities = false; //If true, velocities will be initialized by the module to ensure the desired energy
  //Note that it does not make sense to pass an energy and prevent VerletNVE from initializing velocities to match it.
  auto pg = std::make_shared<ParticleGroup>(sim.pd, sim.sys, "All");
  return std::make_shared<VerletNVE>(sim.pd, pg, sim.sys, par);
}

//Dissipative Particle Dynamics
//DPD is handled by UAMMD as a VerletNVE integrator with a special short range interaction
std::shared_ptr<Integrator> createIntegratorDPD(UAMMD sim){
  using NVE = VerletNVE;
  NVE::Parameters par;
  par.dt = sim.par.dt;
  par.initVelocities = false;
  auto pg = std::make_shared<ParticleGroup>(sim.pd, sim.sys, "All");
  auto verlet = std::make_shared<NVE>(sim.pd, pg, sim.sys, par);
  using DPD = PairForces<Potential::DPD, NeighbourList>;
  Potential::DPD::Parameters dpd_params;
  dpd_params.cutOff = sim.par.cutOff_dpd;
  dpd_params.temperature = sim.par.temperature;
  dpd_params.gamma = sim.par.gamma_dpd;
  dpd_params.A = sim.par.A_dpd;
  dpd_params.dt = par.dt;
  auto pot = std::make_shared<Potential::DPD>(sim.sys, dpd_params);
  DPD::Parameters params;
  params.box = Box(sim.par.L);
  auto pairforces = std::make_shared<DPD>(sim.pd, pg, sim.sys, params, pot);
  verlet->addInteractor(pairforces);
  return verlet;
}

//Smoothed Particle Hydrodynamics
std::shared_ptr<Integrator> createIntegratorSPH(UAMMD sim){
  using NVE = VerletNVE;
  NVE::Parameters par;
  par.dt = sim.par.dt;
  par.initVelocities = false;
  auto pg = std::make_shared<ParticleGroup>(sim.pd, sim.sys, "All");
  auto verlet = std::make_shared<NVE>(sim.pd, pg, sim.sys, par);
  SPH::Parameters params;
  params.box = Box(sim.par.L);
  //Pressure for a given particle "i" in SPH will be computed as gasStiffness·(density_i - restDensity)
  //Where density is computed as a function of the masses of the surroinding particles
  //Particle mass starts as 1, but you can change this in customizations.cuh
  params.support = sim.par.support_sph;   //Cut off distance for the SPH kernel
  params.viscosity = sim.par.viscosity;   //Environment viscosity
  params.gasStiffness = sim.par.gasStiffness_sph;
  params.restDensity = sim.par.restDensity_sph;
  auto sph = std::make_shared<SPH>(sim.pd, pg, sim.sys, params);
  verlet->addInteractor(sph);
  return verlet;
}

//Creates a triply periodic Brownian Dynamics with Hydrodynamic Interactions integration module
std::shared_ptr<Integrator> createIntegratorBDHI(UAMMD sim){
  //There are several hydrodynamics modules, we choose between Positively Split Ewald (PSE) or Force Coupling Method (FCM) here
  // mainly for performance reasons. FCM is faster for small and/or dense systems, but it is limited in the system size by memory.
  // PSE can be slower it temperature>0, but does not have that system size constraints.
  //FCM scales linearly with system size (so doubling the box size in the three dimensions makes it 8 times slower) and number of particles
  //PSE scales linearly with the number of particles, independently of system size. But the "psi" parameter must be tweaked to find the optimal performance for each case.
  //See the wiki for more information about these modules
  real maxL = std::max({sim.par.L.x, sim.par.L.y, sim.par.L.z});
  int maxcells = maxL/sim.par.hydrodynamicRadius;
  //In both modules, particle self diffusion coefficient will be T/(6*pi*viscosity*hydrodynamicRadius) or close to it
  if(maxcells >= 128){
    using Scheme = BDHI::PSE;
    Scheme::Parameters par;
    par.box = Box(sim.par.L);
    par.temperature = sim.par.temperature;
    par.viscosity = sim.par.viscosity;
    par.dt = sim.par.dt;
    par.hydrodynamicRadius = sim.par.hydrodynamicRadius;
    par.tolerance = 1e-4;
    //Balances the load of the algorithm, low values work best for dilute and/or big systems.
    // Higher values will work best for dense and/or small systems.
    par.psi = 1.0/par.hydrodynamicRadius;
    auto bdhi = std::make_shared<BDHI::EulerMaruyama<Scheme>>(sim.pd, sim.sys, par);
    return bdhi;
  }
  else{
    using Scheme = BDHI::FCM;
    Scheme::Parameters par;
    par.box = Box(sim.par.L);
    par.temperature = sim.par.temperature;
    par.viscosity = sim.par.viscosity;
    par.dt = sim.par.dt;
    par.hydrodynamicRadius = sim.par.hydrodynamicRadius;
    par.tolerance = 1e-4;
    auto bdhi = std::make_shared<BDHI::EulerMaruyama<Scheme>>(sim.pd, sim.sys, par);
    return bdhi;
  }
}

//Fluctuating Immersed Boundary
std::shared_ptr<Integrator> createIntegratorFIB(UAMMD sim){
  BDHI::FIB::Parameters par;
  par.temperature = sim.par.temperature;
  par.viscosity = sim.par.viscosity;
  par.hydrodynamicRadius = sim.par.hydrodynamicRadius;
  par.dt = sim.par.dt;
  //par.scheme = BDHI::FIB::IMPROVED_MIDPOINT;
  par.scheme = BDHI::FIB::MIDPOINT;
  par.box = Box(sim.par.L);
  auto pg = std::make_shared<ParticleGroup>(sim.pd, sim.sys, "All");
  return std::make_shared<BDHI::FIB>(sim.pd, pg, sim.sys, par);
}

//Inertial Coupling Method
std::shared_ptr<Integrator> createIntegratorICM(UAMMD sim){
  Hydro::ICM::Parameters par;
  par.temperature = sim.par.temperature;
  par.viscosity = sim.par.viscosity;
  par.density = 1;
  par.hydrodynamicRadius = sim.par.hydrodynamicRadius;
  par.dt = sim.par.dt;
  par.box = Box(sim.par.L);
  return std::make_shared<Hydro::ICM>(sim.pd, sim.sys, par);
}

//Create the integrator as selected via data.main.
//Notice that all Integrators can be passed around as shared_ptr<Integrator>
std::shared_ptr<Integrator> createIntegrator(UAMMD sim){
  if(sim.par.integrator.compare("BD") == 0){
    return createIntegratorBD(sim);
  }
  else if(sim.par.integrator.compare("VerletNVT") == 0){
    return createIntegratorVerletNVT(sim);
  }
  else if(sim.par.integrator.compare("VerletNVE") == 0){
    return createIntegratorVerletNVE(sim);
  }
  else if(sim.par.integrator.compare("DPD") == 0){
    return createIntegratorDPD(sim);
  }
  else if(sim.par.integrator.compare("SPH") == 0){
    return createIntegratorSPH(sim);
  }
  else if(sim.par.integrator.compare("BDHI") == 0){
    return createIntegratorBDHI(sim);
  }
  else if(sim.par.integrator.compare("FIB") == 0){
    return createIntegratorFIB(sim);
  }
  else if(sim.par.integrator.compare("ICM") == 0){
    return createIntegratorICM(sim);
  }

  else{
    sim.sys->log<System::CRITICAL>("Invalid integrator. Choose BD, BDHI, FIB, ICM, VerletNVT, VerletNVE, SPH or DPD");
  }
  return nullptr;
}

std::shared_ptr<Interactor> createShortRangeInteractor(UAMMD sim){
  //The potential used is defined in customizations.cuh and describes an LJ interaction by default
  auto pot = std::make_shared<ShortRangePotential>(sim.par);
  using SR = PairForces<ShortRangePotential, NeighbourList>;
  typename SR::Parameters params;
  params.box = Box(sim.par.L);
  auto pairForces = std::make_shared<SR>(sim.pd, sim.sys, params, pot);
  return pairForces;
}

std::shared_ptr<Interactor> createBondInteractor(UAMMD sim){
  using Bond = HarmonicBond;
  using BF = BondedForces<Bond>;
  typename BF::Parameters params;
  params.file = sim.par.bondFile;
  auto bf = std::make_shared<BF>(sim.pd, sim.sys, params, Bond(sim.par));
  return bf;
}

std::shared_ptr<Interactor> createAngularBondInteractor(UAMMD sim){
  using Bond = Angular;
  using BF = AngularBondedForces<Bond>;
  typename BF::Parameters params;
  params.file = sim.par.angularBondFile;
  auto bf = std::make_shared<BF>(sim.pd, sim.sys, params, Bond(sim.par));
  return bf;
}

std::shared_ptr<Interactor> createTorsionalBondInteractor(UAMMD sim){
  using Bond = Torsional;
  using BF = TorsionalBondedForces<Bond>;
  typename BF::Parameters params;
  params.file = sim.par.torsionalBondFile;
  auto bf = std::make_shared<BF>(sim.pd, sim.sys, params, Bond(sim.par));
  return bf;
}

std::shared_ptr<Interactor> createElectrostaticInteractor(UAMMD sim){
  //Similar to the hydrodynamics modules, electrostatics are available in two algorithms, in this case exposed under the same module
  //One works best for dilute/big systems (Ewald splitting) and the other for dense/small systems (with no Ewald splitting).
  //Ewald splitting is automatically selected according to some system size heuristic.
  //As with BDHI, Ewald splitting is slower but provides an algorithm that scales linearly with the number of particles independently of system size.
  //Non-Ewald splitted version scales linearly with number of particles but also scales linear with system size.
  using Electro = Poisson;
  Electro::Parameters par;
  par.box = Box(sim.par.L);
  par.epsilon = sim.par.permittivity;
  par.tolerance = 1e-4;
  par.gw = sim.par.gaussianWidth;
  real maxL = std::max({sim.par.L.x, sim.par.L.y, sim.par.L.z});
  int maxcells = maxL/par.gw;
  if(maxcells >= 128){
    //Controls Ewald splitting, if the parameter is nor present no Ewald splitting is used.
    par.split = 0.07/par.gw;
  }
  auto elec = std::make_shared<Electro>(sim.pd, sim.sys, par);
  return elec;
}

std::shared_ptr<Interactor> createExternalPotentialInteractor(UAMMD sim){
  //Uses the external potential defined in customizations.cuh
  auto gr = GravityAndWall(sim.par);
  auto ext = std::make_shared<ExternalForces<GravityAndWall>>(sim.pd, sim.sys, gr);
  return ext;  
}


double sumTotalEnergy(std::shared_ptr<Integrator> integrator, std::shared_ptr<ParticleData> particles){
  {
    auto energy = particles->getEnergy(access::location::gpu, access::mode::write);
    thrust::fill(thrust::cuda::par, energy.begin(), energy.end(), real(0.0));
  }
  integrator->sumEnergy();
  for(auto interactor: integrator->getInteractors()){
    interactor->sumEnergy();
  }
  auto energy = particles->getEnergy(access::location::gpu, access::mode::read);
  double totalEnergy = thrust::reduce(thrust::cuda::par, energy.begin(), energy.end(), 0.0);
  return totalEnergy;
}

//Write the current simulation status to a file.
//Only writes positions, as a list of lines containing X Y Z for each particle.
//Each snapshot is separated by a line containing only a "#".
//Particle positions are folded into the primary box via MIC.
void writeSimulation(UAMMD sim){
  //Particles might change ordering once in a while, this array always points to the index of a particle given its name or id.
  //This means that a particle that started at index "i" (and thus has the id or name "i") will always be accesible via
  //id2index[i]
  auto id2index = sim.pd->getIdOrderedIndices(access::cpu);
  static std::ofstream out(sim.par.outfile);
  static std::ofstream eout(sim.par.outfileEnergy);
  static std::ofstream vout(sim.par.outfileVelocities);
  if(eout.good())
    eout<<std::setprecision(2*sizeof(real))<<"# Total energy: "<<sumTotalEnergy(sim.integrator, sim.pd)<<std::endl;
  if(sim.pd->isVelAllocated() and vout.good())
    vout<<"#"<<std::endl;
  Box box(sim.par.L);
  real3 L = box.boxSize;
  out<<"#Lx="<<L.x*0.5<<";Ly="<<L.y*0.5<<";Lz="<<L.z*0.5<<";\n";
  auto pos = sim.pd->getPos(access::location::cpu, access::mode::read);
  auto vel = sim.pd->getVelIfAllocated(access::location::cpu, access::mode::read);
  auto energy = sim.pd->getEnergyIfAllocated(access::location::cpu, access::mode::read);
  fori(0, sim.par.numberParticles){
    //real3 p = box.apply_pbc(make_real3(pos[id2index[i]]));
    real3 p = make_real3(pos[id2index[i]]);
    out<<std::setprecision(2*sizeof(real))<<p<<"\n";
    if(eout.good()) eout<<energy[id2index[i]]<<"\n";
    if(sim.pd->isVelAllocated() and vout.good()) vout<<vel[id2index[i]]<<"\n";
  }
  out<<std::flush;
  vout<<std::flush;
  eout<<std::flush;
}

//Create and run the simulation
int main(int argc, char *argv[]){  
  auto sim = initialize(argc, argv);
  auto bd = createIntegrator(sim);
  sim.integrator = bd;
  if(! sim.par.bondFile.empty()){
    bd->addInteractor(createBondInteractor(sim));
  }
  if(! sim.par.angularBondFile.empty()){
    bd->addInteractor(createAngularBondInteractor(sim));
  }
  if(! sim.par.torsionalBondFile.empty()){
    bd->addInteractor(createTorsionalBondInteractor(sim));
  }
  if(sim.par.useElectrostatics){
    bd->addInteractor(createElectrostaticInteractor(sim));
  }
  if(sim.par.enableGravity){
    bd->addInteractor(createExternalPotentialInteractor(sim));
  }
  if(sim.par.epsilon > 0){
    bd->addInteractor(createShortRangeInteractor(sim));
  }
  forj(0, sim.par.relaxSteps){
    bd->forwardTime();
  }
  Timer tim;
  tim.tic();
  forj(0, sim.par.numberSteps){
    bd->forwardTime();
    if(sim.par.printSteps > 0 and j%sim.par.printSteps==0){
      writeSimulation(sim);
    }
  }
  auto totalTime = tim.toc();
  System::log<System::MESSAGE>("mean FPS: %.2f", sim.par.numberSteps/totalTime);
  return 0;
}

//Writes an example data.main to a file with name given in "datamain"
void writeDefaultDatamain(std::string datamain){
  std::ofstream out(datamain);
  out<<"#About the format of this file:"<<std::endl;
  out<<"#Comments are allowed like this"<<std::endl;
  out<<"#data.main format is [option] [argument1] [argument2]..."<<std::endl;
  out<<"#No arguments can be a valid number of arguments"<<std::endl;
  out<<"#Everything after the expected arguments is ignored"<<std::endl;
  out<<"#The code generic_simulation.cu will read this file and run the simulation it describes"<<std::endl;
  out<<"#Below, you can tweak or turn off and on the different functionalities"<<std::endl;
  out<<"#This example data.main encodes a LJ liquid simulation"<<std::endl;
  out<<""<<std::endl;
  out<<"################################################################"<<std::endl;
  out<<"#Choose the integration scheme, it can be any of the following:"<<std::endl;
  out<<"#  -VerletNVT: Molecular dynamics with a thermostat."<<std::endl;
  out<<"#  -VerletNVE: Molecular dynamics with constant energy (starting with a velocity distribution to match the given temperature)."<<std::endl;
  out<<"#  -BD: Brownian dynamics"<<std::endl;
  out<<"#  -BDHI: Brownian dynamics with Hydrodynamic interactions (via positively split Ewald Rotne-Prager-Yamakawa or Force Coupling method)"<<std::endl;
  out<<"#  -SPH: Smooth Particle Hydrodynamics"<<std::endl;
  out<<"#  -DPD: Dissipative Particle Dynamics"<<std::endl;
  out<<"#  -FIB: Fluctuating Immersed Boundary"<<std::endl;
  out<<"#  -ICM: Inertial Coupling Method"<<std::endl;
  out<<"integrator VerletNVT"<<std::endl;
  out<<""<<std::endl;
  out<<"################################################################"<<std::endl;
  out<<"#Integrator-specific parameters. Only the parameters for the chosen integrator have to be present (they will be ignored for the rest)"<<std::endl;
  out<<"#These parameters will be used by DPD only."<<std::endl;
  out<<"#A_dpd 1   #Conservative force repulsion strength"<<std::endl;
  out<<"#gamma_dpd 4 #Dissipative strength"<<std::endl;
  out<<"#cutOff_dpd 1"<<std::endl;
  out<<"#These parameters will be used by SPH only"<<std::endl;
  out<<"#Pressure for a given particle \"i\" in SPH will be computed as gasStiffness·(density_i - restDensity)"<<std::endl;
  out<<"#Where density is computed as a function of the masses of the surroinding particles"<<std::endl;
  out<<"#Particle mass starts as 1, but you can change this in customizations.cuh"<<std::endl;
  out<<"#gasStiffness_sph 100 "<<std::endl;
  out<<"#support_sph 2.4 #Cut off distance for the SPH kernel"<<std::endl;
  out<<"#restDensity_sph 0.5"<<std::endl;
  out<<"#viscosity 10"<<std::endl;
  out<<"#The hydrodynamic radius of the particles will be used only in BD and BDHI."<<std::endl;
  out<<"hydrodynamicRadius 1.0"<<std::endl;
  out<<"#################################################################"<<std::endl;
  out<<"#Environment parameters"<<std::endl;
  out<<"viscosity 1"<<std::endl;  
  out<<"friction 1 #Used in VerletNVT"<<std::endl;
  out<<"temperature 0.1"<<std::endl;
  out<<"##################################################################"<<std::endl;
  out<<"#Short range potential parameters"<<std::endl;
  out<<"sigma 1 #Diameter of the particles in the LJ potential"<<std::endl;
  out<<"epsilon 1 #Strength of the LJ potential, 0 turns off the interaction"<<std::endl;
  out<<"cutOff  2.5 #CufOff of the LJ potential, 2^(1/6)*sigma will result in an WCA potential"<<std::endl;
  out<<"#################################################################"<<std::endl;
  out<<"#Parameters for other potentials"<<std::endl;
  out<<"#bondFile pair.bonds"<<std::endl;
  out<<"#If not present no bonds will be placed between particles"<<std::endl;
  out<<"#The bond file is a list of bonds (joining particle \"i\" with \"j\" with an spring of constant \"Kspring\" and eq. distance \"r0\")."<<std::endl;
  out<<"#There can also be bonds joining particles with points in space, these are placed after the particle-particle bonds."<<std::endl;
  out<<"#The format of the file is:"<<std::endl;
  out<<"#numberOfBonds"<<std::endl;
  out<<"#i j Kspring r0"<<std::endl;
  out<<"#..."<<std::endl;
  out<<"#numberOfFixedPointBonds"<<std::endl;
  out<<"#x y z Kspring r0"<<std::endl;
  out<<"#..."<<std::endl;
  out<<""<<std::endl;
  out<<"#angularBondFile ang.bonds"<<std::endl;
  out<<"#If not present, no three-particle angular bonds will be placed, similar to the pair-bond file, the format of the angular bonds file is:"<<std::endl;
  out<<"#numberOfBonds"<<std::endl;
  out<<"#i j k Kspring equilibrium_angle"<<std::endl;
  out<<"#..."<<std::endl;
  out<<"#Notice that in this case the bond involves particles i,j and k (being j the middle one). i-j-k"<<std::endl;
  out<<""<<std::endl;
  out<<"#torsionalBondFile tor.bonds"<<std::endl;
  out<<"#If not present, no four-particle torsional bonds will be placed, this time the file format is:"<<std::endl;
  out<<"#numberOfBonds"<<std::endl;
  out<<"#i j k l Kspring equilibrium_angle"<<std::endl;
  out<<"#..."<<std::endl;
  out<<"#Notice that in this case the bond involves particles i,j,k and l. An eq. angle of 0 will result in the 4 particles being in the same plane."<<std::endl;
  out<<""<<std::endl;
  out<<"#Choose if particles interact electrostatically"<<std::endl;
  out<<"#useElectrostatics"<<std::endl;
  out<<"#permittivity 1"<<std::endl;
  out<<"#Charges are modeled as gaussian charge clouds, a point charge will be equivalent to gaussianWidth=0 (a limit which the algorithm cannot handle)."<<std::endl;
  out<<"#gaussianWidth 0.5"<<std::endl;
  out<<"#Optionally charges can be read from a file (a simple list of charges for each particle)."<<std::endl;
  out<<"#Charge defaults to 1 for all particles if this file is not present, again, you may customize this further in customizations.cuh"<<std::endl;
  out<<"#chargeReadFile charges.dat"<<std::endl;
  out<<""<<std::endl;
  out<<"#Enable gravity and wall"<<std::endl;
  out<<"#If present, a wall will be placed at z=-Lz*0.5 and gravity will pull the particles in the z direction"<<std::endl;
  out<<"#enableGravity"<<std::endl;
  out<<"#gravity 1 #gravity strength"<<std::endl;
  out<<"#kwall 20  #Repulsive strength of the wall"<<std::endl;
  out<<"###########################################################################"<<std::endl;
  out<<"#Simulation parameters"<<std::endl;
  out<<"numberParticles 16384  #Number of particles, this number of lines will be read when reading files like \"readFile\""<<std::endl;
  out<<"#readFile init.pos     #If the option is present, initial positions will be read from here. Otherwise particles will start in an FCC lattice"<<std::endl;
  out<<"#The format of this file is simply a sequence of \"numberParticle\" lines with the positions of each particle as: X Y Z"<<std::endl;
  out<<"#Notice that /dev/stdin is a valid file name, in which case initial positions can be piped"<<std::endl;
  out<<""<<std::endl;
  out<<"L 32 32 32             #Domain size in the three dimensions"<<std::endl;
  out<<"outfile /dev/stdout    #Positions will be written to this file"<<std::endl;
  out<<"#outfileVelocities /dev/null #If present velocities will be written to this file (if the integrator uses velocity)"<<std::endl;
  out<<"#outfileEnergy /dev/null    #If present per particle energy will be written to this file for each spanshot"<<std::endl;
  out<<"numberSteps 100000     #Number of simulation steps"<<std::endl;
  out<<"printSteps 500         #Positions will be writen to \"outfile\" every printSteps steps"<<std::endl;
  out<<"relaxSteps  100       #Run this number of steps before writting"<<std::endl;
  out<<"dt 0.001               #Time discretization step"<<std::endl;
  out<<""<<std::endl;
  out<<"#You can use the special shell option to run something in a shell when the data.main is processed"<<std::endl;
  out<<"shell echo hello from data.main > /dev/stderr"<<std::endl;
}

//Reads parameters from file with name given in datamain. If the file is not present, a default will be printed and the execution
//will end.
//This function calls readCustomParamters in customizations.cuh before exiting.
Parameters readParameters(std::string datamain){
  if(not std::ifstream(datamain).good()){
    System::log<System::WARNING>("Parameter file %s not found, autogenerating and exiting", datamain.c_str());
    writeDefaultDatamain(datamain);
    exit(0);
  }
  InputFile in(datamain);
  Parameters par;
  in.getOption("L", InputFile::Required)>>par.L.x>>par.L.y>>par.L.z;
  in.getOption("numberParticles", InputFile::Required)>>par.numberParticles;
  in.getOption("numberSteps", InputFile::Required)>>par.numberSteps;
  in.getOption("printSteps", InputFile::Required)>>par.printSteps;
  in.getOption("relaxSteps", InputFile::Required)>>par.relaxSteps;
  in.getOption("dt", InputFile::Required)>>par.dt;
  in.getOption("temperature", InputFile::Required)>>par.temperature;  
  in.getOption("outfile", InputFile::Required)>>par.outfile;
  in.getOption("outfileVelocities", InputFile::Optional)>>par.outfileVelocities;
  in.getOption("outfileEnergy", InputFile::Optional)>>par.outfileEnergy;
  in.getOption("epsilon", InputFile::Required)>>par.epsilon;
  in.getOption("sigma", InputFile::Required)>>par.sigma;
  in.getOption("cutOff", InputFile::Required)>>par.cutOff;
  in.getOption("readFile", InputFile::Optional)>>par.readFile;  
  in.getOption("bondFile", InputFile::Optional)>>par.bondFile;
  in.getOption("angularBondFile", InputFile::Optional)>>par.angularBondFile;
  in.getOption("torsionalBondFile", InputFile::Optional)>>par.torsionalBondFile;
  in.getOption("integrator", InputFile::Required)>>par.integrator;
  if(par.integrator.compare("SPH")==0){
    in.getOption("support_sph", InputFile::Required)>>par.support_sph;
    in.getOption("viscosity", InputFile::Required)>>par.viscosity;
    in.getOption("gasStiffness_sph", InputFile::Required)>>par.gasStiffness_sph;
    in.getOption("restDensity_sph", InputFile::Required)>>par.restDensity_sph;
  }
  if(par.integrator.compare("DPD")==0){
    in.getOption("A_dpd", InputFile::Required)>>par.A_dpd;
    in.getOption("gamma_dpd", InputFile::Required)>>par.gamma_dpd;
    in.getOption("cutOff_dpd", InputFile::Required)>>par.cutOff_dpd;    
  }
  if(par.integrator.compare("BD")==0 or par.integrator.compare("BDHI")==0 or
     par.integrator.compare("FIB")==0 or par.integrator.compare("ICM")==0){
    in.getOption("hydrodynamicRadius", InputFile::Required)>>par.hydrodynamicRadius;
    in.getOption("viscosity", InputFile::Required)>>par.viscosity;
  }
  if(par.integrator.compare("VerletNVT")==0){
     in.getOption("friction", InputFile::Required)>>par.friction;
  }
  if(in.getOption("useElectrostatics", InputFile::Optional)){
    par.useElectrostatics = true;
    in.getOption("permittivity", InputFile::Required)>>par.permittivity;
    in.getOption("gaussianWidth", InputFile::Required)>>par.gaussianWidth;
    in.getOption("chargeReadFile", InputFile::Optional)>>par.chargeReadFile;
  }
  
  if(in.getOption("enableGravity", InputFile::Optional)){
    par.enableGravity = true;
    in.getOption("gravity", InputFile::Required)>>par.gravity;
    in.getOption("kwall", InputFile::Required)>>par.kwall;
  }  
  readCustomParameters(in, par);
  return par;
}



