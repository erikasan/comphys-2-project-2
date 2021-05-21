#include "../system.h"
#include "../metropolis_langevin.h"
#include "../particle.h"
#include "../WaveFunctions/gaussian_binary.h"
#include "../Hamiltonians/HO_with_coulomb.h"
#include "../InitialStates/initialstate.h"
#include "../InitialStates/randomuniform.h"
#include "../Math/random.h"
#include "../sampler.h"
#include "../GDsampler.h"

#include "../WaveFunctions/simplegaussian.h"

#include <armadillo>
#include <iostream>
#include <string>

using namespace std;
using namespace arma;




int main(int nargs, char **args)
{

  int seed = 2021;

  int numberOfDimensions = 3;
  int numberOfParticles  = 2;
  int numHiddenLayers    = 10;
  int numberOfSteps      = (int) 1e6;
  int equilibration      = (int) 1e5;
  double stepLength      = 0.1;
  double tol             = 1e-8;
  double learningRate    = 0.001;
  int maxIter            = 100;
  double sigma           = 1;
  double omega           = 1./4;

  string filename_blocking = "energies";
  string path = "../../../output/";

  System *system;
  system = new MetropolisLangevin(seed);
  system->setOmega(1);

  // system->m_energyfile = filename_blocking;
  // system->setPath(path);

  system->setSampler(new Sampler(system));
  system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));
  system->setWaveFunction(new Gaussian_Binary(system, numHiddenLayers, sigma));
  system->setHamiltonian(new HO_with_Coulomb(system, omega));
  system->setEquilibrationSteps(equilibration);
  system->setMetropolisSteps(numberOfSteps);
  system->setStepLength(stepLength);
  system->gradientDescent(tol, learningRate, maxIter);

  cout << "Finished!" << endl;

  return 0;
}
