#include "../system.h"
#include "../metropolis_langevin.h"
#include "../particle.h"
#include "../WaveFunctions/gaussian_binary.h"
#include "../Hamiltonians/harmonicoscillator.h"
#include "../InitialStates/initialstate.h"
#include "../InitialStates/randomuniform.h"
#include "../Math/random.h"
#include "../sampler.h"
#include "../GDsampler.h"

#include <armadillo>
#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <omp.h>

using namespace std;
using namespace arma;




int main(int nargs, char **args)
{

  int seed = 2021;

  int numberOfDimensions = 2;
  int numberOfParticles  = 2;
  int numHiddenLayers    = 1;
  int numberOfSteps      = (int) 1e6;
  int equilibration      = (int) 1e5;
  double stepLength      = 0.1;
  double tol             = 1e-5;
  double learningRate    = 0.01;
  int maxIter            = 200;
  double omega           = 1;
  double std             = 0.5;

  string filename_blocking = "no";
  string path = "../../../output/convergence.txt";

  ofstream outfile("../../output/convergence.txt", ofstream::out);
  outfile.close();

  System *system;
  system = new System(seed);
  system->setOmega(1);

  system->setPath(path);
  system->m_energyfile = filename_blocking;

  system->setSampler(new Sampler(system));
  system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));
  system->setWaveFunction(new Gaussian_Binary(system, numHiddenLayers, 1./omega, std));
  system->setHamiltonian(new HarmonicOscillator(system, omega));
  system->setEquilibrationSteps(equilibration);
  system->setMetropolisSteps(numberOfSteps);
  system->setStepLength(stepLength);
  system->gradientDescent(tol, learningRate, maxIter);

  cout << "Finished!" << endl;

  return 0;
}
