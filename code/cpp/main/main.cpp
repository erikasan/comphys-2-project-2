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

  int numberOfDimensions = 1;
  int numberOfParticles  = 1;
  //int numHiddenLayers    = atoi(args[1]);
  int numberOfSteps      = (int) 1e6;
  int equilibration      = (int) 1e5;
  double stepLength      = 0.1;
  double tol             = 1e-5;
  double learningRate    = 0.001;
  int maxIter            = 20;
  double sigma2           = 1;
  //double omega           = atof(args[2]);

  //string filename_blocking = "energy";
  string path= "../";

  System *system;

  double energy;
  ofstream outfile;

  ivec hiddenlayers = regspace<ivec>(1, 10);
  vec omegas = {1};
  #pragma omp parallel for private(system, energy, outfile)
  for (int i = 0; i < hiddenlayers.n_elem; i++){
    int numHiddenLayers = hiddenlayers(i);
    //outfile.open("../../output/energies_SP1D_imp_h=" + to_string(numHiddenLayers));
    system = new MetropolisLangevin(seed + omega);
    system->setOmega(1);

    string filename_blocking = "energy_h=" + to_string(numHiddenLayers);
    system->setPath(path);
    system->m_energyfile = filename_blocking;

    system->setSampler(new Sampler(system));
    system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setWaveFunction(new Gaussian_Binary(system, numHiddenLayers, 1));
    system->setHamiltonian(new HarmonicOscillator(system, 1));
    system->setEquilibrationSteps(equilibration);
    system->setMetropolisSteps(numberOfSteps);
    system->setStepLength(stepLength);
    system->gradientDescent(tol, learningRate, maxIter);
  }

  cout << "Finished!" << endl;
  return 0;
}
