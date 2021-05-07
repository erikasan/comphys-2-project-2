#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/gaussian_binary.h"
#include "../InitialStates/initialstate.h"
#include "../InitialStates/randomuniform.h"
#include "../Math/random.h"

#include "../WaveFunctions/simplegaussian.h"


int main(int nargs, char **args)
{

  int seed = 2021;

  int numberOfDimensions = 2;
  int numberOfParticles  = 2;
  int numHiddenLayers    = 2;

  System *system;
  system = new System(seed);

  system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));

  system->setWaveFunction(new Gaussian_Binary(system, numHiddenLayers));
  system->getWaveFunction()->test_weights_biases();
  return 0;
}
