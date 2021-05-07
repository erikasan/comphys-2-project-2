#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/gaussian_binary.h"
#include "../InitialStates/initialstate.h"
#include "../InitialStates/randomuniform.h"
#include "../Math/random.h"

int main(int nargs, char **args)
{

  int seed = 2021;

  int numberOfDimensions = 2;
  int numberOfParticles  = 2;

  System *system;
  system = new System(seed);

  return 0;
}
