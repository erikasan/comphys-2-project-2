#include "gaussian_binary.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include "../sampler.h"

Gaussian_Binary::Gaussian_Binary(System *system) : WaveFunction(system) {
  // Initialize weights and biases at random
}
