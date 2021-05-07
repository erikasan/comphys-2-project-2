#pragma once
#include "wavefunction.h"

class Gaussian_Binary : public WaveFunction{
public:
  Gaussian_Binary(class System *system);
  double computeDoubleDerivative(std::vector<class Particle*> particles);

private:
  // Weights and biases should be private variables
};
