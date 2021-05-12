#pragma once
#include "wavefunction.h"
#include <armadillo>

using namespace arma;

class Gaussian_Binary : public WaveFunction{
public:
  Gaussian_Binary(class System *system, int N, double sigma);

  double evaluate(std::vector<class Particle*> particles);
  double computeDoubleDerivative(std::vector<class Particle*> particles);
  double gradientTerm(vec x);
  double laplacianTerm(vec x);

  std::vector<double> quantumForce(std::vector<class Particle*> particles);

  // TEMPORARY!!!!
  void test_weights_biases();

private:
  int m_M = 0;
  int m_N = 0;

  double m_sigma = 0;

  mat m_W;
  vec m_a;
  vec m_b;


};
