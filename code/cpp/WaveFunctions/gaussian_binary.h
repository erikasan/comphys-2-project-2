#pragma once
#include "wavefunction.h"

class Gaussian_Binary : public WaveFunction{
public:
  Gaussian_Binary(class System *system, int N);

  double evaluate(std::vector<class Particle*> particles);
  double computeDoubleDerivative(std::vector<class Particle*> particles);
  std::vector<double> quantumForce(std::vector<class Particle*> particles);



  // TEMPORARY!!!!
  void test_weights_biases();

private:
  int m_M = 0;
  int m_N = 0;

  double *m_W = nullptr;
  double *m_a = nullptr;
  double *m_b = nullptr;


};
