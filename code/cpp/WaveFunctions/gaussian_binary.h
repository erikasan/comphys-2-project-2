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

  void sample(std::vector<class Particle*> particles, double localEnergy);
  void computeAverages(double steps);
  void gradientDescent();

  vec grad_a(vec x);
  vec grad_b(vec x);
  mat grad_W(vec x);
  vec localEnergygrad_a(vec x, double localEnergy);
  vec localEnergygrad_b(vec x, double localEnergy);
  mat localEnergygrad_W(vec x, double localEnergy);

  vec convertPositionToArmadillo(std::vector<class Particle*> particles);

  std::vector<double> quantumForce(std::vector<class Particle*> particles);

private:
  int m_M = 0;
  int m_N = 0;

  double m_sigma2 = 0;

  mat m_W;
  vec m_a;
  vec m_b;

  vec m_av_grad_a;
  vec m_av_grad_b;
  mat m_av_grad_W;
  vec m_av_local_energy_grad_a;
  vec m_av_local_energy_grad_b;
  mat m_av_local_energy_grad_W;

};
