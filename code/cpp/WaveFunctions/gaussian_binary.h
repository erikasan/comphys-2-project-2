#pragma once
#include "wavefunction.h"

class Gaussian_Binary : public WaveFunction{
public:
  Gaussian_Binary(class System *system, int N);


  double evaluate(std::vector<class Particle*> particles);
  double evaluate(std::vector<class Particle*> particles, int particle_id);
  double computeDoubleDerivative(std::vector<class Particle*> particles);
  std::vector<double> quantumForce(std::vector<class Particle*> particles, int particle_id);
  std::vector<double> quantumForce(std::vector<class Particle*> particles);
  void updateDistances(std::vector<class Particle*> particles,int particle_id);
  void initiateDistances(std::vector<class Particle*> particles);
  void sample(std::vector<class Particle*> particles, double localEnergy);
  void computeAverages(double steps);
  void gradientDescent();



  // TEMPORARY!!!!
  void test_weights_biases();

private:
  int m_M = 0;
  int m_N = 0;

  double **m_W = nullptr;
  double  *m_a = nullptr;
  double  *m_b = nullptr;


};
