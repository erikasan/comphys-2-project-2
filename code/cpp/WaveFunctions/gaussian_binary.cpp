#include "wavefunction.h"
#include "gaussian_binary.h"
#include "../system.h"
#include "../particle.h"
#include "../sampler.h"

#include <random>
#include <armadillo>
using namespace arma;

// TEMPORARY
#include <iostream>


using namespace std;

Gaussian_Binary::Gaussian_Binary(System *system, int N) : WaveFunction(system)
{
  m_M = (m_system->getNumberOfParticles()) * (m_system->getNumberOfDimensions());
  m_N = N;
  
  m_W.set_size(m_M, m_N);
  m_a.set_size(m_M);
  m_b.set_size(m_N);

  mt19937_64 generator;
  normal_distribution<double> distribution(0, 0.001);

  size_t i, j;

  for (j = 0; j < m_N; j++){    // Armadillo blasphemously uses column-major ordering
    for (i = 0; i < m_M; i++){
      m_W(i, j) = distribution(generator);
    }
  }

  for (i = 0; i < m_M; i++){
    m_a(i) = distribution(generator);
  }

  for (j = 0; j < m_N; j++){
    m_b(j) = distribution(generator);
  }

  return;
}

void Gaussian_Binary::test_weights_biases()
{
  size_t i, j;

  cout << "M = " << m_M << endl;
  cout << "N = " << m_N << endl;
  cout << endl;

  cout << "W" << endl;
  cout << m_W << endl;

  cout << "a" << endl;
  cout << m_a << endl;

  cout << "b" << endl;
  cout << m_b << endl;
  return;
}

double Gaussian_Binary::evaluate(std::vector<class Particle*> particles)
{
  return 0;
}

double Gaussian_Binary::computeDoubleDerivative(std::vector<class Particle*> particles)
{
  return 0;
}

std::vector<double> Gaussian_Binary::quantumForce(std::vector<class Particle*> particles)
{
  std::vector<double> qForce = std::vector<double>();
  return qForce;
}
