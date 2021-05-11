#include "wavefunction.h"
#include "gaussian_binary.h"
#include "../system.h"
#include "../particle.h"
#include "../sampler.h"

#include <random>

// TEMPORARY
#include <iostream>

using namespace std;

Gaussian_Binary::Gaussian_Binary(System *system, int N) : WaveFunction(system)
{
  m_M = (m_system->getNumberOfParticles()) * (m_system->getNumberOfDimensions());
  m_N = N;

  m_W = new double[m_M*m_N];
  m_a = new double[m_M];
  m_b = new double[m_N];

  mt19937_64 generator;
  normal_distribution<double> distribution(0, 0.001);

  size_t i, j;

  for (i = 0; i < m_M; i++){
    m_a[i] = 1;//distribution(generator);
    for (j = 0; j < m_N; j++){
      m_W[i*m_N + j] = 1;//distribution(generator);
    }
  }

  for (j = 0; j < m_N; j++){
    m_b[j] = 1;//distribution(generator);
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
  for (i = 0; i < m_M; i++){
    for (j = 0; j < m_N; j++){
      cout << m_W[i*m_N + j] << " ";
    }
    cout << endl;
  }
  cout << endl;

  cout << "a" << endl;
  for (i = 0; i < m_M; i++){
    cout << m_a[i] << " ";
  }
  cout << endl;
  cout << endl;

  cout << "b" << endl;
  for (i = 0; i < m_N; i++){
    cout << m_b[i] << " ";
  }
  cout << endl;
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
