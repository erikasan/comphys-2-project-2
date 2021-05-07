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
  int M = (m_system->getNumberOfParticles()) * (m_system->getNumberOfDimensions());

  m_M = M;
  m_N = N;

  double W[M][N];
  double a[M];
  double b[N];

  mt19937_64 generator;
  normal_distribution<double> distribution(0, 0.001);

  size_t i, j;

  for (i = 0; i < M; i++){
    a[i] = distribution(generator);
    for (j = 0; j < N; j++){
      W[i][j] = distribution(generator);
    }
  }

  for (j = 0; j < N; j++){
    b[j] = distribution(generator);
  }

  *m_W = &W[0][0];
  m_a  = &a[0];
  m_b  = &b[0];
  return;
}

void Gaussian_Binary::test_weights_biases()
{
  size_t i, j;
  for (i = 0; i < m_M; i++){
    for (j = 0; j < m_N; j++){
      cout << i << " " << j << " " << m_W[i][j];
    }
    cout << endl;
  }

  return;
}

double Gaussian_Binary::evaluate(std::vector<class Particle*> particles)
{
  return 0;
}

double Gaussian_Binary::evaluate(std::vector<class Particle*> particles, int particle_id)
{
  return 0;
}

double Gaussian_Binary::computeDoubleDerivative(std::vector<class Particle*> particles)
{
  return 0;
}

std::vector<double> Gaussian_Binary::quantumForce(std::vector<class Particle*> particles, int particle_id)
{
  std::vector<double> qForce = std::vector<double>();
  return qForce;
}

std::vector<double> Gaussian_Binary::quantumForce(std::vector<class Particle*> particles)
{
  std::vector<double> qForce = std::vector<double>();
  return qForce;
}

void Gaussian_Binary::updateDistances(std::vector<class Particle*> particles,int particle_id)
{
  return;
}

void Gaussian_Binary::initiateDistances(std::vector<class Particle*> particles)
{
  return;
}
void Gaussian_Binary::sample(std::vector<class Particle*> particles, double localEnergy)
{
  return;
}
void Gaussian_Binary::computeAverages(double steps)
{
  return;
}
void Gaussian_Binary::gradientDescent()
{
  return;
}