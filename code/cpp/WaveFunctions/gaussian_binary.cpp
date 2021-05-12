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

Gaussian_Binary::Gaussian_Binary(System *system, int N, double sigma) : WaveFunction(system)
{
  m_M = (m_system->getNumberOfParticles()) * (m_system->getNumberOfDimensions());
  m_N = N;

  m_sigma = sigma;

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

  size_t i, j;

  vec x(m_M);

  int numberOfParticles  = m_system->getNumberOfParticles();
  int numberOfDimensions = m_system->getNumberOfDimensions();

  std::vector<double> position;

  // Convert position vector to armadillo column vector
  for (i = 0; i < numberOfParticles; i++){
    position = particles[i]->getPosition();
    for (j = 0; j < numberOfDimensions; j++){
      x(j) = position[j];
    }
  }

  return -0.5*(Gaussian_Binary::gradientTerm(x) + Gaussian_Binary::laplacianTerm(x));
}

double Gaussian_Binary::gradientTerm(vec x)
{
  double sigma2 = m_sigma*m_sigma;
  vec e(m_N);
  vec gradlnpsi(m_M);

  e = 1/(1 + exp(-(m_b + 1/sigma2*(m_W.t()*x))));

  gradlnpsi = 1/sigma2*(m_a - x + m_W*e);

  return dot(gradlnpsi, gradlnpsi);
}

double Gaussian_Binary::laplacianTerm(vec x)
{
  double sigma2 = m_sigma*m_sigma;

  vec laplacelnpsi(m_M); // For convenience. The actual Laplacian will be the sum over the elements.

  vec g(m_N);
  vec h(m_N);

  g = exp(-(m_b + 1/sigma2*(m_W.t()*x)));

  h = g/square(1 + g);

  laplacelnpsi = -1/sigma2 + 1/(sigma2*sigma2)*(square(m_W)*h);

  return sum(laplacelnpsi);
}


std::vector<double> Gaussian_Binary::quantumForce(std::vector<class Particle*> particles)
{
  std::vector<double> qForce = std::vector<double>();
  return qForce;
}
