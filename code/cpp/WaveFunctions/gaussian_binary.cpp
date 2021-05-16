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

  m_sigma2 = sigma*sigma;

  m_W.set_size(m_M, m_N);
  m_a.set_size(m_M);
  m_b.set_size(m_N);

  m_av_grad_a.zeros(m_M);
  m_av_grad_b.zeros(m_N);
  m_av_grad_W.zeros(m_M, m_N);
  m_av_local_energy_grad_a.zeros(m_M);
  m_av_local_energy_grad_b.zeros(m_N);
  m_av_local_energy_grad_W.zeros(m_M, m_N);

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

}

double Gaussian_Binary::evaluate(std::vector<class Particle*> particles)
{
  vec x = Gaussian_Binary::convertPositionToArmadillo(particles);

  vec x_minus_a = x - m_a;

  vec tmp = 1 + exp(m_b + 1/m_sigma2*(m_W.t()*x));

  return exp(-1/(2*m_sigma2)*dot(x_minus_a, x_minus_a))*prod(tmp);
}

vec Gaussian_Binary::convertPositionToArmadillo(std::vector<class Particle*> particles)
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
      x(i*numberOfDimensions + j) = position[j];
    }
  }
  return x;
}

double Gaussian_Binary::computeDoubleDerivative(std::vector<class Particle*> particles)
{
  vec x = Gaussian_Binary::convertPositionToArmadillo(particles);
  return Gaussian_Binary::gradientTerm(x) + Gaussian_Binary::laplacianTerm(x);
}

double Gaussian_Binary::gradientTerm(vec x)
{
  vec e(m_N);
  vec gradlnpsi(m_M);

  e = 1/(1 + exp(-(m_b + 1/m_sigma2*(m_W.t()*x))));

  gradlnpsi = 1/m_sigma2*(m_a - x + m_W*e);

  return dot(gradlnpsi, gradlnpsi);
}

double Gaussian_Binary::laplacianTerm(vec x)
{
  vec laplacelnpsi(m_M); // For convenience. The actual Laplacian will be the sum over the elements.

  vec g(m_N);
  vec h(m_N);

  g = exp(-(m_b + 1/m_sigma2*(m_W.t()*x)));

  h = g/square(1 + g);

  laplacelnpsi = -1/m_sigma2 + 1/(m_sigma2*m_sigma2)*(square(m_W)*h);

  return sum(laplacelnpsi);
}

vec Gaussian_Binary::grad_a(vec x)
{
  return (x - m_a)/m_sigma2;
}

vec Gaussian_Binary::localEnergygrad_a(vec x, double localEnergy)
{
  return localEnergy*grad_a(x);
}

vec Gaussian_Binary::grad_b(vec x)
{
  return 1/(1 + exp(-(m_b + 1/m_sigma2*(m_W.t()*x))));
}

vec Gaussian_Binary::localEnergygrad_b(vec x, double localEnergy)
{
  return localEnergy*grad_b(x);
}

mat Gaussian_Binary::grad_W(vec x)
{
  return x*(grad_b(x).t())/m_sigma2;
}

mat Gaussian_Binary::localEnergygrad_W(vec x, double localEnergy)
{
  return localEnergy*grad_W(x);
}

void Gaussian_Binary::sample(std::vector<class Particle*> particles, double localEnergy)
{
  vec x = Gaussian_Binary::convertPositionToArmadillo(particles);

  m_av_grad_a += Gaussian_Binary::grad_a(x);
  m_av_grad_b += Gaussian_Binary::grad_b(x);
  m_av_grad_W += Gaussian_Binary::grad_W(x);
  m_av_local_energy_grad_a += Gaussian_Binary::localEnergygrad_a(x, localEnergy);
  m_av_local_energy_grad_b += Gaussian_Binary::localEnergygrad_b(x, localEnergy);
  m_av_local_energy_grad_W += Gaussian_Binary::localEnergygrad_W(x, localEnergy);
  return;
}

void Gaussian_Binary::computeAverages(double steps)
{
  m_av_grad_a /= steps;
  m_av_grad_b /= steps;
  m_av_grad_W /= steps;
  m_av_local_energy_grad_a /= steps;
  m_av_local_energy_grad_b /= steps;
  m_av_local_energy_grad_W /= steps;
  return;
}

void Gaussian_Binary::gradientDescent()
{
  double localEnergy = m_system->getSampler()->getEnergy();
  cout << localEnergy << endl;
  vec a_new, b_new;
  mat W_new;

  a_new = m_a - m_learningRate*2*(m_av_local_energy_grad_a - localEnergy*m_av_grad_a);
  b_new = m_b - m_learningRate*2*(m_av_local_energy_grad_b - localEnergy*m_av_grad_b);
  W_new = m_W - m_learningRate*2*(m_av_local_energy_grad_W - localEnergy*m_av_grad_W);



  if (approx_equal(a_new, m_a, "absdiff", m_tol) & approx_equal(b_new, m_b, "absdiff", m_tol) & approx_equal(W_new, m_W, "absdiff", m_tol)){
    m_system->stopGradientDescent();
  }

  m_a = a_new;
  m_b = b_new;
  m_W = W_new;

  m_av_grad_a.zeros();
  m_av_grad_b.zeros();
  m_av_grad_W.zeros();
  m_av_local_energy_grad_a.zeros();
  m_av_local_energy_grad_b.zeros();
  m_av_local_energy_grad_W.zeros();

  return;
}

std::vector<double> Gaussian_Binary::quantumForce(std::vector<class Particle*> particles)
{
  std::vector<double> qForce = std::vector<double>();
  return qForce;
}
