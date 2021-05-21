#include "HO_with_coulomb.h"
#include "harmonicoscillator.h"
#include <iostream>
#include <cmath>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using namespace std;

HO_with_Coulomb::HO_with_Coulomb(System* system, double omega) : HarmonicOscillator(system, omega) {
}

double HO_with_Coulomb::computeLocalPotentialEnergy(std::vector<Particle*> particles){
  double HO_PE      = HarmonicOscillator::computeLocalPotentialEnergy(particles);
  double Coulomb_PE = 0;

  double term = 0;

  size_t i, j, k;
  int numberOfParticles  = HarmonicOscillator::m_system->getNumberOfParticles();
  int numberOfDimensions = HarmonicOscillator::m_system->getNumberOfDimensions();

  std::vector<double> r_i, r_j;

  for (i = 0; i < numberOfParticles; i++){
    r_i = particles[i]->getPosition();
    for (j = i + 1; j < numberOfParticles; j++){
      r_j = particles[j]->getPosition();
      for (k = 0; k < numberOfDimensions; k++){
        term += (r_i[k] - r_j[k])*(r_i[k] - r_j[k]);
      }
      Coulomb_PE += 1./sqrt(term);
      term = 0;
    }
  }

  return HO_PE + Coulomb_PE;
}
