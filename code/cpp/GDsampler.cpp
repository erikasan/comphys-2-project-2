#include "GDsampler.h"
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"

#include <iostream>
using namespace std;

void GDsampler::gdsampler(std::vector<class Particle*> particles, double localEnergy){
  m_system->getWaveFunction()->sample(particles, localEnergy);
  return;
}
