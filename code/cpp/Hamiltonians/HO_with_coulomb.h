#pragma once
#include "harmonicoscillator.h"
#include <vector>

class HO_with_Coulomb : public HarmonicOscillator {
public:
  HO_with_Coulomb(System* system);
  double computeLocalPotentialEnergy(std::vector<Particle*> particles);
};
