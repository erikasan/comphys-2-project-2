#include "wavefunction.h"
#include <iostream>
using namespace std;

WaveFunction::WaveFunction(System* system) {
  m_system = system;
  cout << "Hello" << endl;
}

void WaveFunction::setTolerance(double tol){
  m_tol = tol;
}

void WaveFunction::setLearningRate(double learningRate){
  m_learningRate = learningRate;
}
