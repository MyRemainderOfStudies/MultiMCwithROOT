/////////////////////////////
/// Generate by ChatGPT
/// Edited 2023/2/16
/////////////////////////////

#include <TFile.h>
#include <TMath.h>
#include <TTree.h>

#include <iostream>
#include <random>

#include "omp.h"

const int nThreads = omp_get_num_threads();
const int nParticles = 1e+5;
const int nParThreads = nParticles / nThreads;
const int nMod = nParticles % nThreads;
const double minTheta = 0.0;
const double maxTheta = TMath::Pi();
const double minPhi = 0.0;
const double maxPhi = TMath::TwoPi();
const double minEta = 0.0;
const double maxEta = 1.0;
std::uniform_real_distribution<double> distTheta(minTheta, maxTheta);
std::uniform_real_distribution<double> distPhi(minPhi, maxPhi);
std::uniform_real_distribution<double> distEta(minEta, maxEta);
std::uniform_real_distribution<double> distEval(0.0, 1.0);
TFile *file = new TFile("output.root", "RECREATE");
// Create a TTree with the desired branch names and types
TTree *tree = new TTree("tree_name", "Tree Title");
// double px, py, pz;
double theta, eta, phi;

void evalFuncEtaTheta(double *result, double _theta, double _eta) {
  *result = 2 * pow(_eta, 2) * ((3 - 2 * _eta) + (2 * _eta - 1) * cos(_theta)) *
            sin(_theta);
}

void montecarlo(int N, int seed) {
  int nInside = 0;
  double theta, eta, phi, eval, res, fNorm = 1. / 4.;
  std::cout << "Seed: " << seed << std::endl;
  std::mt19937 randgen(seed);
  while (nInside < N) {
    theta = distTheta(randgen);
    phi = distPhi(randgen);
    eta = distEta(randgen);
    eval = distEval(randgen);
    evalFuncEtaTheta(&res, theta, eta);
    res *= fNorm;
    if (eval <= res) {
      // px = eta * sin(theta) * cos(phi);
      // py = eta * sin(theta) * sin(phi);
      // pz = eta * cos(theta);
// Atomically write particle data to the TTree to prevent race condition
#pragma omp atomic write
      tree->Fill();
      nInside += 1;
    }
  }
}

int main() {
  // tree->Branch("px", &px, "px/D");
  // tree->Branch("py", &py, "py/D");
  // tree->Branch("pz", &pz, "pz/D");
  tree->Branch("theta", &theta, "theta/D");
  tree->Branch("phi", &phi, "phi/D");
  tree->Branch("eta", &eta, "eta/D");

  std::cout << "threds num: " << nThreads << std::endl;
#pragma omp parallel for
  for (int i = 0; i < nThreads; ++i) {
    if (i == nThreads - 1) {
      montecarlo(nParThreads + nMod, i + 1);
      continue;
    }
    montecarlo(nParThreads, i + 1);
  }

  // Write the TTree and close the file
  file->Write();
  file->Close();

  std::cout << "Successfully wrote TTree to output.root" << std::endl;

  return 0;
}
