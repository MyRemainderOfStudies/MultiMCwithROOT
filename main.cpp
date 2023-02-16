/////////////////////////////
/// Generate by ChatGPT
/////////////////////////////

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "omp.h"

int main() {
  // Create a ROOT file for writing the TTree
  TFile *file = new TFile("output.root", "RECREATE");

  // Create a TTree with the desired branch names and types
  TTree *tree = new TTree("tree_name", "Tree Title");
  float energy, px, py, pz;
  tree->Branch("energy", &energy, "energy/F");
  tree->Branch("px", &px, "px/F");
  tree->Branch("py", &py, "py/F");
  tree->Branch("pz", &pz, "pz/F");

// Generate and write particle data to the TTree in parallel
#pragma omp parallel for
  for (int i = 0; i < 10000; ++i) {
    // Initialize a new random seed for each thread
    unsigned int seed = omp_get_thread_num() + 1;

    // Generate particle energy and momentum with the thread-specific random
    // seed
    energy = /* generate particle energy with random seed */;
    px = /* generate particle momentum in x direction with random seed */;
    py = /* generate particle momentum in y direction with random seed */;
    pz = /* generate particle momentum in z direction with random seed */;

// Atomically write particle data to the TTree to prevent race conditions
#pragma omp atomic write
    tree->Fill();
  }

  // Write the TTree and close the file
  file->Write();
  file->Close();

  std::cout << "Successfully wrote TTree to output.root" << std::endl;

  return 0;
}
