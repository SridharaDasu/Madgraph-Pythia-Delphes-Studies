/*

  A macro to identify generator level taus, their visible part and reconstruct
  taus using EFlow objects in a Delphes file.

  TauLeptons object contains its P4, its highest track P4 and summary of other tau remnants

  Some histograms are made for checking this out in the tauMaker(const char* delphesRootFile) 

  root -l tauMaker.C'("delphes_output.root")'

*/

//------------------------------------------------------------------------------

R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include <string.h>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include "TLorentzVector.h"
#include "Histograms.h"
#include "TauLepton.h"
#include <stdlib.h>

using namespace std;

void WWWAnalysis(const char *inputFile)
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  chain.Print();

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchParticles = treeReader->UseBranch("Particle");
  TClonesArray *branchEFTracks = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFPhotons = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFNHadrons = treeReader->UseBranch("EFlowNeutralHadron");

  // Book histograms

  char outputFile[1024];
  strcpy(outputFile, inputFile);
  *strstr(outputFile, ".root") = '\0';
  strcat(outputFile, "-HtoTauTau.root");

  Histograms histograms(outputFile);
  histograms.Define("GenTauPT", "Generated tau P_{T}", 100, 0.0, 100.0);
  histograms.Define("VisTauPT", "Visible tau P_{T}", 100, 0.0, 100.0);
  histograms.Define("MxTTauPT", "Maximum tau remnant track P_{T}", 100, 0.0, 100.0);
  histograms.Define("TFrTauPT", "Maximum tau remnant track P_{T} / Tau P_{T}", 100, 0.0, 2.0);
  histograms.Define("RecTauPT", "Reconstructed tau P_{T}", 100, 0.0, 100.0);
  histograms.Define("RecTauMxTrkPT", "Reconstructed tau max track P_{T}", 100, 0.0, 100.0);
  histograms.Define("RecTaxPT", "Reconstructed wrong sign tau-like jet P_{T}", 100, 0.0, 100.0);
  histograms.Define("RecTauCh", "Reconstructed tau charge", 41, -20, 20.0);
  histograms.Define("RecTauIR", "Reconstructed tau isolation ratio", 100, 0.0, 10.0);
  histograms.Define("RecTauNP", "Reconstructed tau EFTrack count (n prongs)", 30, 0.0, 30.0);
  histograms.Define("RecTauNG", "Reconstructed tau EFPhoton count", 30, 0.0, 30.0);
  histograms.Define("RecTauNN", "Reconstructed tau EFNHadron count", 30, 0.0, 30.0);
  histograms.Define("RecTauMs", "Reconstructed tau mass", 100, 0.0, 10.0);
  histograms.Define("DeltaR", "DeltaR between max track and other fragments in tau", 100, 0.0, 1.0);
  histograms.Define("MchTauPT", "Visible tau P_{T} for good match", 100, 0.0, 100.0);
  histograms.Define("MchTauCh", "Matched tau charge", 21, -10.0, 10.0);
  histograms.Define("MchTauNPr", "Matched tau nProngs", 30, 0.0, 30.0);
  histograms.Define("MchTauNPh", "Matched tau nPhotons", 30, 0.0, 30.0);
  histograms.Define("MchTauNNH", "Matched tau nNHadrons", 30, 0.0, 30.0);
  histograms.Define("MchTauIso", "Matched tau isolation ratio", 100, 0.0, 1.0);
  histograms.Define("MchPTRes", "Matched tau reco vs vis PT resolution", 100, -1.0, 1.0);
  histograms.Define("MchEtaRes", "Matched tau reco vs vis Eta resolution", 100, -1.0, 1.0);
  histograms.Define("MchPhiRes", "Matched tau reco vs vis Phi resolution", 100, -1.0, 1.0);
  histograms.Define("WrgTauPT", "Visible tau P_{T} for wrong match", 100, 0.0, 100.0);
  histograms.Define("WrgTauCh", "Wrongly matched tau charge", 21, -10.0, 10.0);
  histograms.Define("WrgTauNPr", "Matched tau nProngs", 30, 0.0, 30.0);
  histograms.Define("WrgTauNPh", "Matched tau nPhotons", 30, 0.0, 30.0);
  histograms.Define("WrgTauNNH", "Matched tau nNHadrons", 30, 0.0, 30.0);
  histograms.Define("WrgTauIso", "Wrongly matched tau isolation ratio", 100, 0.0, 1.0);
  histograms.Define("WrgPTRes", "Wrongly matched tau reco vs vis PT resolution", 100, -1.0, 1.0);
  histograms.Define("WrgEtaRes", "Wrongly matched tau reco vs vis Eta resolution", 100, -1.0, 1.0);
  histograms.Define("WrgPhiRes", "Wrongly matched tau reco vs vis Phi resolution", 100, -1.0, 1.0);

  std::array<std::pair<string, string>, 5> histTypes;
  histTypes[0] = {string("GenTau"), string("Generated Tau")};
  histTypes[1] = {string("VisTau"), string("Visible Tau")};
  histTypes[2] = {string("RecTau"), string("Reconstructed Tau")};
  histTypes[3] = {string("MchTau"), string("Correctly Matched Reconstructed Tau")};
  histTypes[4] = {string("WrgTau"), string("Wrongly Matched Reco Tau")};

  size_t NPTPlots = 4;
  for(Int_t i = 0; i < NPTPlots; i++) {
    for (auto& histType : histTypes) {
      string histName = histType.first + std::to_string(i) + std::string("PT");
      string histTitl = histType.second + std::to_string(i) + std::string(" P_{T}");
      histograms.Define(histName, histTitl, 100, 0.0, 100.0);
    }
  }
  
  vector<TauLepton> genTaus;
  vector<TauLepton> visTaus;
  vector<TauLepton> recTaus;
  vector< pair<TauLepton, TauLepton> > mchTaus;

  cout << "Number of entries = " << numberOfEntries << endl;

  // Loop over all events

  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    
    // There should always be a pair of tau+ and tau- in the event from the a-decay
    // Sometimes there may ba a second pair of tau+ and tau- from the decay of a Z-boson
    // Collect all the particles
    
    if (makeGenTaus(branchParticles, genTaus, visTaus)) {
      for (Int_t i = 0; i < genTaus.size(); i++) {
	histograms.Fill("GenTauPT", genTaus[i].P4.Pt());
	histograms.Fill("MxTTauPT", genTaus[i].MaxTrackP4.Pt());
	histograms.Fill("TFrTauPT", genTaus[i].MaxTrackP4.Pt() / genTaus[i].P4.Pt());
	histograms.Fill("VisTauPT", visTaus[i].P4.Pt());
      }
    }

    // Reconstructed taus - these could be there even if there are no real taus in the event
    
    if (makeRecTaus(branchEFTracks, branchEFPhotons, branchEFNHadrons, recTaus)) {
      for (int r = 0; r < recTaus.size(); r++) {
	TauLepton &recTau = recTaus[r];
	double deltaR = recTau.P4.DeltaR(recTau.MaxTrackP4);
	histograms.Fill("DeltaR", deltaR);
	histograms.Fill("RecTauMxTrkPT", recTau.MaxTrackP4.Pt());
	histograms.Fill("RecTauPT", recTau.P4.Pt());
	histograms.Fill("RecTauNP", recTau.nProngs);
	histograms.Fill("RecTauNG", recTau.nPhotons);
	histograms.Fill("RecTauNN", recTau.nNHadrons);
	histograms.Fill("RecTauIR", recTau.isolation / recTau.P4.Pt());
	histograms.Fill("RecTauCh", recTau.charge);
	histograms.Fill("RecTauMs", recTau.P4.M());
      }
    }
    
    // Matched taus - these should be there if the algorithm is any good
    if (makeMchTaus(visTaus, recTaus, mchTaus)) {
      for (int m = 0; m < mchTaus.size(); m++) {
	TauLepton &visTau = mchTaus[m].first;
	TauLepton &recTau = mchTaus[m].second;
	if(visTau.charge == recTau.charge) {
	  histograms.Fill("MchTauPT", visTau.P4.Pt());
	  histograms.Fill("MchTauCh", recTau.charge);
	  histograms.Fill("MchTauIso", recTau.isolation / recTau.P4.Pt());
	  histograms.Fill("MchTauNPr", recTau.nProngs);
	  histograms.Fill("MchTauNPh", recTau.nPhotons);
	  histograms.Fill("MchTauNNH", recTau.nNHadrons);
	  histograms.Fill("MchPTRes", (recTau.P4.Pt() - visTau.P4.Pt()) / visTau.P4.Pt());
	  histograms.Fill("MchEtaRes", (recTau.P4.Eta() - visTau.P4.Eta()) / visTau.P4.Eta());
	  histograms.Fill("MchPhiRes", (recTau.P4.Phi() - visTau.P4.Phi()) / visTau.P4.Phi());
	}
	else {
	  histograms.Fill("WrgTauPT", visTau.P4.Pt());
	  histograms.Fill("WrgTauCh", recTau.charge);
	  histograms.Fill("WrgTauIso", recTau.isolation / recTau.P4.Pt());
	  histograms.Fill("WrgTauNPr", recTau.nProngs);
	  histograms.Fill("WrgTauNPh", recTau.nPhotons);
	  histograms.Fill("WrgTauNNH", recTau.nNHadrons);
	  histograms.Fill("WrgPTRes", (recTau.P4.Pt() - visTau.P4.Pt()) / visTau.P4.Pt());
	  histograms.Fill("WrgEtaRes", (recTau.P4.Eta() - visTau.P4.Eta()) / visTau.P4.Eta());
	  histograms.Fill("WrgPhiRes", (recTau.P4.Phi() - visTau.P4.Phi()) / visTau.P4.Phi());
	}
      }
    }

    for(Int_t i = 0; i < std::min(NPTPlots, genTaus.size()); i++) {
      string histName = "GenTau" + std::to_string(i) + std::string("PT");
      histograms.Fill(histName.c_str(), genTaus[i].P4.Pt());
    }

    for(Int_t i = 0; i < std::min(NPTPlots, visTaus.size()); i++) {
      string histName = "VisTau" + std::to_string(i)  + std::string("PT");
      histograms.Fill(histName.c_str(), visTaus[i].P4.Pt());
    }

    for(Int_t i = 0; i < std::min(NPTPlots, recTaus.size()); i++) {
      string histName = "RecTau" + std::to_string(i)  + std::string("PT");
      histograms.Fill(histName.c_str(), recTaus[i].P4.Pt());
    }

    for(Int_t i = 0; i < std::min(NPTPlots, mchTaus.size()); i++) {
      TauLepton &visTau = mchTaus[i].first;
      TauLepton &recTau = mchTaus[i].second;
      cout << &visTau << &recTau << endl;
      if(visTau.charge == recTau.charge) {
	string histName = "MchTau" + std::to_string(i)  + std::string("PT");
	histograms.Fill(histName.c_str(), recTau.P4.Pt());
	cout << "Did mch taus 2" << endl;
      }
      else {
	string histName = "WrgTau" + std::to_string(i)  + std::string("PT");
	histograms.Fill(histName.c_str(), recTau.P4.Pt());
	cout << "Did wrg taus 2" << endl;
      }
    }

  } // Event loop

}
