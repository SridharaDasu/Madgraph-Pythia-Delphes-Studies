/*

  A macro to identify generator level taus, their visible part and reconstruct
  taus using EFlow objects in a Delphes file.

  Taus object contains its P4, its highest track P4 and summary of other tau remnants

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
#include "TLorentzVector.h"
#include "Histograms.h"
#include "Tau.h"

using namespace std;

void ZTauTauAnalysis(const char *inputFile)
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
  strcat(outputFile, "-histos.root");

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

  histograms.Define("TauPairMass", "Visible oppositely charged tau-pair Mass", 100, 0.0, 200.0);
  histograms.Define("TauPairPT", "Visible oppositely charged tau-pair P_{T}", 100, 0.0, 100.0);
  histograms.Define("TauPairEta", "Visible oppositely charged tau-pair eta", 100, 0.0, 10.0);
  histograms.Define("TauPairPhi", "Visible oppositely charged tau-pair phi", 100, -7.0, 7.0);

  vector<Tau> genTaus;
  vector<Tau> visTaus;
  vector<Tau> recTaus;
  vector< pair<Tau, Tau> > mchTaus;

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
	Tau &recTau = recTaus[r];
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
	if (recTau.nProngs != 1 && recTau.nProngs != 3) continue;
	for (int r2 = r + 1; r2 < recTaus.size(); r2++) {
	  Tau &recTau2 = recTaus[r2];
	  if (recTau2.nProngs != 1 && recTau2.nProngs != 3) continue;
	  if (recTau.charge * recTau2.charge == -1) { // Oppositely charged tau pair
	    TLorentzVector tauPair = recTau.P4 + recTau2.P4;
	    histograms.Fill("TauPairMass", tauPair.M());
	    histograms.Fill("TauPairPT", tauPair.Pt());
	    histograms.Fill("TauPairEta", tauPair.Eta());
	    histograms.Fill("TauPairPhi", tauPair.Phi());
	  }
	}
      }
    }
    
    // Matched taus - these should be there if the algorithm is any good
    if (makeMchTaus(visTaus, recTaus, mchTaus)) {
      for (int m = 0; m < mchTaus.size(); m++) {
	Tau &visTau = mchTaus[m].first;
	Tau &recTau = mchTaus[m].second;
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
      
  } // Event loop

}
