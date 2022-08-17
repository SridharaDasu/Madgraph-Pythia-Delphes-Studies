/*

  A macro to identify generator level taus, their visible part and reconstruct
  taus using EFlow objects in a Delphes file.

  Taus are stored as TLorentzVector objects, charge sorted.

  Reconstructed tau-like objects which are not charge = +/-1 are stored separately.

  Some histograms are made for checking this out.

  To Do: Make a better tau object and store some additional information, 
  e.g., number of charged prongs, electromagnetically decaying neutral hadron 
  fraction, etc. 

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

using namespace std;

class Histograms {
public:
  Histograms(string histFileName) {
    tFile = new TFile(histFileName.c_str(), "recreate");
  }
  virtual ~Histograms() {
    for(map<string, TH1F*>::iterator iter = histograms1D.begin(); iter != histograms1D.end(); ++iter)
      {
	iter->second->Write();
	TCanvas c;
	iter->second->Draw();
	string name = iter->first + ".png";
	c.SaveAs(name.c_str());
      }
    tFile->Close();
  }
  void Define(string name, string title, int nBins, double bLo, double bHi) {
    histograms1D[name] = new TH1F((const char*) name.c_str(), (const char*) title.c_str(), nBins, bLo, bHi);
  }
  void Fill(string name, double value, double weight = 1.0) {
    histograms1D[name]->Fill(value, weight);
  }
private:
  map<string, TH1F*> histograms1D;
  TFile *tFile;
};

class Tau {
public:
  // Constructor
  Tau(TLorentzVector tauP4, TLorentzVector trackP4, int c, int nPr, int nPh, int nNH, double iso = 0.0) {
    P4 = tauP4;
    MaxTrackP4 = trackP4;
    charge = c;
    nProngs = nPr;
    nPhotons = nPh;
    nNHadrons = nNH;
    isolation = iso;
  }
  // Data
  TLorentzVector P4;
  TLorentzVector MaxTrackP4;
  int charge;
  int nProngs;
  int nPhotons;
  int nNHadrons;
  double isolation;
};

void tauMaker(const char *inputFile)
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
  TClonesArray *branchJet      = treeReader->UseBranch("Jet");

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
  histograms.Define("DeltaR", "Average DeltaR between max track and other fragments in tau", 100, 0.0, 1.0);
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

  vector<Tau> genTaus;
  vector<Tau> visTaus;
  vector<Tau> recTaus;
  vector<Tau> mchTaus;

  cout << "Number of entries = " << numberOfEntries << endl;

  // Loop over all events

  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    {
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);

      // There should always be a pair of tau+ and tau- in the event from the a-decay
      // Sometimes there may ba a second pair of tau+ and tau- from the decay of a Z-boson
      // Collect all the particles

      genTaus.clear();
      visTaus.clear();
      recTaus.clear();
      mchTaus.clear();

      for (Int_t p = 0; p < branchParticles->GetEntries(); ++p) {
	GenParticle *particle = (GenParticle*) branchParticles->At(p);
	if (abs(particle->PID) == 16) {
	  GenParticle *tau = 0;
	  GenParticle *mother = (GenParticle*) branchParticles->At(particle->M1);
	  if (mother != 0) {
	    if (abs(mother->PID) == 15) {
	      tau = mother;
	      GenParticle *grandmother = (GenParticle*) branchParticles->At(mother->M1);
	      if (grandmother != 0) {
		if (abs(grandmother->PID) == 15) {
		  tau = grandmother;
		  GenParticle *greatgrandmother = (GenParticle*) branchParticles->At(grandmother->M1);
		  if (greatgrandmother != 0) {
		    if(abs(greatgrandmother->PID) == 15) {
		      tau = greatgrandmother;
		      GenParticle *greatgreatgrandmother = (GenParticle*) branchParticles->At(grandmother->M1);
		      if (greatgreatgrandmother != 0) {
			if(abs(greatgreatgrandmother->PID) == 15) {
			  tau = greatgreatgrandmother;
			}
		      }
		    }
		  }
		}
	      }
	    }
	    else if (mother->PID == 23 || abs(mother->PID) == 24) {
	      // Z or W decay neutrino - just continue
	      continue;
	    }
	    else if (abs(mother->PID) > 410 || abs(mother->PID) < 600) {
	      // D and B meson decay products - just continue
	      continue;
	    }
	    else {
	      cerr << "Mother of a tau neutrino is not a tau! It is " << mother->PID << endl;
	      continue;
	    }
	  }
	  else {
	    cerr << "Mother of a tau neutrino is zero!!" << endl;
	    continue;
	  }
	  GenParticle *theMother = (GenParticle*) branchParticles->At(particle->M1);
	  double chargedDaughterPT = 0;
	  GenParticle *theDaughter = 0;
	  int nProngs = 0;
	  int nPhotons = 0;
	  int nNHadrons = 0;
	  for (Int_t d = theMother->D1; d <= theMother->D2; d++) {
	    GenParticle *daughter = (GenParticle*) branchParticles->At(d);
	    if (theDaughter == 0) theDaughter = daughter;
	    if ((abs(daughter->PID) == 211 || abs(daughter->PID) == 321 || abs(daughter->PID) == 11 || abs(daughter->PID) == 13)) {
	      nProngs++;
	      if(daughter->PT > chargedDaughterPT) {
		theDaughter = daughter;
		chargedDaughterPT = daughter->PT;
	      }
	    }
	    else if(abs(daughter->PID) == 22) {
	      nPhotons++;
	    }
	    else if(abs(daughter->PID) != 12 && abs(daughter->PID) != 14 && abs(daughter->PID) != 16) {
	      nNHadrons++; // Catch all
	    }
	  }
	  int charge = tau->Charge;
	  genTaus.push_back(Tau(tau->P4(), theDaughter->P4(), charge, nProngs, nPhotons, nNHadrons));
	  visTaus.push_back(Tau((tau->P4() - particle->P4()), theDaughter->P4(), charge, nProngs, nPhotons, nNHadrons));
	}
      }

      for (Int_t i = 0; i < genTaus.size(); i++) {
	histograms.Fill("GenTauPT", genTaus[i].P4.Pt());
        histograms.Fill("MxTTauPT", genTaus[i].MaxTrackP4.Pt());
        histograms.Fill("TFrTauPT", genTaus[i].MaxTrackP4.Pt() / genTaus[i].P4.Pt());
	histograms.Fill("VisTauPT", visTaus[i].P4.Pt());
      }

      // Reconstruct taus using EFlow tracks, photons and neutral hadrons
      // Consider every charged track above 5 GeV PT - label at as a potential tau remnant
      // Add in all tracks and clusters within 0.3 of the track
      // Compute sum charge of the tracks
      // Count the number of charged tracks
      // Calculate mass of the tracks
      // Compute isolation as the sum of tracks and clusters in DeltaR = 0.3-0.5 range

      TLorentzVector recoTau;
      TLorentzVector fragmentP4;
      for (Int_t t = 0; t < branchEFTracks->GetEntries(); ++t) {
	Track *track = (Track*) branchEFTracks->At(t);
	if (track->PT > 5.0) {
	  bool selected = true;
	  recoTau.SetPtEtaPhiM(track->PT, track->Eta, track->Phi, track->Mass);
	  double isolation = 0;
	  int charge = track->Charge;
	  int nProngs = 1;
	  double ptSum = 0;
	  double deltaRxPTSum = 0;
	  // if (entry < 10) cout << "Max: EFTrack[" << t << "] = (" << track->PT << ", " << track->Eta << ", " << track->Phi << ", " << track->Mass << ")";
	  for (Int_t a = 0; a < branchEFTracks->GetEntries(); ++a) {
	    if (a != t) {
	      Track *fragment = (Track*) branchEFTracks->At(a);
	      fragmentP4.SetPtEtaPhiM(fragment->PT, fragment->Eta, fragment->Phi, fragment->Mass);
	      double deltaR = recoTau.DeltaR(fragmentP4);
	      if(deltaR < 0.5) {
		if (fragment->PT > track->PT) {
		  // if (entry < 10) cout << " - Deselected" << endl;
		  selected = false;
		  break;
		}
		else {
		  ptSum += fragment->PT;
		  deltaRxPTSum += (deltaR * fragment->PT);
		  if (deltaR < 0.3) {
		    recoTau += fragmentP4;
		    // if (entry < 10) cout << endl << "Frg: EFTrack[" << a << "] = (" << fragment->PT << ", " << fragment->Eta << ", " << fragment->Phi << ", " << fragment->Mass << ")";
		    charge += fragment->Charge;
		    nProngs += 1;
		  }
		  else {
		    isolation += fragment->PT;
		  }
		}
	      }
	    }
	  }
	  // Add in photons (mostly from pizero decays) if within 0.3, and for isolation if between 0.3-0.5
	  int nPhotons = 0;
	  for (Int_t a = t + 1; a < branchEFPhotons->GetEntries(); ++a) {
	    Tower *fragment = (Tower*) branchEFPhotons->At(a);
	    fragmentP4.SetPtEtaPhiM(fragment->ET, fragment->Eta, fragment->Phi, 0.);
	    double deltaR = recoTau.DeltaR(fragmentP4);
	    if(deltaR < 0.5) {
	      ptSum += fragment->ET;
	      deltaRxPTSum += (deltaR * fragment->ET);
	      if (deltaR < 0.3) {
		recoTau += fragmentP4;
		nPhotons++;
	      }
	      else {
		isolation += fragment->ET;
	      }
	    }
	  }
	  // Consider non-pizero remnant neutral hadrons only for isolation
	  int nNHadrons = 0;
	  for (Int_t a = t + 1; a < branchEFNHadrons->GetEntries(); ++a) {
	    Tower *fragment = (Tower*) branchEFNHadrons->At(a);
	    fragmentP4.SetPtEtaPhiM(fragment->ET, fragment->Eta, fragment->Phi, 0.);
	    double deltaR = recoTau.DeltaR(fragmentP4);
	    if(deltaR < 0.5) {
	      ptSum += fragment->ET;
	      deltaRxPTSum += (deltaR * fragment->ET);
	      if (deltaR < 0.3) {
		recoTau += fragmentP4;
		nNHadrons++;
	      }
	      else {
		isolation += fragment->ET;
	      }
	    }
	  }
	  if (selected) {
	    // if (entry < 10) cout << endl << "SelRecoTau = (" << recoTau.Pt() << ", " << recoTau.Eta() << ", " << recoTau.Phi() << ", "<< recoTau.M() << ")" << endl;
	    if (ptSum > 0) {
	      double deltaR = deltaRxPTSum / ptSum;
	      histograms.Fill("DeltaR", deltaR);
	    }
	    histograms.Fill("RecTauMxTrkPT", track->PT);
	    histograms.Fill("RecTauPT", recoTau.Pt());
	    histograms.Fill("RecTauNP", nProngs);
	    histograms.Fill("RecTauNG", nPhotons);
	    histograms.Fill("RecTauNN", nNHadrons);
	    if (nProngs <= 3) {
	      histograms.Fill("RecTauIR", isolation / recoTau.Pt());
	      histograms.Fill("RecTauCh", charge);
	      histograms.Fill("RecTauMs", recoTau.M());
	      recTaus.push_back(Tau(recoTau, track->P4(), charge, nProngs, nPhotons, nNHadrons, isolation));
	    }
	  } // Selected objects
	} // Those above 5 GeV
      } // All tracks

      for (Int_t v = 0; v < visTaus.size(); v++) {
	// Find correct sign match
	for (Int_t r = 0; r < recTaus.size(); r++) {
	  double deltaR = visTaus[v].P4.DeltaR(recTaus[r].P4);
	  if (deltaR < 0.3) {
	    if(visTaus[v].charge == recTaus[r].charge) {
	      mchTaus.push_back(recTaus[r]);
	      histograms.Fill("MchTauPT", visTaus[v].P4.Pt());
	      histograms.Fill("MchTauCh", recTaus[r].charge);
	      histograms.Fill("MchTauIso", recTaus[r].isolation / recTaus[r].P4.Pt());
	      histograms.Fill("MchTauNPr", recTaus[r].nProngs);
	      histograms.Fill("MchTauNPh", recTaus[r].nPhotons);
	      histograms.Fill("MchTauNNH", recTaus[r].nNHadrons);
	      histograms.Fill("MchPTRes", (recTaus[r].P4.Pt() - visTaus[v].P4.Pt()) / visTaus[v].P4.Pt());
	      histograms.Fill("MchEtaRes", (recTaus[r].P4.Eta() - visTaus[v].P4.Eta()) / visTaus[v].P4.Eta());
	      histograms.Fill("MchPhiRes", (recTaus[r].P4.Phi() - visTaus[v].P4.Phi()) / visTaus[v].P4.Phi());
	    }
	    else {
	      histograms.Fill("WrgTauPT", visTaus[v].P4.Pt());
	      histograms.Fill("WrgTauCh", recTaus[r].charge);
	      histograms.Fill("WrgTauIso", recTaus[r].isolation / recTaus[r].P4.Pt());
	      histograms.Fill("WrgTauNPr", recTaus[r].nProngs);
	      histograms.Fill("WrgTauNPh", recTaus[r].nPhotons);
	      histograms.Fill("WrgTauNNH", recTaus[r].nNHadrons);
	      histograms.Fill("WrgPTRes", (recTaus[r].P4.Pt() - visTaus[v].P4.Pt()) / visTaus[v].P4.Pt());
	      histograms.Fill("WrgEtaRes", (recTaus[r].P4.Eta() - visTaus[v].P4.Eta()) / visTaus[v].P4.Eta());
	      histograms.Fill("WrgPhiRes", (recTaus[r].P4.Phi() - visTaus[v].P4.Phi()) / visTaus[v].P4.Phi());
	    }
	    break;
	  }
	}
      }

    } // Event loop
}
