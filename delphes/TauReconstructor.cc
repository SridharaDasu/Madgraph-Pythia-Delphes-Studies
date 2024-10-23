#include <iostream>

#include "modules/TauReconstructor.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

#include "TLorentzVector.h"

TauReconstructor::TauReconstructor() : fMinTauSeedPT(5.), fMaxTauSeedEta(2.5), fMaxTauIsolDeltaR(0.5), fMaxTauCoreDeltaR(0.3) {
}

TauReconstructor::~TauReconstructor() {return;}

void TauReconstructor::Init() {
  fMinTauSeedPT = GetDouble("MinTauSeedPT", 5.0);
  fMaxTauSeedEta = GetDouble("MaxTauSeedEta", 2.5);
  fMaxTauIsolDeltaR = GetDouble("MaxTauIsolDeltaR", 0.5);
  fMaxTauCoreDeltaR = GetDouble("MaxTauCoreDeltaR", 0.3);
  fchadronInputArray = ImportArray(GetString("InputArray", "HCal/eflowTracks"));
  fchadronIterator = fchadronInputArray->MakeIterator();
  fphotonInputArray = ImportArray(GetString("InputArray", "ECal/eflowPhotons"));
  fphotonIterator = fphotonInputArray->MakeIterator();
  fnhadronInputArray = ImportArray(GetString("InputArray", "HCal/eflowNeutralHadrons"));
  fnhadronIterator = fnhadronInputArray->MakeIterator();
  fOutputArray = ExportArray(GetString("OutputArray", "taus"));
}

void TauReconstructor::Process() {
  DelphesFactory *factory = GetFactory();
  std::vector <Candidate *> trackSeeds; 
  fchadronIterator->Reset();
  while(Candidate *track = static_cast<Candidate *>(fchadronIterator->Next())) {
    if(track->Momentum.Pt() < fMinTauSeedPT || fabs(track->Momentum.Eta()) > fMaxTauSeedEta)
      continue;
    else
      trackSeeds.push_back(track);
  }
  std::vector <Tau> recTaus;
  TLorentzVector recoTau;
  for (UInt_t i = 0; i < trackSeeds.size(); i++) {
    Candidate *track = trackSeeds[i];
    fchadronIterator->Reset();
    bool selected = true;
    recoTau.SetPtEtaPhiM(track->Momentum.Pt(), track->Momentum.Eta(), track->Momentum.Phi(), track->Momentum.M());
    double isolation = 0;
    int charge = track->Charge;
    int nProngs = 1;
    while(Candidate* fragment = static_cast<Candidate *>(fchadronIterator->Next())) {
      // if (entry < 10) cout << "Max: EFTrack[" << t << "] = (" << track->PT << ", " << track->Eta << ", " << track->Phi << ", " << track->Mass << ")";
      if (track != fragment) {
	double deltaR = recoTau.DeltaR(fragment->Momentum);
	if(deltaR < fMaxTauIsolDeltaR) {
	  if (fragment->Momentum.Pt() > track->Momentum.Pt()) {
	    // if (entry < 10) cout << " - Deselected" << endl;
	    selected = false;
	    break;
	  }
	  else {
	    if (deltaR < fMaxTauCoreDeltaR) {
	      recoTau += fragment->Momentum;
	      // if (entry < 10) cout << endl << "Frg: EFTrack[" << a << "] = (" << fragment->PT << ", " << fragment->Eta << ", " << fragment->Phi << ", " << fragment->Mass << ")";
	      charge += fragment->Charge;
	      nProngs += 1;
	    }
	    else {
	      isolation += fragment->Momentum.Pt();
	    }
	  }
	}
      }
    }
    // Add in photons (mostly from pizero decays) if within 0.3, and for isolation if between 0.3-0.5
    int nPhotons = 0;
    fphotonIterator->Reset();
    while(Candidate* fragment = static_cast<Candidate *>(fphotonIterator->Next())) {
      double deltaR = recoTau.DeltaR(fragment->Momentum);
      if(deltaR < fMaxTauIsolDeltaR) {
	if (deltaR < fMaxTauCoreDeltaR) {
	  recoTau += fragment->Momentum;
	  nPhotons++;
	}
	else {
	  isolation += fragment->Momentum.Pt();
	}
      }
    }
    // Consider non-pizero remnant neutral hadrons only for isolation
    int nNHadrons = 0;
    fnhadronIterator->Reset();
    while(Candidate* fragment = static_cast<Candidate *>(fphotonIterator->Next())) {
      double deltaR = recoTau.DeltaR(fragment->Momentum);
      if(deltaR < fMaxTauIsolDeltaR) {
	if (deltaR < fMaxTauCoreDeltaR) {
	  recoTau += fragment->Momentum;
	  nNHadrons++;
	}
	else {
	  isolation += fragment->Momentum.Pt();
	}
      }
    }
    if (selected) {
      // if (entry < 10) cout << endl << "SelRecoTau = (" << recoTau.Pt() << ", " << recoTau.Eta() << ", " << recoTau.Phi() << ", "<< recoTau.M() << ")" << endl;
      if (nProngs <= 5) {
	Tau tau;
	tau.PT = recoTau.Pt();
	tau.Eta = recoTau.Eta();
	tau.Phi = recoTau.Phi();
	tau.Mass = recoTau.M();
	// tau.T = track->T;
	tau.Charge = charge;
	tau.nProngs = nProngs;
	tau.nPhotons = nPhotons;
	tau.nNHadrons = nNHadrons;
	tau.isolation = isolation;
	recTaus.push_back(tau);
	Candidate* candidate = factory->NewCandidate();
	candidate->Momentum.SetPtEtaPhiM(recoTau.Pt(), recoTau.Eta(), recoTau.Phi(), recoTau.M());
	candidate->Mass = recoTau.M();
	// candidate->T = track->T;
	candidate->Charge = charge;
	candidate->NCharged = nProngs;
	candidate->NNeutrals = nPhotons;
	candidate->Nclusters = nNHadrons; // Steal this variable
	candidate->IsolationVar = isolation;
	fOutputArray->Add(candidate);
      }
    } // Selected objects
  } // All seeds above 5 GeV
  if (recTaus.size() > 0) {
    std::cout << "N recoTaus = " << recTaus.size() << std::endl;
    for (UInt_t i = 0; i < recTaus.size(); i++) {
      Tau tau = recTaus[i];
      std::cout << tau.PT << tau.Eta << tau.Phi << tau.Mass << tau.Charge << tau.nProngs << tau.nPhotons << tau.nNHadrons << tau.isolation << std::endl;
    }
  }
}

void TauReconstructor::Finish() {
}
