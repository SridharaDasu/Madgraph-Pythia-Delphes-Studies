#ifndef TauLepton__h
#define TauLepton__h
/*

  Taus object contains its P4, its highest track P4 and 
  summary of other tau remnants

  Functions to identify generator level taus, their visible part and reconstruct
  taus using EFlow objects using Delphes objects stored in root format

  #include "TauLepton.h" 
  in your root macro to make use of Tau objects

*/

#include <vector>
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "classes/DelphesClasses.h"

using namespace std;

class TauLepton {
public:
  // Constructor
  TauLepton(TLorentzVector tauP4, TLorentzVector trackP4, int c, int nPr, int nPh, int nNH, double iso = 0.0) {
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

bool makeGenTaus(TClonesArray *branchParticles, vector<TauLepton> &genTaus, vector<TauLepton> &visTaus) {
  genTaus.clear();
  visTaus.clear();
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
      genTaus.push_back(TauLepton(tau->P4(), theDaughter->P4(), charge, nProngs, nPhotons, nNHadrons));
      visTaus.push_back(TauLepton((tau->P4() - particle->P4()), theDaughter->P4(), charge, nProngs, nPhotons, nNHadrons));
    }
  }
  return true;
}

bool makeRecTaus(TClonesArray* branchEFTracks, TClonesArray* branchEFPhotons, TClonesArray* branchEFNHadrons, vector<TauLepton> &recTaus) {
  // Reconstruct taus using EFlow tracks, photons and neutral hadrons
  // Consider every charged track above 5 GeV PT - label at as a potential tau remnant
  // Add in all tracks and clusters within 0.3 of the track
  // Compute sum charge of the tracks
  // Count the number of charged tracks
  // Calculate mass of the tracks
  // Compute isolation as the sum of tracks and clusters in DeltaR = 0.3-0.5 range
  recTaus.clear();
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
      for (Int_t a = 0; a < branchEFPhotons->GetEntries(); ++a) {
	Tower *fragment = (Tower*) branchEFPhotons->At(a);
	fragmentP4.SetPtEtaPhiM(fragment->ET, fragment->Eta, fragment->Phi, 0.);
	double deltaR = recoTau.DeltaR(fragmentP4);
	if(deltaR < 0.5) {
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
      for (Int_t a = 0; a < branchEFNHadrons->GetEntries(); ++a) {
	Tower *fragment = (Tower*) branchEFNHadrons->At(a);
	fragmentP4.SetPtEtaPhiM(fragment->ET, fragment->Eta, fragment->Phi, 0.);
	double deltaR = recoTau.DeltaR(fragmentP4);
	if(deltaR < 0.5) {
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
	if (nProngs <= 5) {
	  recTaus.push_back(TauLepton(recoTau, track->P4(), charge, nProngs, nPhotons, nNHadrons, isolation));
	}
      } // Selected objects
    } // Those above 5 GeV
  } // All tracks
  return true;
}

bool makeRecTausPhotonSeeded(TClonesArray* branchEFTracks, TClonesArray* branchEFPhotons, TClonesArray* branchEFNHadrons, vector<TauLepton> &recTaus, bool clearRecTaus=false) {
  // Reconstruct taus using EFlow tracks, photons and neutral hadrons
  // Consider every charged track above 5 GeV PT - label at as a potential tau remnant
  // Add in all tracks and clusters within 0.3 of the track
  // Compute sum charge of the tracks
  // Count the number of charged tracks
  // Calculate mass of the tracks
  // Compute isolation as the sum of tracks and clusters in DeltaR = 0.3-0.5 range
  if (clearRecTaus) recTaus.clear();
  TLorentzVector recoTau;
  TLorentzVector fragmentP4;
  int nPhotons = 0;
  Track *theTrack = 0;
  for (Int_t p = 0; p < branchEFPhotons->GetEntries(); ++p) {
    Tower *thePhoton = (Tower*) branchEFPhotons->At(p);
    if (thePhoton->ET > 5.0) {
      bool selected = true;
      recoTau.SetPtEtaPhiM(thePhoton->ET, thePhoton->Eta, thePhoton->Phi, 0.0);
      double isolation = 0;
      int charge = 0;
      int nProngs = 0;
      // if (entry < 10) cout << "Max: EFTrack[" << t << "] = (" << track->PT << ", " << track->Eta << ", " << track->Phi << ", " << track->Mass << ")";
      for (Int_t a = 0; a < branchEFTracks->GetEntries(); ++a) {
	Track *fragment = (Track*) branchEFTracks->At(a);
	fragmentP4.SetPtEtaPhiM(fragment->PT, fragment->Eta, fragment->Phi, fragment->Mass);
	double deltaR = recoTau.DeltaR(fragmentP4);
	if(deltaR < 0.5) {
	  if (fragment->PT > thePhoton->ET) {
	    // if (entry < 10) cout << " - Deselected" << endl;
	    selected = false;
	    break;
	  }
	  else {
	    if (deltaR < 0.3) {
	      recoTau += fragmentP4;
	      // if (entry < 10) cout << endl << "Frg: EFTrack[" << a << "] = (" << fragment->PT << ", " << fragment->Eta << ", " << fragment->Phi << ", " << fragment->Mass << ")";
	      charge += fragment->Charge;
	      nProngs += 1;
	      if (theTrack == 0 || fragment->PT > theTrack->PT) {
		theTrack = fragment;
	      }
	    }
	    else {
	      isolation += fragment->PT;
	    }
	  }
	}
	for (Int_t a = 0; a < branchEFPhotons->GetEntries(); ++a) {
	  if (a != p) {
	    Tower *fragment = (Tower*) branchEFPhotons->At(a);
	    fragmentP4.SetPtEtaPhiM(fragment->ET, fragment->Eta, fragment->Phi, 0.);
	    if (fragment->ET > thePhoton->ET) {
	      // if (entry < 10) cout << " - Deselected" << endl;
	      selected = false;
	      break;
	    }
	    else {
	      double deltaR = recoTau.DeltaR(fragmentP4);
	      if(deltaR < 0.5) {
		if (deltaR < 0.3) {
		  recoTau += fragmentP4;
		  nPhotons++;
		}
		else {
		  isolation += fragment->ET;
		}
	      }
	    }
	  }
	}
      }
      // Consider non-pizero remnant neutral hadrons only for isolation
      int nNHadrons = 0;
      for (Int_t a = 0; a < branchEFNHadrons->GetEntries(); ++a) {
	Tower *fragment = (Tower*) branchEFNHadrons->At(a);
	fragmentP4.SetPtEtaPhiM(fragment->ET, fragment->Eta, fragment->Phi, 0.);
	double deltaR = recoTau.DeltaR(fragmentP4);
	if(deltaR < 0.5) {
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
	if (nProngs == 0) {
	  recTaus.push_back(TauLepton(recoTau, TLorentzVector(), charge, nProngs, nPhotons, nNHadrons, isolation));
	}
	else if (nProngs <= 5) {
	  recTaus.push_back(TauLepton(recoTau, theTrack->P4(), charge, nProngs, nPhotons, nNHadrons, isolation));
	}
      } // Selected objects
    } // Those above 5 GeV
  } // All tracks
  return true;
}

bool makeMchTaus(vector<TauLepton> &visTaus, vector<TauLepton> &recTaus, vector< pair<TauLepton, TauLepton> > &mchTaus) {
  mchTaus.clear();
  for (Int_t v = 0; v < visTaus.size(); v++) {
    // Find correct sign match
    for (Int_t r = 0; r < recTaus.size(); r++) {
      double deltaR = visTaus[v].P4.DeltaR(recTaus[r].P4);
      if (deltaR < 0.3) {
	mchTaus.push_back(std::pair(visTaus[v], recTaus[r]));
	break;
      }
    }
  }
  return true;
}

#endif
