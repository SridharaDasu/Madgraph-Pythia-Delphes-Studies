/*

  A macro to identify generator level taus, their visible part and reconstruct
  taus using EFlow objects in a Delphes file.

  Taus object contains its P4, its highest track P4 and summary of other tau remnants

  Some histograms are made for checking this out in the tauMaker(const char* delphesRootFile) 

  root -l tauMaker.C'("delphes_output.root")'

*/

//------------------------------------------------------------------------------

R__LOAD_LIBRARY(libDelphes)
#include "TauLepton.h"
#include "Histograms.h"
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include <string.h>
#include <map>
#include <vector>
#include <string>
#include "TLorentzVector.h"


using namespace std;

void Analysis(const char *inputFile)
//was ZTauTauAnalysis
{
  gSystem->Load("libDelphes");
  gROOT->SetBatch(kTRUE);

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  chain.Print();
  // cout <<"hello"<< endl;
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchParticles = treeReader->UseBranch("Particle");
  TClonesArray *branchJet      = treeReader->UseBranch("Jet"); //reco
  TClonesArray *GenJets = treeReader->UseBranch("GenJet");
  TClonesArray *branchEFTracks = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFPhotons = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFNHadrons = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchMuon     = treeReader->UseBranch("Muon");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMET      = treeReader->UseBranch("MissingET");


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

  histograms.Define("GenTauPairMass", "Generator-level oppositely charged tau-pair Mass", 100, 0.0, 80.0);
  histograms.Define("GenTauPairPT", "Generator-level oppositely charged tau-pair P_{T}", 100, 0.0, 100.0);
  histograms.Define("GenTauPairEta", "Generator-level oppositely charged tau-pair eta", 100, 0.0, 10.0);
  histograms.Define("GenTauPairPhi", "Generator-level oppositely charged tau-pair phi", 100, -7.0, 7.0);

  histograms.Define("VisTauPairMass", "Visible oppositely charged tau-pair Mass", 100, 0.0, 200.0);
  histograms.Define("VisTauPairPT", "Visible oppositely charged tau-pair P_{T}", 100, 0.0, 100.0);
  histograms.Define("VisTauPairEta", "Visible oppositely charged tau-pair eta", 100, 0.0, 10.0);
  histograms.Define("VisTauPairPhi", "Visible oppositely charged tau-pair phi", 100, -7.0, 7.0);

  histograms.Define("RecTauPairMass", "Reconstructed oppositely charged tau-pair Mass", 100, 0.0, 200.0);
  histograms.Define("RecTauPairPT", "Reconstructed oppositely charged tau-pair P_{T}", 100, 0.0, 100.0);
  histograms.Define("RecTauPairEta", "Reconstructed oppositely charged tau-pair eta", 100, 0.0, 10.0);
  histograms.Define("RecTauPairPhi", "Reconstructed oppositely charged tau-pair phi", 100, -7.0, 7.0);

  histograms.Define("MchTauPairMass", "Visible oppositely charged tau-pair Mass", 100, 0.0, 200.0);
  histograms.Define("MchTauPairPT", "Visible oppositely charged tau-pair P_{T}", 100, 0.0, 100.0);
  histograms.Define("MchTauPairEta", "Visible oppositely charged tau-pair eta", 100, 0.0, 10.0);
  histograms.Define("MchTauPairPhi", "Visible oppositely charged tau-pair phi", 100, -7.0, 7.0);

  histograms.Define("WrgTauPairMass", "Visible oppositely charged tau-pair Mass", 100, 0.0, 200.0);
  histograms.Define("WrgTauPairPT", "Visible oppositely charged tau-pair P_{T}", 100, 0.0, 100.0);
  histograms.Define("WrgTauPairEta", "Visible oppositely charged tau-pair eta", 100, 0.0, 10.0);
  histograms.Define("WrgTauPairPhi", "Visible oppositely charged tau-pair phi", 100, -7.0, 7.0);

  histograms.Define("GenJetPairMass", "Generator-level oppositely charged jet-pair Mass", 100, 0.0, 200.0);
  histograms.Define("GenJetPairPT", "Generator-level oppositely charged jet-pair P_{T}", 100, 0.0, 100.0);
  histograms.Define("RecJetPairMass", "Reconstructed oppositely charged jet-pair Mass", 100, 0.0, 200.0);
  histograms.Define("RecJetPairPT", "Reconstructed oppositely charged jet-pair P_{T}", 100, 0.0, 100.0);

  histograms.Define("ZJetPairMass", "Generator-level Z boson jet-pair Mass", 100, 0.0, 200.0);
  histograms.Define("ZJetPairPT", "Generator-level Z boson jet-pair P_{T}", 100, 0.0, 100.0);

  histograms.Define("aJetPairMass", "Generator-level a jet-pair Mass", 100, 0.0, 100.0);
  histograms.Define("NonJetDecayJetPair_aMass", "Class 1-3, 5 a jet-pair Mass", 100, 0.0, 100.0);
  histograms.Define("aTauPairMass", "Generator-level a tau-pair Mass", 100, 0.0, 100.0);
  histograms.Define("Check4Vector", "Jet4Vector Pt", 100, 0.0, 100.0);
  histograms.Define("GoodJetMass", "Cleaned Jet Mass", 100, 0.0, 100.0);
  histograms.Define("deltaR", "deltaR between Jet and Electrons", 100, 0.0, 10.0);

  histograms.Define("Class1Z", "Class 1 Z Mass", 100, 0.0, 100.0);
  histograms.Define("Class2Z", "Class 2 Z Mass", 100, 0.0, 100.0);
  histograms.Define("Class3Z", "Class 3 MET", 100, 0.0, 100.0);
  histograms.Define("Class4Z", "Class 4 Z Mass", 100, 0.0, 100.0);
  histograms.Define("Class5Z", "Class 5 Z Mass", 100, 0.0, 100.0);

  histograms.Define("Class1aGenMass", "Z->Muon Decay Generated a->tautau Recoil (a->bb) Mass", 100, 0.0, 100.0);
  histograms.Define("Class2aGenMass", "Z->Electron Decay Generated a->tautau Recoil (a->bb) Mass", 100, 0.0, 100.0);
  histograms.Define("Class3aGenMass", "Z->Neutrino Decay Generated a->tautau Recoil (a->bb) Mass", 100, 0.0, 100.0);
  histograms.Define("Class4aGenMass", "Z->Jet Decay Generated a->tautau Recoil (a->bb) Mass", 100, 0.0, 100.0);
  histograms.Define("Class5aGenMass", "Z->Tau Decay Generated a->tautau Recoil (a->bb) Mass", 100, 0.0, 100.0);
  
  histograms.Define("Class1aVisMass", "Z->Muon Decay Visible a->tautau Recoil (a->bb) Mass", 100, 0.0, 100.0);
  histograms.Define("Class2aVisMass", "Z->Electron Decay Visible a->tautau Recoil (a->bb) Mass", 100, 0.0, 100.0);
  histograms.Define("Class3aVisMass", "Z->Neutrino Decay Visible a->tautau Recoil (a->bb) Mass", 100, 0.0, 100.0);
  histograms.Define("Class4aVisMass", "Z->Jet Decay Visible a->tautau Recoil (a->bb) Mass", 100, 0.0, 100.0);
  histograms.Define("Class5aVisMass", "Z->Tau Decay Visible a->tautau Recoil (a->bb) Mass", 100, 0.0, 100.0);
  
  histograms.Define("Class1aRecMass", "Z->Muon Decay Reconstructed a->tautau Recoil (a->bb) Mass", 100, 0.0, 100.0);
  histograms.Define("Class2aRecMass", "Z->Electron Decay Reconstructed a->tautau Recoil (a->bb) Mass", 100, 0.0, 100.0);
  histograms.Define("Class3aRecMass", "Z->Neutrino Decay Reconstructed a->tautau Recoil (a->bb) Mass", 100, 0.0, 100.0);
  histograms.Define("Class4aRecMass", "Z->Jet Decay Reconstructed a->tautau Recoil (a->bb) Mass", 100, 0.0, 100.0);
  histograms.Define("Class5aRecMass", "Z->Tau Decay Reconstructed a->tautau Recoil (a->bb) Mass", 100, 0.0, 100.0);
  
  vector<TauLepton> genTaus;
  vector<TauLepton> visTaus;
  vector<TauLepton> recTaus;
  vector< pair<TauLepton, TauLepton> > mchTaus;


  cout << "Number of entries = " << numberOfEntries << endl;

  // Loop over all events

  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    Electron *elec1 = 0;
    Electron *elec2 = 0;
    Electron *electron = 0;
    Muon *muon1 = 0;
    Muon *muon2 = 0;
    Jet *jet = 0;
    Jet *genJet = 0;
    Jet *genJet2 = 0;
    Jet *genJet3 = 0;
    Jet *genJet4 = 0;

    double dimuonMass=0;
    double dielectronMass=0;
    double tauPairMass = 0;
    double jetPairMass = 0;
    double tauPairMass2 = 0;
    double jetPairMass2 = 0;
    MissingET* Met = 0;

    double recTauPairMass = 0;
    double recTauPairMass2 = 0;

    double visTauPairMass = 0;
    double visTauPairMass2 = 0;

    TLorentzVector P4a;
    TLorentzVector P4e(0., 0., 125., 125.);
    TLorentzVector P4p(0., 0., -125., 125.);

    TLorentzVector recTauPair;
    TLorentzVector recTauPair2;

    genJet  = (Jet*) branchJet->At(0);
    genJet2  = (Jet*) branchJet->At(1);
    genJet3  = (Jet*) branchJet->At(2);
    genJet4  = (Jet*) branchJet->At(3);

    bool Class1 = false; //Z -> mu mu
    bool Class2 = false; //Z -> e  e
    bool Class3 = false; //Z -> nu nu
    bool Class4 = false; //Z -> j  j
    bool Class5 = false; //Z -> tau tau
    bool goodJet = true; 

    //Class 1: 2 muons
    if(branchMuon->GetEntries() > 1){
      // Take first two muons
      muon1 = (Muon *) branchMuon->At(0);
      muon2 = (Muon *) branchMuon->At(1);
      dimuonMass=((muon1->P4()) + (muon2->P4())).M();
      if(dimuonMass>80 && dimuonMass<100){
        Class1 = true;

      }
    }

    //Class 2: 2 electrons
    if(branchElectron->GetEntries() > 1) {
      elec1 = (Electron *) branchElectron->At(0);
      elec2 = (Electron *) branchElectron->At(1);
      dielectronMass=((elec1->P4()) + (elec2->P4())).M();
      if(dielectronMass>80 && dielectronMass<100){
        Class2 = true;
      }
    }

    // There should always be a pair of tau+ and tau- in the event from the a-decay
    // Sometimes there may ba a second pair of tau+ and tau- from the decay of a Z-boson
    // Collect all the particles

    TauLepton *genTau = 0;
    TauLepton *genTau2 = 0;
    TauLepton *genTau3 = 0;
    TauLepton *genTau4 = 0;
    bool genTauPairExists = false;
    bool genTauPair2Exists = false;
    TauLepton *visTau = 0;
    TauLepton *visTau2 = 0;
    TauLepton *visTau3 = 0;
    TauLepton *visTau4 = 0;
    bool visTauPairExists = false;
    bool visTauPair2Exists = false;
    if (makeGenTaus(branchParticles, genTaus, visTaus)) {
      for (Int_t i = 0; i < genTaus.size(); i++) {
        histograms.Fill("GenTauPT", genTaus[i].P4.Pt());
        histograms.Fill("MxTTauPT", genTaus[i].MaxTrackP4.Pt());
        histograms.Fill("TFrTauPT", genTaus[i].MaxTrackP4.Pt() / genTaus[i].P4.Pt());
        histograms.Fill("VisTauPT", visTaus[i].P4.Pt());
        for (int i2 = i + 1; i2 < genTaus.size(); i2++) {
          genTau = &genTaus[i];
          genTau2 = &genTaus[i2];
          if (genTau->charge * genTau2->charge == -1) { // Oppositely charged tau pair
            TLorentzVector tauPair = genTau->P4 + genTau2->P4;
	    genTauPairExists = true;
            tauPairMass = tauPair.M();
            histograms.Fill("GenTauPairMass", tauPair.M());
            histograms.Fill("GenTauPairPT", tauPair.Pt());
            histograms.Fill("GenTauPairEta", tauPair.Eta());
            histograms.Fill("GenTauPairPhi", tauPair.Phi());
            for (int i3 = i2 + 1; i3 < genTaus.size(); i3++) {
              for (int i4 = i3 + 1; i4 < genTaus.size(); i4++) {
                genTau3 = &genTaus[i3];
                genTau4 = &genTaus[i4];
                TLorentzVector tauPair2 = genTau3->P4 + genTau4->P4;
		genTauPair2Exists = true;
                tauPairMass2 = tauPair2.M();
              }
            }
          }
          visTau = &visTaus[i];
          visTau2 = &visTaus[i2];
          if (visTau->charge * visTau2->charge == -1) { // Oppositely charged tau pair
            TLorentzVector visTauPair = visTau->P4 + visTau2->P4;
	    visTauPairExists = true;
	    visTauPairMass = visTauPair.M();
            histograms.Fill("VisTauPairMass", visTauPair.M());
            histograms.Fill("VisTauPairPT", visTauPair.Pt());
            histograms.Fill("VisTauPairEta", visTauPair.Eta());
            histograms.Fill("VisTauPairPhi", visTauPair.Phi());
            for (int i3 = i2 + 1; i3 < visTaus.size(); i3++) {
              for (int i4 = i3 + 1; i4 < visTaus.size(); i4++) {
                visTau3 = &visTaus[i3];
                visTau4 = &visTaus[i4];
                TLorentzVector tauPair2 = visTau3->P4 + visTau4->P4;
	        visTauPair2Exists = true;
                visTauPairMass2 = tauPair2.M();
              }
            }
          }
        }
      }
    }

    //Cleaning big mess
    for (Int_t i = 0; i < branchJet->GetEntries(); i++) {
      jet  = (Jet*) branchJet->At(i);
      //Jet = (Jet *) branchJet->At(i);
      TLorentzVector Jet4Vector = jet->P4();
      histograms.Fill("Check4Vector", Jet4Vector.Pt());
      //T.Pt()-make sure i made 4 vector 
      //A.deltaR(B)
      //cout <<"Jet Clean Loop"<< endl;
      if(branchElectron->GetEntries() > 0)
        {
          //cout <<"Electron Clean Loop"<< endl;
          for (Int_t j = 0; j < branchElectron->GetEntries(); j++) {            
            electron  = (Electron*) branchElectron->At(j);
            TLorentzVector Electron4Vector = electron->P4();
            double deltaR = Jet4Vector.DeltaR(Electron4Vector);
            histograms.Fill("deltaR", deltaR);
            if (deltaR<1)
              {
                goodJet = false;
                //cout <<"badjet"<< endl;
              }

            //double deltaRJet = Jet.P4.DeltaR(Jet.P4);
          }
        }
      if (goodJet == true){
        histograms.Fill("GoodJetMass", jet->P4().M());

      }
    }

    if (branchJet->GetEntries() > 1){
      //loop over jet, loop over electrons: if electron and top jet

      TLorentzVector JetPair = (genJet->P4()) + (genJet2->P4());
      jetPairMass = JetPair.M();
      histograms.Fill("GenJetPairMass", JetPair.M());
      histograms.Fill("GenJetPairPT", JetPair.Pt());

      //Class 3: 2 neutrinos
      //QUESTION: should I be doing this with gen Taus?
      //change jets reco->gen
      Met = (MissingET *) branchMET->At(0);
      if(Met->MET > 30 && genTaus.size()==2 && branchJet->GetEntries()==2){
        //gen taus within gen jets
        Class3 = true;
      }
      //Class 4: 2 jets in Z range
      if (branchJet->GetEntries()==4 && JetPair.M()>70 && JetPair.M()<110){
        //need to tag this pair as a Z boson
        Class4 = true;

      }

      //Class 5: 4 taus
      if (genTaus.size()==4){
        Class5 = true;
      }
    }                                            


    // Reconstructed taus - these could be there even if there are no real taus in the event

    TauLepton *recTau = 0;
    TauLepton *recTau2 = 0;
    TauLepton *recTau3 = 0;
    TauLepton *recTau4 = 0;
    bool recTauPairExists = false;
    bool recTauPair2Exists = false;
    if (makeRecTaus(branchEFTracks, branchEFPhotons, branchEFNHadrons, recTaus)) {
      for (int r = 0; r < recTaus.size(); r++) {
        recTau = &recTaus[r];
        double deltaR = recTau->P4.DeltaR(recTau->MaxTrackP4);
        histograms.Fill("DeltaR", deltaR);
        histograms.Fill("RecTauMxTrkPT", recTau->MaxTrackP4.Pt());
        histograms.Fill("RecTauPT", recTau->P4.Pt());
        histograms.Fill("RecTauNP", recTau->nProngs);
        histograms.Fill("RecTauNG", recTau->nPhotons);
        histograms.Fill("RecTauNN", recTau->nNHadrons);
        histograms.Fill("RecTauIR", recTau->isolation / recTau->P4.Pt());
        histograms.Fill("RecTauCh", recTau->charge);
        histograms.Fill("RecTauMs", recTau->P4.M());
        if (recTau->nProngs != 1 && recTau->nProngs != 3) continue;
        for (int r2 = r + 1; r2 < recTaus.size(); r2++) {
          recTau2 = &recTaus[r2];
          if (recTau2->nProngs != 1 && recTau2->nProngs != 3) continue;
          if (recTau->charge * recTau2->charge == -1) { // Oppositely charged tau pair
            recTauPairExists = true;
            TLorentzVector recTauPair = recTau->P4 + recTau2->P4;
            recTauPairMass = recTauPair.M();
            histograms.Fill("RecTauPairMass", recTauPair.M());
            histograms.Fill("RecTauPairPT", recTauPair.Pt());
            histograms.Fill("RecTauPairEta", recTauPair.Eta());
            histograms.Fill("RecTauPairPhi", recTauPair.Phi());
            for (int r3 = r2 + 1; r3 < recTaus.size(); r3++) {
              for (int r4 = r3 + 1; r4 < recTaus.size(); r4++) {
                recTau3 = &recTaus[r3];
                recTau4 = &recTaus[r4];     
                TLorentzVector recTauPair2 = recTau3->P4 + recTau4->P4;
                recTauPairMass2 = recTauPair2.M();
                recTauPair2Exists = true;
              }

            }
          }
        }
      }
    }


    if (Class1 ==true && Class2==false && Class3==false && Class4==false && Class5==false){
      //cout <<"Class1"<< endl;
      if (genTauPairExists) {
	P4a = P4e + P4p - genTau->P4 - genTau2->P4 - muon1->P4() - muon2->P4();
	histograms.Fill("Class1aGenMass", P4a.M());
      }
      if (visTauPairExists) {
	P4a = P4e + P4p - visTau->P4 - visTau2->P4 - muon1->P4() - muon2->P4();
	histograms.Fill("Class1aVisMass", P4a.M());
      }
      if (recTauPairExists) {
	P4a = P4e + P4p - recTau->P4 - recTau2->P4 - muon1->P4() - muon2->P4();
	histograms.Fill("Class1aRecMass", P4a.M());
      }
      // P4a bb code for class 1, use vis taus
      histograms.Fill("Class1Z", dimuonMass);
      histograms.Fill("aJetPairMass", jetPairMass);
      //histograms.Fill("Class1aJetPairMass", jetPairMass);
      histograms.Fill("NonJetDecayJetPair_aMass", jetPairMass);
      if (recTauPairMass<50 && recTauPairMass>10)
        {
          histograms.Fill("aTauPairMass", recTauPairMass);
        }
    }

    if (Class1 ==false && Class2==true && Class3==false && Class4==false && Class5==false){
      //cout <<"Class2"<< endl;
      if (genTauPairExists) {
	P4a = P4e + P4p - genTau->P4 - genTau2->P4 - elec1->P4() - elec2->P4();
	histograms.Fill("Class2aGenMass", P4a.M());
      }
      if (visTauPairExists) {
	P4a = P4e + P4p - visTau->P4 - visTau2->P4 - elec1->P4() - elec2->P4();
	histograms.Fill("Class2aVisMass", P4a.M());
      }
      if (recTauPairExists) {
	P4a = P4e + P4p - recTau->P4 - recTau2->P4 - elec1->P4() - elec2->P4();
	histograms.Fill("Class2aRecMass", P4a.M());
      }
      histograms.Fill("Class2Z", dielectronMass);
      histograms.Fill("aJetPairMass", jetPairMass);
      histograms.Fill("NonJetDecayJetPair_aMass", jetPairMass);
      if (recTauPairMass<50 && recTauPairMass>10)
        {
          histograms.Fill("aTauPairMass", recTauPairMass);
        }
      histograms.Fill("aTauPairMass", recTauPairMass);
    }

    if (Class1 ==false && Class2==false && Class3==true && Class4==false && Class5==false){
      //cout <<"Class3"<< endl;
      //P4a = P4e + P4p - genTau->P4 - genTau2->P4 - muon1->P4 - muon2->P4;
      //QUESTION: IGNORING MET
      Met = (MissingET *) branchMET->At(0);
      histograms.Fill("Class3Z", Met->MET); //just plot as MET
      histograms.Fill("aJetPairMass", jetPairMass);
      //histograms.Fill("Class3aJetPairMass", P4a.M());
      histograms.Fill("NonJetDecayJetPair_aMass", jetPairMass);
      if (recTauPairMass<50 && recTauPairMass>10)
        {
          histograms.Fill("aTauPairMass", recTauPairMass);
        }
    }

    if (Class1 ==false && Class2==false && Class3==false && Class4==true && Class5==false){
      //cout <<"Class4"<< endl;
      if (genTauPairExists) {
	P4a = P4e + P4p - genTau->P4 - genTau2->P4 - genJet->P4() - genJet2->P4();
	histograms.Fill("Class4aGenMass", P4a.M());
      }
      if (visTauPairExists) {
	P4a = P4e + P4p - visTau->P4 - visTau2->P4 - genJet->P4() - genJet2->P4();
	histograms.Fill("Class4aVisMass", P4a.M());
      }
      if (recTauPairExists) {
	P4a = P4e + P4p - recTau->P4 - recTau2->P4 - genJet->P4() - genJet2->P4();
	histograms.Fill("Class4aRecMass", P4a.M());
      }
      histograms.Fill("Class4Z", jetPairMass);
      TLorentzVector JetPair2 = (genJet3->P4()) + (genJet4->P4());
      jetPairMass2 = JetPair2.M();
      histograms.Fill("aJetPairMass", jetPairMass2);
      //histograms.Fill("Class4aJetPairMass", jetPairMass2);
      if (recTauPairMass<50 && recTauPairMass>10)
        {
          histograms.Fill("aTauPairMass", recTauPairMass);
        }
    }

    if (Class1 ==false && Class2==false && Class3==false && Class4==false && Class5==true){
      //cout <<"Class5"<< endl;
      if (genTauPair2Exists) {
	P4a = P4e + P4p - genTau->P4 - genTau2->P4 - genTau3->P4 - genTau4->P4;
	histograms.Fill("Class5aGenMass", P4a.M());
      }
      if (visTauPair2Exists) {
	P4a = P4e + P4p - visTau->P4 - visTau2->P4 - visTau3->P4 - visTau4->P4;
	histograms.Fill("Class5aVisMass", P4a.M());
      }
      if (recTauPair2Exists) {      
	P4a = P4e + P4p - recTau->P4 - recTau2->P4 - recTau3->P4 - recTau4->P4;
	histograms.Fill("Class5aRecMass", P4a.M());
      }
      histograms.Fill("Class5Z", tauPairMass);
      histograms.Fill("aJetPairMass", jetPairMass);
      histograms.Fill("NonJetDecayJetPair_aMass", jetPairMass);
      if (recTauPairMass2<50 && recTauPairMass2>10)
        {
          histograms.Fill("aTauPairMass", recTauPairMass);
        }
      //QUESTION: should leading taus  by default be Z?

    }

    // Matched taus - these should be there if the algorithm is any good
    if (makeMchTaus(visTaus, recTaus, mchTaus)) {
      for (int m = 0; m < mchTaus.size(); m++) {
        TauLepton *visTau = &mchTaus[m].first;
        TauLepton *recTau = &mchTaus[m].second;
        if(visTau->charge == recTau->charge) {
          histograms.Fill("MchTauPT", visTau->P4.Pt());
          histograms.Fill("MchTauCh", recTau->charge);
          histograms.Fill("MchTauIso", recTau->isolation / recTau->P4.Pt());
          histograms.Fill("MchTauNPr", recTau->nProngs);
          histograms.Fill("MchTauNPh", recTau->nPhotons);
          histograms.Fill("MchTauNNH", recTau->nNHadrons);
          histograms.Fill("MchPTRes", (recTau->P4.Pt() - visTau->P4.Pt()) / visTau->P4.Pt());
          histograms.Fill("MchEtaRes", (recTau->P4.Eta() - visTau->P4.Eta()) / visTau->P4.Eta());
          histograms.Fill("MchPhiRes", (recTau->P4.Phi() - visTau->P4.Phi()) / visTau->P4.Phi());
          for (int m2 = m + 1; m2 < mchTaus.size(); m2++) {
            TauLepton *visTau2 = &mchTaus[m2].first;
            TauLepton *recTau2 = &mchTaus[m2].second;
            if (recTau->charge * recTau2->charge == -1) { // Oppositely charged tau pair
              TLorentzVector tauPair = recTau->P4 + recTau2->P4;
              if (visTau2->charge == recTau2->charge) {
                histograms.Fill("MchTauPairMass", tauPair.M());
                histograms.Fill("MchTauPairPT", tauPair.Pt());
                histograms.Fill("MchTauPairEta", tauPair.Eta());
                histograms.Fill("MchTauPairPhi", tauPair.Phi());
              }
              else {
                histograms.Fill("WrgTauPairMass", tauPair.M());
                histograms.Fill("WrgTauPairPT", tauPair.Pt());
                histograms.Fill("WrgTauPairEta", tauPair.Eta());
                histograms.Fill("WrgTauPairPhi", tauPair.Phi());
              }
            }
          }
        }
        else {
          histograms.Fill("WrgTauPT", visTau->P4.Pt());
          histograms.Fill("WrgTauCh", recTau->charge);
          histograms.Fill("WrgTauIso", recTau->isolation / recTau->P4.Pt());
          histograms.Fill("WrgTauNPr", recTau->nProngs);
          histograms.Fill("WrgTauNPh", recTau->nPhotons);
          histograms.Fill("WrgTauNNH", recTau->nNHadrons);
          histograms.Fill("WrgPTRes", (recTau->P4.Pt() - visTau->P4.Pt()) / visTau->P4.Pt());
          histograms.Fill("WrgEtaRes", (recTau->P4.Eta() - visTau->P4.Eta()) / visTau->P4.Eta());
          histograms.Fill("WrgPhiRes", (recTau->P4.Phi() - visTau->P4.Phi()) / visTau->P4.Phi());
          for (int m2 = m + 1; m2 < mchTaus.size(); m2++) {
            TauLepton *recTau2 = &mchTaus[m2].second;
            if (recTau->charge * recTau2->charge == -1) { // Oppositely charged tau pair
              TLorentzVector tauPair = recTau->P4 + recTau2->P4;
              histograms.Fill("WrgTauPairMass", tauPair.M());
              histograms.Fill("WrgTauPairPT", tauPair.Pt());
              histograms.Fill("WrgTauPairEta", tauPair.Eta());
              histograms.Fill("WrgTauPairPhi", tauPair.Phi());
            }
          }
        }
      }
    }
  }
}
