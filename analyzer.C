/*
  Simple macro showing how to access branches from the delphes output root file,
  loop over events, and plot simple quantities such as the jet pt and the di-electron invariant
  mass.
           root -l analyzer.C'("delphes_output.root")'
*/

//------------------------------------------------------------------------------

R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include <string.h>
#include <list>
#include <vector>

void analyzer(const char *inputFile)
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
  TClonesArray *branchJet      = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon     = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton   = treeReader->UseBranch("Photon");
  TClonesArray *branchMET      = treeReader->UseBranch("MissingET");

  // Book histograms

  //change a
  TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 100, 0.0, 250.0);
  TH1 *histDiElecMass = new TH1F("DiElecMass", "M_{inv}(e_{1}, e_{2})", 100, 40.0, 140.0);
  TH1 *histDiPhotMass = new TH1F("DiPhotMass", "M_{inv}(phot_{1}, phot_{2})", 100, 40.0, 140.0); //why does this use a particles?
  TH1 *histDiJetMass = new TH1F("DiJetMass", "M_{inv}(jet_{1}, jet_{2})", 100, 40.0, 140.0);
  TH1 *histMET        = new TH1F("MET", "MET", 100, 0.0, 300.0);
  TH1 *histMuPairMass = new TH1F("MuPairmass", "M_{inv}(mu_{1}, mu_{2})", 100, 40.0, 140.0);
  TH1 *histDiJetMassZeDecay = new TH1F("DiJetMass from Z->e e Decay", "M_{inv}(jet_{1}, jet_{2})", 100, 40.0, 140.0);
  TH1 *histDiJetMassZmuDecay = new TH1F("DiJetMass from Z-> mu mu Decay", "M_{inv}(jet_{1}, jet_{2})", 100, 40.0, 140.0);
  TH1 *histDibJetMassZleptonDecay = new TH1F("DibJetMass from Z->lepton Decay", "M_{inv}(bjet_{1}, bjet_{2})", 100, 40.0, 140.0);
  TH1 *histDiTauMassZleptonDecay = new TH1F("DiTauMass from Z->lepton Decay", "M_{inv}(Tau_{1}, Tau_{2})", 100, 40.0, 140.0);

  TH1 *histNJetsZeDecay = new TH1F("NJets", "Number of Jets", 100, 0, 100);
  TH1 *histJetNCharged = new TH1F("NCharged", "Number of Charged Particles", 100, 0, 100);
  TH1 *histJetDeltaR = new TH1F("JetWidth", "Jet Width - DeltaR", 100, 0, 10.);

  TH1 *histSelectedJet1PT = new TH1F("selected_jet1_pt", "jet P_{T}", 100, 0.0, 250.0);
  TH1 *histSelectedJet2PT = new TH1F("selected_jet2_pt", "jet P_{T}", 100, 0.0, 250.0);
  TH1 *histSelectedTau1PT = new TH1F("selected_tau1_pt", "visible tau1 P_{T}", 100, 0.0, 250.0);
  TH1 *histSelectedTau2PT = new TH1F("selected_tau2_pt", "visible tau2 P_{T}", 100, 0.0, 250.0);


  Double_t weight = 0.0004 * 300 * 1E3 / numberOfEntries;
  Double_t yield = 0;
  Double_t y2 = 0;
  
  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    {
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);

      Jet* jet = 0;
      Jet* jet1  = 0;
      Jet* jet2  = 0;
      Electron *elec1 = 0;
      Electron *elec2 = 0;
      Muon *muon1 = 0;
      Muon *muon2 = 0;
      MissingET* Met = 0;
      Photon* phot1 = 0;
      Photon* phot2 = 0;
      Jet* selectedTau1 = 0;
      Jet* selectedTau2 = 0;
      Jet* selectedJet1 = 0;
      Jet* selectedJet2 = 0; 

      if (branchJet->GetEntries() > 1)
	{
	  jet1  = (Jet*) branchJet->At(0);
	  jet2  = (Jet*) branchJet->At(1);

	  histDiJetMass->Fill(((jet1->P4()) + (jet2->P4())).M()); 
	}
    
      // If event contains at least 2 electrons
      if(branchElectron->GetEntries() > 1) {
	elec1 = (Electron *) branchElectron->At(0);
	elec2 = (Electron *) branchElectron->At(1);
      }

      // If event contains at least 2 muons
      if(branchMuon->GetEntries() > 1)
        {
          // Take first two muons
          muon1 = (Muon *) branchMuon->At(0);
          muon2 = (Muon *) branchMuon->At(1);
          double dimuonMass=((muon1->P4()) + (muon2->P4())).M();

          // Plot their invariant mass
          histMuPairMass->Fill(dimuonMass);

          if (branchJet->GetEntries() > 1)
	    {
	      //plot dijet mass for events with Z boson decay to e+/e-
	      if (dimuonMass>70 && dimuonMass<110)
		{
		  histDiJetMassZmuDecay->Fill(((jet1->P4()) + (jet2->P4())).M());
		}
	    }
          
        }
          
      if(branchMET->GetEntries() > 0)
        {
          MissingET* Met = (MissingET *) branchMET->At(0);
          histMET->Fill(Met->MET, weight);
	}

      if(branchPhoton->GetEntries() > 1)
        {
          Photon* phot1 = (Photon *) branchPhoton->At(0);
          Photon* phot2 = (Photon *) branchPhoton->At(1);
          if(phot1->PT < 20 | phot2->PT < 20 | abs(phot1->Eta) > 2.5 | abs(phot2->Eta) > 2.5) continue;
          histDiPhotMass->Fill(((phot1->P4()) + (phot2->P4())).M(), weight);
        }

      // If event contains at least 1 jet
      // QUESTION: numbers not 2 jets?
      
      
      if (branchJet->GetEntries() >0)
	{
	  // Take first jet
	  jet = (Jet*) branchJet->At(0);
	  // Plot jet transverse momentum
	  histJetPT->Fill(jet->PT);
	}


      // If event contains at least 2 electrons
      if(branchElectron->GetEntries() > 1)
	{
          double dielectronMass=((elec1->P4()) + (elec2->P4())).M();
	  
          // Plot their invariant mass
          histDiElecMass->Fill(dielectronMass);
	  
	  if(!(branchJet->GetEntries() > 1)) continue;
	  
	  if (dielectronMass>70 && dielectronMass<110)
	    {//plot dijet mass for events with Z boson decay to e+/e-
	      histDiJetMassZeDecay->Fill(((jet1->P4()) + (jet2->P4())).M());
	      histNJetsZeDecay->Fill(branchJet->GetEntries());
	      if(branchJet->GetEntries()>3)
                {
                  for (Int_t j= 0; j <= branchJet->GetEntries(); ++j)
                    {
		      
                      jet = (Jet*) branchJet->At(j);
		      
		      histJetNCharged->Fill(jet->NCharged);
                      if (jet->NCharged <= 3)
                        {
                          double deltaR = sqrt(jet->DeltaEta*jet->DeltaEta + jet->DeltaPhi*jet->DeltaPhi);
			  histJetDeltaR->Fill(deltaR);
                          if (deltaR < 0.25) {
			    if (selectedTau1 == 0) {
			      selectedTau1 = jet;
			      histSelectedTau1PT->Fill(selectedTau1->PT);
			    }
			    else if (selectedTau2 == 0) 
			      {
				selectedTau2 = jet;
				histSelectedTau2PT->Fill(selectedTau2->PT);
			      }
			  }
			  
                        }

                      if (selectedJet1 == 0 && jet != selectedTau1 && jet != selectedTau2) {
                        selectedJet1 = jet;
			histSelectedJet1PT->Fill(selectedJet1->PT);
                      }
                      else if (selectedJet2 == 0 && jet != selectedTau1 && jet != selectedTau2) {
                        selectedJet2 = jet;
			histSelectedJet2PT->Fill(selectedJet2->PT);
                      }
		      
                      if (selectedTau1 != 0 && selectedTau2 != 0 && selectedJet1 != 0 && selectedJet2 != 0) {
			cout << "Found all four objects of interest" << endl;
			break;
                      }
                    }
                }
	      if (selectedTau1 == 0 || selectedTau2 == 0 || selectedJet1 == 0 || selectedJet2 == 0) 
		{
		  // We do not have tau-pair and a jet-pair, so give up
		  continue; //jump to beginning of loop
		}
	      
              
	      // We got tau pair and the jet pair
	      // make your pair masses ...
	      histDiTauMassZleptonDecay->Fill((selectedTau1->P4()+selectedTau2->P4()).M());      
	      histDibJetMassZleptonDecay->Fill((selectedJet1->P4()+selectedJet2->P4()).M());      
	      
	    }
	}
      
      
      yield += weight;
      y2 += weight*weight;
      //squared?
      
    }
  
  cout << "Event yield:   " << yield << " +/- " << sqrt(y2) << endl;
  cout << "Selection Eff: " << yield / (weight*numberOfEntries) << endl;

  // Show resulting histograms
  //histJetPT->Draw();
  //histMass->Draw();
  //TCanvas * canvas1 = new TCanvas("canvas1");
  TCanvas * canvas2 = new TCanvas("canvas2");
  //canvas1->cd();
  histJetPT->Draw();
  //canvas2->cd();
  //histMass->Draw();
  //canvas2->SaveAs("la.png");

  char outputFile[1024];
  strcpy(outputFile, inputFile);
  *strstr(outputFile, ".root") = '\0';
  strcat(outputFile, "-histos.root");
  TFile *file1 = new TFile(outputFile, "RECREATE");
  histJetPT->Write();
  histDiElecMass->Write();
  histMuPairMass->Write();
  histMET->Write();
  histDiPhotMass->Write();
  histDiJetMass->Write();
  histDiJetMassZeDecay->Write();
  histDiJetMassZmuDecay->Write();
  histDibJetMassZleptonDecay->Write();
  histDiTauMassZleptonDecay->Write();
  histNJetsZeDecay->Write();
  histJetNCharged->Write();
  histJetDeltaR->Write();
  histSelectedJet1PT->Write();
  histSelectedJet2PT->Write();
  histSelectedTau1PT->Write();
  histSelectedTau2PT->Write();
  file1->Close();

  //distinguish b and tau jets
  //what do we call h2/mass?
  //what does source command do?
  //NCharge <= 3
  //check tau tag's N charge
  //(1) n tracks (2)eta  delR=sqrt(del eta ^2 + del phi^2) require <=0.2: plot delR; plot tau pT
  //jet->NCharged <=    dsqrt
  //questions: goals of plots
  /*
    reference jet inside loop
  */
  //
}
