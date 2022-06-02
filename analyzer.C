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
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton   = treeReader->UseBranch("Photon");
  TClonesArray *branchMET      = treeReader->UseBranch("MissingET");


  // Book histograms
  TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 100, 0.0, 250.0);
  TH1 *histMass = new TH1F("mass", "M_{inv}(e_{1}, e_{2})", 100, 40.0, 140.0);
  TH1 *histDiPhotMass = new TH1F("DiPhotMass", "M_{inv}(a_{1}, a_{2})", 100, 40.0, 140.0);
  TH1 *histMET        = new TH1F("MET", "MET", 100, 0.0, 300.0);

  Double_t weight = 0.0004 * 300 * 1E3 / numberOfEntries;
  Double_t yield = 0;
  Double_t y2 = 0;
  
  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    {
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);

      // If event contains at least 1 jet
      if(branchJet->GetEntries() > 0)
	{
	  // Take first jet
	  Jet *jet = (Jet*) branchJet->At(0);

	  // Plot jet transverse momentum
	  histJetPT->Fill(jet->PT);

	  // Print jet transverse momentum
	  //cout << "Jet pt: "<<jet->PT << endl;
	}

      Electron *elec1, *elec2;

      // If event contains at least 2 electrons
      if(branchElectron->GetEntries() > 1)
	{
	  // Take first two electrons
	  elec1 = (Electron *) branchElectron->At(0);
	  elec2 = (Electron *) branchElectron->At(1);

	  // Plot their invariant mass
	  histMass->Fill(((elec1->P4()) + (elec2->P4())).M());
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
    
      yield += weight;
      y2 += weight*weight;

    }
  
  cout << "Event yield:   " << yield << " +/- " << sqrt(y2) << endl;
  cout << "Selection Eff: " << yield / (weight*numberOfEntries) << endl;

  // Show resulting histograms
  //histJetPT->Draw();
  //histMass->Draw();
  /*TCanvas * canvas1 = new TCanvas("canvas1");
    TCanvas * canvas2 = new TCanvas("canvas2");
    canvas1->cd();
    histJetPT->Draw();
    canvas2->cd();
    histMass->Draw();
  */

  char outputFile[1024];
  strcpy(outputFile, inputFile);
  *strstr(outputFile, ".root") = '\0';
  strcat(outputFile, "-histos.root");
  TFile *file1 = new TFile(outputFile, "UPDATE");
  histJetPT->Write();
  histMass->Write();
  histMET->Write();
  histDiPhotMass->Write();
  file1->Close();

}
