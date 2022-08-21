/*

  A simple class for handling Histograms a bit more efficiently
     1) Chose file name in the constructor
     2) Define histograms (currently 1D only)
     3) Use histogram name to fill anywhere in the code
     4) Destructor automatically saves all the histograms in root file and png files

*/

#include <string.h>
#include <map>
#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"

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
