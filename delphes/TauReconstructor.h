#ifndef TauReconstructor_h
#define TauReconstructor_h

/** \class TauReconstructor
 *
 *  Make taus from tracks
 *
 *  \authors S. Dasu
 *
 */

#include "classes/DelphesModule.h"

#include <map>
#include <string>
#include <vector>

class TObjArray;
class TIterator;

#include "TLorentzVector.h"

class TauReconstructor: public DelphesModule
{
public:
  TauReconstructor();
  ~TauReconstructor();

  void Init();
  void Process();
  void Finish();

private:

  TLorentzVector P4;
  TLorentzVector MaxTrackP4;
  int charge;
  int nProngs;
  int nPhotons;
  int nNHadrons;
  double isolation;

  Double_t fMinTauSeedPT;
  Double_t fMaxTauSeedEta;  
  Double_t fMaxTauIsolDeltaR;
  Double_t fMaxTauCoreDeltaR;

  TObjArray *fchadronInputArray;
  TIterator *fchadronIterator;

  TObjArray *fnhadronInputArray;
  TIterator *fnhadronIterator;

  TObjArray *fphotonInputArray;
  TIterator *fphotonIterator;

  TObjArray *fOutputArray;

  ClassDef(TauReconstructor, 1)
};

#endif

