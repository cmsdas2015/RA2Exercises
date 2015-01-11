#ifndef SELECTION_H
#define SELECTION_H

// Event selection helper functions and cut values.
//
// Author: Matthias Schroeder
//         matthias.schroeder@AT@desy.de
//         November 2013
class Selection {
public:
  // Returns result of delta-phi selection
  static bool deltaPhi(double dPhi1, double dPhi2, double dPhi3) {
    return dPhi1 > 0.5 && dPhi2 > 0.5 && dPhi3 > 0.3;
  }

  // Returns result of HT baseline selection
  static bool ht(double ht) { return ht > 500.; }

  // Returns result of MHT baseline selection
  static bool mht(double mht) { return mht > 200.; }

  // Returns result of N(jets) baseline selection
  static bool nJets(unsigned int nJets) { return nJets >= 3; }

  // Returns result of search-bin selection (to be applied on top
  // of baseline selection)
  static unsigned int searchBin(double ht, double mht, unsigned int nJets);

  // Cut values for HT jets
  static float htJetPtMin() { return 50.; }
  static float htJetEtaMax() { return 2.5; }

  // Cut values for MHT jets
  static float mhtJetPtMin() { return 30.; }
  static float mhtJetEtaMax() { return 5.0; }

};


unsigned int Selection::searchBin(double ht, double mht, unsigned int nJets) {
  unsigned int bin = 0;

  if( nJets >= 3 && nJets <= 6 ) {
    
    if( ht > 500 && ht < 800 ) {

      if(      mht > 200 && mht < 500 ) bin = 1;
      else if( mht > 500 && mht < 750 ) bin = 2;
      else if( mht > 750              ) bin = 3;

    } else if( ht > 800  && ht < 1200 ) {

      if(      mht > 200 && mht < 500 ) bin = 4;
      else if( mht > 500 && mht < 750 ) bin = 5;
      else if( mht > 750              ) bin = 6;

    } else if( ht > 1200 ) {

      if(      mht > 200 && mht < 500 ) bin = 7;
      else if( mht > 500 && mht < 750 ) bin = 8;
      else if( mht > 750              ) bin = 9;

    }
    
  } else if( nJets >= 7 && nJets <= 8 ) {
    
    if( ht > 500 && ht < 800 ) {

      if(      mht > 200 && mht < 500 ) bin = 10;
      else if( mht > 500 && mht < 750 ) bin = 11;
      else if( mht > 750              ) bin = 12;

    } else if( ht > 800  && ht < 1200 ) {

      if(      mht > 200 && mht < 500 ) bin = 13;
      else if( mht > 500 && mht < 750 ) bin = 14;
      else if( mht > 750              ) bin = 15;

    } else if( ht > 1200 ) {

      if(      mht > 200 && mht < 500 ) bin = 16;
      else if( mht > 500 && mht < 750 ) bin = 17;
      else if( mht > 750              ) bin = 18;

    }

  } else if( nJets >= 9 ) {
    
    if( ht > 500 && ht < 800 ) {

      if(      mht > 200 && mht < 500 ) bin = 19;
      else if( mht > 500 && mht < 750 ) bin = 20;
      else if( mht > 750              ) bin = 21;

    } else if( ht > 800  && ht < 1200 ) {

      if(      mht > 200 && mht < 500 ) bin = 22;
      else if( mht > 500 && mht < 750 ) bin = 23;
      else if( mht > 750              ) bin = 24;

    } else if( ht > 1200 ) {

      if(      mht > 200 && mht < 500 ) bin = 25;
      else if( mht > 500 && mht < 750 ) bin = 26;
      else if( mht > 750              ) bin = 27;

    }

  }

  return bin;
}
#endif

