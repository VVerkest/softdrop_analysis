//  sdfunctions.hh
//  Veronica Verkest June 4, 2018

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "TROOT.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TClonesArray.h"
#include "TLatex.h"
#include "TMathText.h"
#include "TProfile.h"

// TStarJetPico
#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTower.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"
#include "TStarJetPicoTriggerInfo.h"
#include "TStarJetPicoUtils.h"

#include <string>
#include <iostream>
#include <sstream>
#include <sstream>
#include <iomanip>
#include <cmath>

#ifndef sdfunctions_hh
#define sdfunctions_hh

namespace sd_analysis {
  
  std::vector<fastjet::PseudoJet> GatherParticles ( TStarJetVectorContainer<TStarJetVector> * container , double etaCutVal, double partMinPtVal, std::vector<fastjet::PseudoJet> & rawParticles );

  void FillJetInfo_WEIGHTED(std::vector<fastjet::PseudoJet> &rawJets, TTree* Tree, int event, double &jPt, double &jEta, double &jPhi, double &jE, int &jEvent,int &jncons,double &wt,double weight);

  void FillSDJetInfo_WEIGHTED ( fastjet::PseudoJet &rawJets, TTree* Tree, int event, double &jPt, double &jEta, double &jPhi, double &jE, int &jEvent,int &jncons,double &wt,double &jRg,double &jZg,double weight,double rg_val,double zg_val );
  
  bool Vz_candidate( TStarJetPicoEventHeader* header, double VzCut );

  void PrintJet( fastjet::PseudoJet jet );

  double LookupXsec(TString::TString & currentfile );

  void InitReaderPythia( TStarJetPicoReader & reader, TChain* chain, int nEvents );
  
}

#endif
