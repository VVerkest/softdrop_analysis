//  sdfunctions.cxx
//  Veronica Verkest June 4, 2018

#include "sdfunctions.hh"
#include "fastjet/contrib/SoftDrop.hh"

namespace sd_analysis {




  std::vector<fastjet::PseudoJet> GatherParticles ( TStarJetVectorContainer<TStarJetVector> * container , double etaCutVal, double partMinPtVal, std::vector<fastjet::PseudoJet> & rawParticles ){
    for ( int i=0; i < container->GetEntries() ; ++i ) {
      TStarJetVector* sv = container->Get(i);
      fastjet::PseudoJet current = fastjet::PseudoJet( *sv );
      current.set_user_index( sv->GetCharge() );

      if ( std::abs(current.eta()) > 1.0 )      { continue; }  // removes particles with eta>|1|
      if ( current.pt() < 0.2 )      { continue; }  // removes particles with pt<0.2GeV

      rawParticles.push_back(current);
    }
    return rawParticles;
  }



  void FillJetInfo_WEIGHTED ( std::vector<fastjet::PseudoJet> &rawJets, TTree* Tree, int event, double &jPt, double &jEta, double &jPhi, double &jE, int &jEvent,int &jncons,double &wt,double weight ) {
    for ( int j=0; j<rawJets.size(); j++ ) {   // FILL JET INFO WEIGHED BY X-SECTION
      jPt = rawJets[j].pt();    jEta = rawJets[j].eta();    jPhi = rawJets[j].phi();
      jE = rawJets[j].e();    jEvent = event;    wt = weight;
      std::vector<fastjet::PseudoJet> Cons = rawJets[j].constituents(); jncons = Cons.size();
      Tree->Fill();
    }
  }



  void FillSDJetInfo_WEIGHTED ( fastjet::PseudoJet &rawJet, TTree* Tree, int event, double &jPt, double &jEta, double &jPhi, double &jE, int &jEvent,int &jncons,double &wt,double &jRg,double &jZg,double weight,double rg_val,double zg_val ) {      // FILL JET INFO WEIGHED BY X-SECTION
      jPt = rawJet.pt();    jEta = rawJet.eta();    jPhi = rawJet.phi();
      jE = rawJet.e();    jEvent = event;    wt = weight;
      std::vector<fastjet::PseudoJet> Cons = rawJet.constituents(); jncons = Cons.size();
      jRg = rg_val;      jZg = zg_val;
      Tree->Fill();
  }

  
  
  bool Vz_candidate( TStarJetPicoEventHeader* header, double VzCut ) {
    double vertexZ = header->GetPrimaryVertexZ();
    if ( vertexZ >= VzCut || vertexZ <= -VzCut )
      return false;      //    |Vz| <= 30.0 cm
    return true;
  }



  void PrintJet( fastjet::PseudoJet jet) {
    if (jet == 0) {
      std::cout << " 0 ";
    } else {
      std::cout << " pt = " << jet.pt() << "    m = " << jet.m() << "    y = " << jet.rap() << "    phi = " << jet.phi() << std::endl;
    }
  }
  


  double LookupXsec(TString::TString & currentfile ) {

    static const Double_t Xsec[12] = {
      1.0,        // Placeholder for 2-3
      1.30E+09, // 3-4
      3.15E+08, // 4-5
      1.37E+08, // 5-7
      2.30E+07, // 7-9
      5.53E+06, // 9-11
      2.22E+06, // 11-15
      3.90E+05, // 15-25
      1.02E+04, // 25-35
      5.01E+02, // 35-45
      2.86E+01, // 45-55
      1.46E+00 // 55-65
    };

    static const Double_t Nmc[12] = {
      1, // 2-3
      672518, // 3-4
      672447, // 4-5
      393498, // 5-7
      417659, // 7-9
      412652, // 9-11
      419030, // 11-15
      396744, // 15-25
      399919, // 25-35
      119995, // 35-45
      117999, // 45-55
      119999 // 55-65
    };

    Double_t w[12];
    for ( int i=0; i<12 ; ++i ){
      w[i] = Xsec[i] / Nmc[i];
      // w[i] = Nmc[i] / Xsec[i];
    }

    if ( currentfile.Contains("picoDst_3_4") ) return w[1];
    if ( currentfile.Contains("picoDst_4_5") ) return w[2];
    if ( currentfile.Contains("picoDst_5_7") ) return w[3];
    if ( currentfile.Contains("picoDst_7_9") ) return w[4];
    if ( currentfile.Contains("picoDst_9_11") ) return w[5];
    if ( currentfile.Contains("picoDst_11_15") ) return w[6];
    if ( currentfile.Contains("picoDst_15_25") ) return w[7];
    if ( currentfile.Contains("picoDst_25_35") ) return w[8];
    if ( currentfile.Contains("picoDst_35_45") ) return w[9];
    if ( currentfile.Contains("picoDst_45_55") ) return w[10];
    if ( currentfile.Contains("picoDst_55_65") ) return w[11];
    return 1;
  }


  

  
  //  INITIATE PYTHIA READER
  void InitReaderPythia( TStarJetPicoReader & reader, TChain* chain, int nEvents ) {

    std::string collisionType = "pp";
      
    // First tolower() on the analysisType
    // shouldnt be necessary....
    std::transform(collisionType.begin(), collisionType.end(), collisionType.begin(), ::tolower);
    
    // set the chain
    reader.SetInputChain( chain );
    // apply hadronic correction - subtract 100% of charged track energy from towers
    reader.SetApplyFractionHadronicCorrection( true );
    reader.SetFractionHadronicCorrection( 0.9999 );
    reader.SetRejectTowerElectrons( kFALSE );
    
    // Event and track selection
    // -------------------------
    
    TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
    evCuts->SetVertexZCut ( 1000 );
    evCuts->SetRefMultCut( 0 );
    evCuts->SetMaxEventPtCut( 1000 );
    evCuts->SetMaxEventEtCut( 1000 );

    // Tracks cuts
    TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
    trackCuts->SetDCACut( 100 );
    trackCuts->SetMinNFitPointsCut( -1 );
    trackCuts->SetFitOverMaxPointsCut( -1 );    
    
    std::cout << "Using these track cuts:" << std::endl;
    std::cout << " dca : " << trackCuts->GetDCACut(  ) << std::endl;
    std::cout << " nfit : " <<   trackCuts->GetMinNFitPointsCut( ) << std::endl;
    std::cout << " nfitratio : " <<   trackCuts->GetFitOverMaxPointsCut( ) << std::endl;
    
    // Towers
    TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
    towerCuts->SetMaxEtCut( 9999.0 );
    towerCuts->AddBadTowers( "src/dummy_tower_list.txt" );

    std::cout << "Using these tower cuts:" << std::endl;
    std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
    std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;
    
    // V0s: Turn off
    reader.SetProcessV0s(false);
    
    // Initialize the reader
    reader.Init( nEvents ); //runs through all events with -1
  }

  
}
