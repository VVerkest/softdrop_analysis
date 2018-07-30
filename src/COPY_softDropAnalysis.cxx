//  SOFT-DROP ANALYSIS
//  Veronica Verkest    June 4, 2018

//  "Explore the parameter space" of soft-drop grooming
//  Parameters:  Zg, Beta, centrality    ???
//  Observables:  Aj, Rg    ???

#include "sdfunctions.hh"

using namespace fastjet;
using namespace std;
using namespace sd_analysis;

int main () {

  const double numEvents = -1;       //  NUMBER OF EVENTS  (-1 runs all)
  const double pi = 3.141592653;
  const double R = 0.4;
  const double absMaxVz = 30.0;   // |Vz|<=30 cm
  const double etaCut = 1.0;
  const double partMinPt = 0.2;      //  particle Pt >= 0.2 GeV
  const double jetMinPt = 2.0;      //  Jet Pt >= 2.0 GeV

  const int zVals = 2;
  const int betaVals = 7;
  
  double z_cut[zVals] = {0.10, 0.50};   // Use a symmetry cut z > z_cut R^beta
  TString zString[zVals] = { "z_01_", "z_05_" };
  double beta[betaVals]  = { -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0 };
  TString betaString[betaVals] = { "_b_n20", "_b_n10", "_b_n05", "_b_00", "_b_05", "_b_10", "_b_20", };
  
  TH1::SetDefaultSumw2( );  // Histograms will calculate gaussian errors
  TH2::SetDefaultSumw2( );
  TH3::SetDefaultSumw2( );

  for ( int z=0; z<zVals; ++z ) {
    
    // for ( int b=0; b<betaVals; ++b ) {

    int b=3;
    
      TString outFileName = "out/softDropAnalysis_" + zString[z] + betaString[b] + ".root";
      contrib::SoftDrop sd( beta[b], z_cut[z], R );   // give the soft drop groomer a short name
      cout << endl << endl << "SoftDrop groomer is: " << endl << sd.description() << endl << endl;

      TChain* Chain = new TChain( "JetTreeMc" );     Chain->Add( "AddedGeantPythia/picoDst*" );
      TStarJetPicoReader Reader;                               InitReaderPythia( Reader, Chain, numEvents );
      TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVector* sv;
      TStarJetVectorContainer<TStarJetVector> * container;

      JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
      TString fileName;  //  For x-sec
      vector<PseudoJet> rawParticles, Jets, rawJets, jetCons, sdJets, sdCons;

      double jetPt, jetEta, jetPhi, jetE, wt, Rg, Zg;
      int EventID, nCons;

      //  ALL JETS
      TTree * jetTree = new TTree( "jetTree", "Jet Tree" );
      TBranch *j_Pt;  TBranch *j_Eta;  TBranch *j_Phi;  TBranch *j_E;  TBranch *j_EventID;  TBranch *j_NCons;  TBranch *j_Weight;
      j_Pt = jetTree->Branch("jetPt", &jetPt);  j_Eta = jetTree->Branch("jetEta", &jetEta);  j_Phi = jetTree->Branch("jetPhi", &jetPhi);
      j_E = jetTree->Branch("jetE", &jetE);  j_EventID = jetTree->Branch("EventID", &EventID);  j_NCons = jetTree->Branch("nCons", &nCons);
      j_Weight = jetTree->Branch("wt", &wt);

      //  SOFT-DROP GROOMED JETS
      TTree * sdTree = new TTree( "sdTree", "Soft-Drop Jet Tree" );
      TBranch *sd_Pt;  TBranch *sd_Eta;  TBranch *sd_Phi;  TBranch *sd_E;  TBranch *sd_EventID;  TBranch *sd_NCons;  TBranch *sd_Weight;  TBranch *sd_Rg;  TBranch *sd_Zg;
      sd_Pt = sdTree->Branch("jetPt", &jetPt);  sd_Eta = sdTree->Branch("jetEta", &jetEta);  sd_Phi = sdTree->Branch("jetPhi", &jetPhi);
      sd_E = sdTree->Branch("jetE", &jetE);  sd_EventID = sdTree->Branch("EventID", &EventID);  sd_NCons = sdTree->Branch("nCons", &nCons);
      sd_Weight = sdTree->Branch("wt", &wt);  sd_Rg = sdTree->Branch("Rg", &Rg);  sd_Zg = sdTree->Branch("Zg", &Zg);
  
      int nEvents = 0;    int NJets = 0;    int ID;    int NsdJets = 0;
      double weight, Rval, Zval;
  
      // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      while ( Reader.NextEvent() ) {

	ID = Reader.GetNOfCurrentEvent();
	nEvents++;    Reader.PrintStatus(20); 

	rawParticles.clear();  Jets.clear();  rawJets.clear();  jetCons.clear();  sdJets.clear();  sdCons.clear();  //  clear containers

	event = Reader.GetEvent();    header = event->GetHeader();

	container = Reader.GetOutputContainer();

	if ( Vz_candidate( header, absMaxVz ) == false ) { continue; }
	ID = Reader.GetNOfCurrentEvent();
	fileName =  Reader.GetInputChain()->GetCurrentFile()->GetName();
	weight = LookupXsec( fileName );
    
	GatherParticles ( container, etaCut, partMinPt, rawParticles);

	//  CREATE JET SELECTOR
	Selector etaSelector = SelectorAbsEtaMax( 1.0-R );    Selector ptSelector = SelectorPtMin(jetMinPt);    Selector etaPtSelector = etaSelector && ptSelector;

	ClusterSequence jetCluster( rawParticles, jet_def );           //  CLUSTER ALL JETS
	vector<PseudoJet> rawJets = sorted_by_pt( etaPtSelector( jetCluster.inclusive_jets() ) );     // EXTRACT SELECTED JETS
	FillJetInfo_WEIGHTED( rawJets, jetTree, ID, jetPt, jetEta, jetPhi, jetE, EventID, nCons, wt, weight );
	NJets += rawJets.size();

	// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ PERFORM SOFT-DROP GROOMING ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
	for (unsigned i = 0; i < rawJets.size(); i++) {
	  sdJets.push_back( sd(rawJets[i]) );
    
	  Rval = sdJets[i].structure_of<contrib::SoftDrop>().delta_R();
	  Zval = sdJets[i].structure_of<contrib::SoftDrop>().symmetry();
	  // if ( Rval == -1 ) Rval = 0.0;
	  // if ( Zval == -1 ) Zval = 0.0;
      
	  assert(sdJets[i] != 0); //because soft drop is a groomer (not a tagger), it should always return a soft-dropped jet
	  // cout << endl << " Original  jet:  " ;      PrintJet( rawJets[i] );
	  // cout << "SoftDropped jet: ";      PrintJet( sdJets[i] );
	  // cout << "Rg between subjets: " << Rval << "          " ;
	  // cout << "  symmetry measure (z): " << Zval << endl;
      
	  FillSDJetInfo_WEIGHTED ( sdJets[i], sdTree, ID, jetPt, jetEta, jetPhi, jetE, EventID, nCons, wt, Rg, Zg, weight, Rval, Zval );
	}

	NsdJets += sdJets.size();
      }
      // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      TFile * outFile = new TFile( outFileName,"RECREATE" );
      cout << endl << endl<< "Writing file to:  " << outFileName << endl << endl;
      jetTree->Write();
      sdTree->Write();
      outFile->Close();
      //   }
  }
  
  return 0;
}
