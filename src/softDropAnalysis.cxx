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
  const double jetMaxPt = 60.0;      //  Jet Pt <= 60.0 GeV

  const int zVals = 2;
  const int betaVals = 7;

  TString histoName, histoTitle, NAME;
  
  double z_cut[zVals] = {0.10, 0.50};   // Use a symmetry cut z > z_cut R^beta
  TString zString[zVals] = { "z_01_", "z_05_" };
  TString zTitleString[zVals] = { "0.10", "0.50" };
  double beta[betaVals]  = { -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0 };
  TString betaString[betaVals] = { "_b_n20", "_b_n10", "_b_n05", "_b_00", "_b_05", "_b_10", "_b_20", };
  TString betaTitleString[betaVals]  = { "-2.0", "-1.0", "-0.5", "0.0", "0.5", "1.0", "2.0" };

  int nEvents = 0;    int NJets = 0;    int ID;    int NsdJets[betaVals][zVals];
  double weight, Rval, Zval;
  TString fileName;  //  For x-sec
  vector<PseudoJet> rawParticles, Jets, rawJets, sdJets;
  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION

  TChain* Chain = new TChain( "JetTreeMc" );     Chain->Add( "AddedGeantPythia/picoDst*" );
  TStarJetPicoReader Reader;                               InitReaderPythia( Reader, Chain, numEvents );
  TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVector* sv;
  TStarJetVectorContainer<TStarJetVector> * container;
      
  TH1::SetDefaultSumw2( );  // Histograms will calculate gaussian errors
  TH2::SetDefaultSumw2( );
  TH3::SetDefaultSumw2( );

  TH1D *h_nJets = new TH1D( "nJets", "Number of Accepted Jets", 4,0,4);
  const int zbCombos = 4;
  const char *binLabelChar[zbCombos];
  string binLabelString[zbCombos] = { "z: 0.1, #beta: -0.5", "z: 0.1, #beta: 0.0", "z: 0.1, #beta: 2.0", "z: 0.5, #beta: 2.0" };
  for ( int i=0; i<zbCombos; ++i ) {
    binLabelChar[i] = binLabelString[i].c_str();
    h_nJets->GetXaxis()->SetBinLabel( i+1, binLabelChar[i]);
  }
  
  TH2D *h_RgVSpt[betaVals][zVals];
  TH2D *h_ZgVSpt[betaVals][zVals];
  TFile * outFile = new TFile ("out/softDropAnalysis.root","RECREATE");
  
  TTree * sdTree[betaVals][zVals];
  double jetPt[betaVals][zVals], jetEta[betaVals][zVals], jetPhi[betaVals][zVals], jetE[betaVals][zVals], sdPt[betaVals][zVals], sdEta[betaVals][zVals], sdPhi[betaVals][zVals];
  double sdE[betaVals][zVals], wt[betaVals][zVals], Rg[betaVals][zVals], Zg[betaVals][zVals];  int EventID[betaVals][zVals], nJetCons[betaVals][zVals], nSDCons[betaVals][zVals];

  int nZgBins;  double minZgBin, maxZgBin;
  
  for ( int z=0; z<zVals; ++z ) {
    for ( int b=0; b<betaVals; ++b ) {
      
      histoName = "h_Rg__" + zString[z] + betaString[b];
      histoTitle = "R_{g}:  z_{cut} = " + zTitleString[z] + ",  #beta = " + betaTitleString[b] + ";p_{T} (GeV);R_{g}";
      h_RgVSpt[b][z] = new TH2D( histoName, histoTitle, 70,0,70, 30,0,0.6);

      nZgBins = 60;  minZgBin = 0.0;  maxZgBin = 0.6;
      if ( z == 1 ) { nZgBins = 60;  minZgBin = 0.0;  maxZgBin = 0.6; }
      
      histoName = "h_Zg__" + zString[z] + betaString[b];
      histoTitle = "Z_{g}:  z_{cut} = " + zTitleString[z] + ",  #beta = " + betaTitleString[b] + ";p_{T} (GeV);z_{g}";
      h_ZgVSpt[b][z] = new TH2D( histoName, histoTitle, 70,0,70, nZgBins, minZgBin, maxZgBin);

      NAME = "sdTree__" + zString[z] + betaString[b];
      sdTree[b][z] = new TTree( NAME, "Soft-Drop Jet Tree" );

      sdTree[b][z]->Branch("jetPt", &jetPt[b][z]);  sdTree[b][z]->Branch("jetEta", &jetEta[b][z]);  sdTree[b][z]->Branch("jetPhi", &jetPhi[b][z]);
      sdTree[b][z]->Branch("jetE", &jetE[b][z]);  sdTree[b][z]->Branch("EventID", &EventID[b][z]);  sdTree[b][z]->Branch("nJetCons", &nJetCons[b][z]);
      sdTree[b][z]->Branch("wt", &wt[b][z]);  sdTree[b][z]->Branch("sdPt", &sdPt[b][z]);  sdTree[b][z]->Branch("sdEta", &sdEta[b][z]);
      sdTree[b][z]->Branch("sdPhi", &sdPhi[b][z]);  sdTree[b][z]->Branch("sdE", &sdE[b][z]);  sdTree[b][z]->Branch("nSDCons", &nSDCons[b][z]);
      sdTree[b][z]->Branch("Rg", &Rg[b][z]);  sdTree[b][z]->Branch("Zg", &Zg[b][z]);

      NsdJets[b][z] =0;

    }
  }



  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    ID = Reader.GetNOfCurrentEvent();
    Reader.PrintStatus(20); 

    rawParticles.clear();  rawJets.clear();  //  clear containers

    event = Reader.GetEvent();    header = event->GetHeader();

    container = Reader.GetOutputContainer();

    if ( Vz_candidate( header, absMaxVz ) == false ) { continue; }
    ID = Reader.GetNOfCurrentEvent();
    fileName =  Reader.GetInputChain()->GetCurrentFile()->GetName();
    weight = LookupXsec( fileName );
    
    GatherParticles ( container, etaCut, partMinPt, rawParticles);

    //  CREATE JET SELECTOR
    Selector etaSelector = SelectorAbsEtaMax( 1.0-R );    Selector ptMinSelector = SelectorPtMin(jetMinPt);    Selector ptMaxSelector = SelectorPtMax(jetMaxPt);
    Selector etaPtSelector = etaSelector && ptMinSelector && ptMaxSelector;

    ClusterSequence jetCluster( rawParticles, jet_def );           //  CLUSTER ALL JETS
    vector<PseudoJet> rawJets = sorted_by_pt( etaPtSelector( jetCluster.inclusive_jets() ) );     // EXTRACT SELECTED JETS
    NJets += rawJets.size();
    

    for ( int z=0; z<zVals; ++z ) {
      for ( int b=0; b<betaVals; ++b ) {

	sdJets.clear();
	contrib::SoftDrop sd( beta[b], z_cut[z], R );   // give the soft drop groomer a short name
	
	// if ( nEvents == 0 ) { cout << endl << endl << "SoftDrop groomer is: " << endl << sd.description() << endl << endl; }

	for (unsigned i = 0; i < rawJets.size(); i++) {
      
	  sdJets.push_back( sd(rawJets[i]) );            // ~ ~ ~ ~ ~ ~ ~ ~ PERFORM SOFT-DROP GROOMING ~ ~ ~ ~ ~ ~ ~ ~ 
    
	  Rval = sdJets[i].structure_of<contrib::SoftDrop>().delta_R();
	  Zval = sdJets[i].structure_of<contrib::SoftDrop>().symmetry();
	  // if ( Rval == -1 ) Rval = 0.0;	  if ( Zval == -1 ) Zval = 0.0;
	  assert(sdJets[i] != 0); //because soft drop is a groomer (not a tagger), it should always return a soft-dropped jet

	  jetPt[b][z] = rawJets[i].pt();    jetEta[b][z] = rawJets[i].eta();    jetPhi[b][z] = rawJets[i].phi();
	  jetE[b][z] = rawJets[i].e();    EventID[b][z] = ID;    wt[b][z] = weight;
	  std::vector<fastjet::PseudoJet> Cons= rawJets[i].constituents();
	  nJetCons[b][z] = Cons.size();

	  sdPt[b][z] = sdJets[i].pt();    sdEta[b][z] = sdJets[i].eta();    sdPhi[b][z] = sdJets[i].phi();
	  sdE[b][z] = sdJets[i].e();
	  std::vector<fastjet::PseudoJet> sdCons = sdJets[i].constituents();
	  nSDCons[b][z] = sdCons.size();

	  Rg[b][z] = Rval;    Zg[b][z] = Zval;

	  h_RgVSpt[b][z]->Fill( rawJets[i].pt(), Rval, weight );
	  h_ZgVSpt[b][z]->Fill( rawJets[i].pt(), Zval, weight );
	    
	  sdTree[b][z]->Fill();
	  
	  if ( Rval != -1 && Zval != -1) { NsdJets[b][z] += 1; }

	}
      }
    }
    nEvents++;

  }
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  // cout << NsdJets[2][0] << "   " << NsdJets[3][0] << "   " << NsdJets[6][0] << "   " << NsdJets[6][1] << endl;

  h_nJets->Fill( 0.0, NsdJets[2][0] );
  h_nJets->Fill( 1, NsdJets[3][0] );
  h_nJets->Fill( 2, NsdJets[6][0] );
  h_nJets->Fill( 3, NsdJets[6][1] );
  
  for ( int z=0; z<zVals; ++z ) {
    for ( int b=0; b<betaVals; ++b ) {
      sdTree[b][z]->Write();
      h_RgVSpt[b][z]->Write();
      h_ZgVSpt[b][z]->Write();
    }
  }
  
  h_nJets->Write();
  outFile->Close();
  
  return 0;
}
  // double z_cut[zVals] = {0.10, 0.50};   // Use a symmetry cut z > z_cut R^beta
  // double beta[betaVals]  = { -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0 };
