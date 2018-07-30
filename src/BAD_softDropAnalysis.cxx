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

  TString histoName, histoTitle, NAME;
  
  double z_cut[zVals] = {0.10, 0.50};   // Use a symmetry cut z > z_cut R^beta
  TString zString[zVals] = { "z_01_", "z_05_" };
  TString zTitleString[zVals] = { "0.10", "0.50" };
  double beta[betaVals]  = { -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0 };
  TString betaString[betaVals] = { "_b_n20", "_b_n10", "_b_n05", "_b_00", "_b_05", "_b_10", "_b_20", };
  TString betaTitleString[betaVals]  = { "-2.0", "-1.0", "-0.5", "0.0", "0.5", "1.0", "2.0" };

  TH1::SetDefaultSumw2( );  // Histograms will calculate gaussian errors
  TH2::SetDefaultSumw2( );
  TH3::SetDefaultSumw2( );
  
  TH2D *h_RgVSpt[betaVals][zVals];
  TH2D *h_ZgVSpt[betaVals][zVals];
  TFile * outFile = new TFile ("out/softDropAnalysis.root","RECREATE");
  
  TString outFileName;
  TTree * sdTree[betaVals][zVals];

  for ( int z=0; z<zVals; ++z ) {
    
    for ( int b=0; b<betaVals; ++b ) {

      contrib::SoftDrop sd( beta[b], z_cut[z], R );   // give the soft drop groomer a short name
      cout << endl << endl << "SoftDrop groomer is: " << endl << sd.description() << endl << endl;

      TChain* Chain = new TChain( "JetTreeMc" );     Chain->Add( "AddedGeantPythia/picoDst*" );
      TStarJetPicoReader Reader;                               InitReaderPythia( Reader, Chain, numEvents );
      TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVector* sv;
      TStarJetVectorContainer<TStarJetVector> * container;

      histoName = "h_Rg__" + zString[z] + betaString[b];
      histoTitle = "R_{g}:  z_{cut} = " + zTitleString[z] + ",  #beta = " + betaTitleString[b] + ";p_{T} (GeV);R_{g}";
      h_RgVSpt[b][z] = new TH2D( histoName, histoTitle, 70,0,70, 40,0,1);

      histoName = "h_Zg__" + zString[z] + betaString[b];
      histoTitle = "Z_{g}:  z_{cut} = " + zTitleString[z] + ",  #beta = " + betaTitleString[b] + ";p_{T} (GeV);z_{g}";
      h_ZgVSpt[b][z] = new TH2D( histoName, histoTitle, 70,0,70, 50,0,1);

      JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
      TString fileName;  //  For x-sec
      vector<PseudoJet> rawParticles, Jets, rawJets, sdJets;

      double jetPt, jetEta, jetPhi, jetE, sdPt, sdEta, sdPhi, sdE, wt, Rg, Zg;
      int EventID, nJetCons, nSDCons;


      //      T H I S      I S      B A D      A N D      Y O U      S H O U L D      F E E L       B A D
      //  ... okay I think it's better now...
      NAME = "sdTree__" + zString[z] + betaString[b];
      sdTree[b][z] = new TTree( NAME, "Soft-Drop Jet Tree" );

      sdTree[b][z]->Branch("jetPt", &jetPt[b][z]);  sdTree[b][z]->Branch("jetEta", &jetEta[b][z]);  sdTree[b][z]->Branch("jetPhi", &jetPhi[b][z]);
      sdTree[b][z]->Branch("jetE", &jetE[b][z]);  sdTree[b][z]->Branch("EventID", &EventID[b][z]);  sdTree[b][z]->Branch("nJetCons", &nJetCons[b][z]);
      sdTree[b][z]->Branch("wt", &wt[b][z]);  sdTree[b][z]->Branch("sdPt", &sdPt[b][z]);  sdTree[b][z]->Branch("sdEta", &sdEta[b][z]);
      sdTree[b][z]->Branch("sdPhi", &sdPhi[b][z]);  sdTree[b][z]->Branch("sdE", &sdE[b][z]);  sdTree[b][z]->Branch("nSDCons", &nSDCons[b][z]);
      sdTree[b][z]->Branch("Rg", &Rg[b][z]);  sdTree[b][z]->Branch("Zg", &Zg[b][z]);
  
      int nEvents = 0;    int NJets = 0;    int ID;    int NsdJets = 0;
      double weight, Rval, Zval;
  
      // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      while ( Reader.NextEvent() ) {

	ID = Reader.GetNOfCurrentEvent();
	nEvents++;    Reader.PrintStatus(20); 

	rawParticles.clear();  Jets.clear();  rawJets.clear();  sdJets.clear();  //  clear containers

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
	NJets += rawJets.size();



	for (unsigned i = 0; i < rawJets.size(); i++) {
	  
	  sdJets.push_back( sd(rawJets[i]) );            // ~ ~ ~ ~ ~ ~ ~ ~ PERFORM SOFT-DROP GROOMING ~ ~ ~ ~ ~ ~ ~ ~ 
    
	  Rval = sdJets[i].structure_of<contrib::SoftDrop>().delta_R();
	  Zval = sdJets[i].structure_of<contrib::SoftDrop>().symmetry();
	  // if ( Rval == -1 ) Rval = 0.0;	  if ( Zval == -1 ) Zval = 0.0;
	  assert(sdJets[i] != 0); //because soft drop is a groomer (not a tagger), it should always return a soft-dropped jet

	  jetPt = rawJets[i].pt();    jetEta = rawJets[i].eta();    jetPhi = rawJets[i].phi();
	  jetE = rawJets[i].e();    EventID = ID;    wt = weight;
	  std::vector<fastjet::PseudoJet> Cons = rawJets[i].constituents();
	  nJetCons = Cons.size();

	  sdPt = sdJets[i].pt();    sdEta = sdJets[i].eta();    sdPhi = sdJets[i].phi();
	  sdE = sdJets[i].e();
	  std::vector<fastjet::PseudoJet> sdCons = sdJets[i].constituents();
	  nSDCons = sdCons.size();

	  Rg = Rval;    Zg = Zval;

	  h_RgVSpt[b][z]->Fill( rawJets[i].pt(), Rval, weight );
	  h_ZgVSpt[b][z]->Fill( rawJets[i].pt(), Zval, weight );
	    
	  // sdTree[b][z]->Fill();
	}

	NsdJets += sdJets.size();
      }
      // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      //h_RgVSpt[b][z]->Scale( 1./h_RgVSpt[b][z]->Integral() );
      h_RgVSpt[b][z]->Write();
      //h_ZgVSpt[b][z]->Scale( 1./h_ZgVSpt[b][z]->Integral() );
      h_ZgVSpt[b][z]->Write();
      
      // cout << endl << endl<< "Writing file to:  " << outFileName << endl << endl;
      // sdTree[b][z]->Write();
    }
  }

  // for ( int z=0; z<zVals; ++z ) {
  //   for ( int b=0; b<betaVals; ++b ) {
  //     outFile[betaVals][zVals]->Close();
  //   }
  // }

  outFile->Close();
  
  return 0;
}
