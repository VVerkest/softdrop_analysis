
void RgZg2DPlot() {
  double R = 0.4;
  TString NAME, TITLE;
  TString saveName, histoName, histoTitle, canvasName, canvasTitle, treeName;
  const int zVals = 2;        const int betaVals = 7;        const int nPtBins = 5;
  double jetPt[betaVals][zVals], jetEta[betaVals][zVals], jetPhi[betaVals][zVals], jetE[betaVals][zVals], sdPt[betaVals][zVals], sdEta[betaVals][zVals], sdPhi[betaVals][zVals];
  double sdE[betaVals][zVals], wt[betaVals][zVals], Rg[betaVals][zVals], Zg[betaVals][zVals];  int EventID[betaVals][zVals], nJetCons[betaVals][zVals], nSDCons[betaVals][zVals];
  int nEntries;
    
  double ptBinLo[nPtBins] = { 15, 20, 25, 30, 40 };
  double ptBinHi[nPtBins] = { 20, 25, 30, 40, 70 };
  TString ptBinString[nPtBins] = { "15-20 GeV", "20-25 GeV", "25-30 GeV", "30-40 GeV", "40+ GeV" };
  TString ptBinName[nPtBins] = { "_15_20", "_20_25", "_25_30", "_30_40", "_40" };
  
  double z_cut[zVals] = { 0.10, 0.50 };   // Use a symmetry cut z > z_cut R^beta
  TString zString[zVals] = { "z_01_", "z_05_" };
  TString zTitleString[zVals] = { "0.10", "0.50" };
  int marker[zVals] = { 8, 22 };
  double beta[betaVals]  = { -2.0, -1.0, -0.5, 0.00, 0.5, 1.0, 2.0 };
  TString betaString[betaVals] = { "_b_n20", "_b_n10", "_b_n05", "_b_00", "_b_05", "_b_10", "_b_20", };
  TString betaTitleString[betaVals]  = { "-2.0", "-1.0", "-0.5", "0.0", "0.5", "1.0", "2.0" };
  int color[betaVals] = { 100, 93, 81, 64, 58, 51, 6 };

  TH2D * RgVSpt[betaVals][zVals];
  TH2D * ZgVSpt[betaVals][zVals];
  
  TString inFileName = "out/softDropAnalysis.root";
  TFile* inFile = new TFile( inFileName, "READ" );
  TTree * SDtree[betaVals][zVals];
  TString outDir = "pdfs/";
  
  for ( int z=0; z<zVals; ++z ) {
    for ( int b=0; b<betaVals; ++b ) {

      treeName = "sdTree__" + zString[z] + betaString[b];
      SDtree[b][z] = (TTree*) inFile->Get( treeName );

      SDtree[b][z]->SetBranchAddress("jetPt",&jetPt[b][z]);  SDtree[b][z]->SetBranchAddress("jetEta",&jetEta[b][z]);  SDtree[b][z]->SetBranchAddress("jetPhi",&jetPhi[b][z]);
      SDtree[b][z]->SetBranchAddress("jetE",&jetE[b][z]);  SDtree[b][z]->SetBranchAddress("wt",&wt[b][z]);  SDtree[b][z]->SetBranchAddress("EventID",&EventID[b][z]);
      SDtree[b][z]->SetBranchAddress("sdPt",&sdPt[b][z]);  SDtree[b][z]->SetBranchAddress("sdEta",&sdEta[b][z]);  SDtree[b][z]->SetBranchAddress("sdPhi",&sdPhi[b][z]);
      SDtree[b][z]->SetBranchAddress("sdE",&sdE[b][z]);  SDtree[b][z]->SetBranchAddress("EventID",&EventID[b][z]);  SDtree[b][z]->SetBranchAddress("nJetCons",&nJetCons[b][z]);
      SDtree[b][z]->SetBranchAddress("nSDCons",&nSDCons[b][z]);  SDtree[b][z]->SetBranchAddress("Rg",&Rg[b][z]);  SDtree[b][z]->SetBranchAddress("Zg",&Zg[b][z]);

      histoName = "h_Rg__" + zString[z] + betaString[b];
      RgVSpt[b][z] = (TH2D*) inFile->Get( histoName );
      
      histoName = "h_Zg__" + zString[z] + betaString[b];
      ZgVSpt[b][z] = (TH2D*) inFile->Get( histoName );
      
    }
  }

  

  TFile* outFile = new TFile( "out/RgZg2DPlots.root", "RECREATE" );
      
  for ( int i=0; i<nPtBins; ++i ) {

  
    for ( int z=0; z<zVals; ++z ) {
      for ( int b=0; b<betaVals; ++b ) {

  	if ( z==0 ) {    if ( b==0|| b==1 || b==4 || b==5  ) { continue; }    }
  	if ( z==1 && b!=6) { continue; }
	
	// histoName = "h_Rg__" + zString[z] + betaString[b] + "_" + ptBinName[i];
	// RgVSpt[b][z]->SetName( histoName );
  	// RgVSpt[b][z]->SetAxisRange( ptBinLo[i], ptBinHi[i], "x");
	// histoName = "h_Zg__" + zString[z] + betaString[b] + "_" + ptBinName[i];
	// ZgVSpt[b][z]->SetName( histoName );
  	// ZgVSpt[b][z]->SetAxisRange( ptBinLo[i], ptBinHi[i], "x");
	
  	RgVSpt[b][z]->Write();
  	ZgVSpt[b][z]->Write();
      }
    }



  }


  
  inFile->Close();    //outFile->Close();
      
}
