
void RgZgPlot() {
  double R = 0.4;
  TString NAME, TITLE;
  TString lpf = "lpf";
  TString saveName, histoName, histoTitle, treeName;
  const int zVals = 2;        const int betaVals = 7;        const int nPtBins = 5;
  double jetPt[betaVals][zVals], jetEta[betaVals][zVals], jetPhi[betaVals][zVals], jetE[betaVals][zVals], sdPt[betaVals][zVals], sdEta[betaVals][zVals], sdPhi[betaVals][zVals];
  double sdE[betaVals][zVals], wt[betaVals][zVals], Rg[betaVals][zVals], Zg[betaVals][zVals];  int EventID[betaVals][zVals], nJetCons[betaVals][zVals], nSDCons[betaVals][zVals];
  int nEntries;

  TCanvas *r1[betaVals][zVals], *z1[betaVals][zVals];
  TLegend *rleg, *zleg;
  TLegendEntry *rentry, *zentry;
    
  double ptBinLo[nPtBins] = { 15, 20, 25, 30, 40 };
  double ptBinHi[nPtBins] = { 20, 25, 30, 40, 70 };
  TString ptBinString[nPtBins] = { "15-20 GeV", "20-25 GeV", "25-30 GeV", "30-40 GeV", "40+ GeV" };
  TString ptBinName[nPtBins] = { "_15_20", "_20_25", "_25_30", "_30_40", "_40" };
  int color[nPtBins] = { 100, 93, 81, 63, 51 };
  
  double z_cut[zVals] = { 0.10, 0.50 };   // Use a symmetry cut z > z_cut R^beta
  TString zString[zVals] = { "z_01_", "z_05_" };
  TString zTitleString[zVals] = { "0.10", "0.50" };
  double beta[betaVals]  = { -2.0, -1.0, -0.5, 0.00, 0.5, 1.0, 2.0 };
  TString betaString[betaVals] = { "_b_n20", "_b_n10", "_b_n05", "_b_00", "_b_05", "_b_10", "_b_20", };
  TString betaTitleString[betaVals]  = { "-2.0", "-1.0", "-0.5", "0.0", "0.5", "1.0", "2.0" };

  TH1D * h_Rg[nPtBins];
  TH1D * h_Zg[nPtBins];
  
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

  TFile* outFile = new TFile( "out/ZgRgPtAnalysis.root", "RECREATE" );
  
      for ( int z=0; z<zVals; ++z ) {
        for ( int b=0; b<betaVals; ++b ) {
      
      nEntries = SDtree[b][z]->GetEntries();

      NAME = "r1_" + zString[z] + betaString[b];
      r1[b][z] = new TCanvas( NAME , "" ,0 ,23 ,1280 ,700 );       gStyle->SetOptStat(0);   // CANVAS

      rleg = new TLegend(0.75, 0.68, 0.9, 0.9,NULL,"brNDC");    // LEGEND
      rleg->SetBorderSize(1);   rleg->SetLineColor(1);   rleg->SetLineStyle(1);   rleg->SetLineWidth(1);   rleg->SetFillColor(0);   rleg->SetFillStyle(1001);


      for ( int i=0; i<nPtBins; ++i ) {
	gStyle->SetOptTitle(1);
	NAME = "h_Rg_" + zString[z] + betaString[b] + "__pt" + ptBinName[i];
	TITLE = ptBinString[i];
	h_Rg[i] = RgVSpt[b][z]->ProjectionY( NAME , ptBinLo[i], ptBinHi[i] );       // PROJECT
	h_Rg[i]->Scale( 1./h_Rg[i]->Integral() );                     // NORMALIZE
	//	gStyle->SetOptTitle(0);
	h_Rg[i]->Write();
	h_Rg[i]->SetLineWidth(2);
	h_Rg[i]->SetLineColor( color[i] );
	h_Rg[i]->SetMarkerColor( color[i] );
	h_Rg[i]->Draw("SAME");                                                    // DRAW
	rentry=rleg->AddEntry( NAME, TITLE, lpf );                            // ADD TO LEGEND
	rentry->SetLineColor( color[i] );   rentry->SetMarkerColor( color[i] );
	rentry->SetFillStyle(1001);   rentry->SetTextFont(42);
	rentry->SetLineStyle(1);   rentry->SetLineWidth(2);
	rentry->SetMarkerStyle(1);   rentry->SetMarkerSize(1);
	lpf += "lpf";
      }

      rleg->Draw();
      r1[b][z]->Modified();
      r1[b][z]->cd();
      r1[b][z]->SetSelected(r1[b][z]);
      saveName = outDir + "RgPlot_" + zString[z] + betaString[b] + ".pdf";
      r1[b][z]->SaveAs( saveName ,"PDF");
      r1[b][z]->Write();
      r1[b][z]->Close();

      NAME = "z1_" + zString[z] + betaString[b];
      z1[b][z] = new TCanvas( NAME , "" ,0 ,23 ,1280 ,700 );       gStyle->SetOptStat(0);   // CANVAS
  
      zleg = new TLegend(0.75, 0.68, 0.9, 0.9,NULL,"brNDC");    // LEGEND
      zleg->SetBorderSize(1);   zleg->SetLineColor(1);   zleg->SetLineStyle(1);   zleg->SetLineWidth(1);   zleg->SetFillColor(0);   zleg->SetFillStyle(1001);

      
      for ( int i=0; i<nPtBins; ++i ) {
	gStyle->SetOptTitle(1);
	NAME = "h_Zg_" + zString[z] + betaString[b] + "__pt" + ptBinName[i];
	TITLE = ptBinString[i];
	h_Zg[i] = ZgVSpt[b][z]->ProjectionY( NAME , ptBinLo[i], ptBinHi[i] );       // PROJECT
	h_Zg[i]->Scale( 1./h_Zg[i]->Integral() );                     // NORMALIZE
	//	gStyle->SetOptTitle(0);
	h_Zg[i]->Write();
	h_Zg[i]->SetLineWidth(2);
	h_Zg[i]->SetLineColor( color[i] );
	h_Zg[i]->SetMarkerColor( color[i] );
	h_Zg[i]->Draw("SAME");                                                    // DRAW
	zentry=zleg->AddEntry( NAME, TITLE, lpf );                            // ADD TO LEGEND
	zentry->SetLineColor( color[i] );   zentry->SetMarkerColor( color[i] );
	zentry->SetFillStyle(1001);   zentry->SetTextFont(42);
	zentry->SetLineStyle(1);   zentry->SetLineWidth(2);
	zentry->SetMarkerStyle(1);   zentry->SetMarkerSize(1);
	lpf += "lpf";
      }

      zleg->Draw();
      z1[b][z]->Modified();
      z1[b][z]->cd();
      z1[b][z]->SetSelected(z1[b][z]);
      // saveName = outDir + "ZgPlot.pdf";
      saveName = outDir + "ZgPlot_" + zString[z] + betaString[b] + ".pdf";
      z1[b][z]->SaveAs( saveName ,"PDF");
      z1[b][z]->Write();
      z1[b][z]->Close();

      for ( int i=0; i<nPtBins; ++i ) {      h_Rg[i]->Clear();      h_Zg[i]->Clear();      }
      
    }
  }

      inFile->Close();    //outFile->Close();
      
}
