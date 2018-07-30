
void sdplot() {
  double R = 0.4;
  TString NAME, TITLE;
  TString saveName, histoName, histoTitle, canvasName, canvasTitle, treeName;
  const int zVals = 2;        const int betaVals = 7;        const int nPtBins = 5;
  double jetPt[betaVals][zVals], jetEta[betaVals][zVals], jetPhi[betaVals][zVals], jetE[betaVals][zVals], sdPt[betaVals][zVals], sdEta[betaVals][zVals], sdPhi[betaVals][zVals];
  double sdE[betaVals][zVals], wt[betaVals][zVals], Rg[betaVals][zVals], Zg[betaVals][zVals];  int EventID[betaVals][zVals], nJetCons[betaVals][zVals], nSDCons[betaVals][zVals];
  int nEntries;  int nZgBins;  double minZgBin;  double maxZgBin;
  double baseRg, groomRg, baseZg, groomZg, ptVal, dPt;

  TCanvas *r1[nPtBins][zVals], *z1[nPtBins][zVals];
  TLegend *rleg[nPtBins][zVals], *zleg[nPtBins][zVals], *dPtLeg[betaVals][zVals];
  TLegendEntry *rentry, *zentry, *dPtEntry;
    
  double ptBinLo[nPtBins] = { 15, 20, 25, 30, 40 };
  double ptBinHi[nPtBins] = { 20, 25, 30, 40, 70 };
  TString ptBinString[nPtBins] = { "15-20 GeV", "20-25 GeV", "25-30 GeV", "30-40 GeV", "40+ GeV" };
  TString ptBinName[nPtBins] = { "_15_20", "_20_25", "_25_30", "_30_40", "_40" };
  int ptColor[nPtBins] = { 100, 93, 81, 63, 51 };

  
  double z_cut[zVals] = { 0.10, 0.50 };   // Use a symmetry cut z > z_cut R^beta
  TString zString[zVals] = { "z_01_", "z_05_" };
  TString zTitleString[zVals] = { "0.10", "0.50" };
  int marker[zVals] = { 8, 22 };
  double beta[betaVals]  = { -2.0, -1.0, -0.5, 0.00, 0.5, 1.0, 2.0 };
  TString betaString[betaVals] = { "_b_n20", "_b_n10", "_b_n05", "_b_00", "_b_05", "_b_10", "_b_20", };
  TString betaTitleString[betaVals]  = { "-2.0", "-1.0", "-0.5", "0.0", "0.5", "1.0", "2.0" };
  int color[betaVals] = { 100, 93, 81, 64, 58, 51, 6 };

  TH2D * h_dPt_vs_pt[betaVals][zVals];
  
  TH1D * h_Rg[nPtBins][betaVals][zVals];
  TH1D * h_Zg[nPtBins][betaVals][zVals];

  TH2D * RgVSpt[betaVals][zVals];
  TH2D * ZgVSpt[betaVals][zVals];

  TH3D * Rg_base_pt[betaVals][zVals];     // Baseline:  z_cut=0.1, beta = 0
  TH3D * Zg_base_pt[betaVals][zVals];     // b=3, z=0

  const int nCorrPlots = 4;
  TH2D * RgCorr[nCorrPlots][nPtBins];
  TH2D * ZgCorr[nCorrPlots][nPtBins];
  
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

      histoName = "baseRg__" + zString[z] + betaString[b];
      // histoTitle = "R_{g} (#beta="+betaTitleString[b]+", z="+zTitleString[z]+") vs. Base R_{g}";
      // histoTitle += ";R_{g}^{base};R_{g} (#beta="+betaTitleString[b]+", z="+zTitleString[z]+");p_{T}^{Pythia} (GeV)";
      Rg_base_pt[b][z] = new TH3D( histoName, "", 30,0,0.6, 30,0,0.6, 60,0,60 );
      
      nZgBins = 60;  minZgBin = 0.0;  maxZgBin = 0.6;
      if ( z == 1 ) { nZgBins = 60;  minZgBin = 0.0;  maxZgBin = 0.6; }
      histoName = "baseZg__" + zString[z] + betaString[b];
      // histoTitle = "Z_{g} (#beta="+betaTitleString[b]+", z="+zTitleString[z]+") vs. Base Z_{g}";
      // histoTitle += ";Z_{g}^{base};Z_{g} (#beta="+betaTitleString[b]+", z="+zTitleString[z]+");p_{T}^{Pythia} (GeV)";
      Zg_base_pt[b][z] = new TH3D( histoName, "", nZgBins,minZgBin,maxZgBin, 30,0,0.6, 60,0,60 );

      histoName = "h_dPt_vs_pt__" + zString[z] + betaString[b];
      histoTitle = "#Delta p_{T} vs. p_{T}^{Pythia};p_{T}^{Pythia} (GeV);#Delta p_{T} (GeV)";
      h_dPt_vs_pt[b][z] = new TH2D( histoName, histoTitle, 60,0,60, 50,-40,10);

    }
  }

  TFile* outFile = new TFile( "out/sdParameterAnalysis.root", "RECREATE" );
  
  for ( int z=0; z<zVals; ++z ) {
    for ( int b=0; b<betaVals; ++b ) {
      
      if ( z==0 ) {    if ( b==0|| b==1 || b==4 || b==5  ) { continue; }    }
      if ( z==1 && b!=6) { continue; }
      
      nEntries = SDtree[b][z]->GetEntries();

      for ( int i=0; i<nEntries; ++i ) {
	SDtree[b][z]->GetEntry(i);
	ptVal = jetPt[b][z];
	groomZg = Zg[b][z];
	SDtree[3][0]->GetEntry(i);
	baseZg = Zg[3][0];
	// TLeaf* zgLeaf = SDtree[3][0]->GetLeaf("Zg"); cout << zgLeaf->GetValue() << "   " << Zg[3][0] <<endl;
	Zg_base_pt[b][z]->Fill( baseZg, groomZg, ptVal, wt[b][z] );

	dPt = sdPt[b][z] - jetPt[b][z];
	h_dPt_vs_pt[b][z]->Fill( ptVal, dPt, wt[b][z] );
	
	groomRg = Rg[b][z];
	SDtree[3][0]->GetEntry(i);
	baseRg = Rg[3][0];
	// TLeaf* rgLeaf = SDtree[3][0]->GetLeaf("Rg"); cout << rgLeaf->GetValue() << "   " << Rg[3][0] <<endl;
	Rg_base_pt[b][z]->Fill( baseRg, groomRg, ptVal, wt[b][z] );
      }

      Rg_base_pt[b][z]->Write();
      Zg_base_pt[b][z]->Write();
      h_dPt_vs_pt[b][z]->Write();
      
    }
  }


  TCanvas *dPtCanvas = new TCanvas( "dPtCanvas" , "" ,0 ,23 ,1280 ,700 );
  TH1D * dPtHisto = new TH1D( "dPtHisto", "", 25,-20,5 );
  for ( int z=0; z<zVals; ++z ) {
    for ( int b=0; b<betaVals; ++b ) {

      if ( z==0 ) {    if ( b==0|| b==1 || b==4 || b==5  ) { continue; }    }
      if ( z==1 && b!=6) { continue; }
      dPtLeg[b][z] = new TLegend(0.14,0.665,0.293,0.886,NULL,"brNDC");    // LEGEND
      dPtLeg[b][z]->SetBorderSize(1);   dPtLeg[b][z]->SetLineColor(1);   dPtLeg[b][z]->SetLineStyle(1);
      dPtLeg[b][z]->SetLineWidth(1);   dPtLeg[b][z]->SetFillColor(0);   dPtLeg[b][z]->SetFillStyle(1001);
      
      gStyle->SetOptStat(0);
      dPtCanvas->SetLogy();
      histoTitle = "#Delta p_{T}: #beta="+betaTitleString[b]+", z_{cut}="+zTitleString[z] + ";#Delta p_{T} (GeV)";
      dPtHisto->SetTitle( histoTitle );
      dPtHisto->SetMaximum(1.0);
      dPtHisto->Draw();
      
      TH1D* h_dPt[nPtBins];
      
      TString lpf = "lpf";
	
      for ( int i=0; i<nPtBins; ++i ) {

	histoName = "h_dPt__" + zString[z] + betaString[b] + "_" + ptBinName[i];

	h_dPt[i] = h_dPt_vs_pt[b][z]->ProjectionY( histoName, ptBinLo[i], ptBinHi[i] );
	h_dPt[i]->Scale( 1./h_dPt[i]->Integral("width") );
	h_dPt[i]->Write();

	h_dPt[i]->SetLineColor( ptColor[i] );
	h_dPt[i]->SetMarkerColor( ptColor[i] );
	h_dPt[i]->SetMarkerStyle( marker[z] );

	h_dPt[i]->Draw("SAME");

	NAME = "h_dPt__" + zString[z] + betaString[b] + "_" + ptBinName[i];
	TITLE = ptBinString[i];
	
	dPtEntry=dPtLeg[b][z]->AddEntry( NAME, TITLE, lpf );                            // ADD TO LEGEND
	dPtEntry->SetLineColor( ptColor[i] );   dPtEntry->SetMarkerColor( ptColor[i] );
	dPtEntry->SetFillStyle(1001);   dPtEntry->SetTextFont(42);
	dPtEntry->SetLineStyle(1);   dPtEntry->SetLineWidth(2);
	dPtEntry->SetMarkerStyle( marker[z] );   dPtEntry->SetMarkerSize(1);
	lpf += "lpf";
      }

      dPtLeg[b][z]->Draw();
      saveName = "pdfs/dPt__" + zString[z] + betaString[b] + ".pdf";
      dPtCanvas->SaveAs( saveName, "PDF" );
      //dPtCanvas->Clear(); 
    }
  }

  // int cindex = 0;
  // TCanvas* baseCanvas = new TCanvas( "bc" , "" ,0 ,23 ,1280 ,700 );
  
  // for ( int z=0; z<zVals; ++z ) {
  //   for ( int b=0; b<betaVals; ++b ) {

  //     if ( z==0 ) {    if ( b==0|| b==1 || b==4 || b==5  ) { continue; }    }
  //     if ( z==1 && b!=6) { continue; }

  //     for ( int i=0; i<nPtBins; ++i ) {

  // 	Rg_base_pt[b][z]->SetAxisRange( ptBinLo[i], ptBinHi[i], "z" );
  // 	RgCorr[cindex][i] = (TH2D*) Rg_base_pt[b][z]->Project3D("YX");
  // 	histoName = "baseRg__" + zString[z] + betaString[b] + "_" + ptBinName[i];
  // 	histoTitle = "R_{g} (#beta="+betaTitleString[b]+", z="+zTitleString[z]+") vs. Base R_{g} :  " + ptBinString[i];
  // 	histoTitle += ";R_{g}^{base};R_{g} (#beta="+betaTitleString[b]+", z="+zTitleString[z]+");p_{T}^{Pythia} (GeV)";
  // 	RgCorr[cindex][i]->SetName( histoName );
  // 	RgCorr[cindex][i]->SetTitle( histoTitle );
  // 	RgCorr[cindex][i]->Write();
  // 	RgCorr[cindex][i]->Draw("COLZ");
  // 	saveName = "pdfs/" + histoName + ".pdf";
  // 	baseCanvas->SaveAs( saveName, "PDF");

  // 	Zg_base_pt[b][z]->SetAxisRange( ptBinLo[i], ptBinHi[i], "z" );
  // 	ZgCorr[cindex][i] = (TH2D*) Zg_base_pt[b][z]->Project3D("YX");
  // 	histoName = "baseZg__" + zString[z] + betaString[b] + "_" + ptBinName[i];
  // 	histoTitle = "Z_{g} (#beta="+betaTitleString[b]+", z="+zTitleString[z]+") vs. Base Z_{g} :  " + ptBinString[i];
  // 	histoTitle += ";Z_{g}^{base};Z_{g} (#beta="+betaTitleString[b]+", z="+zTitleString[z]+");p_{T}^{Pythia} (GeV)";
  // 	ZgCorr[cindex][i]->SetName( histoName );
  // 	ZgCorr[cindex][i]->SetTitle( histoTitle );
  // 	ZgCorr[cindex][i]->Write();
  // 	ZgCorr[cindex][i]->Draw("COLZ");
  // 	saveName = "pdfs/" + histoName + ".pdf";
  // 	baseCanvas->SaveAs( saveName, "PDF");
	
  //     }
  //     cindex++;
  //   }
  // }

  // baseCanvas->Close();

  /*     
	 for ( int i=0; i<nPtBins; ++i ) {


	 //      Rg PLOTS     
	 for ( int z=0; z<zVals; ++z ) {

	 rleg[i][z] = new TLegend(0.75, 0.68, 0.9, 0.9,NULL,"brNDC");    // LEGEND
	 rleg[i][z]->SetBorderSize(1);   rleg[i][z]->SetLineColor(1);   rleg[i][z]->SetLineStyle(1);   rleg[i][z]->SetLineWidth(1);   rleg[i][z]->SetFillColor(0);   rleg[i][z]->SetFillStyle(1001);
      
	 TString lpf = "lpf";

	 canvasName = "Rg_" + zString[z] + ptBinName[i] +"GeV";
	 canvasTitle = "Rg: " + ptBinString[i] + ", z_{cut} = " + zTitleString[z];
	 r1[i][z] = new TCanvas( canvasName , canvasTitle ,0 ,23 ,1280 ,709 );
      
	 for ( int b=0; b<betaVals; ++b ) {

	 gStyle->SetOptTitle(1);
	 TITLE = "Rg: " + ptBinString[i] + ",  z_{cut} = " + zTitleString[z] + ";R_{g}";
	 h_Rg[i][b][z]->SetTitle( TITLE );
	
	 h_Rg[i][b][z]->SetLineColor( color[b] );
	 h_Rg[i][b][z]->SetMarkerColor( color[b] );	
	 h_Rg[i][b][z]->SetMarkerStyle( marker[z] );
	 // h_Rg[i][b][z]->Set();
	 // h_Rg[i][b][z]->Set();
	 h_Rg[i][b][z]->Draw("SAME");
	 NAME = "Rg_" + zString[z] + betaString[b] + ptBinName[i] +"GeV";
	 TITLE = "#beta = " + betaTitleString[b];
	 rentry=rleg[i][z]->AddEntry( NAME, TITLE, lpf );                            // ADD TO LEGEND
	 rentry->SetLineColor( color[b] );   rentry->SetMarkerColor( color[b] );
	 rentry->SetFillStyle(1001);   rentry->SetTextFont(42);
	 rentry->SetLineStyle(1);   rentry->SetLineWidth(2);
	 rentry->SetMarkerStyle( marker[z] );   rentry->SetMarkerSize(1);
	 lpf += "lpf";
	 }
      
	 rleg[i][z]->Draw();
	 r1[i][z]->Modified();
	 r1[i][z]->cd();
	 r1[i][z]->SetSelected( r1[i][z] );
	 saveName = "pdfs/" + NAME + ".pdf";
	 r1[i][z]->SaveAs( saveName, "PDF" );
	 }
	 // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~



	 //      Zg PLOTS     
	 for ( int z=0; z<zVals; ++z ) {

	 zleg[i][z] = new TLegend(0.75, 0.68, 0.9, 0.9,NULL,"brNDC");    // LEGEND
	 zleg[i][z]->SetBorderSize(1);   zleg[i][z]->SetLineColor(1);   zleg[i][z]->SetLineStyle(1);   zleg[i][z]->SetLineWidth(1);   zleg[i][z]->SetFillColor(0);   zleg[i][z]->SetFillStyle(1001);
      
	 TString lpf = "lpf";

	 canvasName = "Zg_" + zString[z] + ptBinName[i] +"GeV";
	 canvasTitle = "Zg: " + ptBinString[i] + ", z_{cut} = " + zTitleString[z];
	 z1[i][z] = new TCanvas( canvasName , canvasTitle ,0 ,23 ,1280 ,709 );
      
	 for ( int b=0; b<betaVals; ++b ) {

	 gStyle->SetOptTitle(1);
	 TITLE = "Zg: " + ptBinString[i] + ",  z_{cut} = " + zTitleString[z] + ";Z_{g}";
	 h_Zg[i][b][z]->SetTitle( TITLE );
	
	 h_Zg[i][b][z]->SetLineColor( color[b] );
	 h_Zg[i][b][z]->SetMarkerColor( color[b] );	
	 h_Zg[i][b][z]->SetMarkerStyle( marker[z] );
	 // h_Zg[i][b][z]->Set();
	 // h_Zg[i][b][z]->Set();
	 h_Zg[i][b][z]->Draw("SAME");
	 NAME = "Zg_" + zString[z] + betaString[b] + ptBinName[i] +"GeV";
	 TITLE = "#beta = " + betaTitleString[b];
	 zentry=zleg[i][z]->AddEntry( NAME, TITLE, lpf );                            // ADD TO LEGEND
	 zentry->SetLineColor( color[b] );   zentry->SetMarkerColor( color[b] );
	 zentry->SetFillStyle(1001);   zentry->SetTextFont(42);
	 zentry->SetLineStyle(1);   zentry->SetLineWidth(2);
	 zentry->SetMarkerStyle( marker[z] );   zentry->SetMarkerSize(1);
	 lpf += "lpf";
	 }
      
	 zleg[i][z]->Draw();
	 z1[i][z]->Modified();
	 z1[i][z]->cd();
	 z1[i][z]->SetSelected( z1[i][z] );
	 saveName = "pdfs/" + NAME + ".pdf";
	 z1[i][z]->SaveAs( saveName, "PDF" );
	 }
	 // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
	 }
  */

  
  inFile->Close();    //outFile->Close();
      
}
