
void RgPlot() {
  const float pi = 3.141592;
  double R = 0.4;
  const int nPtBins = 5;
  TString NAME, TITLE;
  int nEntries, index, EventID, nCons;
  double jetPt, jetEta, jetPhi, jetE, wt, Rg, Zg, div, scale;
  TString lpf = "lpf";
  
  TFile* inFile = new TFile( "out/softDropAnalysis_z_01__b_00.root", "READ" );
  TTree * Jtree = (TTree*) inFile->Get("jetTree");
  TTree * SDtree = (TTree*) inFile->Get("sdTree");

  Jtree->SetBranchAddress("jetPt",&jetPt);  Jtree->SetBranchAddress("jetEta",&jetEta);  Jtree->SetBranchAddress("jetPhi",&jetPhi);  Jtree->SetBranchAddress("jetE",&jetE);
  Jtree->SetBranchAddress("wt",&wt);  Jtree->SetBranchAddress("EventID",&EventID);  Jtree->SetBranchAddress("nCons",&nCons);  SDtree->SetBranchAddress("jetPt",&jetPt);
  SDtree->SetBranchAddress("jetEta",&jetEta);  SDtree->SetBranchAddress("jetPhi",&jetPhi);  SDtree->SetBranchAddress("jetE",&jetE);  SDtree->SetBranchAddress("wt",&wt);
  SDtree->SetBranchAddress("EventID",&EventID);  SDtree->SetBranchAddress("nCons",&nCons);  SDtree->SetBranchAddress("Rg",&Rg);  SDtree->SetBranchAddress("Zg",&Zg);
  
  double ptBinLo[nPtBins] = { 10, 15, 20, 30, 40 };
  double ptBinHi[nPtBins] = { 15, 20, 30, 40, 70 };
  TString ptBinString[nPtBins] = { "10-15 GeV", "15-20 GeV", "20-30 GeV", "30-40 GeV", "40+ GeV" };
  TString ptBinName[nPtBins] = { "_10_15", "_15_20", "_20_30", "_30_40", "_40" };
  int color[nPtBins] = { 100, 93, 81, 63, 51 };

  TH2D * h_RgVSpt = new TH2D( "h_RgVSpt", "R_{g};p_{T} (GeV);R_{g}",70,0,70,20,0,0.5);
  TH1D * h_Rg[nPtBins];

  nEntries = SDtree->GetEntries();
  
  for ( int i=0; i<nEntries; ++i ) {     // FILL 2D HISTOGRAM
    SDtree->GetEntry(i);
    Jtree->GetEntry(i);
    h_RgVSpt->Fill( jetPt, Rg, wt );
  }
  
  TCanvas * c1 = new TCanvas( "c1" , "" ,0 ,23 ,1280 ,700 );       gStyle->SetOptStat(0);   // CANVAS

  TLegend *leg = new TLegend(0.75, 0.68, 0.9, 0.9,NULL,"brNDC");    // LEGEND
  leg->SetBorderSize(1);   leg->SetLineColor(1);   leg->SetLineStyle(1);   leg->SetLineWidth(1);   leg->SetFillColor(0);   leg->SetFillStyle(1001);
  TLegendEntry *entry;

  for ( int i=0; i<nPtBins; ++i ) {
    NAME = "h_Rg" + ptBinName[i];
    TITLE = ptBinString[i];
    h_Rg[i] = h_RgVSpt->ProjectionY( NAME , ptBinLo[i], ptBinHi[i] );       // PROJECT
    h_Rg[i]->SetAxisRange(0,2,"x");
    h_Rg[i]->Scale( 1./h_Rg[i]->Integral() );                     // NORMALIZE
    h_Rg[i]->SetLineWidth(2);
    h_Rg[i]->SetLineColor( color[i] );
    h_Rg[i]->SetMarkerColor( color[i] );
    h_Rg[i]->Draw("SAME");                                                    // DRAW
    entry=leg->AddEntry( NAME, TITLE, lpf );                            // ADD TO LEGEND
    entry->SetLineColor( color[i] );   entry->SetMarkerColor( color[i] );
    entry->SetFillStyle(1001);   entry->SetTextFont(42);
    entry->SetLineStyle(1);   entry->SetLineWidth(2);
    entry->SetMarkerStyle(1);   entry->SetMarkerSize(1);
    lpf += "lpf";
      }

  leg->Draw();
  c1->Modified();
  c1->cd();
  c1->SetSelected(c1);
  c1->SaveAs("pdfs/RgPlot.pdf","PDF");

}
