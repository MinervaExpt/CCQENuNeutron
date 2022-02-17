{
  TFile* f = new TFile("../test_muon_ela_SF_neutronW.root","read");
  TH2D* h_qelike = (TH2D*) f->Get("h_nHits_totalE_qelike");
  TH2D* h_qelikenot = (TH2D*) f->Get("h_nHits_totalE_qelikenot");
  h_qelike->SetMarkerSize(0.2);
  h_qelike->SetMarkerColor(kRed);
  h_qelikenot->SetMarkerSize(0.1);
  h_qelikenot->SetMarkerColor(kBlue);

  h_qelikenot->GetXaxis()->SetTitle("nHits");
  h_qelikenot->GetYaxis()->SetTitle("TotalE");
  h_qelike->GetXaxis()->SetTitle("nHits");
  h_qelike->GetYaxis()->SetTitle("TotalE");
  h_qelikenot->Draw("p");
  h_qelike->Draw("psame");

  TLine* line = new TLine(0,0,100,1000);
  line->Draw();

  TLegend* leg = new TLegend(.7,.7,.9,.9);
  leg->AddEntry(h_qelike,"QELike");
  leg->AddEntry(h_qelikenot,"Not QELike");
  c1->SetGrid();
  c1->Print("hHits_TotalE.pdf");


  h_qelikenot->Rebin2D(2,2);
  h_qelike->Rebin2D(2,2);
  h_qelike->Divide(h_qelikenot);
  h_qelike->Draw("colz");
  //h_qelike->GetZaxis()->SetRangeUser(0,2);
  gStyle->SetPalette(53);
  PlotUtils::MnvPlotter::SetROOT6Palette( 70 );
  line->SetLineColor(kRed);
  TLatex* latex = new TLatex();
  latex->DrawLatexNDC(.3,.95,"QELike/QELikeNot");
  line->Draw();
  c1->SetLogz();
  c1->Print("hHits_TotalE_ratio.pdf");
}
