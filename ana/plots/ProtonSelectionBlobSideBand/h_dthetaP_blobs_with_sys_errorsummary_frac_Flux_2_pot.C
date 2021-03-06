{
//=========Macro generated from canvas: cNeutA/Neutron Angulars
//=========  (Thu Jan 30 11:35:13 2020) by ROOT version5.34/36
   TCanvas *cNeutA = new TCanvas("cNeutA", "Neutron Angulars",1,1,900,726);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cNeutA->SetHighLightColor(2);
   cNeutA->Range(-257.1429,-0.0789354,257.1429,0.4473672);
   cNeutA->SetFillColor(0);
   cNeutA->SetBorderMode(0);
   cNeutA->SetBorderSize(2);
   cNeutA->SetLeftMargin(0.15);
   cNeutA->SetRightMargin(0.15);
   cNeutA->SetBottomMargin(0.15);
   cNeutA->SetFrameLineWidth(2);
   cNeutA->SetFrameBorderMode(0);
   cNeutA->SetFrameLineWidth(2);
   cNeutA->SetFrameBorderMode(0);
   Double_t xAxis368[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *h_tmp_err_errSum_4271__360 = new TH1D("h_tmp_err_errSum_4271__360","dthetaP:ptmu",18, xAxis368);
   h_tmp_err_errSum_4271__360->SetMinimum(1e-05);
   h_tmp_err_errSum_4271__360->SetMaximum(0.4);
   h_tmp_err_errSum_4271__360->SetDirectory(0);
   h_tmp_err_errSum_4271__360->SetStats(0);
   h_tmp_err_errSum_4271__360->SetLineWidth(3);
   h_tmp_err_errSum_4271__360->GetXaxis()->CenterTitle(true);
   h_tmp_err_errSum_4271__360->GetXaxis()->SetNdivisions(509);
   h_tmp_err_errSum_4271__360->GetXaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__360->GetXaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__360->GetXaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__360->GetXaxis()->SetTitleOffset(1.15);
   h_tmp_err_errSum_4271__360->GetYaxis()->SetTitle("Fractional Uncertainty");
   h_tmp_err_errSum_4271__360->GetYaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__360->GetYaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__360->GetYaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__360->GetYaxis()->SetTitleOffset(1.2);
   h_tmp_err_errSum_4271__360->GetZaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__360->GetZaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__360->GetZaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__360->GetZaxis()->SetTitleOffset(0.75);
   h_tmp_err_errSum_4271__360->Draw("HIST");
   Double_t xAxis369[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *tmp_vertError_Flux__361 = new TH1D("tmp_vertError_Flux__361","h_dthetaP_ptmu_data_2track_weighted_nobck_blobs_projX_Flux",18, xAxis369);
   tmp_vertError_Flux__361->SetBinContent(1,0.01328932);
   tmp_vertError_Flux__361->SetBinContent(2,0.02499593);
   tmp_vertError_Flux__361->SetBinContent(3,0.03895124);
   tmp_vertError_Flux__361->SetBinContent(4,0.04592862);
   tmp_vertError_Flux__361->SetBinContent(5,0.03762732);
   tmp_vertError_Flux__361->SetBinContent(6,0.06406374);
   tmp_vertError_Flux__361->SetBinContent(7,2.727068);
   tmp_vertError_Flux__361->SetBinContent(8,0.1333004);
   tmp_vertError_Flux__361->SetBinContent(9,0.116224);
   tmp_vertError_Flux__361->SetBinContent(10,0.04003518);
   tmp_vertError_Flux__361->SetBinContent(11,0.1554535);
   tmp_vertError_Flux__361->SetBinContent(12,0.08450886);
   tmp_vertError_Flux__361->SetBinContent(13,0.0657859);
   tmp_vertError_Flux__361->SetBinContent(14,0.05264905);
   tmp_vertError_Flux__361->SetBinContent(15,0.0401827);
   tmp_vertError_Flux__361->SetBinContent(16,0.0205636);
   tmp_vertError_Flux__361->SetBinContent(17,0.0233965);
   tmp_vertError_Flux__361->SetBinContent(18,0.02111652);
   tmp_vertError_Flux__361->SetEntries(20);
   tmp_vertError_Flux__361->SetDirectory(0);
   tmp_vertError_Flux__361->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#9966ff");
   tmp_vertError_Flux__361->SetLineColor(ci);
   tmp_vertError_Flux__361->SetLineWidth(3);
   tmp_vertError_Flux__361->GetXaxis()->CenterTitle(true);
   tmp_vertError_Flux__361->GetXaxis()->SetNdivisions(509);
   tmp_vertError_Flux__361->GetXaxis()->SetLabelFont(42);
   tmp_vertError_Flux__361->GetXaxis()->SetLabelSize(0.05);
   tmp_vertError_Flux__361->GetXaxis()->SetTitleSize(0.06);
   tmp_vertError_Flux__361->GetXaxis()->SetTitleOffset(1.15);
   tmp_vertError_Flux__361->GetYaxis()->SetLabelFont(42);
   tmp_vertError_Flux__361->GetYaxis()->SetLabelSize(0.05);
   tmp_vertError_Flux__361->GetYaxis()->SetTitleSize(0.06);
   tmp_vertError_Flux__361->GetYaxis()->SetTitleOffset(1.2);
   tmp_vertError_Flux__361->GetZaxis()->SetLabelFont(42);
   tmp_vertError_Flux__361->GetZaxis()->SetLabelSize(0.05);
   tmp_vertError_Flux__361->GetZaxis()->SetTitleSize(0.06);
   tmp_vertError_Flux__361->GetZaxis()->SetTitleOffset(0.75);
   tmp_vertError_Flux__361->Draw("HIST SAME");
   
   TLegend *leg = new TLegend(0.18,0.85,0.38,0.89,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("tmp_vertError_Flux","Flux","l");

   ci = TColor::GetColor("#9966ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   leg->Draw();
   Double_t xAxis370[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *h_tmp_err_errSum_4271__362 = new TH1D("h_tmp_err_errSum_4271__362","dthetaP:ptmu",18, xAxis370);
   h_tmp_err_errSum_4271__362->SetMinimum(1e-05);
   h_tmp_err_errSum_4271__362->SetMaximum(0.4);
   h_tmp_err_errSum_4271__362->SetDirectory(0);
   h_tmp_err_errSum_4271__362->SetStats(0);
   h_tmp_err_errSum_4271__362->SetLineWidth(3);
   h_tmp_err_errSum_4271__362->GetXaxis()->CenterTitle(true);
   h_tmp_err_errSum_4271__362->GetXaxis()->SetNdivisions(509);
   h_tmp_err_errSum_4271__362->GetXaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__362->GetXaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__362->GetXaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__362->GetXaxis()->SetTitleOffset(1.15);
   h_tmp_err_errSum_4271__362->GetYaxis()->SetTitle("Fractional Uncertainty");
   h_tmp_err_errSum_4271__362->GetYaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__362->GetYaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__362->GetYaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__362->GetYaxis()->SetTitleOffset(1.2);
   h_tmp_err_errSum_4271__362->GetZaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__362->GetZaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__362->GetZaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__362->GetZaxis()->SetTitleOffset(0.75);
   h_tmp_err_errSum_4271__362->Draw("sameaxis");
   cNeutA->Modified();
   cNeutA->cd();
   cNeutA->SetSelected(cNeutA);
}
