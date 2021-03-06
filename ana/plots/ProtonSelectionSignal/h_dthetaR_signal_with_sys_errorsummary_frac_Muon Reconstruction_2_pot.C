{
//=========Macro generated from canvas: cNeutA/Neutron Angulars
//=========  (Thu Jan 30 11:35:14 2020) by ROOT version5.34/36
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
   Double_t xAxis437[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *h_tmp_err_errSum_4271__429 = new TH1D("h_tmp_err_errSum_4271__429","dthetaP:ptmu",18, xAxis437);
   h_tmp_err_errSum_4271__429->SetMinimum(1e-05);
   h_tmp_err_errSum_4271__429->SetMaximum(0.4);
   h_tmp_err_errSum_4271__429->SetDirectory(0);
   h_tmp_err_errSum_4271__429->SetStats(0);
   h_tmp_err_errSum_4271__429->SetLineWidth(3);
   h_tmp_err_errSum_4271__429->GetXaxis()->CenterTitle(true);
   h_tmp_err_errSum_4271__429->GetXaxis()->SetNdivisions(509);
   h_tmp_err_errSum_4271__429->GetXaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__429->GetXaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__429->GetXaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__429->GetXaxis()->SetTitleOffset(1.15);
   h_tmp_err_errSum_4271__429->GetYaxis()->SetTitle("Fractional Uncertainty");
   h_tmp_err_errSum_4271__429->GetYaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__429->GetYaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__429->GetYaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__429->GetYaxis()->SetTitleOffset(1.2);
   h_tmp_err_errSum_4271__429->GetZaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__429->GetZaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__429->GetZaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__429->GetZaxis()->SetTitleOffset(0.75);
   h_tmp_err_errSum_4271__429->Draw("HIST");
   Double_t xAxis438[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *tmp_vertError_Flux__430 = new TH1D("tmp_vertError_Flux__430","h_dthetaP_ptmu_data_2track_weighted_nobck_micblob_projX_Flux",18, xAxis438);
   tmp_vertError_Flux__430->SetBinContent(1,0.01823313);
   tmp_vertError_Flux__430->SetBinContent(2,0.009534842);
   tmp_vertError_Flux__430->SetBinContent(3,0.01005605);
   tmp_vertError_Flux__430->SetBinContent(4,0.00948148);
   tmp_vertError_Flux__430->SetBinContent(5,0.03690299);
   tmp_vertError_Flux__430->SetBinContent(6,0.01604169);
   tmp_vertError_Flux__430->SetBinContent(7,0.05619378);
   tmp_vertError_Flux__430->SetBinContent(8,0.08110101);
   tmp_vertError_Flux__430->SetBinContent(9,0.01066144);
   tmp_vertError_Flux__430->SetBinContent(10,0.00915663);
   tmp_vertError_Flux__430->SetBinContent(11,0.01682011);
   tmp_vertError_Flux__430->SetBinContent(12,0.01777313);
   tmp_vertError_Flux__430->SetBinContent(13,0.0336798);
   tmp_vertError_Flux__430->SetBinContent(14,0.02372045);
   tmp_vertError_Flux__430->SetBinContent(15,0.009516493);
   tmp_vertError_Flux__430->SetBinContent(16,0.009864455);
   tmp_vertError_Flux__430->SetBinContent(17,0.004337403);
   tmp_vertError_Flux__430->SetBinContent(18,0.008003298);
   tmp_vertError_Flux__430->SetEntries(20);
   tmp_vertError_Flux__430->SetDirectory(0);
   tmp_vertError_Flux__430->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#9966ff");
   tmp_vertError_Flux__430->SetLineColor(ci);
   tmp_vertError_Flux__430->SetLineWidth(3);
   tmp_vertError_Flux__430->GetXaxis()->CenterTitle(true);
   tmp_vertError_Flux__430->GetXaxis()->SetNdivisions(509);
   tmp_vertError_Flux__430->GetXaxis()->SetLabelFont(42);
   tmp_vertError_Flux__430->GetXaxis()->SetLabelSize(0.05);
   tmp_vertError_Flux__430->GetXaxis()->SetTitleSize(0.06);
   tmp_vertError_Flux__430->GetXaxis()->SetTitleOffset(1.15);
   tmp_vertError_Flux__430->GetYaxis()->SetLabelFont(42);
   tmp_vertError_Flux__430->GetYaxis()->SetLabelSize(0.05);
   tmp_vertError_Flux__430->GetYaxis()->SetTitleSize(0.06);
   tmp_vertError_Flux__430->GetYaxis()->SetTitleOffset(1.2);
   tmp_vertError_Flux__430->GetZaxis()->SetLabelFont(42);
   tmp_vertError_Flux__430->GetZaxis()->SetLabelSize(0.05);
   tmp_vertError_Flux__430->GetZaxis()->SetTitleSize(0.06);
   tmp_vertError_Flux__430->GetZaxis()->SetTitleOffset(0.75);
   tmp_vertError_Flux__430->Draw("HIST SAME");
   
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
   Double_t xAxis439[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *h_tmp_err_errSum_4271__431 = new TH1D("h_tmp_err_errSum_4271__431","dthetaP:ptmu",18, xAxis439);
   h_tmp_err_errSum_4271__431->SetMinimum(1e-05);
   h_tmp_err_errSum_4271__431->SetMaximum(0.4);
   h_tmp_err_errSum_4271__431->SetDirectory(0);
   h_tmp_err_errSum_4271__431->SetStats(0);
   h_tmp_err_errSum_4271__431->SetLineWidth(3);
   h_tmp_err_errSum_4271__431->GetXaxis()->CenterTitle(true);
   h_tmp_err_errSum_4271__431->GetXaxis()->SetNdivisions(509);
   h_tmp_err_errSum_4271__431->GetXaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__431->GetXaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__431->GetXaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__431->GetXaxis()->SetTitleOffset(1.15);
   h_tmp_err_errSum_4271__431->GetYaxis()->SetTitle("Fractional Uncertainty");
   h_tmp_err_errSum_4271__431->GetYaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__431->GetYaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__431->GetYaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__431->GetYaxis()->SetTitleOffset(1.2);
   h_tmp_err_errSum_4271__431->GetZaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__431->GetZaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__431->GetZaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__431->GetZaxis()->SetTitleOffset(0.75);
   h_tmp_err_errSum_4271__431->Draw("sameaxis");
   Double_t xAxis440[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *h_tmp_err_errSum_4271__432 = new TH1D("h_tmp_err_errSum_4271__432","dthetaP:ptmu",18, xAxis440);
   h_tmp_err_errSum_4271__432->SetMinimum(1e-05);
   h_tmp_err_errSum_4271__432->SetMaximum(0.4);
   h_tmp_err_errSum_4271__432->SetDirectory(0);
   h_tmp_err_errSum_4271__432->SetStats(0);
   h_tmp_err_errSum_4271__432->SetLineWidth(3);
   h_tmp_err_errSum_4271__432->GetXaxis()->CenterTitle(true);
   h_tmp_err_errSum_4271__432->GetXaxis()->SetNdivisions(509);
   h_tmp_err_errSum_4271__432->GetXaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__432->GetXaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__432->GetXaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__432->GetXaxis()->SetTitleOffset(1.15);
   h_tmp_err_errSum_4271__432->GetYaxis()->SetTitle("Fractional Uncertainty");
   h_tmp_err_errSum_4271__432->GetYaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__432->GetYaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__432->GetYaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__432->GetYaxis()->SetTitleOffset(1.2);
   h_tmp_err_errSum_4271__432->GetZaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__432->GetZaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__432->GetZaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__432->GetZaxis()->SetTitleOffset(0.75);
   h_tmp_err_errSum_4271__432->Draw("sameaxis");
   Double_t xAxis441[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *h_tmp_err_errSum_4271__433 = new TH1D("h_tmp_err_errSum_4271__433","dthetaP:ptmu",18, xAxis441);
   h_tmp_err_errSum_4271__433->SetMinimum(1e-05);
   h_tmp_err_errSum_4271__433->SetMaximum(0.4);
   h_tmp_err_errSum_4271__433->SetDirectory(0);
   h_tmp_err_errSum_4271__433->SetStats(0);
   h_tmp_err_errSum_4271__433->SetLineWidth(3);
   h_tmp_err_errSum_4271__433->GetXaxis()->CenterTitle(true);
   h_tmp_err_errSum_4271__433->GetXaxis()->SetNdivisions(509);
   h_tmp_err_errSum_4271__433->GetXaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__433->GetXaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__433->GetXaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__433->GetXaxis()->SetTitleOffset(1.15);
   h_tmp_err_errSum_4271__433->GetYaxis()->SetTitle("Fractional Uncertainty");
   h_tmp_err_errSum_4271__433->GetYaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__433->GetYaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__433->GetYaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__433->GetYaxis()->SetTitleOffset(1.2);
   h_tmp_err_errSum_4271__433->GetZaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__433->GetZaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__433->GetZaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__433->GetZaxis()->SetTitleOffset(0.75);
   h_tmp_err_errSum_4271__433->Draw("sameaxis");
   Double_t xAxis442[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *h_tmp_err_errSum_4271__434 = new TH1D("h_tmp_err_errSum_4271__434","dthetaP:ptmu",18, xAxis442);
   h_tmp_err_errSum_4271__434->SetMinimum(1e-05);
   h_tmp_err_errSum_4271__434->SetMaximum(0.4);
   h_tmp_err_errSum_4271__434->SetDirectory(0);
   h_tmp_err_errSum_4271__434->SetStats(0);
   h_tmp_err_errSum_4271__434->SetLineWidth(3);
   h_tmp_err_errSum_4271__434->GetXaxis()->CenterTitle(true);
   h_tmp_err_errSum_4271__434->GetXaxis()->SetNdivisions(509);
   h_tmp_err_errSum_4271__434->GetXaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__434->GetXaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__434->GetXaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__434->GetXaxis()->SetTitleOffset(1.15);
   h_tmp_err_errSum_4271__434->GetYaxis()->SetTitle("Fractional Uncertainty");
   h_tmp_err_errSum_4271__434->GetYaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__434->GetYaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__434->GetYaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__434->GetYaxis()->SetTitleOffset(1.2);
   h_tmp_err_errSum_4271__434->GetZaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__434->GetZaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__434->GetZaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__434->GetZaxis()->SetTitleOffset(0.75);
   h_tmp_err_errSum_4271__434->Draw("sameaxis");
   Double_t xAxis443[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *h_tmp_err_errSum_4271__435 = new TH1D("h_tmp_err_errSum_4271__435","dthetaP:ptmu",18, xAxis443);
   h_tmp_err_errSum_4271__435->SetMinimum(1e-05);
   h_tmp_err_errSum_4271__435->SetMaximum(0.4);
   h_tmp_err_errSum_4271__435->SetDirectory(0);
   h_tmp_err_errSum_4271__435->SetStats(0);
   h_tmp_err_errSum_4271__435->SetLineWidth(3);
   h_tmp_err_errSum_4271__435->GetXaxis()->CenterTitle(true);
   h_tmp_err_errSum_4271__435->GetXaxis()->SetNdivisions(509);
   h_tmp_err_errSum_4271__435->GetXaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__435->GetXaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__435->GetXaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__435->GetXaxis()->SetTitleOffset(1.15);
   h_tmp_err_errSum_4271__435->GetYaxis()->SetTitle("Fractional Uncertainty");
   h_tmp_err_errSum_4271__435->GetYaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__435->GetYaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__435->GetYaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__435->GetYaxis()->SetTitleOffset(1.2);
   h_tmp_err_errSum_4271__435->GetZaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__435->GetZaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__435->GetZaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__435->GetZaxis()->SetTitleOffset(0.75);
   h_tmp_err_errSum_4271__435->Draw("sameaxis");
   Double_t xAxis444[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *h_tmp_err_errSum_4271__436 = new TH1D("h_tmp_err_errSum_4271__436","dthetaP:ptmu",18, xAxis444);
   h_tmp_err_errSum_4271__436->SetMinimum(1e-05);
   h_tmp_err_errSum_4271__436->SetMaximum(0.4);
   h_tmp_err_errSum_4271__436->SetDirectory(0);
   h_tmp_err_errSum_4271__436->SetStats(0);
   h_tmp_err_errSum_4271__436->SetLineWidth(3);
   h_tmp_err_errSum_4271__436->GetXaxis()->CenterTitle(true);
   h_tmp_err_errSum_4271__436->GetXaxis()->SetNdivisions(509);
   h_tmp_err_errSum_4271__436->GetXaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__436->GetXaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__436->GetXaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__436->GetXaxis()->SetTitleOffset(1.15);
   h_tmp_err_errSum_4271__436->GetYaxis()->SetTitle("Fractional Uncertainty");
   h_tmp_err_errSum_4271__436->GetYaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__436->GetYaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__436->GetYaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__436->GetYaxis()->SetTitleOffset(1.2);
   h_tmp_err_errSum_4271__436->GetZaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__436->GetZaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__436->GetZaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__436->GetZaxis()->SetTitleOffset(0.75);
   h_tmp_err_errSum_4271__436->Draw("sameaxis");
   Double_t xAxis445[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *h_tmp_err_errSum_4271__437 = new TH1D("h_tmp_err_errSum_4271__437","dthetaP:ptmu",18, xAxis445);
   h_tmp_err_errSum_4271__437->SetMinimum(1e-05);
   h_tmp_err_errSum_4271__437->SetMaximum(0.4);
   h_tmp_err_errSum_4271__437->SetDirectory(0);
   h_tmp_err_errSum_4271__437->SetStats(0);
   h_tmp_err_errSum_4271__437->SetLineWidth(3);
   h_tmp_err_errSum_4271__437->GetXaxis()->CenterTitle(true);
   h_tmp_err_errSum_4271__437->GetXaxis()->SetNdivisions(509);
   h_tmp_err_errSum_4271__437->GetXaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__437->GetXaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__437->GetXaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__437->GetXaxis()->SetTitleOffset(1.15);
   h_tmp_err_errSum_4271__437->GetYaxis()->SetTitle("Fractional Uncertainty");
   h_tmp_err_errSum_4271__437->GetYaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__437->GetYaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__437->GetYaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__437->GetYaxis()->SetTitleOffset(1.2);
   h_tmp_err_errSum_4271__437->GetZaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__437->GetZaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__437->GetZaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__437->GetZaxis()->SetTitleOffset(0.75);
   h_tmp_err_errSum_4271__437->Draw("sameaxis");
   Double_t xAxis446[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *h_tmp_err_errSum_4271__438 = new TH1D("h_tmp_err_errSum_4271__438","dthetaP:ptmu",18, xAxis446);
   h_tmp_err_errSum_4271__438->SetMinimum(1e-05);
   h_tmp_err_errSum_4271__438->SetMaximum(0.4);
   h_tmp_err_errSum_4271__438->SetDirectory(0);
   h_tmp_err_errSum_4271__438->SetStats(0);
   h_tmp_err_errSum_4271__438->SetLineWidth(3);
   h_tmp_err_errSum_4271__438->GetXaxis()->CenterTitle(true);
   h_tmp_err_errSum_4271__438->GetXaxis()->SetNdivisions(509);
   h_tmp_err_errSum_4271__438->GetXaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__438->GetXaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__438->GetXaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__438->GetXaxis()->SetTitleOffset(1.15);
   h_tmp_err_errSum_4271__438->GetYaxis()->SetTitle("Fractional Uncertainty");
   h_tmp_err_errSum_4271__438->GetYaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__438->GetYaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__438->GetYaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__438->GetYaxis()->SetTitleOffset(1.2);
   h_tmp_err_errSum_4271__438->GetZaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__438->GetZaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__438->GetZaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__438->GetZaxis()->SetTitleOffset(0.75);
   h_tmp_err_errSum_4271__438->Draw("sameaxis");
   Double_t xAxis447[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *h_tmp_err_errSum_4271__439 = new TH1D("h_tmp_err_errSum_4271__439","dthetaP:ptmu",18, xAxis447);
   h_tmp_err_errSum_4271__439->SetMinimum(1e-05);
   h_tmp_err_errSum_4271__439->SetMaximum(0.4);
   h_tmp_err_errSum_4271__439->SetDirectory(0);
   h_tmp_err_errSum_4271__439->SetStats(0);
   h_tmp_err_errSum_4271__439->SetLineWidth(3);
   h_tmp_err_errSum_4271__439->GetXaxis()->CenterTitle(true);
   h_tmp_err_errSum_4271__439->GetXaxis()->SetNdivisions(509);
   h_tmp_err_errSum_4271__439->GetXaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__439->GetXaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__439->GetXaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__439->GetXaxis()->SetTitleOffset(1.15);
   h_tmp_err_errSum_4271__439->GetYaxis()->SetTitle("Fractional Uncertainty");
   h_tmp_err_errSum_4271__439->GetYaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__439->GetYaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__439->GetYaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__439->GetYaxis()->SetTitleOffset(1.2);
   h_tmp_err_errSum_4271__439->GetZaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__439->GetZaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__439->GetZaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__439->GetZaxis()->SetTitleOffset(0.75);
   h_tmp_err_errSum_4271__439->Draw("sameaxis");
   Double_t xAxis448[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *h_tmp_err_errSum_4271__440 = new TH1D("h_tmp_err_errSum_4271__440","dthetaP:ptmu",18, xAxis448);
   h_tmp_err_errSum_4271__440->SetMinimum(1e-05);
   h_tmp_err_errSum_4271__440->SetMaximum(0.4);
   h_tmp_err_errSum_4271__440->SetDirectory(0);
   h_tmp_err_errSum_4271__440->SetStats(0);
   h_tmp_err_errSum_4271__440->SetLineWidth(3);
   h_tmp_err_errSum_4271__440->GetXaxis()->CenterTitle(true);
   h_tmp_err_errSum_4271__440->GetXaxis()->SetNdivisions(509);
   h_tmp_err_errSum_4271__440->GetXaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__440->GetXaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__440->GetXaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__440->GetXaxis()->SetTitleOffset(1.15);
   h_tmp_err_errSum_4271__440->GetYaxis()->SetTitle("Fractional Uncertainty");
   h_tmp_err_errSum_4271__440->GetYaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__440->GetYaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__440->GetYaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__440->GetYaxis()->SetTitleOffset(1.2);
   h_tmp_err_errSum_4271__440->GetZaxis()->SetLabelFont(42);
   h_tmp_err_errSum_4271__440->GetZaxis()->SetLabelSize(0.05);
   h_tmp_err_errSum_4271__440->GetZaxis()->SetTitleSize(0.06);
   h_tmp_err_errSum_4271__440->GetZaxis()->SetTitleOffset(0.75);
   h_tmp_err_errSum_4271__440->Draw("sameaxis");
   cNeutA->Modified();
   cNeutA->cd();
   cNeutA->SetSelected(cNeutA);
}
