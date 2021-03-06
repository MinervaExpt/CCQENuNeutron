{
//=========Macro generated from canvas: cNeutA/Neutron Angulars
//=========  (Thu Jan 30 11:35:12 2020) by ROOT version5.34/36
   TCanvas *cNeutA = new TCanvas("cNeutA", "Neutron Angulars",1,1,900,726);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cNeutA->SetHighLightColor(2);
   cNeutA->Range(-257.1429,-556.4506,257.1429,3153.22);
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
   Double_t xAxis146[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *tmpMC__138 = new TH1D("tmpMC__138","dthetaP:ptmu",18, xAxis146);
   tmpMC__138->SetBinContent(1,66.48904);
   tmpMC__138->SetBinContent(2,82.19329);
   tmpMC__138->SetBinContent(3,130.6603);
   tmpMC__138->SetBinContent(4,211.2014);
   tmpMC__138->SetBinContent(5,339.661);
   tmpMC__138->SetBinContent(6,528.3071);
   tmpMC__138->SetBinContent(7,846.8248);
   tmpMC__138->SetBinContent(8,1143.128);
   tmpMC__138->SetBinContent(9,1645.034);
   tmpMC__138->SetBinContent(10,1688.734);
   tmpMC__138->SetBinContent(11,1168.084);
   tmpMC__138->SetBinContent(12,684.7981);
   tmpMC__138->SetBinContent(13,467.1951);
   tmpMC__138->SetBinContent(14,372.0023);
   tmpMC__138->SetBinContent(15,225.3816);
   tmpMC__138->SetBinContent(16,138.628);
   tmpMC__138->SetBinContent(17,79.92888);
   tmpMC__138->SetBinContent(18,64.80535);
   tmpMC__138->SetBinError(1,20.09539);
   tmpMC__138->SetBinError(2,25.46587);
   tmpMC__138->SetBinError(3,44.2094);
   tmpMC__138->SetBinError(4,76.24735);
   tmpMC__138->SetBinError(5,138.27);
   tmpMC__138->SetBinError(6,197.6981);
   tmpMC__138->SetBinError(7,311.7293);
   tmpMC__138->SetBinError(8,365.3049);
   tmpMC__138->SetBinError(9,445.3938);
   tmpMC__138->SetBinError(10,458.9698);
   tmpMC__138->SetBinError(11,361.4284);
   tmpMC__138->SetBinError(12,283.8341);
   tmpMC__138->SetBinError(13,178.1318);
   tmpMC__138->SetBinError(14,152.6432);
   tmpMC__138->SetBinError(15,81.05067);
   tmpMC__138->SetBinError(16,47.61911);
   tmpMC__138->SetBinError(17,24.93476);
   tmpMC__138->SetBinError(18,19.12344);
   tmpMC__138->SetMinimum(0);
   tmpMC__138->SetMaximum(2819.349);
   tmpMC__138->SetEntries(398.3572);
   tmpMC__138->SetDirectory(0);
   tmpMC__138->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ffcccc");
   tmpMC__138->SetFillColor(ci);

   ci = TColor::GetColor("#000099");
   tmpMC__138->SetLineColor(ci);
   tmpMC__138->SetLineWidth(2);
   tmpMC__138->SetMarkerStyle(0);
   tmpMC__138->GetXaxis()->CenterTitle(true);
   tmpMC__138->GetXaxis()->SetNdivisions(509);
   tmpMC__138->GetXaxis()->SetLabelFont(42);
   tmpMC__138->GetXaxis()->SetLabelSize(0.05);
   tmpMC__138->GetXaxis()->SetTitleSize(0.06);
   tmpMC__138->GetXaxis()->SetTitleOffset(1.15);
   tmpMC__138->GetYaxis()->SetLabelFont(42);
   tmpMC__138->GetYaxis()->SetLabelSize(0.05);
   tmpMC__138->GetYaxis()->SetTitleSize(0.06);
   tmpMC__138->GetYaxis()->SetTitleOffset(1.2);
   tmpMC__138->GetZaxis()->SetLabelFont(42);
   tmpMC__138->GetZaxis()->SetLabelSize(0.05);
   tmpMC__138->GetZaxis()->SetTitleSize(0.06);
   tmpMC__138->GetZaxis()->SetTitleOffset(0.75);
   tmpMC__138->Draw("E2");
   Double_t xAxis147[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *tmpMC__139 = new TH1D("tmpMC__139","dthetaP:ptmu",18, xAxis147);
   tmpMC__139->SetBinContent(1,66.48904);
   tmpMC__139->SetBinContent(2,82.19329);
   tmpMC__139->SetBinContent(3,130.6603);
   tmpMC__139->SetBinContent(4,211.2014);
   tmpMC__139->SetBinContent(5,339.661);
   tmpMC__139->SetBinContent(6,528.3071);
   tmpMC__139->SetBinContent(7,846.8248);
   tmpMC__139->SetBinContent(8,1143.128);
   tmpMC__139->SetBinContent(9,1645.034);
   tmpMC__139->SetBinContent(10,1688.734);
   tmpMC__139->SetBinContent(11,1168.084);
   tmpMC__139->SetBinContent(12,684.7981);
   tmpMC__139->SetBinContent(13,467.1951);
   tmpMC__139->SetBinContent(14,372.0023);
   tmpMC__139->SetBinContent(15,225.3816);
   tmpMC__139->SetBinContent(16,138.628);
   tmpMC__139->SetBinContent(17,79.92888);
   tmpMC__139->SetBinContent(18,64.80535);
   tmpMC__139->SetBinError(1,20.09539);
   tmpMC__139->SetBinError(2,25.46587);
   tmpMC__139->SetBinError(3,44.2094);
   tmpMC__139->SetBinError(4,76.24735);
   tmpMC__139->SetBinError(5,138.27);
   tmpMC__139->SetBinError(6,197.6981);
   tmpMC__139->SetBinError(7,311.7293);
   tmpMC__139->SetBinError(8,365.3049);
   tmpMC__139->SetBinError(9,445.3938);
   tmpMC__139->SetBinError(10,458.9698);
   tmpMC__139->SetBinError(11,361.4284);
   tmpMC__139->SetBinError(12,283.8341);
   tmpMC__139->SetBinError(13,178.1318);
   tmpMC__139->SetBinError(14,152.6432);
   tmpMC__139->SetBinError(15,81.05067);
   tmpMC__139->SetBinError(16,47.61911);
   tmpMC__139->SetBinError(17,24.93476);
   tmpMC__139->SetBinError(18,19.12344);
   tmpMC__139->SetMinimum(0);
   tmpMC__139->SetMaximum(2819.349);
   tmpMC__139->SetEntries(398.3572);
   tmpMC__139->SetDirectory(0);
   tmpMC__139->SetStats(0);
   tmpMC__139->SetLineColor(2);
   tmpMC__139->SetLineWidth(3);
   tmpMC__139->SetMarkerStyle(0);
   tmpMC__139->GetXaxis()->CenterTitle(true);
   tmpMC__139->GetXaxis()->SetNdivisions(509);
   tmpMC__139->GetXaxis()->SetLabelFont(42);
   tmpMC__139->GetXaxis()->SetLabelSize(0.05);
   tmpMC__139->GetXaxis()->SetTitleSize(0.06);
   tmpMC__139->GetXaxis()->SetTitleOffset(1.15);
   tmpMC__139->GetYaxis()->SetLabelFont(42);
   tmpMC__139->GetYaxis()->SetLabelSize(0.05);
   tmpMC__139->GetYaxis()->SetTitleSize(0.06);
   tmpMC__139->GetYaxis()->SetTitleOffset(1.2);
   tmpMC__139->GetZaxis()->SetLabelFont(42);
   tmpMC__139->GetZaxis()->SetLabelSize(0.05);
   tmpMC__139->GetZaxis()->SetTitleSize(0.06);
   tmpMC__139->GetZaxis()->SetTitleOffset(0.75);
   tmpMC__139->Draw("SAME HIST");
   Double_t xAxis148[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *tmp_data__140 = new TH1D("tmp_data__140","dthetaP:ptmu",18, xAxis148);
   tmp_data__140->SetBinContent(1,95.16374);
   tmp_data__140->SetBinContent(2,84.36831);
   tmp_data__140->SetBinContent(3,141.0335);
   tmp_data__140->SetBinContent(4,282.6689);
   tmp_data__140->SetBinContent(5,664.0867);
   tmp_data__140->SetBinContent(6,640.2204);
   tmp_data__140->SetBinContent(7,21.78605);
   tmp_data__140->SetBinContent(8,520.1758);
   tmp_data__140->SetBinContent(9,670.6489);
   tmp_data__140->SetBinContent(10,1879.566);
   tmp_data__140->SetBinContent(11,435.653);
   tmp_data__140->SetBinContent(12,679.9636);
   tmp_data__140->SetBinContent(13,629.3979);
   tmp_data__140->SetBinContent(14,476.081);
   tmp_data__140->SetBinContent(15,318.4664);
   tmp_data__140->SetBinContent(16,245.1307);
   tmp_data__140->SetBinContent(17,84.76965);
   tmp_data__140->SetBinContent(18,59.61659);
   tmp_data__140->SetBinError(1,22.2797);
   tmp_data__140->SetBinError(2,31.66483);
   tmp_data__140->SetBinError(3,69.77777);
   tmp_data__140->SetBinError(4,129.596);
   tmp_data__140->SetBinError(5,237.1786);
   tmp_data__140->SetBinError(6,317.3246);
   tmp_data__140->SetBinError(7,486.2364);
   tmp_data__140->SetBinError(8,570.9452);
   tmp_data__140->SetBinError(9,685.8927);
   tmp_data__140->SetBinError(10,700.156);
   tmp_data__140->SetBinError(11,582.0047);
   tmp_data__140->SetBinError(12,496.1075);
   tmp_data__140->SetBinError(13,322.0441);
   tmp_data__140->SetBinError(14,240.278);
   tmp_data__140->SetBinError(15,135.4775);
   tmp_data__140->SetBinError(16,66.22195);
   tmp_data__140->SetBinError(17,30.89073);
   tmp_data__140->SetBinError(18,21.78352);
   tmp_data__140->SetEntries(112.8412);
   tmp_data__140->SetDirectory(0);
   tmp_data__140->SetStats(0);
   tmp_data__140->SetMarkerStyle(20);
   tmp_data__140->GetXaxis()->CenterTitle(true);
   tmp_data__140->GetXaxis()->SetNdivisions(509);
   tmp_data__140->GetXaxis()->SetLabelFont(42);
   tmp_data__140->GetXaxis()->SetLabelSize(0.05);
   tmp_data__140->GetXaxis()->SetTitleSize(0.06);
   tmp_data__140->GetXaxis()->SetTitleOffset(1.15);
   tmp_data__140->GetYaxis()->SetLabelFont(42);
   tmp_data__140->GetYaxis()->SetLabelSize(0.05);
   tmp_data__140->GetYaxis()->SetTitleSize(0.06);
   tmp_data__140->GetYaxis()->SetTitleOffset(1.2);
   tmp_data__140->GetZaxis()->SetLabelFont(42);
   tmp_data__140->GetZaxis()->SetLabelSize(0.05);
   tmp_data__140->GetZaxis()->SetTitleSize(0.06);
   tmp_data__140->GetZaxis()->SetTitleOffset(0.75);
   tmp_data__140->Draw("SAME E1 X0");
   
   TLegend *leg = new TLegend(0.565,0.73,0.825,0.89,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("tmp_data","Data","lep");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("tmpMC","Simulation","fl");

   ci = TColor::GetColor("#ffcccc");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineColor(2);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   leg->Draw();
   Double_t xAxis149[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *tmpMC__141 = new TH1D("tmpMC__141","dthetaP:ptmu",18, xAxis149);
   tmpMC__141->SetBinContent(1,66.48904);
   tmpMC__141->SetBinContent(2,82.19329);
   tmpMC__141->SetBinContent(3,130.6603);
   tmpMC__141->SetBinContent(4,211.2014);
   tmpMC__141->SetBinContent(5,339.661);
   tmpMC__141->SetBinContent(6,528.3071);
   tmpMC__141->SetBinContent(7,846.8248);
   tmpMC__141->SetBinContent(8,1143.128);
   tmpMC__141->SetBinContent(9,1645.034);
   tmpMC__141->SetBinContent(10,1688.734);
   tmpMC__141->SetBinContent(11,1168.084);
   tmpMC__141->SetBinContent(12,684.7981);
   tmpMC__141->SetBinContent(13,467.1951);
   tmpMC__141->SetBinContent(14,372.0023);
   tmpMC__141->SetBinContent(15,225.3816);
   tmpMC__141->SetBinContent(16,138.628);
   tmpMC__141->SetBinContent(17,79.92888);
   tmpMC__141->SetBinContent(18,64.80535);
   tmpMC__141->SetBinError(1,20.09539);
   tmpMC__141->SetBinError(2,25.46587);
   tmpMC__141->SetBinError(3,44.2094);
   tmpMC__141->SetBinError(4,76.24735);
   tmpMC__141->SetBinError(5,138.27);
   tmpMC__141->SetBinError(6,197.6981);
   tmpMC__141->SetBinError(7,311.7293);
   tmpMC__141->SetBinError(8,365.3049);
   tmpMC__141->SetBinError(9,445.3938);
   tmpMC__141->SetBinError(10,458.9698);
   tmpMC__141->SetBinError(11,361.4284);
   tmpMC__141->SetBinError(12,283.8341);
   tmpMC__141->SetBinError(13,178.1318);
   tmpMC__141->SetBinError(14,152.6432);
   tmpMC__141->SetBinError(15,81.05067);
   tmpMC__141->SetBinError(16,47.61911);
   tmpMC__141->SetBinError(17,24.93476);
   tmpMC__141->SetBinError(18,19.12344);
   tmpMC__141->SetMinimum(0);
   tmpMC__141->SetMaximum(2819.349);
   tmpMC__141->SetEntries(398.3572);
   tmpMC__141->SetDirectory(0);
   tmpMC__141->SetStats(0);

   ci = TColor::GetColor("#ffcccc");
   tmpMC__141->SetFillColor(ci);

   ci = TColor::GetColor("#000099");
   tmpMC__141->SetLineColor(ci);
   tmpMC__141->SetLineWidth(2);
   tmpMC__141->SetMarkerStyle(0);
   tmpMC__141->GetXaxis()->CenterTitle(true);
   tmpMC__141->GetXaxis()->SetNdivisions(509);
   tmpMC__141->GetXaxis()->SetLabelFont(42);
   tmpMC__141->GetXaxis()->SetLabelSize(0.05);
   tmpMC__141->GetXaxis()->SetTitleSize(0.06);
   tmpMC__141->GetXaxis()->SetTitleOffset(1.15);
   tmpMC__141->GetYaxis()->SetLabelFont(42);
   tmpMC__141->GetYaxis()->SetLabelSize(0.05);
   tmpMC__141->GetYaxis()->SetTitleSize(0.06);
   tmpMC__141->GetYaxis()->SetTitleOffset(1.2);
   tmpMC__141->GetZaxis()->SetLabelFont(42);
   tmpMC__141->GetZaxis()->SetLabelSize(0.05);
   tmpMC__141->GetZaxis()->SetTitleSize(0.06);
   tmpMC__141->GetZaxis()->SetTitleOffset(0.75);
   tmpMC__141->Draw("sameaxis");
   Double_t xAxis150[19] = {-180, -125, -80, -55, -35, -25, -15, -10, -5, 0, 5, 10, 15, 25, 35, 55, 80, 125, 180}; 
   
   TH1D *mnv_tmp_data_statonly__142 = new TH1D("mnv_tmp_data_statonly__142","dthetaP:ptmu",18, xAxis150);
   mnv_tmp_data_statonly__142->SetBinContent(1,95.16374);
   mnv_tmp_data_statonly__142->SetBinContent(2,84.36831);
   mnv_tmp_data_statonly__142->SetBinContent(3,141.0335);
   mnv_tmp_data_statonly__142->SetBinContent(4,282.6689);
   mnv_tmp_data_statonly__142->SetBinContent(5,664.0867);
   mnv_tmp_data_statonly__142->SetBinContent(6,640.2204);
   mnv_tmp_data_statonly__142->SetBinContent(7,21.78605);
   mnv_tmp_data_statonly__142->SetBinContent(8,520.1758);
   mnv_tmp_data_statonly__142->SetBinContent(9,670.6489);
   mnv_tmp_data_statonly__142->SetBinContent(10,1879.566);
   mnv_tmp_data_statonly__142->SetBinContent(11,435.653);
   mnv_tmp_data_statonly__142->SetBinContent(12,679.9636);
   mnv_tmp_data_statonly__142->SetBinContent(13,629.3979);
   mnv_tmp_data_statonly__142->SetBinContent(14,476.081);
   mnv_tmp_data_statonly__142->SetBinContent(15,318.4664);
   mnv_tmp_data_statonly__142->SetBinContent(16,245.1307);
   mnv_tmp_data_statonly__142->SetBinContent(17,84.76965);
   mnv_tmp_data_statonly__142->SetBinContent(18,59.61659);
   mnv_tmp_data_statonly__142->SetBinError(1,18.81258);
   mnv_tmp_data_statonly__142->SetBinError(2,25.18773);
   mnv_tmp_data_statonly__142->SetBinError(3,54.15073);
   mnv_tmp_data_statonly__142->SetBinError(4,95.92023);
   mnv_tmp_data_statonly__142->SetBinError(5,189.0873);
   mnv_tmp_data_statonly__142->SetBinError(6,238.9351);
   mnv_tmp_data_statonly__142->SetBinError(7,386.5452);
   mnv_tmp_data_statonly__142->SetBinError(8,424.1903);
   mnv_tmp_data_statonly__142->SetBinError(9,453.1438);
   mnv_tmp_data_statonly__142->SetBinError(10,458.5533);
   mnv_tmp_data_statonly__142->SetBinError(11,420.1796);
   mnv_tmp_data_statonly__142->SetBinError(12,388.1987);
   mnv_tmp_data_statonly__142->SetBinError(13,240.0707);
   mnv_tmp_data_statonly__142->SetBinError(14,185.0476);
   mnv_tmp_data_statonly__142->SetBinError(15,95.72835);
   mnv_tmp_data_statonly__142->SetBinError(16,55.06851);
   mnv_tmp_data_statonly__142->SetBinError(17,25.06326);
   mnv_tmp_data_statonly__142->SetBinError(18,18.00181);
   mnv_tmp_data_statonly__142->SetEntries(112.8412);
   mnv_tmp_data_statonly__142->SetDirectory(0);
   mnv_tmp_data_statonly__142->SetStats(0);
   mnv_tmp_data_statonly__142->SetMarkerStyle(20);
   mnv_tmp_data_statonly__142->GetXaxis()->CenterTitle(true);
   mnv_tmp_data_statonly__142->GetXaxis()->SetNdivisions(509);
   mnv_tmp_data_statonly__142->GetXaxis()->SetLabelFont(42);
   mnv_tmp_data_statonly__142->GetXaxis()->SetLabelSize(0.05);
   mnv_tmp_data_statonly__142->GetXaxis()->SetTitleSize(0.06);
   mnv_tmp_data_statonly__142->GetXaxis()->SetTitleOffset(1.15);
   mnv_tmp_data_statonly__142->GetYaxis()->SetLabelFont(42);
   mnv_tmp_data_statonly__142->GetYaxis()->SetLabelSize(0.05);
   mnv_tmp_data_statonly__142->GetYaxis()->SetTitleSize(0.06);
   mnv_tmp_data_statonly__142->GetYaxis()->SetTitleOffset(1.2);
   mnv_tmp_data_statonly__142->GetZaxis()->SetLabelFont(42);
   mnv_tmp_data_statonly__142->GetZaxis()->SetLabelSize(0.05);
   mnv_tmp_data_statonly__142->GetZaxis()->SetTitleSize(0.06);
   mnv_tmp_data_statonly__142->GetZaxis()->SetTitleOffset(0.75);
   mnv_tmp_data_statonly__142->Draw("SAME E1 X0");
   TLatex *   tex = new TLatex(0.18,0.8925,"MINER#nuA Preliminary");
tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextColor(4);
   tex->SetTextFont(112);
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.25,0.85,"POT-Normalized");
tex->SetNDC();
   tex->SetTextAlign(22);
   tex->SetTextColor(9);
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.25,0.82,"Data POT: 1.44E+20");
tex->SetNDC();
   tex->SetTextAlign(22);
   tex->SetTextColor(9);
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.25,0.79,"MC POT: 6.08E+20");
tex->SetNDC();
   tex->SetTextAlign(22);
   tex->SetTextColor(9);
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
   cNeutA->Modified();
   cNeutA->cd();
   cNeutA->SetSelected(cNeutA);
}
