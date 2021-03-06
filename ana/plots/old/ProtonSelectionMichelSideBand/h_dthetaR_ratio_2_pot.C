{
//=========Macro generated from canvas: cNeutA/Neutron Angulars
//=========  (Sun Nov 17 14:42:28 2019) by ROOT version5.34/36
   TCanvas *cNeutA = new TCanvas("cNeutA", "Neutron Angulars",1,1,900,726);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cNeutA->SetHighLightColor(2);
   cNeutA->Range(-257.0125,-2.092105,257.0125,1.855263);
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
   Double_t xAxis86[41] = {-179.9087, -170.9133, -161.9179, -152.9224, -143.927, -134.9316, -125.9361, -116.9407, -107.9452, -98.94981, -89.95437, -80.95894, -71.9635, -62.96806, -53.97262, -44.97719, -35.98175, -26.98631, -17.99087, -8.995437, 0, 8.995437, 17.99087, 26.98631, 35.98175, 44.97719, 53.97262, 62.96806, 71.9635, 80.95894, 89.95437, 98.94981, 107.9452, 116.9407, 125.9361, 134.9316, 143.927, 152.9224, 161.9179, 170.9133, 179.9087}; 
   
   PlotUtils::MnvH1D *tmp_ratio__80 = new PlotUtils::MnvH1D("tmp_ratio__80","dthetaR:ptmu",40, xAxis86);
   tmp_ratio__80->SetBinContent(1,1.270594);
   tmp_ratio__80->SetBinContent(2,-13.61178);
   tmp_ratio__80->SetBinContent(3,-2.180014);
   tmp_ratio__80->SetBinContent(4,-1.673273);
   tmp_ratio__80->SetBinContent(5,-1.209761);
   tmp_ratio__80->SetBinContent(6,-0.8270913);
   tmp_ratio__80->SetBinContent(7,-1.464188);
   tmp_ratio__80->SetBinContent(8,2.97348);
   tmp_ratio__80->SetBinContent(9,1.7438);
   tmp_ratio__80->SetBinContent(10,2.163257);
   tmp_ratio__80->SetBinContent(11,1.183012);
   tmp_ratio__80->SetBinContent(12,1.991301);
   tmp_ratio__80->SetBinContent(13,0.5931825);
   tmp_ratio__80->SetBinContent(14,0.8085469);
   tmp_ratio__80->SetBinContent(15,0.4157553);
   tmp_ratio__80->SetBinContent(16,-0.4494238);
   tmp_ratio__80->SetBinContent(17,0.7455156);
   tmp_ratio__80->SetBinContent(18,-0.8802684);
   tmp_ratio__80->SetBinContent(19,0.2487657);
   tmp_ratio__80->SetBinContent(20,1.128743);
   tmp_ratio__80->SetBinContent(21,1.206454);
   tmp_ratio__80->SetBinContent(22,5.050504);
   tmp_ratio__80->SetBinContent(23,-1.07176);
   tmp_ratio__80->SetBinContent(24,-0.9707679);
   tmp_ratio__80->SetBinContent(25,9.762419);
   tmp_ratio__80->SetBinContent(26,-7.732307);
   tmp_ratio__80->SetBinContent(27,3.517555);
   tmp_ratio__80->SetBinContent(28,-6.038777);
   tmp_ratio__80->SetBinContent(29,-15.48655);
   tmp_ratio__80->SetBinContent(30,6.983151);
   tmp_ratio__80->SetBinContent(31,1.780579);
   tmp_ratio__80->SetBinContent(32,1.591815);
   tmp_ratio__80->SetBinContent(33,-2.730195);
   tmp_ratio__80->SetBinContent(34,2.543315);
   tmp_ratio__80->SetBinContent(35,0.7629915);
   tmp_ratio__80->SetBinContent(36,5.449637);
   tmp_ratio__80->SetBinContent(37,-2.155728);
   tmp_ratio__80->SetBinContent(38,0.3347788);
   tmp_ratio__80->SetBinContent(39,-0.7125434);
   tmp_ratio__80->SetBinContent(40,-9.67911);
   tmp_ratio__80->SetBinError(1,20.12773);
   tmp_ratio__80->SetBinError(2,99.34624);
   tmp_ratio__80->SetBinError(3,14.49021);
   tmp_ratio__80->SetBinError(4,3.342392);
   tmp_ratio__80->SetBinError(5,1.530232);
   tmp_ratio__80->SetBinError(6,6.141803);
   tmp_ratio__80->SetBinError(7,1.864015);
   tmp_ratio__80->SetBinError(8,4.518184);
   tmp_ratio__80->SetBinError(9,2.046946);
   tmp_ratio__80->SetBinError(10,2.184006);
   tmp_ratio__80->SetBinError(11,1.434161);
   tmp_ratio__80->SetBinError(12,1.793937);
   tmp_ratio__80->SetBinError(13,0.8289395);
   tmp_ratio__80->SetBinError(14,0.9868224);
   tmp_ratio__80->SetBinError(15,0.8316798);
   tmp_ratio__80->SetBinError(16,0.9510906);
   tmp_ratio__80->SetBinError(17,1.033282);
   tmp_ratio__80->SetBinError(18,0.8588423);
   tmp_ratio__80->SetBinError(19,0.5573504);
   tmp_ratio__80->SetBinError(20,0.6082448);
   tmp_ratio__80->SetBinError(21,0.7320499);
   tmp_ratio__80->SetBinError(22,7.235918);
   tmp_ratio__80->SetBinError(23,4.37845);
   tmp_ratio__80->SetBinError(24,5.029808);
   tmp_ratio__80->SetBinError(25,87.13841);
   tmp_ratio__80->SetBinError(26,59.20781);
   tmp_ratio__80->SetBinError(27,83.33004);
   tmp_ratio__80->SetBinError(28,36.32996);
   tmp_ratio__80->SetBinError(29,192.455);
   tmp_ratio__80->SetBinError(30,43.14385);
   tmp_ratio__80->SetBinError(31,7.716359);
   tmp_ratio__80->SetBinError(32,5.50027);
   tmp_ratio__80->SetBinError(33,10.01152);
   tmp_ratio__80->SetBinError(34,8.587015);
   tmp_ratio__80->SetBinError(35,4.186146);
   tmp_ratio__80->SetBinError(36,66.28867);
   tmp_ratio__80->SetBinError(37,6.759026);
   tmp_ratio__80->SetBinError(38,1.956145);
   tmp_ratio__80->SetBinError(39,5.431251);
   tmp_ratio__80->SetBinError(40,52.30006);
   tmp_ratio__80->SetMinimum(-1.5);
   tmp_ratio__80->SetMaximum(1.5);
   tmp_ratio__80->SetEntries(14.62439);
   tmp_ratio__80->SetDirectory(0);
   tmp_ratio__80->SetStats(0);
   tmp_ratio__80->SetLineWidth(3);
   tmp_ratio__80->SetMarkerStyle(20);
   tmp_ratio__80->GetXaxis()->CenterTitle(true);
   tmp_ratio__80->GetXaxis()->SetNdivisions(509);
   tmp_ratio__80->GetXaxis()->SetLabelFont(42);
   tmp_ratio__80->GetXaxis()->SetLabelSize(0.05);
   tmp_ratio__80->GetXaxis()->SetTitleSize(0.06);
   tmp_ratio__80->GetXaxis()->SetTitleOffset(1.15);
   tmp_ratio__80->GetYaxis()->SetTitle("Data / MC");
   tmp_ratio__80->GetYaxis()->SetLabelFont(42);
   tmp_ratio__80->GetYaxis()->SetLabelSize(0.05);
   tmp_ratio__80->GetYaxis()->SetTitleSize(0.06);
   tmp_ratio__80->GetYaxis()->SetTitleOffset(1.2);
   tmp_ratio__80->GetZaxis()->SetLabelFont(42);
   tmp_ratio__80->GetZaxis()->SetLabelSize(0.05);
   tmp_ratio__80->GetZaxis()->SetTitleSize(0.06);
   tmp_ratio__80->GetZaxis()->SetTitleOffset(0.75);
   tmp_ratio__80->Draw("X0");
   Double_t xAxis87[41] = {-179.9087, -170.9133, -161.9179, -152.9224, -143.927, -134.9316, -125.9361, -116.9407, -107.9452, -98.94981, -89.95437, -80.95894, -71.9635, -62.96806, -53.97262, -44.97719, -35.98175, -26.98631, -17.99087, -8.995437, 0, 8.995437, 17.99087, 26.98631, 35.98175, 44.97719, 53.97262, 62.96806, 71.9635, 80.95894, 89.95437, 98.94981, 107.9452, 116.9407, 125.9361, 134.9316, 143.927, 152.9224, 161.9179, 170.9133, 179.9087}; 
   
   TH1D *h_dthetaR_ptmu_mc_projX_TotalError__81 = new TH1D("h_dthetaR_ptmu_mc_projX_TotalError__81","dthetaR:ptmu",40, xAxis87);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(1,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(2,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(3,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(4,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(5,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(6,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(7,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(8,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(9,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(10,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(11,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(12,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(13,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(14,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(15,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(16,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(17,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(18,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(19,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(20,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(21,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(22,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(23,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(24,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(25,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(26,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(27,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(28,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(29,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(30,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(31,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(32,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(33,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(34,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(35,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(36,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(37,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(38,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(39,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinContent(40,1);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(1,7.253662);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(2,2.984055);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(3,1.845761);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(4,0.6316987);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(5,0.368895);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(6,1.611082);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(7,0.7708192);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(8,1.103525);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(9,1.10062);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(10,1.219601);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(11,1.385047);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(12,1.366537);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(13,1.278368);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(14,1.726924);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(15,1.900666);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(16,2.579046);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(17,2.614081);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(18,2.235248);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(19,1.47452);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(20,0.9492769);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(21,0.6363984);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(22,1.140876);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(23,1.267822);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(24,1.326585);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(25,1.816667);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(26,1.868929);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(27,5.871048);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(28,2.810855);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(29,5.047772);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(30,2.838866);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(31,2.189517);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(32,1.809451);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(33,2.098567);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(34,1.831281);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(35,1.445549);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(36,6.809079);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(37,1.888157);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(38,1.212539);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(39,2.35888);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetBinError(40,3.399901);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetMaximum(14.64363);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetEntries(82);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetDirectory(0);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ffcccc");
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetFillColor(ci);

   ci = TColor::GetColor("#ffcccc");
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetLineColor(ci);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetLineWidth(3);
   h_dthetaR_ptmu_mc_projX_TotalError__81->SetMarkerStyle(0);
   h_dthetaR_ptmu_mc_projX_TotalError__81->GetXaxis()->SetRange(1,40);
   h_dthetaR_ptmu_mc_projX_TotalError__81->GetXaxis()->CenterTitle(true);
   h_dthetaR_ptmu_mc_projX_TotalError__81->GetXaxis()->SetNdivisions(509);
   h_dthetaR_ptmu_mc_projX_TotalError__81->GetXaxis()->SetLabelFont(42);
   h_dthetaR_ptmu_mc_projX_TotalError__81->GetXaxis()->SetLabelSize(0.05);
   h_dthetaR_ptmu_mc_projX_TotalError__81->GetXaxis()->SetTitleSize(0.06);
   h_dthetaR_ptmu_mc_projX_TotalError__81->GetXaxis()->SetTitleOffset(1.15);
   h_dthetaR_ptmu_mc_projX_TotalError__81->GetYaxis()->SetLabelFont(42);
   h_dthetaR_ptmu_mc_projX_TotalError__81->GetYaxis()->SetLabelSize(0.05);
   h_dthetaR_ptmu_mc_projX_TotalError__81->GetYaxis()->SetTitleSize(0.06);
   h_dthetaR_ptmu_mc_projX_TotalError__81->GetYaxis()->SetTitleOffset(1.2);
   h_dthetaR_ptmu_mc_projX_TotalError__81->GetZaxis()->SetLabelFont(42);
   h_dthetaR_ptmu_mc_projX_TotalError__81->GetZaxis()->SetLabelSize(0.05);
   h_dthetaR_ptmu_mc_projX_TotalError__81->GetZaxis()->SetTitleSize(0.06);
   h_dthetaR_ptmu_mc_projX_TotalError__81->GetZaxis()->SetTitleOffset(0.75);
   h_dthetaR_ptmu_mc_projX_TotalError__81->Draw("E2 same ][");
   Double_t xAxis88[41] = {-179.9087, -170.9133, -161.9179, -152.9224, -143.927, -134.9316, -125.9361, -116.9407, -107.9452, -98.94981, -89.95437, -80.95894, -71.9635, -62.96806, -53.97262, -44.97719, -35.98175, -26.98631, -17.99087, -8.995437, 0, 8.995437, 17.99087, 26.98631, 35.98175, 44.97719, 53.97262, 62.96806, 71.9635, 80.95894, 89.95437, 98.94981, 107.9452, 116.9407, 125.9361, 134.9316, 143.927, 152.9224, 161.9179, 170.9133, 179.9087}; 
   
   PlotUtils::MnvH1D *tmp_ratio__82 = new PlotUtils::MnvH1D("tmp_ratio__82","dthetaR:ptmu",40, xAxis88);
   tmp_ratio__82->SetBinContent(1,1.270594);
   tmp_ratio__82->SetBinContent(2,-13.61178);
   tmp_ratio__82->SetBinContent(3,-2.180014);
   tmp_ratio__82->SetBinContent(4,-1.673273);
   tmp_ratio__82->SetBinContent(5,-1.209761);
   tmp_ratio__82->SetBinContent(6,-0.8270913);
   tmp_ratio__82->SetBinContent(7,-1.464188);
   tmp_ratio__82->SetBinContent(8,2.97348);
   tmp_ratio__82->SetBinContent(9,1.7438);
   tmp_ratio__82->SetBinContent(10,2.163257);
   tmp_ratio__82->SetBinContent(11,1.183012);
   tmp_ratio__82->SetBinContent(12,1.991301);
   tmp_ratio__82->SetBinContent(13,0.5931825);
   tmp_ratio__82->SetBinContent(14,0.8085469);
   tmp_ratio__82->SetBinContent(15,0.4157553);
   tmp_ratio__82->SetBinContent(16,-0.4494238);
   tmp_ratio__82->SetBinContent(17,0.7455156);
   tmp_ratio__82->SetBinContent(18,-0.8802684);
   tmp_ratio__82->SetBinContent(19,0.2487657);
   tmp_ratio__82->SetBinContent(20,1.128743);
   tmp_ratio__82->SetBinContent(21,1.206454);
   tmp_ratio__82->SetBinContent(22,5.050504);
   tmp_ratio__82->SetBinContent(23,-1.07176);
   tmp_ratio__82->SetBinContent(24,-0.9707679);
   tmp_ratio__82->SetBinContent(25,9.762419);
   tmp_ratio__82->SetBinContent(26,-7.732307);
   tmp_ratio__82->SetBinContent(27,3.517555);
   tmp_ratio__82->SetBinContent(28,-6.038777);
   tmp_ratio__82->SetBinContent(29,-15.48655);
   tmp_ratio__82->SetBinContent(30,6.983151);
   tmp_ratio__82->SetBinContent(31,1.780579);
   tmp_ratio__82->SetBinContent(32,1.591815);
   tmp_ratio__82->SetBinContent(33,-2.730195);
   tmp_ratio__82->SetBinContent(34,2.543315);
   tmp_ratio__82->SetBinContent(35,0.7629915);
   tmp_ratio__82->SetBinContent(36,5.449637);
   tmp_ratio__82->SetBinContent(37,-2.155728);
   tmp_ratio__82->SetBinContent(38,0.3347788);
   tmp_ratio__82->SetBinContent(39,-0.7125434);
   tmp_ratio__82->SetBinContent(40,-9.67911);
   tmp_ratio__82->SetBinError(1,20.12773);
   tmp_ratio__82->SetBinError(2,99.34624);
   tmp_ratio__82->SetBinError(3,14.49021);
   tmp_ratio__82->SetBinError(4,3.342392);
   tmp_ratio__82->SetBinError(5,1.530232);
   tmp_ratio__82->SetBinError(6,6.141803);
   tmp_ratio__82->SetBinError(7,1.864015);
   tmp_ratio__82->SetBinError(8,4.518184);
   tmp_ratio__82->SetBinError(9,2.046946);
   tmp_ratio__82->SetBinError(10,2.184006);
   tmp_ratio__82->SetBinError(11,1.434161);
   tmp_ratio__82->SetBinError(12,1.793937);
   tmp_ratio__82->SetBinError(13,0.8289395);
   tmp_ratio__82->SetBinError(14,0.9868224);
   tmp_ratio__82->SetBinError(15,0.8316798);
   tmp_ratio__82->SetBinError(16,0.9510906);
   tmp_ratio__82->SetBinError(17,1.033282);
   tmp_ratio__82->SetBinError(18,0.8588423);
   tmp_ratio__82->SetBinError(19,0.5573504);
   tmp_ratio__82->SetBinError(20,0.6082448);
   tmp_ratio__82->SetBinError(21,0.7320499);
   tmp_ratio__82->SetBinError(22,7.235918);
   tmp_ratio__82->SetBinError(23,4.37845);
   tmp_ratio__82->SetBinError(24,5.029808);
   tmp_ratio__82->SetBinError(25,87.13841);
   tmp_ratio__82->SetBinError(26,59.20781);
   tmp_ratio__82->SetBinError(27,83.33004);
   tmp_ratio__82->SetBinError(28,36.32996);
   tmp_ratio__82->SetBinError(29,192.455);
   tmp_ratio__82->SetBinError(30,43.14385);
   tmp_ratio__82->SetBinError(31,7.716359);
   tmp_ratio__82->SetBinError(32,5.50027);
   tmp_ratio__82->SetBinError(33,10.01152);
   tmp_ratio__82->SetBinError(34,8.587015);
   tmp_ratio__82->SetBinError(35,4.186146);
   tmp_ratio__82->SetBinError(36,66.28867);
   tmp_ratio__82->SetBinError(37,6.759026);
   tmp_ratio__82->SetBinError(38,1.956145);
   tmp_ratio__82->SetBinError(39,5.431251);
   tmp_ratio__82->SetBinError(40,52.30006);
   tmp_ratio__82->SetMinimum(-1.5);
   tmp_ratio__82->SetMaximum(1.5);
   tmp_ratio__82->SetEntries(14.62439);
   tmp_ratio__82->SetDirectory(0);
   tmp_ratio__82->SetStats(0);
   tmp_ratio__82->SetLineWidth(3);
   tmp_ratio__82->SetMarkerStyle(20);
   tmp_ratio__82->GetXaxis()->CenterTitle(true);
   tmp_ratio__82->GetXaxis()->SetNdivisions(509);
   tmp_ratio__82->GetXaxis()->SetLabelFont(42);
   tmp_ratio__82->GetXaxis()->SetLabelSize(0.05);
   tmp_ratio__82->GetXaxis()->SetTitleSize(0.06);
   tmp_ratio__82->GetXaxis()->SetTitleOffset(1.15);
   tmp_ratio__82->GetYaxis()->SetTitle("Data / MC");
   tmp_ratio__82->GetYaxis()->SetLabelFont(42);
   tmp_ratio__82->GetYaxis()->SetLabelSize(0.05);
   tmp_ratio__82->GetYaxis()->SetTitleSize(0.06);
   tmp_ratio__82->GetYaxis()->SetTitleOffset(1.2);
   tmp_ratio__82->GetZaxis()->SetLabelFont(42);
   tmp_ratio__82->GetZaxis()->SetLabelSize(0.05);
   tmp_ratio__82->GetZaxis()->SetTitleSize(0.06);
   tmp_ratio__82->GetZaxis()->SetTitleOffset(0.75);
   tmp_ratio__82->Draw("SAME AXIS");
   Double_t xAxis89[41] = {-179.9087, -170.9133, -161.9179, -152.9224, -143.927, -134.9316, -125.9361, -116.9407, -107.9452, -98.94981, -89.95437, -80.95894, -71.9635, -62.96806, -53.97262, -44.97719, -35.98175, -26.98631, -17.99087, -8.995437, 0, 8.995437, 17.99087, 26.98631, 35.98175, 44.97719, 53.97262, 62.96806, 71.9635, 80.95894, 89.95437, 98.94981, 107.9452, 116.9407, 125.9361, 134.9316, 143.927, 152.9224, 161.9179, 170.9133, 179.9087}; 
   
   PlotUtils::MnvH1D *tmp_ratio__83 = new PlotUtils::MnvH1D("tmp_ratio__83","dthetaR:ptmu",40, xAxis89);
   tmp_ratio__83->SetBinContent(1,1.270594);
   tmp_ratio__83->SetBinContent(2,-13.61178);
   tmp_ratio__83->SetBinContent(3,-2.180014);
   tmp_ratio__83->SetBinContent(4,-1.673273);
   tmp_ratio__83->SetBinContent(5,-1.209761);
   tmp_ratio__83->SetBinContent(6,-0.8270913);
   tmp_ratio__83->SetBinContent(7,-1.464188);
   tmp_ratio__83->SetBinContent(8,2.97348);
   tmp_ratio__83->SetBinContent(9,1.7438);
   tmp_ratio__83->SetBinContent(10,2.163257);
   tmp_ratio__83->SetBinContent(11,1.183012);
   tmp_ratio__83->SetBinContent(12,1.991301);
   tmp_ratio__83->SetBinContent(13,0.5931825);
   tmp_ratio__83->SetBinContent(14,0.8085469);
   tmp_ratio__83->SetBinContent(15,0.4157553);
   tmp_ratio__83->SetBinContent(16,-0.4494238);
   tmp_ratio__83->SetBinContent(17,0.7455156);
   tmp_ratio__83->SetBinContent(18,-0.8802684);
   tmp_ratio__83->SetBinContent(19,0.2487657);
   tmp_ratio__83->SetBinContent(20,1.128743);
   tmp_ratio__83->SetBinContent(21,1.206454);
   tmp_ratio__83->SetBinContent(22,5.050504);
   tmp_ratio__83->SetBinContent(23,-1.07176);
   tmp_ratio__83->SetBinContent(24,-0.9707679);
   tmp_ratio__83->SetBinContent(25,9.762419);
   tmp_ratio__83->SetBinContent(26,-7.732307);
   tmp_ratio__83->SetBinContent(27,3.517555);
   tmp_ratio__83->SetBinContent(28,-6.038777);
   tmp_ratio__83->SetBinContent(29,-15.48655);
   tmp_ratio__83->SetBinContent(30,6.983151);
   tmp_ratio__83->SetBinContent(31,1.780579);
   tmp_ratio__83->SetBinContent(32,1.591815);
   tmp_ratio__83->SetBinContent(33,-2.730195);
   tmp_ratio__83->SetBinContent(34,2.543315);
   tmp_ratio__83->SetBinContent(35,0.7629915);
   tmp_ratio__83->SetBinContent(36,5.449637);
   tmp_ratio__83->SetBinContent(37,-2.155728);
   tmp_ratio__83->SetBinContent(38,0.3347788);
   tmp_ratio__83->SetBinContent(39,-0.7125434);
   tmp_ratio__83->SetBinContent(40,-9.67911);
   tmp_ratio__83->SetBinError(1,20.12773);
   tmp_ratio__83->SetBinError(2,99.34624);
   tmp_ratio__83->SetBinError(3,14.49021);
   tmp_ratio__83->SetBinError(4,3.342392);
   tmp_ratio__83->SetBinError(5,1.530232);
   tmp_ratio__83->SetBinError(6,6.141803);
   tmp_ratio__83->SetBinError(7,1.864015);
   tmp_ratio__83->SetBinError(8,4.518184);
   tmp_ratio__83->SetBinError(9,2.046946);
   tmp_ratio__83->SetBinError(10,2.184006);
   tmp_ratio__83->SetBinError(11,1.434161);
   tmp_ratio__83->SetBinError(12,1.793937);
   tmp_ratio__83->SetBinError(13,0.8289395);
   tmp_ratio__83->SetBinError(14,0.9868224);
   tmp_ratio__83->SetBinError(15,0.8316798);
   tmp_ratio__83->SetBinError(16,0.9510906);
   tmp_ratio__83->SetBinError(17,1.033282);
   tmp_ratio__83->SetBinError(18,0.8588423);
   tmp_ratio__83->SetBinError(19,0.5573504);
   tmp_ratio__83->SetBinError(20,0.6082448);
   tmp_ratio__83->SetBinError(21,0.7320499);
   tmp_ratio__83->SetBinError(22,7.235918);
   tmp_ratio__83->SetBinError(23,4.37845);
   tmp_ratio__83->SetBinError(24,5.029808);
   tmp_ratio__83->SetBinError(25,87.13841);
   tmp_ratio__83->SetBinError(26,59.20781);
   tmp_ratio__83->SetBinError(27,83.33004);
   tmp_ratio__83->SetBinError(28,36.32996);
   tmp_ratio__83->SetBinError(29,192.455);
   tmp_ratio__83->SetBinError(30,43.14385);
   tmp_ratio__83->SetBinError(31,7.716359);
   tmp_ratio__83->SetBinError(32,5.50027);
   tmp_ratio__83->SetBinError(33,10.01152);
   tmp_ratio__83->SetBinError(34,8.587015);
   tmp_ratio__83->SetBinError(35,4.186146);
   tmp_ratio__83->SetBinError(36,66.28867);
   tmp_ratio__83->SetBinError(37,6.759026);
   tmp_ratio__83->SetBinError(38,1.956145);
   tmp_ratio__83->SetBinError(39,5.431251);
   tmp_ratio__83->SetBinError(40,52.30006);
   tmp_ratio__83->SetMinimum(-1.5);
   tmp_ratio__83->SetMaximum(1.5);
   tmp_ratio__83->SetEntries(14.62439);
   tmp_ratio__83->SetDirectory(0);
   tmp_ratio__83->SetStats(0);
   tmp_ratio__83->SetLineWidth(3);
   tmp_ratio__83->SetMarkerStyle(20);
   tmp_ratio__83->GetXaxis()->CenterTitle(true);
   tmp_ratio__83->GetXaxis()->SetNdivisions(509);
   tmp_ratio__83->GetXaxis()->SetLabelFont(42);
   tmp_ratio__83->GetXaxis()->SetLabelSize(0.05);
   tmp_ratio__83->GetXaxis()->SetTitleSize(0.06);
   tmp_ratio__83->GetXaxis()->SetTitleOffset(1.15);
   tmp_ratio__83->GetYaxis()->SetTitle("Data / MC");
   tmp_ratio__83->GetYaxis()->SetLabelFont(42);
   tmp_ratio__83->GetYaxis()->SetLabelSize(0.05);
   tmp_ratio__83->GetYaxis()->SetTitleSize(0.06);
   tmp_ratio__83->GetYaxis()->SetTitleOffset(1.2);
   tmp_ratio__83->GetZaxis()->SetLabelFont(42);
   tmp_ratio__83->GetZaxis()->SetLabelSize(0.05);
   tmp_ratio__83->GetZaxis()->SetTitleSize(0.06);
   tmp_ratio__83->GetZaxis()->SetTitleOffset(0.75);
   tmp_ratio__83->Draw("SAME X0");
   TLine *line = new TLine(-179.9087,1,179.9087,1);
   line->SetLineColor(36);
   line->SetLineStyle(2);
   line->SetLineWidth(3);
   line->Draw();
   Double_t xAxis90[41] = {-179.9087, -170.9133, -161.9179, -152.9224, -143.927, -134.9316, -125.9361, -116.9407, -107.9452, -98.94981, -89.95437, -80.95894, -71.9635, -62.96806, -53.97262, -44.97719, -35.98175, -26.98631, -17.99087, -8.995437, 0, 8.995437, 17.99087, 26.98631, 35.98175, 44.97719, 53.97262, 62.96806, 71.9635, 80.95894, 89.95437, 98.94981, 107.9452, 116.9407, 125.9361, 134.9316, 143.927, 152.9224, 161.9179, 170.9133, 179.9087}; 
   
   PlotUtils::MnvH1D *tmp_ratio__84 = new PlotUtils::MnvH1D("tmp_ratio__84","dthetaR:ptmu",40, xAxis90);
   tmp_ratio__84->SetBinContent(1,1.270594);
   tmp_ratio__84->SetBinContent(2,-13.61178);
   tmp_ratio__84->SetBinContent(3,-2.180014);
   tmp_ratio__84->SetBinContent(4,-1.673273);
   tmp_ratio__84->SetBinContent(5,-1.209761);
   tmp_ratio__84->SetBinContent(6,-0.8270913);
   tmp_ratio__84->SetBinContent(7,-1.464188);
   tmp_ratio__84->SetBinContent(8,2.97348);
   tmp_ratio__84->SetBinContent(9,1.7438);
   tmp_ratio__84->SetBinContent(10,2.163257);
   tmp_ratio__84->SetBinContent(11,1.183012);
   tmp_ratio__84->SetBinContent(12,1.991301);
   tmp_ratio__84->SetBinContent(13,0.5931825);
   tmp_ratio__84->SetBinContent(14,0.8085469);
   tmp_ratio__84->SetBinContent(15,0.4157553);
   tmp_ratio__84->SetBinContent(16,-0.4494238);
   tmp_ratio__84->SetBinContent(17,0.7455156);
   tmp_ratio__84->SetBinContent(18,-0.8802684);
   tmp_ratio__84->SetBinContent(19,0.2487657);
   tmp_ratio__84->SetBinContent(20,1.128743);
   tmp_ratio__84->SetBinContent(21,1.206454);
   tmp_ratio__84->SetBinContent(22,5.050504);
   tmp_ratio__84->SetBinContent(23,-1.07176);
   tmp_ratio__84->SetBinContent(24,-0.9707679);
   tmp_ratio__84->SetBinContent(25,9.762419);
   tmp_ratio__84->SetBinContent(26,-7.732307);
   tmp_ratio__84->SetBinContent(27,3.517555);
   tmp_ratio__84->SetBinContent(28,-6.038777);
   tmp_ratio__84->SetBinContent(29,-15.48655);
   tmp_ratio__84->SetBinContent(30,6.983151);
   tmp_ratio__84->SetBinContent(31,1.780579);
   tmp_ratio__84->SetBinContent(32,1.591815);
   tmp_ratio__84->SetBinContent(33,-2.730195);
   tmp_ratio__84->SetBinContent(34,2.543315);
   tmp_ratio__84->SetBinContent(35,0.7629915);
   tmp_ratio__84->SetBinContent(36,5.449637);
   tmp_ratio__84->SetBinContent(37,-2.155728);
   tmp_ratio__84->SetBinContent(38,0.3347788);
   tmp_ratio__84->SetBinContent(39,-0.7125434);
   tmp_ratio__84->SetBinContent(40,-9.67911);
   tmp_ratio__84->SetBinError(1,20.12773);
   tmp_ratio__84->SetBinError(2,99.34624);
   tmp_ratio__84->SetBinError(3,14.49021);
   tmp_ratio__84->SetBinError(4,3.342392);
   tmp_ratio__84->SetBinError(5,1.530232);
   tmp_ratio__84->SetBinError(6,6.141803);
   tmp_ratio__84->SetBinError(7,1.864015);
   tmp_ratio__84->SetBinError(8,4.518184);
   tmp_ratio__84->SetBinError(9,2.046946);
   tmp_ratio__84->SetBinError(10,2.184006);
   tmp_ratio__84->SetBinError(11,1.434161);
   tmp_ratio__84->SetBinError(12,1.793937);
   tmp_ratio__84->SetBinError(13,0.8289395);
   tmp_ratio__84->SetBinError(14,0.9868224);
   tmp_ratio__84->SetBinError(15,0.8316798);
   tmp_ratio__84->SetBinError(16,0.9510906);
   tmp_ratio__84->SetBinError(17,1.033282);
   tmp_ratio__84->SetBinError(18,0.8588423);
   tmp_ratio__84->SetBinError(19,0.5573504);
   tmp_ratio__84->SetBinError(20,0.6082448);
   tmp_ratio__84->SetBinError(21,0.7320499);
   tmp_ratio__84->SetBinError(22,7.235918);
   tmp_ratio__84->SetBinError(23,4.37845);
   tmp_ratio__84->SetBinError(24,5.029808);
   tmp_ratio__84->SetBinError(25,87.13841);
   tmp_ratio__84->SetBinError(26,59.20781);
   tmp_ratio__84->SetBinError(27,83.33004);
   tmp_ratio__84->SetBinError(28,36.32996);
   tmp_ratio__84->SetBinError(29,192.455);
   tmp_ratio__84->SetBinError(30,43.14385);
   tmp_ratio__84->SetBinError(31,7.716359);
   tmp_ratio__84->SetBinError(32,5.50027);
   tmp_ratio__84->SetBinError(33,10.01152);
   tmp_ratio__84->SetBinError(34,8.587015);
   tmp_ratio__84->SetBinError(35,4.186146);
   tmp_ratio__84->SetBinError(36,66.28867);
   tmp_ratio__84->SetBinError(37,6.759026);
   tmp_ratio__84->SetBinError(38,1.956145);
   tmp_ratio__84->SetBinError(39,5.431251);
   tmp_ratio__84->SetBinError(40,52.30006);
   tmp_ratio__84->SetMinimum(-1.5);
   tmp_ratio__84->SetMaximum(1.5);
   tmp_ratio__84->SetEntries(14.62439);
   tmp_ratio__84->SetDirectory(0);
   tmp_ratio__84->SetStats(0);
   tmp_ratio__84->SetLineWidth(3);
   tmp_ratio__84->SetMarkerStyle(20);
   tmp_ratio__84->GetXaxis()->CenterTitle(true);
   tmp_ratio__84->GetXaxis()->SetNdivisions(509);
   tmp_ratio__84->GetXaxis()->SetLabelFont(42);
   tmp_ratio__84->GetXaxis()->SetLabelSize(0.05);
   tmp_ratio__84->GetXaxis()->SetTitleSize(0.06);
   tmp_ratio__84->GetXaxis()->SetTitleOffset(1.15);
   tmp_ratio__84->GetYaxis()->SetTitle("Data / MC");
   tmp_ratio__84->GetYaxis()->SetLabelFont(42);
   tmp_ratio__84->GetYaxis()->SetLabelSize(0.05);
   tmp_ratio__84->GetYaxis()->SetTitleSize(0.06);
   tmp_ratio__84->GetYaxis()->SetTitleOffset(1.2);
   tmp_ratio__84->GetZaxis()->SetLabelFont(42);
   tmp_ratio__84->GetZaxis()->SetLabelSize(0.05);
   tmp_ratio__84->GetZaxis()->SetTitleSize(0.06);
   tmp_ratio__84->GetZaxis()->SetTitleOffset(0.75);
   tmp_ratio__84->Draw("sameaxis");
   cNeutA->Modified();
   cNeutA->cd();
   cNeutA->SetSelected(cNeutA);
}
