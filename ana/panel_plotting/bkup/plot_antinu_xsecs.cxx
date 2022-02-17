// #include "../mec_common.h"
// #include "../syst_common.h"
// #include "../plot_common.h"

//  #include "../util/plot/myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
// #include "../util/plot/plot.h"

#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "TEnv.h"
#include "TChain.h"
#include "TF2.h"
#include "Math/DistFunc.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TFile.h"

#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

using namespace PlotUtils;

//======================================================================
TString uniq()
{
  static int i=0;
  return TString::Format("uniq%d", i++);
}


//======================================================================
TGraphErrors* histToGraph(TH1* h, bool includeZeros=true)
{
  TGraphErrors* grE=new TGraphErrors;
  grE->SetLineColor(h->GetLineColor());
  grE->SetLineStyle(h->GetLineStyle());
  grE->SetLineWidth(h->GetLineWidth());

  for(int i=0; i<h->GetNbinsX(); ++i){
    int bin=i+1;
    double binVal=h->GetBinContent(bin);
    if(includeZeros || binVal!=0){
      grE->SetPoint(i, h->GetBinCenter(bin), binVal);
      grE->SetPointError(i, 0, h->GetBinError(bin));
    }
  }
  return grE;
}

//======================================================================
int main()
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  TFile f("/pnfs/minerva/persistent/users/schellma/May2017/bigrun_more_v25_mec1_phil1_rpa1_2017-05-13_1117_qelikelo/cross_sections/eroica/cross_sections_muonpz_muonpt_lowangleqelike_minerva.root");
  MnvH2D* dataMnv=(MnvH2D*)f.Get("cross_sections_muonpt_muonpz_data");
  MnvH2D* mcMnv=(MnvH2D*)f.Get("cross_sections_muonpt_muonpz_mc");

  dataMnv->Scale(1e39);
  mcMnv->Scale(1e39);

  TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* data=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());

  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(0.7);
  data->SetLineColor(kBlack);

  dataStat->SetMarkerStyle(1);
  dataStat->SetLineColor(kBlack);

  GridCanvas* gc=new GridCanvas(uniq(), 4, 3, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  for(int i=0; i<11; ++i){
    TPad* pad=(TPad*)gc->cd(i+1);
    // pad->SetGridy();
    TH1* mcproj=mc->ProjectionY(uniq(), i+1, i+1);
    mcproj->SetLineColor(kRed);
    mcproj->SetMarkerColor(kRed);
    mcproj->Draw("hist");

    TGraphErrors* grData=histToGraph(data->ProjectionY(uniq(), i+1, i+1));
    TGraphErrors* grDataStat=histToGraph(dataStat->ProjectionY(uniq(), i+1, i+1));
    grDataStat->Draw("e");
    grData->Draw("ep");
    double binmin=mc->GetXaxis()->GetBinLowEdge(i+1);
    double binmax=binmin+mc->GetXaxis()->GetBinWidth(i+1);
    TLatex* la=new TLatex(1-pad->GetRightMargin()-0.01,
                          1-pad->GetTopMargin()-0.02,
                          TString::Format("%.1f < #it{p_{z}} < %.1f GeV", binmin, binmax));
    la->SetTextAlign(33); // top right
    la->SetNDC();
    la->SetTextFont(42);
    la->SetTextSize(0.035);
    la->Draw();
    mcproj->GetXaxis()->SetNdivisions(4);
    mcproj->GetYaxis()->SetNdivisions(4);
  }
  gc->SetYLimits(1e-3, 0.7);

  gc->SetXTitle("Muon transverse momentum (GeV)");
  gc->SetYTitle("Cross section");
  gc->ResetPads();
  gc->Draw();

  gc->Print("antinu-2d-xsec-pt.eps");

  GridCanvas* gc2=new GridCanvas(uniq(), 3, 2, 800, 500);
  gc2->SetRightMargin(0.01);
  gc2->SetLeftMargin(0.1);
  gc2->ResetPads();

  for(int i=0; i<6; ++i){
    TPad* pad=(TPad*)gc2->cd(i+1);
    TH1* mcproj=mc->ProjectionX(uniq(), i+1, i+1);
    mcproj->SetLineColor(kRed);
    mcproj->SetMarkerColor(kRed);
    mcproj->Draw("hist");

    TGraphErrors* grData=histToGraph(data->ProjectionX(uniq(), i+1, i+1));
    TGraphErrors* grDataStat=histToGraph(dataStat->ProjectionX(uniq(), i+1, i+1));
    grDataStat->Draw("e");
    grData->Draw("ep");
    double binmin=mc->GetYaxis()->GetBinLowEdge(i+1);
    double binmax=binmin+mc->GetYaxis()->GetBinWidth(i+1);
    TLatex* la=new TLatex(1-pad->GetRightMargin()-0.01,
                          1-pad->GetTopMargin()-0.02,
                          TString::Format("%.1f < #it{p_{T}} < %.1f GeV", binmin, binmax));
    la->SetTextAlign(33); // top right
    la->SetNDC();
    la->SetTextFont(42);
    la->SetTextSize(0.035);
    la->Draw();
    mcproj->GetXaxis()->SetNdivisions(4);
    mcproj->GetYaxis()->SetNdivisions(4);
  }
  gc2->SetYLimits(1e-3, 0.7);

  gc2->SetXTitle("Muon longitudinal momentum (GeV)");
  gc2->SetYTitle("Cross section");
  gc2->ResetPads();
  gc2->Draw();

  gc2->Print("antinu-2d-xsec-pz.eps");

}
