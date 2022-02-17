//#include "myPlotStyle.h"
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
#include "localColor.h"

#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"

using namespace PlotUtils;

void makePlots(bool doMultipliers)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  TFile f("/pnfs/minerva/persistent/users/schellma/May2017/bigrun_more_v25_mec1_phil1_rpa1_2017-05-13_1117_qelikelo/cross_sections/eroica/cross_sections_muonpz_muonpt_lowangleqelike_minerva.root");
  TFile f2("/pnfs/minerva/persistent/users/schellma/May2017/AntiNuModels_Paper_qelike_muonvars_anglecut_paper.root");
  MnvH2D* dataMnv=(MnvH2D*)f.Get("cross_sections_muonpt_muonpz_data");
  MnvH2D* mcMnv=(MnvH2D*)f.Get("cross_sections_muonpt_muonpz_mc");
  TH2D* nomGenieMnv = (MnvH2D*)f2.Get("Untuned_GENIE");
  TH2D* bestGenieMnv = (TH2D*)f2.Get("MnvGENIE_RPA_MEC_notune");
  dataMnv->GetXaxis()->SetTitle("p_{||} (GeV)");

  dataMnv->Scale(1e39, "width");
  mcMnv->Scale(1e39, "width");
  nomGenieMnv->Scale(1e39, "width");
  bestGenieMnv->Scale(1e39, "width");
  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* data=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mc_tmp=new TH2D(mcMnv->GetCVHistoWithStatError());
  TH2* nomGenie_tmp=nomGenieMnv;
  TH2* bestGenie_tmp=bestGenieMnv;
  
  //clone the 2D so I can change 1.5 to 1.49999... stupid hack;
  const int nBinsX = mc_tmp->GetNbinsX();
  const int nBinsY = mc_tmp->GetNbinsY();
  double xbinEdges[nBinsX+1];
  double ybinEdges[nBinsY+1];
  mc_tmp->GetXaxis()->GetLowEdge(xbinEdges);
  xbinEdges[nBinsX]=15;
  mc_tmp->GetYaxis()->GetLowEdge(ybinEdges);
  ybinEdges[nBinsY]=1.5;
  ybinEdges[0]=0.0001;
  TH2* mc = new TH2D("mc","mc",nBinsX,xbinEdges,nBinsY,ybinEdges);
  TH2* nomGenie = new TH2D("nomgine","nomgenie",nBinsX,xbinEdges,nBinsY,ybinEdges);
  TH2* bestGenie = new TH2D("bestgine","bestgenie",nBinsX,xbinEdges,nBinsY,ybinEdges);
  for(int i=0;i<nBinsX;i++){
    for(int j=0;j<nBinsY;j++){
      mc->SetBinContent(i+1,j+1,mc_tmp->GetBinContent(i+1,j+1));
      mc->SetBinError(i+1,j+1,mc_tmp->GetBinError(i+1,j+1));
      nomGenie->SetBinContent(i+1,j+1,nomGenie_tmp->GetBinContent(i+1,j+1));
      nomGenie->SetBinError(i+1,j+1,nomGenie_tmp->GetBinError(i+1,j+1));
      bestGenie->SetBinContent(i+1,j+1,bestGenie_tmp->GetBinContent(i+1,j+1));
      bestGenie->SetBinError(i+1,j+1,bestGenie_tmp->GetBinError(i+1,j+1));
      cout << nomGenie->GetBinContent(i+1,j+1) << "\t" << bestGenie->GetBinContent(i+1,j+1) << endl;
    }
  }
  
  vector<int> mycolors = getColors(1);
  // These line and marker styles will be propagated to the 1D plots
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);
  
  nomGenie->SetLineColor(kBlue);
  nomGenie->SetLineWidth(2);

  // These line and marker styles will be propagated to the 1D plots
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(0.7);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2);

  dataStat->SetMarkerStyle(1);
  dataStat->SetLineColor(kBlack);
  dataStat->SetLineWidth(2);

  bestGenie->SetLineColor(mycolors[10]);
  bestGenie->SetLineWidth(2);

  // Make a list of the histograms we want to draw, along with the
  // draw options we want to use for them. You can add "graph" to the
  // draw options if you want the histogram to be converted to a graph
  // and then drawn. In that case the draw options are interpreted as
  // options to TGraphErrors::Draw().
  //
  // I don't know what happens if you put a "graph" first in the list,
  // so don't do that. Make sure the first item doesn't have "graph"
  // in its options
  std::vector<std::pair<TH2*, const char*> > histAndOpts;
  histAndOpts.push_back(std::make_pair(mc,       "hist"));
  histAndOpts.push_back(std::make_pair(nomGenie, "hist"));
  histAndOpts.push_back(std::make_pair(bestGenie, "hist"));
  histAndOpts.push_back(std::make_pair(dataStat, "graph e"));
  histAndOpts.push_back(std::make_pair(data,     "graph ep"));



  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range
  double multipliers[]={1, 1, 1, 1,
                        2, 4, 4, 4,
                        10, 20, 30, 1};

  GridCanvas* gc=plotpT1DAntiNu(histAndOpts, doMultipliers ? multipliers : NULL);
  // Set the y range manually. Can also use gc->Remax() to guess automatically
  gc->SetYLimits(0, 3.9);

  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TLegend* leg=new TLegend(0.8, 0.17, 1, 0.4);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(data, "MINERvA data", "lpe");
  leg->AddEntry(mc, "MINERvA Tune", "l");
  leg->AddEntry(nomGenie, "GENIE 2.8.4", "l");
  leg->AddEntry(bestGenie, "2p2h and RPA", "l");
  leg->Draw();

  gc->Print(doMultipliers ? "antinu-2d-xsec-pt-multiplier.eps" : "antinu-2d-xsec-pt.eps");
  gc->Print(doMultipliers ? "antinu-2d-xsec-pt-multiplier.png" : "antinu-2d-xsec-pt.png");
  gc->Print(doMultipliers ? "antinu-2d-xsec-pt-multiplier.C" : "antinu-2d-xsec-pt.C");

  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  // Values to multiply each bin by to get them on a similar range
  double  multipliers2[]={5, 2, 1,
                         1, 2, 50 };

  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  GridCanvas* gc2=plotpz1DAntiNu(histAndOpts, doMultipliers ? multipliers2 : NULL);
  gc2->SetYLimits(0, 3.9);
  gc2->Print(doMultipliers ? "antinu-2d-xsec-pz-multiplier.eps" : "antinu-2d-xsec-pz.eps");
  gc2->Print(doMultipliers ? "antinu-2d-xsec-pz-multiplier.png" : "antinu-2d-xsec-pz.png");
  gc2->Print(doMultipliers ? "antinu-2d-xsec-pz-multiplier.C" : "antinu-2d-xsec-pz.C");

}

int main()
{
  makePlots(true);
  makePlots(false);

  return 0;
}
