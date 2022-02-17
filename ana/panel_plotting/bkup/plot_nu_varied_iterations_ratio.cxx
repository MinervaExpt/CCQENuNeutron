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

  TFile f("/minerva/data/users/drut1186/ana_hists/test_iterations_CV/CrossSection_iterations_1_SpecialSampleIncluded.root");
  TFile f2("/minerva/data/users/drut1186/ana_hists/test_iterations_CV/CrossSection_iterations_2_SpecialSampleIncluded.root");
  TFile f3("/minerva/data/users/drut1186/ana_hists/test_iterations_CV/CrossSection_iterations_3_SpecialSampleIncluded.root");
  TFile f4("/minerva/data/users/drut1186/ana_hists/test_iterations_CV/CrossSection_iterations_5_SpecialSampleIncluded.root");
  TFile f5("/minerva/data/users/drut1186/ana_hists/test_iterations_CV/CrossSection_iterations_10_SpecialSampleIncluded.root");
  TFile f6("/minerva/data/users/drut1186/ana_hists/test_iterations_CV/CrossSection_iterations_20_SpecialSampleIncluded.root");

  MnvH2D* dataMnv=(MnvH2D*)f.Get("h_pzmu_ptmu_data_nobck_unfold_effcor_cross_section");
  MnvH2D* dataMnv2=(MnvH2D*)f2.Get("h_pzmu_ptmu_data_nobck_unfold_effcor_cross_section");
  MnvH2D* dataMnv3=(MnvH2D*)f3.Get("h_pzmu_ptmu_data_nobck_unfold_effcor_cross_section");
  MnvH2D* dataMnv4=(MnvH2D*)f4.Get("h_pzmu_ptmu_data_nobck_unfold_effcor_cross_section");
  MnvH2D* dataMnv5=(MnvH2D*)f5.Get("h_pzmu_ptmu_data_nobck_unfold_effcor_cross_section");
  MnvH2D* dataMnv6=(MnvH2D*)f6.Get("h_pzmu_ptmu_data_nobck_unfold_effcor_cross_section");


  dataMnv->GetXaxis()->SetTitle("p_{||} (GeV)");

  dataMnv->Scale(1e39, "width");
  dataMnv2->Scale(1e39, "width");
  dataMnv3->Scale(1e39, "width");
  dataMnv4->Scale(1e39, "width");
  dataMnv5->Scale(1e39, "width");
  dataMnv6->Scale(1e39, "width");


  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* data=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* datasys=new TH2D(dataMnv6->GetCVHistoWithError());
  TH2* data2=new TH2D(dataMnv2->GetCVHistoWithError());
  TH2* data3=new TH2D(dataMnv3->GetCVHistoWithError());
  TH2* data4=new TH2D(dataMnv4->GetCVHistoWithError());
  TH2* data5=new TH2D(dataMnv5->GetCVHistoWithError());
  TH2* data6=new TH2D(dataMnv6->GetCVHistoWithError());
  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = getColors(1);

  // These line and marker styles will be propagated to the 1D plots

  data->SetLineColor(mycolors[4]);
  data2->SetLineColor(mycolors[5]);
  data3->SetLineColor(mycolors[6]);
  data4->SetLineColor(mycolors[7]);
  data5->SetLineColor(mycolors[8]);
  data6->SetLineColor(mycolors[9]);

  data->SetLineWidth(2);
  data2->SetLineWidth(2);
  data3->SetLineWidth(2);
  data4->SetLineWidth(2);
  data5->SetLineWidth(2);
  data6->SetLineWidth(2);

  data->Divide(data6);
  data3->Divide(data6);
  data4->Divide(data6);
  data5->Divide(data6);

  datasys->Divide(data6);
  data2->Divide(data6);
  data6->Divide(data6);


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
  histAndOpts.push_back(std::make_pair(data,     "hist"));
  histAndOpts.push_back(std::make_pair(data2,     "hist"));
  histAndOpts.push_back(std::make_pair(data3,     "hist"));
  histAndOpts.push_back(std::make_pair(data4,     "hist"));
  histAndOpts.push_back(std::make_pair(data5,     "hist"));
  histAndOpts.push_back(std::make_pair(data6,     "hist"));
  histAndOpts.push_back(std::make_pair(datasys,     "graph e"));
  //  histAndOpts.push_back(std::make_pair(dataStat, "hist"));

  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range
  double multipliers[]={1, 1, 1, 1,
                        2, 2, 2, 2,
                        10, 20, 30, 50};

  GridCanvas* gc=plotpT1D(histAndOpts, doMultipliers ? multipliers : NULL);
  // Set the y range manually. Can also use gc->Remax() to guess automatically
  gc->SetYLimits(0, 2);
  gc->SetYTitle("ratio to nominal iterations");
  gc->Modified();

  gc->Print(doMultipliers ? "nu-varied-iteration-ratio-pt-multiplier.eps" : "nu-varied-iteration-ratio-pt.eps");
  gc->Print(doMultipliers ? "nu-varied-iteration-ratio-pt-multiplier.png" : "nu-varied-iteration-ratio-pt.png");
  gc->Print(doMultipliers ? "nu-varied-iteration-ratio-pt-multiplier.C" : "nu-varied-iteration-ratio-pt.C");

  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TCanvas *c1 = new TCanvas();
  TLegend* leg=new TLegend(0, 0, 1, 1);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(data, "MINERvA data 1 iteration", "lpe");
  leg->AddEntry(data2, "MINERvA data 2 iteration", "lpe");
  leg->AddEntry(data3, "MINERvA data 3 iteration", "lpe");
  leg->AddEntry(data4, "MINERvA data 5 iteration", "lpe");
  leg->AddEntry(data5, "MINERvA data 10 iteration", "lpe");
  leg->AddEntry(data6, "MINERvA data 20 iteration", "lpe");
  leg->Draw("");
  c1->Print("nu-varied-iteration-legend.png");
  c1->Print("nu-varied-iteration-legend.C");
  c1->Print("nu-varied-iteration-legend.eps");

  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  // Values to multiply each bin by to get them on a similar range
  double  multipliers2[]={15, 5, 2, 1, 1,
                         1, 1, 1, 1,2 };

  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  GridCanvas* gc2=plotpz1D(histAndOpts, doMultipliers ? multipliers2 : NULL);
  gc2->SetYLimits(0, 2);
  gc2->SetYTitle("ratio to nominal iterations");
  gc2->Modified();
  gc2->Print(doMultipliers ? "nu-varied-iteration-ratio-pz-multiplier.eps" : "nu-varied-iteration-ratio-pz.eps");
  gc2->Print(doMultipliers ? "nu-varied-iteration-ratio-pz-multiplier.png" : "nu-varied-iteration-ratio-pz.png");
  gc2->Print(doMultipliers ? "nu-varied-iteration-ratio-pz-multiplier.C" : "nu-varied-iteration-ratio-pz.C");

}

int main()
{
  //  makePlots(true);
  makePlots(false);

  return 0;
}
