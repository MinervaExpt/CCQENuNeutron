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
#include "localColor.h"

#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"

using namespace PlotUtils;

void makePlots(bool doMultipliers, string location)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  TFile f(Form("%s_CV/CrossSection_per_nucleon_iterations_4_SpecialSampleIncluded.root_pzmu.root",location.c_str()));
  TFile f2("/minerva/app/users/drut1186/cmtuser/Minerva_v10r8p9_CCQENu2DLE_hadronv2/Ana/CCQENu2DLE/ana/plot_macros_pub/CrossSections/GENIEModels.root");
  MnvH2D* dataMnv=(MnvH2D*)f.Get("h_pzmu_ptmu_data_nobck_unfold_effcor_cross_section");
  MnvH2D* mcMnv=(MnvH2D*)f.Get("h_pzmu_ptmu_data_nobck_unfold_effcor_cross_section");

  dataMnv->Scale(1e39, "width");
  mcMnv->Scale(1e39, "width");

  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* data=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());
  TH2* mod1 = (TH2D*)f2.Get("GENIE");
  TH2* mod2 = (TH2D*)f2.Get("PionTuned");
  TH2* mod3 = (TH2D*)f2.Get("GENIE+RPA");
  TH2* mod4 = (TH2D*)f2.Get("GENIE+2p2h");
  TH2* mod5 = (TH2D*)f2.Get("GENIE+RPA+2p2h");
  TH2* mod6 = (TH2D*)f2.Get("GENIE+RPA+2p2h+Recoil+NievesRPARes");
  TH2* mod7 = (TH2D*)f2.Get("GENIE+RPA+2p2h+Recoil+MINOSRPARes");

  mod1->Scale(1e39, "width");
  mod2->Scale(1e39, "width");
  mod3->Scale(1e39, "width");
  mod4->Scale(1e39, "width");
  mod5->Scale(1e39, "width");
  mod6->Scale(1e39, "width");
  mod7->Scale(1e39, "width");

  vector<int> mycolors = getColors(1);

  // These line and marker styles will be propagated to the 1D plots
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);
  
  mod1->SetLineColor(kBlue);
  mod1->SetLineWidth(2);

  mod2->SetLineColor(mycolors[7]);
  mod2->SetLineWidth(2);
  mod3->SetLineColor(mycolors[8]);
  mod3->SetLineWidth(2);
  mod4->SetLineColor(mycolors[9]);
  mod4->SetLineWidth(2);
  mod5->SetLineColor(mycolors[10]);
  mod5->SetLineWidth(2);
  mod6->SetLineColor(mycolors[11]);
  mod6->SetLineWidth(2);
  mod7->SetLineColor(mycolors[12]);
  mod7->SetLineWidth(2);


  data->Divide(mc);
  dataStat->Divide(mc);
  mod2->Divide(mc);
  mod3->Divide(mc);
  mod4->Divide(mc);
  mod5->Divide(mc);
  mod6->Divide(mc);
  mod7->Divide(mc);
  mod1->Divide(mc);
  mc->Divide(mc);


  // These line and marker styles will be propagated to the 1D plots
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(0.7);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2);

  dataStat->SetMarkerStyle(1);
  dataStat->SetLineColor(kBlack);
  dataStat->SetLineWidth(2);

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
  histAndOpts.push_back(std::make_pair(mc,       "hist"));//0
  histAndOpts.push_back(std::make_pair(mod1,     "hist"));
  histAndOpts.push_back(std::make_pair(mod2,     "hist"));
  histAndOpts.push_back(std::make_pair(mod3,     "hist"));
  histAndOpts.push_back(std::make_pair(mod4,     "hist"));
  histAndOpts.push_back(std::make_pair(mod5,     "hist"));//5
  histAndOpts.push_back(std::make_pair(mod6,     "hist"));
  histAndOpts.push_back(std::make_pair(mod7,     "hist"));
  histAndOpts.push_back(std::make_pair(dataStat, "hist"));
  histAndOpts.push_back(std::make_pair(data,     "hist"));//9

  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range
  double multipliers[]={1, 1, 1, 1,
                        2, 2, 2, 2,
                        10, 20, 30, 50};

  
  GridCanvas* gc=plotpT1D(histAndOpts, doMultipliers ? multipliers : NULL);
  // Set the y range manually. Can also use gc->Remax() to guess automatically
  //  gc->SetYLimits(0, 3.9);
  gc->SetYLimits(0.8,1.2);
  gc->SetYTitle("Ratio to Default GENIE");
  gc->Modified();
  gc->Print(doMultipliers ? "nu-2d-xsec-models-dataxsec-pt-multiplier.eps" : "nu-2d-xsec-models-dataxsec-pt.eps");
  gc->Print(doMultipliers ? "nu-2d-xsec-models-dataxsec-pt-multiplier.png" : "nu-2d-xsec-models-dataxsec-pt.png");
  gc->Print(doMultipliers ? "nu-2d-xsec-models-dataxsec-pt-multiplier.C" : "nu-2d-xsec-models-dataxsec-pt.C");

  for(int i=0;i<histAndOpts.size()-4;i++){
    std::vector<std::pair<TH2*, const char*> > histAndOpts_subset;
    histAndOpts_subset.push_back(histAndOpts[0]);
    histAndOpts_subset.push_back(histAndOpts[1]);
    histAndOpts_subset.push_back(histAndOpts[i+2]);
    histAndOpts_subset.push_back(histAndOpts[8]);
    histAndOpts_subset.push_back(histAndOpts[9]);
    GridCanvas* gc_sub=plotpT1D(histAndOpts_subset, doMultipliers ? multipliers : NULL);
    // Set the y range manually. Can also use gc->Remax() to guess automatically
    //  gc->SetYLimits(0, 3.9);
    gc_sub->SetYLimits(0.8,1.2);
    gc_sub->SetYTitle("Ratio to Default GENIE");
    gc_sub->Modified();
    gc_sub->Print(doMultipliers ? Form("nu-2d-xsec-models-dataxsec-pt-multiplier-group-%d.eps",i) : Form("nu-2d-xsec-models-dataxsec-pt-group-%d.eps",i));
    gc_sub->Print(doMultipliers ? Form("nu-2d-xsec-models-dataxsec-pt-multiplier-group-%d.PNG",i) : Form("nu-2d-xsec-models-dataxsec-pt-group-%d.png",i));
    gc_sub->Print(doMultipliers ? Form("nu-2d-xsec-models-dataxsec-pt-multiplier-group-%d.C",i) : Form("nu-2d-xsec-models-dataxsec-pt-group-%d.C",i));


  }
  

  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  // Values to multiply each bin by to get them on a similar range
  double  multipliers2[]={5, 1, 1, 1, 1,
                         1, 3, 15, 100, 500};

  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  GridCanvas* gc2=plotpz1D(histAndOpts, doMultipliers ? multipliers2 : NULL);
  //  gc2->SetYLimits(0, 3.9);
  gc2->SetYLimits(0.8,1.2);
  gc2->SetYTitle("Ratio to Default GENIE");
  gc2->Modified();
  gc2->Print(doMultipliers ? "nu-2d-xsec-models-pz-multiplier.eps" : "nu-2d-xsec-models-pz.eps");
  gc2->Print(doMultipliers ? "nu-2d-xsec-models-pz-multiplier.png" : "nu-2d-xsec-models-pz.png");
  gc2->Print(doMultipliers ? "nu-2d-xsec-models-pz-multiplier.C" : "nu-2d-xsec-models-pz.C");

  for(int i=0;i<histAndOpts.size()-4;i++){
    std::vector<std::pair<TH2*, const char*> > histAndOpts_subset;
    histAndOpts_subset.push_back(histAndOpts[0]);
    histAndOpts_subset.push_back(histAndOpts[1]);
    histAndOpts_subset.push_back(histAndOpts[i+2]);
    histAndOpts_subset.push_back(histAndOpts[8]);
    histAndOpts_subset.push_back(histAndOpts[9]);
    GridCanvas* gc2_sub=plotpz1D(histAndOpts_subset, doMultipliers ? multipliers : NULL);
    // Set the y range manually. Can also use gc->Remax() to guess automatically
    //  gc->SetYLimits(0, 3.9);
    gc2_sub->SetYLimits(0.8,1.2);
    gc2_sub->SetYTitle("Ratio to Default GENIE");
    gc2_sub->Modified();

    gc2_sub->Print(doMultipliers ? Form("nu-2d-xsec-models-dataxsec-pz-multiplier-group-%d.eps",i) : Form("nu-2d-xsec-models-dataxsec-pz-group-%d.eps",i));
    gc2_sub->Print(doMultipliers ? Form("nu-2d-xsec-models-dataxsec-pz-multiplier-group-%d.PNG",i) : Form("nu-2d-xsec-models-dataxsec-pz-group-%d.png",i));
    gc2_sub->Print(doMultipliers ? Form("nu-2d-xsec-models-dataxsec-pz-multiplier-group-%d.C",i) : Form("nu-2d-xsec-models-dataxsec-pz-group-%d.C",i));


  }



  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TLegend* leg=new TLegend(0,0,1,1);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.07);
  leg->AddEntry(mc, "MINERvA Tune v1", "l");
  leg->AddEntry(data, "MINERvA data", "lpe");
  leg->AddEntry(mod1, "GENIE 2.8.4", "l");
  leg->AddEntry(mod2, "Pion Tune","l");
  leg->AddEntry(mod3, "RPA", "l");
  leg->AddEntry(mod4, "2p2h", "l");
  leg->AddEntry(mod5, "RPA+2p2h", "l");
  leg->AddEntry(mod6, "MINERvA Tune v1 + Res. RPA (Nieves)", "l");
  leg->AddEntry(mod7, "MINERvA Tune v1 + Res. RPA (MINOS)", "l");

  TCanvas c1("Legend","Legend",10,10,800,500);
  leg->Draw();
  c1.Print("nu-2d-xsec-models-legend.eps");
  c1.Print("nu-2d-xsec-models-legend.png");
  c1.Print("nu-2d-xsec-models-legend.C");
}

int main(int argc, char* argv[])
{
  //  makePlots(true);
  makePlots(false,argv[1]);

  return 0;
}
