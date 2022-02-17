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


#include "plot.h"

using namespace PlotUtils;


//======================================================================
void makePlots(bool doMultiplier)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  TFile f("/minerva/data/users/drut1186/ana_hists/wc_ver2.0/CrossSection_SpecialSampleIncluded.root");
  MnvH2D* dataMnv=(MnvH2D*)f.Get("h_pzmu_ptmu_data_nobck_unfold_effcor_cross_section");
  MnvH2D* mcMnv=(MnvH2D*)f.Get("h_pzmu_ptmu_mc_nobck_unfold_effcor_cross_section");

  dataMnv->Scale(1e39, "width");
  mcMnv->Scale(1e39, "width");

  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* data=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());

  // These line and marker styles will be propagated to the 1D plots
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(0.7);
  data->SetLineColor(kBlack);

  dataStat->SetMarkerStyle(1);
  dataStat->SetLineColor(kBlack);

  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range
  double multipliers[]={1, 1, 1, 1,
                        2, 2, 2, 2,
                        10, 10, 10, 10};

  if(!doMultiplier){
    for(int i=0; i<12; ++i) multipliers[i]=1;
  }

  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), 4, 3, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<12; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    // Get the 1D projection for MC and set its colour
    TH1* mcproj=mc->ProjectionY(uniq(), i+1, i+1);
    mcproj->Scale(multipliers[i]);
    mcproj->SetLineColor(kRed);
    mcproj->SetMarkerColor(kRed);
    mcproj->Draw("hist");

    // We'll get the data 1D projections and convert them to
    // histograms. The advantage of the conversion is that we can make
    // ROOT not draw the error bar in x, which I don't like
    TGraphErrors* grData=histToGraph(data->ProjectionY(uniq(), i+1, i+1), multipliers[i]);
    TGraphErrors* grDataStat=histToGraph(dataStat->ProjectionY(uniq(), i+1, i+1), multipliers[i]);
    grDataStat->Draw("e");
    grData->Draw("ep");

    // Draw the label that says what pz bin this is.
    // First get the bin edge values
    double binmin=mc->GetXaxis()->GetBinLowEdge(i+1);
    double binmax=binmin+mc->GetXaxis()->GetBinWidth(i+1);

    // Now create the TLatex. The calculation for where to put it
    // requires a bit of explanation: each of the pads is actually the
    // same size as the whole GridCanvas. They get to be in different
    // positions because I set the pads' margins differently (which
    // puts the frames in different places). I want to put the label
    // near the top right of the frame, so I calculate relative to the
    // right margin and top margin
    TLatex* la=new TLatex(1-pad->GetRightMargin()-0.01,
                          1-pad->GetTopMargin()-0.02,
                          TString::Format("%.1f < #it{p_{z}} < %.1f GeV", binmin, binmax));
    la->SetTextAlign(33); // top right
    la->SetNDC();
    la->SetTextFont(42);
    la->SetTextSize(0.035);
    la->Draw();

    // Do the same for the multiplier text, if necessary
    if(doMultiplier && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.0f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();

    }

    // Fewer tick marks: less visual noise
    mcproj->GetXaxis()->SetNdivisions(4);
    mcproj->GetYaxis()->SetNdivisions(4);
  }
  // Manually set the y range, but we could also do gc->Remax() to attempt to set it automatically
  gc->SetYLimits(0, 3.9);

  gc->SetXTitle("Muon transverse momentum (GeV)");
  gc->SetYTitle("Cross section");
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  gc->Print(doMultiplier ? "nu-2d-xsec-pt-multiplier.eps" : "nu-2d-xsec-pt.eps");

  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  // Values to multiply each bin by to get them on a similar range
  double multipliers2[]={1, 1, 1, 1, 1,
                         1, 3, 15, 100, 500};

  if(!doMultiplier){
    for(int i=0; i<10; ++i) multipliers2[i]=1;
  }

  GridCanvas* gc2=new GridCanvas(uniq(), 5, 2, 800, 500);
  gc2->SetRightMargin(0.01);
  gc2->SetLeftMargin(0.1);
  gc2->ResetPads();

  for(int i=0; i<10; ++i){
    TPad* pad=(TPad*)gc2->cd(i+1);
    TH1* mcproj=mc->ProjectionX(uniq(), i+1, i+1);
    mcproj->Scale(multipliers2[i]);
    mcproj->SetLineColor(kRed);
    mcproj->SetMarkerColor(kRed);
    mcproj->Draw("hist");

    TGraphErrors* grData=histToGraph(data->ProjectionX(uniq(), i+1, i+1), multipliers2[i]);
    TGraphErrors* grDataStat=histToGraph(dataStat->ProjectionX(uniq(), i+1, i+1), multipliers2[i]);
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

    // Do the same for the multiplier text, if necessary
    if(doMultiplier && multipliers2[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.0f", multipliers2[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();

    }

    mcproj->GetXaxis()->SetNdivisions(4);
    mcproj->GetYaxis()->SetNdivisions(4);
  }
  gc2->SetYLimits(0, 3.9);
  // gc2->Remax();

  gc2->SetXTitle("Muon longitudinal momentum (GeV)");
  gc2->SetYTitle("Cross section");
  gc2->ResetPads();
  gc2->Draw();

  gc2->Print(doMultiplier ? "nu-2d-xsec-pz-multiplier.eps" : "nu-2d-xsec-pz.eps");

}

int main()
{
  makePlots(true);
  makePlots(false);
}
