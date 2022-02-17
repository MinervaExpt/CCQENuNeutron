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

void makePlots(bool doMultipliers,bool doGenies,string location)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  TFile f1(Form("%s_CV/CrossSection_per_nucleon_iterations_4_SpecialSampleIncluded.root_pzmu.root",location.c_str()));//Final result
  TFile f2(Form("%s_default/CrossSection_per_nucleon_iterations_4_CombinedPlaylists.root_pzmu.root",location.c_str()));//Default GENIE
  TFile f3(Form("%s_pion_rpa/CrossSection_per_nucleon_iterations_4_SpecialSampleIncluded.root_pzmu.root",location.c_str()));//GENIE+2p2h+RPA (aka no tune)
  MnvH2D* dataMnv=(MnvH2D*)f1.Get("h_pzmu_ptmu_data_nobck_unfold_effcor_cross_section");
  MnvH2D* mcMnv=(MnvH2D*)f1.Get("h_pzmu_ptmu_mc_nobck_unfold_effcor_cross_section");
  MnvH2D* nomGenieMnv = (MnvH2D*)f2.Get("h_pzmu_ptmu_mc_nobck_unfold_effcor_cross_section");
  MnvH2D* bestGenieMnv = (MnvH2D*)f3.Get("h_pzmu_ptmu_mc_nobck_unfold_effcor_cross_section");


  MnvH2D* mcMnv_qelike_qe = (MnvH2D*)f1.Get("h_pzmu_ptmu_cross_section_qelike_qe");//Get from N track
  MnvH2D* mcMnv_qelike_res = (MnvH2D*)f1.Get("h_pzmu_ptmu_cross_section_qelike_res");//Get from N track
  MnvH2D* mcMnv_qelike_dis = (MnvH2D*)f1.Get("h_pzmu_ptmu_cross_section_qelike_dis");//Get from N track
  MnvH2D* mcMnv_qelike_2p2h = (MnvH2D*)f1.Get("h_pzmu_ptmu_cross_section_qelike_2p2h");//Get from N track
  MnvH2D* mcMnv_qelike_2p2h_no_lowrec = (MnvH2D*)f3.Get("h_pzmu_ptmu_cross_section_qelike_2p2h");//Get from N track

  dataMnv->GetXaxis()->SetTitle("p_{||} (GeV)");

  dataMnv->Scale(1e39, "width");
  mcMnv->Scale(1e39, "width");
  nomGenieMnv->Scale(1e39, "width");
  bestGenieMnv->Scale(1e39, "width");

  mcMnv_qelike_qe->Scale(1e39,"width");
  mcMnv_qelike_res->Scale(1e39,"width");
  mcMnv_qelike_dis->Scale(1e39,"width");
  mcMnv_qelike_2p2h->Scale(1e39,"width");
  mcMnv_qelike_2p2h_no_lowrec->Scale(1e39,"width");

  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* data=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());
  TH2* nomGenie=new TH2D(nomGenieMnv->GetCVHistoWithStatError());
  TH2* bestGenie = new TH2D(bestGenieMnv->GetCVHistoWithStatError());

  TH2* mc_qelike_qe = new TH2D(mcMnv_qelike_qe->GetCVHistoWithStatError());
  TH2* mc_qelike_res = new TH2D(mcMnv_qelike_res->GetCVHistoWithStatError());
  TH2* mc_qelike_dis = new TH2D(mcMnv_qelike_dis->GetCVHistoWithStatError());
  TH2* mc_qelike_2p2h = new TH2D(mcMnv_qelike_2p2h->GetCVHistoWithStatError());
  TH2* mc_qelike_2p2h_no_lowrec = new TH2D(mcMnv_qelike_2p2h_no_lowrec->GetCVHistoWithStatError());

  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = getColors(2);
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);
  
  nomGenie->SetLineColor(kBlue);
  nomGenie->SetLineWidth(2);

  bestGenie->SetLineColor(mycolors[10]);
  bestGenie->SetLineWidth(2);

  mc_qelike_qe->SetLineColor(mycolors[3]);
  mc_qelike_res->SetLineColor(mycolors[4]);
  mc_qelike_dis->SetLineColor(mycolors[5]);
  mc_qelike_2p2h->SetLineColor(mycolors[6]);
  mc_qelike_2p2h_no_lowrec->SetLineColor(mycolors[6]);
  mc_qelike_2p2h_no_lowrec->SetLineStyle(3);

  // These line and marker styles will be propagated to the 1D plots
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(0.5);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2);

  dataStat->SetLineColor(kBlack);


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
  histAndOpts.push_back(std::make_pair(dataStat, "histpe1"));
  histAndOpts.push_back(std::make_pair(mc,       "graph0 l"));
  if(doGenies){
    histAndOpts.push_back(std::make_pair(nomGenie, "hist l"));
    histAndOpts.push_back(std::make_pair(bestGenie, "hist l"));
  }
  else{
    histAndOpts.push_back(std::make_pair(mc_qelike_qe,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_res,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_dis,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_2p2h,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_2p2h_no_lowrec,       "hist l"));
  }
  histAndOpts.push_back(std::make_pair(data,     "histpe1"));



  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range
  double multipliers[]={1, 1, 1, 1,
                        2, 4, 6, 8,
                        10, 20, 30, 50};

  GridCanvas* gc=plotpT1D(histAndOpts, doMultipliers ? multipliers : NULL);
  // Set the y range manually. Can also use gc->Remax() to guess automatically
  gc->SetYLimits(0, 4.59);
  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TLegend* leg=new TLegend(0.17, 0.7, 0.31, 0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(data, "MINERvA data", "lpe");
  leg->AddEntry(mc, "MnvGENIE", "l");
  if(doGenies){
    leg->AddEntry(nomGenie, "GENIE 2.8.4", "l");
    leg->AddEntry(bestGenie, "2p2h and RPA", "l");
  }
  else{
    leg->AddEntry(mc_qelike_qe,"QE","l");
    leg->AddEntry(mc_qelike_res,"Resonant","l");
    leg->AddEntry(mc_qelike_dis,"DIS","l");
    leg->AddEntry(mc_qelike_2p2h,"2p2h","l");
    leg->AddEntry(mc_qelike_2p2h_no_lowrec,"2p2h without fit","l");
  }
  leg->Draw("SAME");
  if(doGenies){
    gc->Print(doMultipliers ? "nu-2d-xsec-genies-pt-multiplier.eps" : "nu-2d-xsec-genies-pt.eps");
    gc->Print(doMultipliers ? "nu-2d-xsec-genies-pt-multiplier.png" : "nu-2d-xsec-genies-pt.png");
    gc->Print(doMultipliers ? "nu-2d-xsec-genies-pt-multiplier.C" : "nu-2d-xsec-genies-pt.C");
  }
  else{
    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier.eps" : "nu-2d-xsec-comps-pt.eps");
    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier.png" : "nu-2d-xsec-comps-pt.png");
    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier.C" : "nu-2d-xsec-comps-pt.C");
  }

  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  // Values to multiply each bin by to get them on a similar range
  double  multipliers2[]={30, 8, 3, 2, 1,
			  1, 1, 1, 2, 4,
                          20, 100, 1000};

  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  GridCanvas* gc2=plotpz1D(histAndOpts, doMultipliers ? multipliers2 : NULL);
  gc2->SetYLimits(0, 5.49);
  gc2->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  gc2->Modified();
  if(doGenies){
    gc2->Print(doMultipliers ? "nu-2d-xsec-genies-pz-multiplier.eps" : "nu-2d-xsec-genies-pz.eps");
    gc2->Print(doMultipliers ? "nu-2d-xsec-genies-pz-multiplier.png" : "nu-2d-xsec-genies-pz.png");
    gc2->Print(doMultipliers ? "nu-2d-xsec-genies-pz-multiplier.C" : "nu-2d-xsec-genies-pz.C");
  }
  else{
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.eps" : "nu-2d-xsec-comps-pz.eps");
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.png" : "nu-2d-xsec-comps-pz.png");
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.C" : "nu-2d-xsec-comps-pz.C");
  }

}

int main(int argc, char* argv[])
{
  makePlots(true,true,argv[1]);
  makePlots(true,false,argv[1]);
  makePlots(false,true,argv[1]);
  makePlots(false,false,argv[1]);

  return 0;
}
