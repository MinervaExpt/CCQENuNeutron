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

void makePlots(bool doMultipliers, bool doTracks, bool doBkgSub, string location)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  //three files 1-track, 2+track, N-track
  //CV
  TFile f1(Form("%s_CV/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal_SpecialSampleIncluded.root",location.c_str()));//1track
  TFile f2(Form("%s_CV/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_SpecialSampleIncluded.root",location.c_str()));//2track
  TFile f3(Form("%s_CV/MuonEventSelection_MakeFlux-1_Multiplicity-0_Sample-Signal_SpecialSampleIncluded.root",location.c_str()));//Ntrack
  //CV without low recoil tune
  TFile f4(Form("%s_pion_rpa/MuonEventSelection_MakeFlux-1_Multiplicity-0_Sample-Signal_SpecialSampleIncluded.root",location.c_str()));//NTrack
  //Constraint File
  TFile f5(Form("%s_CV/SideBandFit_SpecialSampleIncluded.root",location.c_str()));

  //need pzmuptmu
  MnvH2D* mcMnv=(MnvH2D*)f3.Get("h_pzmu_ptmu_mc");//Get from N track
  MnvH2D* dataMnv = (MnvH2D*)f3.Get("h_pzmu_ptmu_data");//Get from N track
  MnvH2D* SingleTrackMnv = (MnvH2D*)f1.Get("h_pzmu_ptmu_data");//Get from 1 track
  MnvH2D* MultiTrackMnv = (MnvH2D*)f2.Get("h_pzmu_ptmu_data");//Get from 2+ track
  MnvH2D* SingleTrackMCMnv = (MnvH2D*)f1.Get("h_pzmu_ptmu_mc");//Get from 1 track
  MnvH2D* MultiTrackMCMnv = (MnvH2D*)f2.Get("h_pzmu_ptmu_mc");//Get from 2+ track

  MnvH2D* mcMnv_qelike_qe = (MnvH2D*)f3.Get("h_pzmu_ptmu_qelike_qe");//Get from N track
  MnvH2D* mcMnv_qelike_res = (MnvH2D*)f3.Get("h_pzmu_ptmu_qelike_res");//Get from N track
  MnvH2D* mcMnv_qelike_dis = (MnvH2D*)f3.Get("h_pzmu_ptmu_qelike_dis");//Get from N track
  MnvH2D* mcMnv_qelike_2p2h = (MnvH2D*)f3.Get("h_pzmu_ptmu_qelike_oth");//Get from N track
  MnvH2D* mcMnv_qelike_2p2h_no_lowrec = (MnvH2D*)f4.Get("h_pzmu_ptmu_qelike_oth");//Get from N track
  //bkgs
  MnvH2D* mcMnv_1track_bkg = (MnvH2D*)f1.Get("h_pzmu_ptmu_qelikenot");//Get from 1 track
  MnvH2D* mcMnv_2track_bkg = (MnvH2D*)f2.Get("h_pzmu_ptmu_qelikenot");//Get from 2+ track
  //the constraint
  MnvH2D* bkgConstraint_1track = (MnvH2D*)f5.Get("h_weights_1track_pzptbins_qelikenot");//get from 1 track
  MnvH2D* bkgConstraint_2track = (MnvH2D*)f5.Get("h_weights_2track_pzptbins_qelikenot");//get from 2 track
  //apply constraint
  mcMnv_1track_bkg->Multiply(mcMnv_1track_bkg,bkgConstraint_1track);
  mcMnv_2track_bkg->Multiply(mcMnv_2track_bkg,bkgConstraint_2track);
  
  MnvH2D* mcMnv_bkg = (MnvH2D*)mcMnv_1track_bkg->Clone("tmp");
  mcMnv_bkg->Add(mcMnv_2track_bkg);
  
  dataMnv->GetXaxis()->SetTitle("p_{||} (GeV)");

  if(doBkgSub){
    dataMnv->Add(mcMnv_bkg,-1);
    mcMnv->Add(mcMnv_bkg,-1);
    SingleTrackMnv->Add(mcMnv_1track_bkg,-1);
    MultiTrackMnv->Add(mcMnv_2track_bkg,-1);
    SingleTrackMCMnv->Add(mcMnv_1track_bkg,-1);
    MultiTrackMCMnv->Add(mcMnv_2track_bkg,-1);
  }
  
  dataMnv->Scale(1e-3, "width");
  mcMnv->Scale(1e-3, "width");
  SingleTrackMnv->Scale(1e-3, "width");
  MultiTrackMnv->Scale(1e-3, "width");
  SingleTrackMCMnv->Scale(1e-3, "width");
  MultiTrackMCMnv->Scale(1e-3, "width");
  mcMnv_qelike_qe->Scale(1e-3, "width");
  mcMnv_qelike_res->Scale(1e-3, "width");
  mcMnv_qelike_dis->Scale(1e-3, "width");
  mcMnv_qelike_2p2h->Scale(1e-3, "width");
  mcMnv_qelike_2p2h_no_lowrec->Scale(1e-3, "width");
  mcMnv_bkg->Scale(1e-3, "width");

  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* data=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());
  TH2* mc_qelike_qe = new TH2D(mcMnv_qelike_qe->GetCVHistoWithStatError());
  TH2* mc_qelike_res = new TH2D(mcMnv_qelike_res->GetCVHistoWithStatError());
  TH2* mc_qelike_dis = new TH2D(mcMnv_qelike_dis->GetCVHistoWithStatError());
  TH2* mc_qelike_2p2h = new TH2D(mcMnv_qelike_2p2h->GetCVHistoWithStatError());
  TH2* mc_qelike_2p2h_no_lowrec = new TH2D(mcMnv_qelike_2p2h_no_lowrec->GetCVHistoWithStatError());
  TH2* mc_bkg = new TH2D(mcMnv_bkg->GetCVHistoWithStatError());
  TH2* data_1track = new TH2D(SingleTrackMnv->GetCVHistoWithStatError());
  TH2* data_2track = new TH2D(MultiTrackMnv->GetCVHistoWithStatError());
  TH2* mc_1track = new TH2D(SingleTrackMCMnv->GetCVHistoWithStatError());
  TH2* mc_2track = new TH2D(MultiTrackMCMnv->GetCVHistoWithStatError());



  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = getColors(2);
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);

  mc_1track->SetLineColor(kRed);
  mc_1track->SetLineStyle(2);
  mc_2track->SetLineColor(kRed);
  mc_2track->SetLineStyle(9);
  //need to add signal and bkg colors
  mc_qelike_qe->SetLineColor(mycolors[3]);
  mc_qelike_res->SetLineColor(mycolors[4]);
  mc_qelike_dis->SetLineColor(mycolors[5]);
  mc_qelike_2p2h->SetLineColor(mycolors[6]);
  mc_qelike_2p2h_no_lowrec->SetLineColor(mycolors[6]);
  mc_qelike_2p2h_no_lowrec->SetLineStyle(2);
  mc_bkg->SetLineColor(mycolors[10]);
  //need to add 1track and 2 track
  
  // These line and marker styles will be propagated to the 1D plots
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(0.7);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2);

  dataStat->SetMarkerStyle(1);
  dataStat->SetLineColor(kBlack);
  dataStat->SetLineWidth(2);

  data_1track->SetLineColor(kBlack);
  data_1track->SetLineStyle(2);
  data_2track->SetLineColor(kBlack);
  data_2track->SetLineStyle(9);


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

  //Do by track breakdown
  if(doTracks){
    histAndOpts.push_back(std::make_pair(mc_1track,       "histl"));
    histAndOpts.push_back(std::make_pair(mc_2track,       "histl"));
    histAndOpts.push_back(std::make_pair(data_1track,       "histl"));
    histAndOpts.push_back(std::make_pair(data_2track,       "histl"));
  }
  
  //do by physics process
  else{
    histAndOpts.push_back(std::make_pair(mc_qelike_qe,       "histl"));
    histAndOpts.push_back(std::make_pair(mc_qelike_res,       "histl"));
    histAndOpts.push_back(std::make_pair(mc_qelike_dis,       "histl"));
    histAndOpts.push_back(std::make_pair(mc_qelike_2p2h,       "histl"));
    histAndOpts.push_back(std::make_pair(mc_qelike_2p2h_no_lowrec,       "histl"));
    histAndOpts.push_back(std::make_pair(mc_bkg,       "histl"));
  }

  histAndOpts.push_back(std::make_pair(dataStat, "graph e"));
  histAndOpts.push_back(std::make_pair(data,     "graph ep"));



  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range
  double multipliers[]={1, 1, 1, 1,
                        1, 2, 2, 2,
                        10, 20, 30, 50};

  GridCanvas* gc=plotpT1D(histAndOpts, doMultipliers ? multipliers : NULL);
  // Set the y range manually. Can also use gc->Remax() to guess automatically
  gc->SetYLimits(0, 75);
  gc->SetYTitle("Event Rate(x10^{-3}) per GeV^{2}");
  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TLegend* leg=new TLegend(0.17, 0.7, 0.31, 0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(data, "MINERvA data", "lpe");
  leg->AddEntry(mc, "MINERvA Tune", "l");
  if(doTracks){
    leg->AddEntry(mc_1track, "MC 1-track","l");
    leg->AddEntry(mc_2track, "MC Multi-track","l");
    leg->AddEntry(data_1track, "data 1-track","l");
    leg->AddEntry(data_2track, "data Multi-track","l");
  }
  else{
    leg->AddEntry(mc_qelike_qe,"QE","l");
    leg->AddEntry(mc_qelike_res,"Resonant","l");
    leg->AddEntry(mc_qelike_dis,"DIS","l");
    leg->AddEntry(mc_qelike_2p2h,"2p2h","l");
    leg->AddEntry(mc_qelike_2p2h_no_lowrec,"2p2h without fit","l");
    leg->AddEntry(mc_bkg,"Background","l");

  }

  leg->Draw("SAME");
  if(doTracks){
    if(doBkgSub){
      gc->Print(doMultipliers ? "nu-2d-evtrate-bkgsub-trackmult-pt-multiplier.eps" : "nu-2d-evtrate-bkgsub-trackmult-pt.eps");
      gc->Print(doMultipliers ? "nu-2d-evtrate-bkgsub-trackmult-pt-multiplier.png" : "nu-2d-evtrate-bkgsub-trackmult-pt.png");
      gc->Print(doMultipliers ? "nu-2d-evtrate-bkgsub-trackmult-pt-multiplier.C" : "nu-2d-evtrate-bkgsub-trackmult-pt.C");
    }
    else{
      gc->Print(doMultipliers ? "nu-2d-evtrate-trackmult-pt-multiplier.eps" : "nu-2d-evtrate-trackmult-pt.eps");
      gc->Print(doMultipliers ? "nu-2d-evtrate-trackmult-pt-multiplier.png" : "nu-2d-evtrate-trackmult-pt.png");
      gc->Print(doMultipliers ? "nu-2d-evtrate-trackmult-pt-multiplier.C" : "nu-2d-evtrate-trackmult-pt.C");
    }
  }
  else{
    gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier.eps" : "nu-2d-evtrate-model-pt.eps");
    gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier.png" : "nu-2d-evtrate-model-pt.png");
    gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier.C" : "nu-2d-evtrate-model-pt.C");
  }

  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  // Values to multiply each bin by to get them on a similar range
  double  multipliers2[]={15, 3, 2, 1, 1,
			  1, 1, 1, 2, 4,
			  15, 50, 200};

  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  GridCanvas* gc2=plotpz1D(histAndOpts, doMultipliers ? multipliers2 : NULL);
  gc2->SetYLimits(0, 75);
  gc2->SetYTitle("Event Rate(x10^{-3}) per GeV^{2}");
  gc2->Modified();
  if(doTracks){
    if(doBkgSub){
      gc2->Print(doMultipliers ? "nu-2d-evtrate-bkgsub-trackmult-pz-multiplier.eps" : "nu-2d-evtrate-bkgsub-trackmult-pz.eps");
      gc2->Print(doMultipliers ? "nu-2d-evtrate-bkgsub-trackmult-pz-multiplier.png" : "nu-2d-evtrate-bkgsub-trackmult-pz.png");
      gc2->Print(doMultipliers ? "nu-2d-evtrate-bkgsub-trackmult-pz-multiplier.C" : "nu-2d-evtrate-bkgsub-trackmult-pz.C");
    }
    else{
      gc2->Print(doMultipliers ? "nu-2d-evtrate-trackmult-pz-multiplier.eps" : "nu-2d-evtrate-trackmult-pz.eps");
      gc2->Print(doMultipliers ? "nu-2d-evtrate-trackmult-pz-multiplier.png" : "nu-2d-evtrate-trackmult-pz.png");
      gc2->Print(doMultipliers ? "nu-2d-evtrate-trackmult-pz-multiplier.C" : "nu-2d-evtrate-trackmult-pz.C");
    }
  }
  else{
    gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier.eps" : "nu-2d-evtrate-model-pz.eps");
    gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier.png" : "nu-2d-evtrate-model-pz.png");
    gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier.C" : "nu-2d-evtrate-model-pz.C");
  }

}

int main(int argc, char* argv[])
{

  string location = argv[1];
  //multipliers
    makePlots(true,true,false,location);
  makePlots(true,true,true,location);
  makePlots(true,false,false,location);
  //standard
  makePlots(false,true,false,location);
  makePlots(false,true,true,location);
  makePlots(false,false,false,location);

  return 0;
}
