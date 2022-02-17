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

TH2D* ChangeBins(TH2*input){

  //clone the 2D so I can change 1.5 to 1.49999... stupid hack;
  const int nBinsX = input->GetNbinsX();
  const int nBinsY = input->GetNbinsY();
  double xbinEdges[nBinsX+1];
  double ybinEdges[nBinsY+1];
  input->GetXaxis()->GetLowEdge(xbinEdges);
  xbinEdges[nBinsX]=15;
  input->GetYaxis()->GetLowEdge(ybinEdges);
  ybinEdges[nBinsY]=1.5;
  ybinEdges[0]=0.0001;
  string name = input->GetName();
  name+="_squash";
  TH2D* myhist = new TH2D(name.c_str(),name.c_str(),nBinsX,xbinEdges,nBinsY,ybinEdges);
  for(int i=0;i<nBinsX;i++){
    for(int j=0;j<nBinsY;j++){
      myhist->SetBinContent(i+1,j+1,input->GetBinContent(i+1,j+1));
      myhist->SetBinError(i+1,j+1,input->GetBinError(i+1,j+1));
    }
  }

  return myhist;
}


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

  dataMnv->Scale(1e39, "width");
  mcMnv->Scale(1e39, "width");

  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat_tmp=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* data_tmp=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mc_tmp=new TH2D(mcMnv->GetCVHistoWithStatError());
  TH2* mod_tmp1 = (TH2D*)f2.Get("gfg_norpa_Nieves_ma099");
  TH2* mod_tmp2 = (TH2D*)f2.Get("gfg_norpa_nomec_ma099");
  TH2* mod_tmp3 = (TH2D*)f2.Get("gfg_norpa_tem_ma099");
  TH2* mod_tmp4 = (TH2D*)f2.Get("gfg_rpa_Nieves_ma099");
  TH2* mod_tmp5 = (TH2D*)f2.Get("gfg_rpa_tem_ma099");
  TH2* mod_tmp6 = (TH2D*)f2.Get("lfg_norpa_tem_ma099");
  TH2* mod_tmp7 = (TH2D*)f2.Get("sf_norpa_nomec_ma099");
  TH2* mod_tmp8 = (TH2D*)f2.Get("Untuned_GENIE");
  TH2* mod_tmp9 = (TH2D*)f2.Get("MnvGENIE_noRPA_noMEC");
  TH2* mod_tmp10 = (TH2D*)f2.Get("MnvGENIE_RPA_noMEC");
  TH2* mod_tmp11 = (TH2D*)f2.Get("MnvGENIE_noRPA_MEC");
  TH2* mod_tmp12 = (TH2D*)f2.Get("MnvGENIE_RPA_MEC_notune");


  //do width before I squash the first bin to zero suppress
  mod_tmp1->Scale(1e39, "width");
  mod_tmp2->Scale(1e39, "width");
  mod_tmp3->Scale(1e39, "width");
  mod_tmp4->Scale(1e39, "width");
  mod_tmp5->Scale(1e39, "width");
  mod_tmp6->Scale(1e39, "width");
  mod_tmp7->Scale(1e39, "width");
  mod_tmp8->Scale(1e39, "width");
  mod_tmp9->Scale(1e39, "width");
  mod_tmp10->Scale(1e39, "width");
  mod_tmp11->Scale(1e39, "width");
  mod_tmp12->Scale(1e39, "width");
  //Now squash!!
  TH2* dataStat =ChangeBins(dataStat_tmp);
  TH2* data =ChangeBins(data_tmp);
  TH2* mc =ChangeBins(mc_tmp);
  TH2* mod1 =ChangeBins(mod_tmp1);
  TH2* mod2 =ChangeBins(mod_tmp2);
  TH2* mod3 =ChangeBins(mod_tmp3);
  TH2* mod4 =ChangeBins(mod_tmp4);
  TH2* mod5 =ChangeBins(mod_tmp5);
  TH2* mod6 =ChangeBins(mod_tmp6);
  TH2* mod7 =ChangeBins(mod_tmp7);
  TH2* mod8 =ChangeBins(mod_tmp8);
  TH2* mod9 =ChangeBins(mod_tmp9);
  TH2* mod10 =ChangeBins(mod_tmp10);
  TH2* mod11 =ChangeBins(mod_tmp11);
  TH2* mod12 =ChangeBins(mod_tmp12);
  
  
  vector<int> mycolors = getColors(1);


  // These line and marker styles will be propagated to the 1D plots
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);
  
  mod8->SetLineColor(kBlue);
  mod8->SetLineWidth(2);

  mod1->SetLineColor(mycolors[0]);
  mod1->SetLineWidth(2);
  mod2->SetLineColor(mycolors[1]);
  mod2->SetLineWidth(2);
  mod3->SetLineColor(mycolors[2]);
  mod3->SetLineWidth(2);
  mod4->SetLineColor(mycolors[3]);
  mod4->SetLineWidth(2);
  mod5->SetLineColor(mycolors[4]);
  mod5->SetLineWidth(2);
  mod6->SetLineColor(mycolors[5]);
  mod6->SetLineWidth(2);
  mod7->SetLineColor(mycolors[6]);
  mod7->SetLineWidth(2);
  mod9->SetLineColor(mycolors[7]);
  mod9->SetLineWidth(2);
  mod10->SetLineColor(mycolors[8]);
  mod10->SetLineWidth(2);
  mod11->SetLineColor(mycolors[9]);
  mod11->SetLineWidth(2);
  mod12->SetLineColor(mycolors[10]);
  mod12->SetLineWidth(2);



 

  //all ratio to untuned genie (mod8)
  data->Divide(mod8);
  dataStat->Divide(mod8);
  mod1->Divide(mod8);
  mod2->Divide(mod8);
  mod3->Divide(mod8);
  mod4->Divide(mod8);
  mod5->Divide(mod8);
  mod6->Divide(mod8);
  mod7->Divide(mod8);
  mod9->Divide(mod8);
  mod10->Divide(mod8);
  mod11->Divide(mod8);
  mod12->Divide(mod8);
  mc->Divide(mod8);
  mod8->Divide(mod8);



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
  // histAndOpts.push_back(std::make_pair(mod1,     "hist"));
  // histAndOpts.push_back(std::make_pair(mod2,     "hist"));
  // histAndOpts.push_back(std::make_pair(mod3,     "hist"));
  // histAndOpts.push_back(std::make_pair(mod4,     "hist"));//4
  // histAndOpts.push_back(std::make_pair(mod5,     "hist"));
  // histAndOpts.push_back(std::make_pair(mod6,     "hist"));
  // histAndOpts.push_back(std::make_pair(mod7,     "hist"));
  histAndOpts.push_back(std::make_pair(mod8,     "hist"));
  // histAndOpts.push_back(std::make_pair(mod9,     "hist"));//9
  // histAndOpts.push_back(std::make_pair(mod10,     "hist"));
  // histAndOpts.push_back(std::make_pair(mod11,     "hist"));
  histAndOpts.push_back(std::make_pair(mod12,     "hist"));
  histAndOpts.push_back(std::make_pair(dataStat, "graph e"));
  histAndOpts.push_back(std::make_pair(data,     "graph ep"));//14

  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range
  double multipliers[]={1, 1, 1, 1,
                        2, 2, 2, 2,
                        10, 20, 30};
  
  // for(int i=0;i<histAndOpts.size()-4;i++){
  //   std::vector<std::pair<TH2*, const char*> > histAndOpts_subset;
  //   histAndOpts_subset.push_back(histAndOpts[0]);
  //   histAndOpts_subset.push_back(histAndOpts[8]);
  //   histAndOpts_subset.push_back(histAndOpts[i+1]);
  //   histAndOpts_subset.push_back(histAndOpts[13]);
  //   histAndOpts_subset.push_back(histAndOpts[14]);
  int i = 0;
    GridCanvas* gc=plotpT1DAntiNu(histAndOpts, doMultipliers ? multipliers : NULL);
    // Set the y range manually. Can also use gc->Remax() to guess automatically
    //  gc->SetYLimits(0, 3.9);
    gc->SetYLimits(0, 2.49);
    
    //  leg->Draw();
    gc->SetYTitle("Ratio to Default GENIE");
    gc->Modified();
    gc->Print(doMultipliers ? Form("anti-nu2d-xsec-models-pt-multiplier-model-%d.eps",i) : Form("anti-nu2d-xsec-models-pt-model-%d.eps",i));
    gc->Print(doMultipliers ? Form("anti-nu2d-xsec-models-pt-multiplier-model-%d.png",i) : Form("anti-nu2d-xsec-models-pt-model-%d.png",i));
    gc->Print(doMultipliers ? Form("anti-nu2d-xsec-models-pt-multiplier-model-%d.C",i) : Form("anti-nu2d-xsec-models-pt-model-%d.C",i));
    //}
  
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TLegend* leg=new TLegend(0, 0, 1, 1);
  //  leg->SetNColumns(2);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.07);
  leg->AddEntry(mc, "MINERvA Tune v1", "l");
  leg->AddEntry(data, "MINERvA data", "lpe");
  leg->AddEntry(mod8, "GENIE 2.8.4", "l");
  leg->AddEntry(mod9, "Pion Tune","l");
  leg->AddEntry(mod10, "RPA","l");
  leg->AddEntry(mod11, "2p2h","l");
  leg->AddEntry(mod12, "RPA+2p2h","l");
  leg->AddEntry(mod1, "gfg_norpa_Nieves_ma099","l");
  leg->AddEntry(mod2, "gfg_norpa_nomec_ma099","l");
  leg->AddEntry(mod3, "gfg_norpa_tem_ma099","l");
  leg->AddEntry(mod4, "gfg_rpa_Nieves_ma099","l");
  leg->AddEntry(mod5, "gfg_rpa_tem_ma099","l");
  leg->AddEntry(mod6, "lfg_norpa_tem_ma099","l");
  leg->AddEntry(mod7, "sf_norpa_nomec_ma099","l");



  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  // Values to multiply each bin by to get them on a similar range
  double  multipliers2[]={5, 1, 1, 1, 1,
                         1, 3, 15, 100, 500};

  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  for(int i=0;i<histAndOpts.size()-4;i++){
    std::vector<std::pair<TH2*, const char*> > histAndOpts_subset;
    histAndOpts_subset.push_back(histAndOpts[0]);
    histAndOpts_subset.push_back(histAndOpts[8]);
    histAndOpts_subset.push_back(histAndOpts[i+1]);
    histAndOpts_subset.push_back(histAndOpts[13]);
    histAndOpts_subset.push_back(histAndOpts[14]);
    GridCanvas* gc2=plotpz1DAntiNu(histAndOpts, doMultipliers ? multipliers2 : NULL);
    //  gc2->SetYLimits(0, 3.9);
    gc2->SetYLimits(0, 2.49);
    gc2->SetYTitle("Ratio to Default GENIE");
    gc2->Modified();
    gc2->Print(doMultipliers ? Form("anti-nu2d-xsec-models-pz-multiplier-model-%d.eps",i) : Form("anti-nu2d-xsec-models-pz-model-%d.eps",i));
    gc2->Print(doMultipliers ? Form("anti-nu2d-xsec-models-pz-multiplier-model-%d.png",i) : Form("anti-nu2d-xsec-models-pz-model-%d.png",i));
    gc2->Print(doMultipliers ? Form("anti-nu2d-xsec-models-pz-multiplier-model-%d.C",i) : Form("anti-nu2d-xsec-models-pz-model-%d.C",i));
  }
  TCanvas c1("Legend","Legend",10,10,800,500);
  leg->Draw();
  c1.Print("anti-nu2d-xsec-models-Legend.eps");
  c1.Print("anti-nu2d-xsec-models-Legend.png");
  c1.Print("anti-nu2d-xsec-models-Legend.C");
}

int main()
{
  //  makePlots(true);
  makePlots(false);

  return 0;
}
