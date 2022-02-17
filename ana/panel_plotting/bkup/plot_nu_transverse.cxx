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

void makePlots(bool doMultipliers,bool doGenies,string location,string var)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  TFile f1(Form("%s_CV/CrossSection_CombinedPlaylists.root_%s.root",location.c_str(),var.c_str()));//Final result
  TFile f2(Form("%s_default/CrossSection_CombinedPlaylists.root_%s.root",location.c_str(),var.c_str()));//Default GENIE
  TFile f3(Form("%s_pion_rpa/CrossSection_CombinedPlaylists.root_%s.root",location.c_str(),var.c_str()));//GENIE+2p2h+RPA (aka no tune)
  MnvH2D* dataMnv=(MnvH2D*)f1.Get(Form("h_%s_ptmu_data_nobck_unfold_effcor_cross_section",var.c_str()));
  MnvH2D* mcMnv=(MnvH2D*)f1.Get(Form("h_%s_ptmu_mc_nobck_unfold_effcor_cross_section",var.c_str()));
  MnvH2D* nomGenieMnv = (MnvH2D*)f1.Get(Form("h_%s_ptmu_mc_nobck_unfold_effcor_cross_section",var.c_str()));//2
  MnvH2D* bestGenieMnv = (MnvH2D*)f1.Get(Form("h_%s_ptmu_mc_nobck_unfold_effcor_cross_section",var.c_str()));//3


  MnvH2D* mcMnv_qelike_qe = (MnvH2D*)f1.Get(Form("h_%s_ptmu_cross_section_qelike_qe",var.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_res = (MnvH2D*)f1.Get(Form("h_%s_ptmu_cross_section_qelike_res",var.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_dis = (MnvH2D*)f1.Get(Form("h_%s_ptmu_cross_section_qelike_dis",var.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_2p2h = (MnvH2D*)f1.Get(Form("h_%s_ptmu_cross_section_qelike_2p2h",var.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_2p2h_no_lowrec = (MnvH2D*)f1.Get(Form("h_%s_ptmu_cross_section_qelike_2p2h",var.c_str()));//Get from N track//3

  //  dataMnv->GetXaxis()->SetTitle("p_{||} (GeV)");

  dataMnv->Scale(1e40, "width");
  mcMnv->Scale(1e40, "width");
  nomGenieMnv->Scale(1e40, "width");
  bestGenieMnv->Scale(1e40, "width");

  mcMnv_qelike_qe->Scale(1e40,"width");
  mcMnv_qelike_res->Scale(1e40,"width");
  mcMnv_qelike_dis->Scale(1e40,"width");
  mcMnv_qelike_2p2h->Scale(1e40,"width");
  mcMnv_qelike_2p2h_no_lowrec->Scale(1e40,"width");

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
  histAndOpts.push_back(std::make_pair(mc,       "hist l"));
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
  double multipliers[]={3, 5, 5, 5, 5, 
			3, 3, 3, 3, 3,
			3, 3, 2, 2, 2, 
			2, 2, 2};

  double multipliersdphit[]={0.5,0.75,1,1,2,2,
			     2,2,2,2,4,4,
			     8,8,8,8,8,8,
			     8,8,8,8,8,8,
			     15,15,15,15,15,15};
  double multiplierspn[]={5,1,0.75,0.75,0.9,1,
			  2,2,2,3,3,3,
			  5,10,10,10,10,10,
			  100,100,100,100,100,100,
			  1000};


  GridCanvas* gc=NULL;
  if(var=="dalphat")gc= plotYAxis1D(histAndOpts, "p_{t}" ,"#delta#alpha_{t}" ,doMultipliers ? multipliers : NULL);
  if(var=="dphit")gc= plotYAxis1D(histAndOpts, "p_{t}" ,"#delta#phi_{t}" ,doMultipliers ? multipliersdphit : NULL);
  if(var=="pn")gc= plotYAxis1D(histAndOpts, "p_{t}" ,"p_{n}" ,doMultipliers ? multiplierspn : NULL);
  if(var=="dpt")gc= plotYAxis1D(histAndOpts, "p_{t}" ,"#delta p_{t}" ,doMultipliers ? multiplierspn : NULL);
  if(var=="dptx")gc= plotYAxis1D(histAndOpts, "p_{t}" ,"#delta p_{tx}" ,doMultipliers ? multiplierspn : NULL);
  if(var=="dpty")gc= plotYAxis1D(histAndOpts, "p_{t}" ,"#delta p_{ty}" ,doMultipliers ? multiplierspn : NULL);
  if(var=="signed")gc= plotYAxis1D(histAndOpts, "p_{t}" ,"sign" ,doMultipliers ? multiplierspn : NULL);
  if(var=="signeddalphat")gc= plotYAxis1D(histAndOpts, "p_{t}" ,"signed #delta#alpha_{t}" ,doMultipliers ? multiplierspn : NULL);
  if(var=="signeddphit")gc= plotYAxis1D(histAndOpts, "p_{t}" ,"signed #delta#phi_{t}" ,doMultipliers ? multiplierspn : NULL);

  // Set the y range manually. Can also use gc->Remax() to guess automatically
  if(var=="dalphat")gc->SetYLimits(0, 1.49);
  if(var=="dphit")gc->SetYLimits(0, 1.49);
  if(var=="pn")gc->SetYLimits(0,99.99);
  if(var=="dpt")gc->SetYLimits(0,99.99);
  if(var=="dptx")gc->SetYLimits(0,99.99);
  if(var=="dpty")gc->SetYLimits(0,99.99);
  if(var=="signed")gc->SetYLimits(0,99.99);
  if(var=="signeddalphat")gc->SetYLimits(0,99.99);
  if(var=="signeddphit")gc->SetYLimits(0,99.99);
  
  //Label thy axis
  if(var=="dalphat")gc->SetYTitle("d^{2}#sigma/dp_{T}d#delta#alpha_{t} (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  if(var=="dphit")gc->SetYTitle("d^{2}#sigma/dp_{T}d#delta#phi_{t} (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  if(var=="pn")gc->SetYTitle("d^{2}#sigma/dp_{T}d#p_{n} (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  if(var=="dpt")gc->SetYTitle("d^{2}#sigma/dp_{T}d#delta p_{t} (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  if(var=="dptx")gc->SetYTitle("d^{2}#sigma/dp_{T}d#delta p_{tx} (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  if(var=="dpty")gc->SetYTitle("d^{2}#sigma/dp_{T}d#delta p_{ty} (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  if(var=="signed")gc->SetYTitle("d^{2}#sigma/dp_{T}dsign (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  if(var=="signeddalphat")gc->SetYTitle("d^{2}#sigma/dp_{T}dsigned#delta#alpha_{t} (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  if(var=="signeddphit")gc->SetYTitle("d^{2}#sigma/dp_{T}dsigned#delta#phi{t} (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  
    
  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  //  TLegend* leg=new TLegend(0.17, 0.7, 0.31, 0.9);
  TLegend* leg=new TLegend(0.7, 0.1, 0.9, 0.3);
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
  if(var!="dphit")leg->Draw("SAME");

  if(doGenies){
    gc->Print(doMultipliers ? Form("nu-2d-xsec-genies-pt-%s-multiplier.eps",var.c_str()) : Form("nu-2d-xsec-genies-pt-%s.eps",var.c_str()));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-genies-pt-%s-multiplier.png",var.c_str()) : Form("nu-2d-xsec-genies-pt-%s.png",var.c_str()));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-genies-pt-%s-multiplier.C",var.c_str()) : Form("nu-2d-xsec-genies-pt-%s.C",var.c_str()));
  }
  else{
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-pt-%s-multiplier.eps",var.c_str()) : Form("nu-2d-xsec-comps-pt-%s.eps",var.c_str()));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-pt-%s-multiplier.png",var.c_str()) : Form("nu-2d-xsec-comps-pt-%s.png",var.c_str()));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-pt-%s-multiplier.C",var.c_str()) : Form("nu-2d-xsec-comps-pt-%s.C",var.c_str()));
  }

  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  // Values to multiply each bin by to get them on a similar range
  double  multipliers2[]={60, 15, 10, 5,
			  2,  2,  2 , 2,
			  2,  2,  10 , 40,
			  400};
  double multipliers2dphit[]={50,10,7,3,
			      1,0.75,0.75,0.75,
			      0.75,1,10,100,
			      3000};
  double multipliers2pn[]={15,7,3,2,
			   1,0.8,0.8,0.8,
			   0.8,2,10,50,
			   500};
			   


  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  GridCanvas* gc2=NULL;
  if(var=="dalphat")gc2 =plotXAxis1D(histAndOpts, "#delta#alpha_{t}" , "p_{t}" ,doMultipliers ? multipliers2 : NULL);
  if(var=="dphit")gc2 =plotXAxis1D(histAndOpts, "#delta#phi_{t}" , "p_{t}" ,doMultipliers ? multipliers2dphit : NULL);
  if(var=="pn")gc2=plotXAxis1D(histAndOpts, "p_{n}" , "p_{t}" ,doMultipliers ? multipliers2pn : NULL);
  if(var=="dpt")gc2=plotXAxis1D(histAndOpts, "#delta p_{t}" , "p_{t}" ,doMultipliers ? multipliers2pn : NULL);
  if(var=="dptx")gc2=plotXAxis1D(histAndOpts, "#delta p_{tx}" , "p_{t}" ,doMultipliers ? multipliers2pn : NULL);
  if(var=="dpty")gc2=plotXAxis1D(histAndOpts, "#delta p_{ty}" , "p_{t}" ,doMultipliers ? multipliers2pn : NULL);
  if(var=="signed")gc2=plotXAxis1D(histAndOpts, "signed" , "p_{t}" ,doMultipliers ? multipliers2pn : NULL);
  if(var=="signeddalphat")gc2=plotXAxis1D(histAndOpts, "signed #delta#alpha_{t}" , "p_{t}" ,doMultipliers ? multipliers2pn : NULL);
  if(var=="signeddphit")gc2=plotXAxis1D(histAndOpts, "signed #delta#phi_{t}" , "p_{t}" ,doMultipliers ? multipliers2pn : NULL);



  if(var=="dalphat") gc2->SetYLimits(0, 1.49);
  if(var=="dphit") gc2->SetYLimits(0, 1.49);
  if(var=="pn") gc2->SetYLimits(0, 99.99);
  if(var=="dpt") gc2->SetYLimits(0, 99.99);
  if(var=="dptx") gc2->SetYLimits(0, 99.99);
  if(var=="dpty") gc2->SetYLimits(0, 99.99);
  if(var=="signed") gc2->SetYLimits(0, 99.99);
  if(var=="signeddalphat") gc2->SetYLimits(0, 99.99);
  if(var=="signeddphit") gc2->SetYLimits(0, 99.99);

  if(var=="dalphat")gc2->SetYTitle("d^{2}#sigma/dp_{T}d#delta#alpha_{t} (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  if(var=="dphit")gc2->SetYTitle("d^{2}#sigma/dp_{T}d#delta#phi_{t} (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  if(var=="pn")gc2->SetYTitle("d^{2}#sigma/dp_{T}dp_{n} (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  if(var=="dpt")gc2->SetYTitle("d^{2}#sigma/dp_{T}d#delta p_{t} (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  if(var=="dptx")gc2->SetYTitle("d^{2}#sigma/dp_{T}d#delta p_{tx} (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  if(var=="dpty")gc2->SetYTitle("d^{2}#sigma/dp_{T}d#delta p_{ty} (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  if(var=="signed")gc2->SetYTitle("d^{2}#sigma/dp_{T}d#signed (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  if(var=="signeddalphat")gc2->SetYTitle("d^{2}#sigma/dp_{T}dsigned #delta#alpha_{t} (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  if(var=="signeddphitt")gc2->SetYTitle("d^{2}#sigma/dp_{T}dsigned #delta#phi_{t} (x10^{-40} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  gc2->Modified();
  if(doGenies){
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-genies-pz-%s-multiplier.eps",var.c_str()) : Form("nu-2d-xsec-genies-pz-%s.eps",var.c_str()));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-genies-pz-%s-multiplier.png",var.c_str()) : Form("nu-2d-xsec-genies-pz-%s.png",var.c_str()));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-genies-pz-%s-multiplier.C",var.c_str()) : Form("nu-2d-xsec-genies-pz-%s.C",var.c_str()));
  }
  else{
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-pz-%s-multiplier.eps",var.c_str()) : Form("nu-2d-xsec-comps-pz-%s.eps",var.c_str()));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-pz-%s-multiplier.png",var.c_str()) : Form("nu-2d-xsec-comps-pz-%s.png",var.c_str()));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-pz-%s-multiplier.C",var.c_str()) : Form("nu-2d-xsec-comps-pz-%s.C",var.c_str()));
  }

}

int main(int argc, char* argv[])
{
  makePlots(true,true,argv[1],argv[2]);
  makePlots(true,false,argv[1],argv[2]);
  makePlots(false,true,argv[1],argv[2]);
  makePlots(false,false,argv[1],argv[2]);

  return 0;
}
