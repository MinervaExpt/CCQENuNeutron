#include "myPlotStyle.h"
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

void makePlots(bool doMultipliers,bool track_2,bool notRatio)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  gStyle->SetNdivisions(555, "XYZ");


  TFile *f = NULL;
  if(track_2) f = new TFile("/minerva/app/users/drut1186/cmtuser/Minerva_v10r8p9_CCQENu2DLE_hadronv2/Ana/CCQENu2DLE/ana/plot_macros_pub/Vertex/VtxEnergy_2track.root");
  else f = new TFile("/minerva/app/users/drut1186/cmtuser/Minerva_v10r8p9_CCQENu2DLE_hadronv2/Ana/CCQENu2DLE/ana/plot_macros_pub/Vertex/VtxEnergy_1track.root");

  MnvH2D* dataMnv=(MnvH2D*)f->Get("h_vtx_new_data"); // data
  MnvH2D* m1=(MnvH2D*)f->Get("h_vtx_new_qelike_qe");
  MnvH2D* m2=(MnvH2D*)f->Get("h_vtx_new_qelike_res");
  MnvH2D* m3=(MnvH2D*)f->Get("h_vtx_new_qelike_dis");
  MnvH2D* m4=(MnvH2D*)f->Get("h_vtx_new_qelike_oth");
  MnvH2D* m5=(MnvH2D*)f->Get("h_vtx_new_qelikenot");
  MnvH2D* m6=(MnvH2D*)f->Get("h_vtx_new_qelike");
  //Get the total
  cout << "Total" << endl;
  MnvH2D* mtotal=(MnvH2D*)m6->Clone("h_vtx_new_mc");
  mtotal->Add(m5);

  cout << m5->GetNbinsX() << "\t" << m6->GetNbinsX() << endl;
  //norm them
  dataMnv->Scale(1.0,"width");
  m1->Scale(1.0,"width");
  m2->Scale(1.0,"width");
  m3->Scale(1.0,"width");
  m4->Scale(1.0,"width");
  m5->Scale(1.0,"width");
  mtotal->Scale(1.0,"width");

  dataMnv->GetXaxis()->SetTitle("Vertex Energy (MeV)");

  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  cout << "Getting the TH2" << endl;

  TH2* data=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* datasys=new TH2D(dataMnv->GetCVHistoWithError());

  TH2* mcstatsys = new TH2D(mtotal->GetCVHistoWithError()); 
  TH2* mcqe = new TH2D(m1->GetCVHistoWithStatError());
  TH2* mcres = new TH2D(m2->GetCVHistoWithStatError());
  TH2* mcdis = new TH2D(m3->GetCVHistoWithStatError());
  TH2* mc2p2h = new TH2D(m4->GetCVHistoWithStatError());
  TH2* mcnot = new TH2D(m5->GetCVHistoWithStatError());



  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = getColors(2);

  // These line and marker styles will be propagated to the 1D plots

  data->SetLineColor(kBlack);
  data->SetLineWidth(4);
  data->SetMarkerColor(kBlack);
  data->SetMarkerStyle(8);


  datasys->SetLineColor(kBlack);
  datasys->SetLineWidth(4);
  datasys->SetMarkerColor(kBlack);
  datasys->SetMarkerStyle(8);

  mcstatsys->SetLineColor(kRed);
  mcstatsys->SetLineWidth(4);
  mcstatsys->SetMarkerColor(kRed);
  
  mcqe->SetLineColor(mycolors[2]);
  mcqe->SetLineWidth(4);
  mcqe->SetMarkerColor(mycolors[2]);

  mcres->SetLineColor(mycolors[3]);
  mcres->SetLineWidth(4);
  mcres->SetMarkerColor(mycolors[3]);

  mcdis->SetLineColor(mycolors[4]);
  mcdis->SetLineWidth(4);
  mcdis->SetMarkerColor(mycolors[4]);

  mc2p2h->SetLineColor(mycolors[5]);
  mc2p2h->SetLineWidth(4);
  mc2p2h->SetMarkerColor(mycolors[5]);

  mcnot->SetLineColor(mycolors[6]);
  mcnot->SetLineWidth(4);
  mcnot->SetMarkerColor(mycolors[6]);
  

  



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
  histAndOpts.push_back(std::make_pair(mcstatsys,    "histe"));
  histAndOpts.push_back(std::make_pair(data,     "p"));
  histAndOpts.push_back(std::make_pair(mcqe,     "histl"));
  histAndOpts.push_back(std::make_pair(mcres,     "histl"));
  histAndOpts.push_back(std::make_pair(mcdis,     "histl"));
  histAndOpts.push_back(std::make_pair(mc2p2h,     "histl"));
  histAndOpts.push_back(std::make_pair(mcnot,     "histl"));

  histAndOpts.push_back(std::make_pair(datasys,     "pe"));


  //  histAndOpts.push_back(std::make_pair(dataStat, "hist"));

  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz


 // Values to multiply each bin by to get them on a similar range
  double  multipliers[]={1, 1, 1, 1, 1,
                         1, 1, 1, 1, 1,
			 1, 1, 1, 1, 1};

  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  GridCanvas* gc=NULL;
  if(notRatio) gc=plotvtxrate1D(histAndOpts, doMultipliers ? multipliers : NULL,track_2);
  else gc=plotvtx1D(histAndOpts, doMultipliers ? multipliers : NULL,mcstatsys,track_2);
  // Set the y range manually. Can also use gc->Remax() to guess automatically
  //if(track_2)gc->SetYLimits(0.4, 1.75);
  if(notRatio)gc->SetYLimits(0.01,3000);
  else gc->SetYLimits(0.00001,1.4999);
  if(notRatio)  gc->SetLogy(1);

  gc->SetYTitle("Ratio to Nominal GENIE");

  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right

  double x1 = track_2? 0.7:0.75;
  double x2 = track_2? 1:1;
  double y1 = track_2? 0.3:0;
  double y2 = track_2? 0.7:0.4;

  TLegend* leg=new TLegend(x1,y1,x2,y2);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(data, "MINERvA data", "lpe");
  leg->AddEntry(mcstatsys, "MINERvA Tune", "lpe");
  leg->AddEntry(mcqe, "QE-like, QE", "lpe");
  leg->AddEntry(mcres, "QE-like, Resonant", "lpe");
  leg->AddEntry(mcdis, "QE-like, DIS", "lpe");
  leg->AddEntry(mc2p2h, "QE-like, 2p2h", "lpe");
  leg->AddEntry(mcnot, "Background", "lpe");


  
  leg->Draw("same");
  gc->Modified();

  if(!doMultipliers){
    if(notRatio){
      gc->Print(track_2 ? "nu-vtx-rate-pt-2track.eps" : "nu-vtx-rate-pt.eps");
      gc->Print(track_2 ? "nu-vtx-rate-pt-2track.png" : "nu-vtx-rate-pt.png");
      gc->Print(track_2 ? "nu-vtx-rate-pt-2track.C" : "nu-vtx-rate-pt.C");
    }
    else{
      gc->Print(track_2 ? "nu-vtx-rate-ratio-pt-2track.eps" : "nu-vtx-rate-ratio-pt.eps");
      gc->Print(track_2 ? "nu-vtx-rate-ratio-pt-2track.png" : "nu-vtx-rate-ratio-pt.png");
      gc->Print(track_2 ? "nu-vtx-rate-ratio-pt-2track.C" : "nu-vtx-rate-ratio-pt.C");
    }
  }
  else{
    if(notRatio){
      gc->Print(track_2 ? "nu-vtx-rate-pt-multipliers-2track.eps" : "nu-vtx-rate-pt-multipliers.eps");
      gc->Print(track_2 ? "nu-vtx-rate-pt-multipliers-2track.png" : "nu-vtx-rate-pt-multipliers.png");
      gc->Print(track_2 ? "nu-vtx-rate-pt-multipliers-2track.C" : "nu-vtx-rate-pt-multipliers.C");
    }
    else{
      gc->Print(track_2 ? "nu-vtx-rate-ratio-pt-multipliers-2track.eps" : "nu-vtx-rate-ratio-pt-multipliers.eps");
      gc->Print(track_2 ? "nu-vtx-rate-ratio-pt-multipliers-2track.png" : "nu-vtx-rate-ratio-pt-multipliers.png");
      gc->Print(track_2 ? "nu-vtx-rate-ratio-pt-multipliers-2track.C" : "nu-vtx-rate-ratio-pt-multipliers.C");
    }
  }



}

int main()
{
  cout << "Creating the 1 track sample" << endl;
  makePlots(false,false,false);
  makePlots(false,false,true);
  cout << "Creating the 2 track sample" << endl;
  makePlots(false,true,false);
  makePlots(false,true,true);



  return 0;
}
