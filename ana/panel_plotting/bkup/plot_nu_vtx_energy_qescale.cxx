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

void makePlots(bool doMultipliers,bool track_2)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  gStyle->SetNdivisions(555, "XYZ");


  TFile *f = NULL;
  if(track_2) f = new TFile("/minerva/app/users/drut1186/cmtuser/Minerva_v10r8p9_CCQENu2DLE_hadronv2/Ana/CCQENu2DLE/ana/plot_macros_pub/Vertex/Vertex_Constrained_2DPlots_Mult_2.root");
  else f = new TFile("/minerva/app/users/drut1186/cmtuser/Minerva_v10r8p9_CCQENu2DLE_hadronv2/Ana/CCQENu2DLE/ana/plot_macros_pub/Vertex/Vertex_Constrained_2DPlots_Mult_1.root");

  MnvH2D* dataMnv=(MnvH2D*)f->Get("h_vtx_new_data"); // data
  MnvH2D* lowrec = (MnvH2D*)f->Get("lowrec"); // CV MC also where I get the low recoil models 
  MnvH2D* mnn=(MnvH2D*)lowrec->GetVertErrorBand("Low_recoil_fits")->GetHist(0)->Clone("modelnn");
  MnvH2D* mnp=(MnvH2D*)lowrec->GetVertErrorBand("Low_recoil_fits")->GetHist(1)->Clone("modelnp");
  MnvH2D* mqe=(MnvH2D*)lowrec->GetVertErrorBand("Low_recoil_fits")->GetHist(2)->Clone("modelqe");
  MnvH2D* mqe90=(MnvH2D*)lowrec->GetVertErrorBand("Low_recoil_fits")->GetHist(3)->Clone("modelqe90");
  MnvH2D* mqe80=(MnvH2D*)lowrec->GetVertErrorBand("Low_recoil_fits")->GetHist(4)->Clone("modelqe80");
  MnvH2D* mqe70=(MnvH2D*)lowrec->GetVertErrorBand("Low_recoil_fits")->GetHist(5)->Clone("modelqe70");
  MnvH2D* mqe60=(MnvH2D*)lowrec->GetVertErrorBand("Low_recoil_fits")->GetHist(6)->Clone("modelqe60");
  MnvH2D* mqe50=(MnvH2D*)lowrec->GetVertErrorBand("Low_recoil_fits")->GetHist(7)->Clone("modelqe50");
  MnvH2D* mqe40=(MnvH2D*)lowrec->GetVertErrorBand("Low_recoil_fits")->GetHist(8)->Clone("modelqe40");
  MnvH2D* mqe30=(MnvH2D*)lowrec->GetVertErrorBand("Low_recoil_fits")->GetHist(9)->Clone("modelqe30");
  MnvH2D* mqe20=(MnvH2D*)lowrec->GetVertErrorBand("Low_recoil_fits")->GetHist(10)->Clone("modelqe20");
  MnvH2D* mqe10=(MnvH2D*)lowrec->GetVertErrorBand("Low_recoil_fits")->GetHist(11)->Clone("modelqe10");


  lowrec->PopVertErrorBand("Low_recoil_fits");//remove the low recoil error from the MC error since that is something we'd like to measure??
  
  dataMnv->GetXaxis()->SetTitle("Vertex Energy (MeV)");


  TH2* basePlus90 = new TH2D(mqe90->GetCVHistoWithStatError());
  TH2* basePlus80 = new TH2D(mqe80->GetCVHistoWithStatError());
  TH2* basePlus70 = new TH2D(mqe70->GetCVHistoWithStatError());
  TH2* basePlus60 = new TH2D(mqe60->GetCVHistoWithStatError());
  TH2* basePlus50 = new TH2D(mqe50->GetCVHistoWithStatError());
  TH2* basePlus40 = new TH2D(mqe40->GetCVHistoWithStatError());
  TH2* basePlus30 = new TH2D(mqe30->GetCVHistoWithStatError());
  TH2* basePlus20 = new TH2D(mqe20->GetCVHistoWithStatError());
  TH2* basePlus10 = new TH2D(mqe10->GetCVHistoWithStatError());








  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  cout << "Getting the TH2" << endl;

  TH2* data=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* datasys=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mcstatsys = new TH2D(lowrec->GetCVHistoWithError()); // This is the full tune MC CV with model errors (excludes signal and low recoil fit mods

  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = getColors(0);

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
  mcstatsys->SetFillColor(kRed);
  mcstatsys->SetFillStyle(3444);

  //  basePlus100->SetLineColor(mycolors[2]);
  basePlus90->SetLineColor(mycolors[3]);
  basePlus80->SetLineColor(mycolors[4]);
  basePlus70->SetLineColor(mycolors[5]);
  basePlus60->SetLineColor(mycolors[6]);
  basePlus50->SetLineColor(mycolors[7]);
  basePlus40->SetLineColor(mycolors[8]);
  basePlus30->SetLineColor(mycolors[9]);
  basePlus20->SetLineColor(mycolors[10]);
  basePlus10->SetLineColor(mycolors[11]);
  //  basePlus0->SetLineColor(mycolors[12]);

  //  basePlus100->SetLineWidth(3);
  basePlus90->SetLineWidth(3);
  basePlus80->SetLineWidth(3);
  basePlus70->SetLineWidth(3);
  basePlus60->SetLineWidth(3);
  basePlus50->SetLineWidth(5);
  basePlus40->SetLineWidth(3);
  basePlus30->SetLineWidth(3);
  basePlus20->SetLineWidth(3);
  basePlus10->SetLineWidth(5);
  //  basePlus0->SetLineWidth(3);


  //  basePlus100->SetLineStyle(1);
  basePlus90->SetLineStyle(1);
  basePlus80->SetLineStyle(1);
  basePlus70->SetLineStyle(1);
  basePlus60->SetLineStyle(1);
  basePlus50->SetLineStyle(1);
  basePlus40->SetLineStyle(1);
  basePlus30->SetLineStyle(1);
  basePlus20->SetLineStyle(1);
  basePlus10->SetLineStyle(1);
  //  basePlus0->SetLineStyle(1);

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
  //  histAndOpts.push_back(std::make_pair(basePlus100,     "histl"));
  histAndOpts.push_back(std::make_pair(basePlus90,     "histl"));
  histAndOpts.push_back(std::make_pair(basePlus80,     "histl"));
  histAndOpts.push_back(std::make_pair(basePlus70,     "histl"));
  histAndOpts.push_back(std::make_pair(basePlus60,     "histl"));
  histAndOpts.push_back(std::make_pair(basePlus50,     "histl"));
  histAndOpts.push_back(std::make_pair(basePlus40,     "histl"));
  histAndOpts.push_back(std::make_pair(basePlus30,     "histl"));
  histAndOpts.push_back(std::make_pair(basePlus20,     "histl"));
  histAndOpts.push_back(std::make_pair(basePlus10,     "histl"));
  //  histAndOpts.push_back(std::make_pair(basePlus0,     "histl"));
  //  histAndOpts.push_back(std::make_pair(m3tmp,     "hist"));
  //  histAndOpts.push_back(std::make_pair(m1tmp,     "histl"));
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
  GridCanvas* gc=plotvtx1D(histAndOpts, false ,mcstatsys,track_2);
  // Set the y range manually. Can also use gc->Remax() to guess automatically
  if(track_2)gc->SetYLimits(0.5, 1.5);
  else gc->SetYLimits(0.75, 1.2);
  gc->SetYTitle("Ratio to MnvGENIE");

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
  leg->AddEntry(basePlus90, "90% QE Variant", "lpe");
  leg->AddEntry(basePlus80, "80% QE Variant", "lpe");
  leg->AddEntry(basePlus70, "70% QE Variant", "lpe");
  leg->AddEntry(basePlus60, "60% QE Variant", "lpe");
  leg->AddEntry(basePlus50, "50% QE Variant", "lpe");
  leg->AddEntry(basePlus40, "40% QE Variant", "lpe");
  leg->AddEntry(basePlus30, "30% QE Variant", "lpe");
  leg->AddEntry(basePlus20, "20% QE Variant", "lpe");
  leg->AddEntry(basePlus10, "10% QE Variant", "lpe");

  
  leg->Draw("same");
  gc->Modified();

  gc->Print(track_2 ? "nu-vtx-ratio-pt-2track.eps" : "nu-vtx-ratio-pt.eps");
  gc->Print(track_2 ? "nu-vtx-ratio-pt-2track.png" : "nu-vtx-ratio-pt.png");
  gc->Print(track_2 ? "nu-vtx-ratio-pt-2track.C" : "nu-vtx-ratio-pt.C");
  



}

int main()
{
  //  makePlots(true);
  cout << "Creating the 1 track sample" << endl;
  makePlots(false,false);


  //  makePlots(true,false);
  cout << "Creating the 2 track sample" << endl;
  makePlots(false,true);
  //  makePlots(true,true);
  return 0;
}
