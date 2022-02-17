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

void makePlots(bool doMultipliers,bool track_2, bool doGENIE,bool doRPA)
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
  MnvH2D* m1=(MnvH2D*)f->Get("GENIE");
  MnvH2D* m2=(MnvH2D*)f->Get("GENIE_RPA");
  MnvH2D* m3=(MnvH2D*)f->Get("GENIE_RPA_2p2h");
  MnvH2D* m4=(MnvH2D*)f->Get("GENIE_2p2h");
  MnvH2D* m5=(MnvH2D*)f->Get("GENIE_PionTune");
  MnvH2D* m6=(MnvH2D*)f->Get("GENIE_RPA_2p2h_recoil_NievesRPARes");
  MnvH2D* m7=(MnvH2D*)f->Get("GENIE_RPA_2p2h_recoil_MINOSRPARes");

  
  MnvH2D* mnn=(MnvH2D*)lowrec->GetVertErrorBand("Low_recoil_fits")->GetHist(0)->Clone("modelnn");
  MnvH2D* mnp=(MnvH2D*)lowrec->GetVertErrorBand("Low_recoil_fits")->GetHist(1)->Clone("modelnp");
  MnvH2D* mqe=(MnvH2D*)lowrec->GetVertErrorBand("Low_recoil_fits")->GetHist(2)->Clone("modelqe");
  lowrec->PopVertErrorBand("Low_recoil_fits");//remove the low recoil error from the MC error since that is something we'd like to measure??
  
  dataMnv->GetXaxis()->SetTitle("Vertex Energy (MeV)");

  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  cout << "Getting the TH2" << endl;

  TH2* data=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* datasys=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mcstatsys = new TH2D(lowrec->GetCVHistoWithError()); // This is the full tune MC CV with model errors (excludes signal and low recoil fit mods

  TH2* model1 = new TH2D(m1->GetCVHistoWithStatError());
  TH2* model2 = new TH2D(m2->GetCVHistoWithStatError());
  TH2* model3 = new TH2D(m3->GetCVHistoWithStatError());
  TH2* model4 = new TH2D(m4->GetCVHistoWithStatError());
  TH2* model5 = new TH2D(m5->GetCVHistoWithStatError());
  TH2* model6 = new TH2D(m6->GetCVHistoWithStatError());
  TH2* model7 = new TH2D(m7->GetCVHistoWithStatError());
  TH2* modelnn = new TH2D(mnn->GetCVHistoWithStatError());
  TH2* modelnp = new TH2D(mnp->GetCVHistoWithStatError());
  TH2* modelqe = new TH2D(mqe->GetCVHistoWithStatError());

  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = getColors(1);

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


  model1->SetLineColor(kBlue);
  model1->SetLineWidth(4);
  model1->SetMarkerColor(kBlue);


  model2->SetLineColor(mycolors[3]);
  model2->SetLineWidth(4);
  model2->SetMarkerColor(mycolors[3]);


  model3->SetLineColor(mycolors[6]);
  model3->SetLineWidth(4);
  model3->SetMarkerColor(mycolors[6]);


  model4->SetLineColor(mycolors[8]);
  model4->SetLineWidth(4);
  model4->SetMarkerColor(mycolors[8]);


  model5->SetLineColor(mycolors[4]);
  model5->SetLineWidth(4);
  model5->SetMarkerColor(mycolors[4]);


  model6->SetLineColor(mycolors[7]);
  model6->SetLineWidth(4);
  model6->SetMarkerColor(mycolors[7]);


  model7->SetLineColor(mycolors[16]);
  model7->SetLineWidth(4);
  model7->SetMarkerColor(mycolors[16]);


  modelnn->SetLineColor(mycolors[16]);
  modelnn->SetLineWidth(4);
  modelnn->SetLineStyle(0);
  modelnn->SetMarkerColor(mycolors[16]);


  modelnp->SetLineColor(mycolors[10]);
  modelnp->SetLineWidth(4);
  modelnp->SetLineStyle(0);
  modelnp->SetMarkerColor(mycolors[10]);


  modelqe->SetLineColor(mycolors[12]);
  modelqe->SetLineWidth(4);
  modelqe->SetLineStyle(0);
  modelqe->SetMarkerColor(mycolors[12]);






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
  if(doGENIE){
    histAndOpts.push_back(std::make_pair(model2,     "histc"));
    histAndOpts.push_back(std::make_pair(model3,     "histc"));
    histAndOpts.push_back(std::make_pair(model4,     "histc"));
    histAndOpts.push_back(std::make_pair(model5,     "histc"));
  }
  else if(doRPA){
    histAndOpts.push_back(std::make_pair(model2,     "histc"));
    histAndOpts.push_back(std::make_pair(model3,     "histc"));
    histAndOpts.push_back(std::make_pair(model6,     "histc"));
    histAndOpts.push_back(std::make_pair(model7,     "histc"));

  }
  else{
    histAndOpts.push_back(std::make_pair(modelnn,     "histc"));
    histAndOpts.push_back(std::make_pair(modelnp,     "histc"));
    histAndOpts.push_back(std::make_pair(modelqe,     "histc"));

  }
  histAndOpts.push_back(std::make_pair(model1,     "histc"));
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
  GridCanvas* gc=plotvtx1D(histAndOpts, doMultipliers ? multipliers : NULL,model1,track_2);
  // Set the y range manually. Can also use gc->Remax() to guess automatically
  if(track_2)gc->SetYLimits(0.4, 1.75);
  else gc->SetYLimits(0.4, 2.21);
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
  leg->AddEntry(model1, "GENIE 2.8.4","lpe");
  if(doGENIE){
    leg->AddEntry(model2, "GENIE+RPA","lpe");
    leg->AddEntry(model3, "GENIE+RPA+2p2h","lpe");
    leg->AddEntry(model4, "GENIE+2p2h","lpe");
    leg->AddEntry(model5, "GENIE+PionTune","lpe");
  }
  else if(doRPA){
    leg->AddEntry(model2, "GENIE+RPA","lpe");
    leg->AddEntry(model3, "GENIE+RPA+2p2h","lpe");
    leg->AddEntry(model6, "MINERvA Tune+NievesRPAOnPion","lpe");
    leg->AddEntry(model7, "MINERvA Tune+MINOSRPAOnPion","lpe");
  }
  else{
    leg->AddEntry(modelnn, "MINERvA Tune nn-only","lpe");
    leg->AddEntry(modelnp, "MINERvA Tune np-only","lpe");
    leg->AddEntry(modelqe, "MINERvA Tune QE-only","lpe");
  }
  
  leg->Draw("same");
  gc->Modified();

  if(doMultipliers){
    if(doGENIE){
      gc->Print(track_2 ? "nu-vtx-ratio-pt-2track.eps" : "nu-vtx-ratio-pt.eps");
      gc->Print(track_2 ? "nu-vtx-ratio-pt-2track.png" : "nu-vtx-ratio-pt.png");
      gc->Print(track_2 ? "nu-vtx-ratio-pt-2track.C" : "nu-vtx-ratio-pt.C");
    }
    else if(doRPA){
      gc->Print(track_2 ? "nu-vtx-ratio-pt-pionrpa-2track.eps" : "nu-vtx-ratio-pt-pionrpa.eps");
      gc->Print(track_2 ? "nu-vtx-ratio-pt-pionrpa-2track.png" : "nu-vtx-ratio-pt-pionrpa.png");
      gc->Print(track_2 ? "nu-vtx-ratio-pt-pionrpa-2track.C" : "nu-vtx-ratio-pt-pionrpa.C");

    }
    else{
      gc->Print(track_2 ? "nu-vtx-ratio-pt-minervatune-2track.eps" : "nu-vtx-ratio-pt-minervatune.eps");
      gc->Print(track_2 ? "nu-vtx-ratio-pt-minervatune-2track.png" : "nu-vtx-ratio-pt-minervatune.png");
      gc->Print(track_2 ? "nu-vtx-ratio-pt-minervatune-2track.C" : "nu-vtx-ratio-pt-minervatune.C");
    }
  }
  else{
    if(doGENIE){
      gc->Print(track_2 ? "nu-vtx-ratio-pt-multipliers-2track.eps" : "nu-vtx-ratio-pt-multipliers.eps");
      gc->Print(track_2 ? "nu-vtx-ratio-pt-multipliers-2track.png" : "nu-vtx-ratio-pt-multipliers.png");
      gc->Print(track_2 ? "nu-vtx-ratio-pt-multipliers-2track.C" : "nu-vtx-ratio-pt-multipliers.C");
    }
    else if(doRPA){
      gc->Print(track_2 ? "nu-vtx-ratio-pt-pionrpa-multipliers-2track.eps" : "nu-vtx-ratio-pt-pionrpa-multipliers.eps");
      gc->Print(track_2 ? "nu-vtx-ratio-pt-pionrpa-multipliers-2track.png" : "nu-vtx-ratio-pt-pionrpa-multipliers.png");
      gc->Print(track_2 ? "nu-vtx-ratio-pt-pionrpa-multipliers-2track.C" : "nu-vtx-ratio-pt-pionrpa-multipliers.C");

    }
    else{
      gc->Print(track_2 ? "nu-vtx-ratio-pt-minervatune-multipliers-2track.eps" : "nu-vtx-ratio-pt-minervatune-multipliers.eps");
      gc->Print(track_2 ? "nu-vtx-ratio-pt-minervatune-multipliers-2track.png" : "nu-vtx-ratio-pt-minervatune-multipliers.png");
      gc->Print(track_2 ? "nu-vtx-ratio-pt-minervatune-multipliers-2track.C" : "nu-vtx-ratio-pt-minervatune-multipliers.C");
    }
  }



}

int main()
{
  //  makePlots(true);
  cout << "Creating the 1 track sample" << endl;
  makePlots(false,false,true,false);
  makePlots(false,false,false,false);
  makePlots(false,false,false,true);

  //  makePlots(true,false);
  cout << "Creating the 2 track sample" << endl;
  makePlots(false,true,true,false);
  makePlots(false,true,false,false);
  makePlots(false,true,false,true);
  //  makePlots(true,true);
  return 0;
}
