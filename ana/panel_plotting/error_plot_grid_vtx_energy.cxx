//  #include "../util/plot/myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvPlotter.h"

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

#include <iostream>

using namespace std;
using namespace PlotUtils;

// =====================================================================
// axis==1: x axis. axis==2: y axis

// Take a vector of 1D histograms which all have the same binning and
// stack them together into a 2D histogram. The axis argument says
// which axis to stack them on
TH2* concatenateHists(vector<TH1*> hists1D, int axis = 1, bool track2=false)
{
  //  assert(hists1D.size());

  const int nbins = track2 ? 4:8;  
  double track2bins[] = {0,0.15,0.4,0.7,2.5};
  double track1bins[] = {0,0.15,0.25,0.33,0.4,0.47,0.55,0.7,2.5};

  cout << "concatenateHists with " << hists1D.size() << " hists, axis=" << axis << " 2track sample " << track2 <<endl;
  cout << hists1D[0] ->GetName()<< endl;
  TH2* ret=NULL;
  if(track2)ret = new TH2D(uniq(), TString::Format(";;%s", hists1D[0]->GetXaxis()->GetTitle()),hists1D[0]->GetXaxis()->GetNbins(), hists1D[0]->GetXaxis()->GetXbins()->GetArray(), nbins,track2bins);
  else     ret = new TH2D(uniq(), TString::Format(";;%s", hists1D[0]->GetXaxis()->GetTitle()),hists1D[0]->GetXaxis()->GetNbins(), hists1D[0]->GetXaxis()->GetXbins()->GetArray(), nbins,track1bins);
  

  cout << "Return has " << ret->GetYaxis()->GetNbins() << endl;

  ret->SetLineColor(hists1D[0]->GetLineColor());
  ret->SetLineStyle(hists1D[0]->GetLineStyle());
  ret->SetLineWidth(3);
  for(unsigned int iHist=0; iHist<hists1D.size(); ++iHist){
    cout << "Analyzing hist " << hists1D[iHist]->GetName() << endl;
    for(int j=0; j<hists1D[0]->GetXaxis()->GetNbins()+2; ++j){
      int ixBin= j; // Want this to be vtx energy
      int iyBin= iHist+1; // This is the nth bin
      double content=hists1D[iHist]->GetBinContent(j);
      ret->SetBinContent(ixBin, iyBin, content);
    }
  }

  return ret;
}

// =====================================================================
vector<std::pair<TH2*, const char*> > getSystHistsAndOpts(MnvH2D* data, bool track2)
{
  MnvPlotter plotter;
  plotter.ApplyStyle(kCCQENuStyle);

  const int nBins = track2 ? 4:8;  
  int minb2[4]={1,3,6,9};
  int maxb2[4]={2,5,8,13};

  int minb[8]={1,3,4,5,6,7,8,9};
  int maxb[8]={2,3,4,5,6,7,8,13};
  
  // For each bin in the other variable, make a vector of the
  // systematic histograms
  vector<vector<TH1*> > hists;
  hists.resize(nBins);

  cout << "We have " << nBins << " bins to process" << endl;
  // Get MnvPlotter to plot all the histograms, and slurp them into histsPT
  for(int i=0; i<nBins; ++i){
    // First plot the histograms in the dummy canvas...
    TCanvas c;
    track2 ? cout << i << "\t" << minb2[i] << "\t" << maxb2[i] << endl : cout << i << "\t" << minb[i] << "\t" << maxb[i] << endl;
    MnvH1D* proj=track2 ? data->ProjectionX(uniq(), minb2[i], maxb2[i]) : data->ProjectionX(uniq(), minb[i], maxb[i]);
    plotter.DrawErrorSummary(proj, "N", "FSI Models", true, -1, false, "", true);
    //c.SaveAs(Form("Projattempt_%d.C",i)); // Enable this if you want to debug the plots!!
    std::vector<TH1*> padHists=getPadHists(&c);
    hists[i]=padHists;
  }
  


  // concatenateHists wants a vector of hists for each of the bins of
  // a given systematic. But histsPT is the other way round (the inner
  // vector loops over systematics).  So we have this fiddly loop to
  // make a transposed version of the vector-of-vector

  // It would have been easier to just pass the original
  // vector-of-vector into concatenateHists, and tell it which
  // systematic we wanted, but I've written and debugged this now, so
  // not changing it

  //  First index is systematic, second
  // index is bin
  vector<vector<TH1*> > hists_transpose;
  int nSyst=hists[0].size();
  cout << "There are " << nSyst << " systematics" << endl;
  hists_transpose.resize(nSyst); 

  //each hist[x] is a set of systematics for projected pt bin N
  //hists[x][syst]
  //x = number of pt bins (4 or 8)
  //syst = 9
  for(int iSyst=0; iSyst<nSyst; ++iSyst){//loop o 8
    for(unsigned int iBin=0; iBin<hists.size(); ++iBin){  //loob of? 4?
      hists_transpose[iSyst].push_back(hists[iBin][iSyst]);
    }
  }

  vector<std::pair<TH2*, const char*> > hists2D;
  cout << "I will plot "<<hists_transpose.size() << " histograms per bin "<< endl;
  for(int iSyst=0; iSyst<hists_transpose.size(); ++iSyst){
    cout << "On systematic " << iSyst << endl;
    TH2* h2d=concatenateHists(hists_transpose[iSyst],1, track2);
    // We want to draw all of these histograms as graphs, and exclude
    // the zero bins, to get rid of ROOT artifacts. The "graph0" draw
    // option does that (and I made it safe to pass all graphs)
    hists2D.push_back(std::make_pair(h2d, "hist"));
  }

  cout << "I've created a set of histograms to plot" << hists2D.size() << endl;
  return hists2D;
}

// =====================================================================
void makePlots(bool track2)
{
  // This turns out to be complicated (well, I could have made it less
  // bad, but this way is complicated and generalizable):
  //
  // MnvPlotter knows how to make the histograms we need, including
  // fancy stuff like grouping systematics together and sticking to
  // colour schemes. But it doesn't know about GridCanvas, and it
  // goes off and draws the histograms in its own way
  //
  // So we're going to take projections, and ask MnvPlotter to plot
  // the projection in a dummy canvas. Then we'll grab all the
  // histograms from that canvas and hold onto them.
  //
  // We could just grab all those 1D histograms and plot them straight
  // into the appropriate panel of our GridCanvas, but I want to reuse
  // the plotpz1D and plotpT1D functions in plot.h, because they know
  // how to do fanciness like squashing the tail in pz. Those
  // functions take a vector of 2D histograms, so we have to take our
  // 1D histograms from MnvPlotter, and stack them back together into
  // 2D histograms. That's what the concatenateHists function does

  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  TFile *f = NULL;
  if(track2) f = new TFile("/minerva/app/users/drut1186/cmtuser/Minerva_v10r8p9_CCQENu2DLE_hadronv2/Ana/CCQENu2DLE/ana/plot_macros_pub/Vertex/Vertex_Constrained_2DPlots_Mult_2.root");
  else f = new TFile("/minerva/app/users/drut1186/cmtuser/Minerva_v10r8p9_CCQENu2DLE_hadronv2/Ana/CCQENu2DLE/ana/plot_macros_pub/Vertex/Vertex_Constrained_2DPlots_Mult_1.root");
  MnvH2D* data=(MnvH2D*)f->Get("lowrec");

  vector<std::pair<TH2*, const char*> > hists2D=getSystHistsAndOpts(data,  track2);
  cout <<   "Number of hists" << hists2D.size() << endl;;
  if(track2){
    GridCanvas* gcPT=plotvtxrate1D(hists2D, NULL, true,true);
    double x1 = track2? 0.7:0.75;
    double x2 = track2? 1:1;
    double y1 = track2? 0.3:0;
    double y2 = track2? 0.7:0.4;
    TLegend* leg=new TLegend(x1,y1,x2,y2);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(hists2D[0].first, "Total Uncertainty", "l");
    leg->AddEntry(hists2D[1].first, "Statistical", "l");
    leg->AddEntry(hists2D[2].first, "Low Recoil Fits", "l");
    leg->AddEntry(hists2D[3].first, "FSI", "l");
    leg->AddEntry(hists2D[4].first, "Flux", "l");
    leg->AddEntry(hists2D[5].first, "Muon Reconstruction", "l");
    leg->AddEntry(hists2D[6].first, "Other", "l");
    leg->AddEntry(hists2D[7].first, "Model", "l");

    gcPT->SetYTitle("Fractional uncertainty");
    gcPT->SetYLimits(0, 0.49);
    gcPT->Modified();
    leg->Draw("SAME");
    gcPT->Print("errors-2track-vtx.eps");
    gcPT->Print("errors-2track-vtx.png");
    gcPT->Print("errors-2track-vtx.C");
  }
  else{
    GridCanvas* gcPT=plotvtxrate1D(hists2D, NULL, false,true);
    double x1 = track2? 0.7:0.75;
    double x2 = track2? 1:1;
    double y1 = track2? 0.3:0;
    double y2 = track2? 0.7:0.4;
    TLegend* leg=new TLegend(x1,y1,x2,y2);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(hists2D[0].first, "Total Uncertainty", "l");
    leg->AddEntry(hists2D[1].first, "Statistical", "l");
    leg->AddEntry(hists2D[2].first, "Low Recoil Fit", "l");
    leg->AddEntry(hists2D[3].first, "FSI Models", "l");
    leg->AddEntry(hists2D[4].first, "Flux", "l");
    leg->AddEntry(hists2D[5].first, "Muon Reconstruction", "l");
    leg->AddEntry(hists2D[6].first, "Others", "l");
    leg->AddEntry(hists2D[7].first, "Interaction Models", "l");

    gcPT->SetYTitle("Fractional uncertainty");
    gcPT->SetYLimits(0, 0.49);
    gcPT->Modified();
    leg->Draw("SAME");
    gcPT->Print("errors-vtx.eps");
    gcPT->Print("errors-vtx.png");
    gcPT->Print("errors-vtx.C");
  }
}

int main()
{
  makePlots(true);
  makePlots(false);
  return 0;
}
