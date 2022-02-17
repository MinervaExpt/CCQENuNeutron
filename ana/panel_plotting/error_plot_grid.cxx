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
TH2* concatenateHists(vector<TH1*>& hists1D, int axis)
{
  assert(hists1D.size());

  // cout << "concatenateHists with " << hists1D.size() << " hists, axis=" << axis << endl;
  const int nyBins=13;
  const double yBins[nyBins+1]={0, 0.07, 0.15, 0.25, 0.33, 0.4, 0.47, 0.55, 0.7, 0.85, 1, 1.25, 1.5, 2.5};

  const int nxBins=12;
  const double xBins[nxBins+1]={1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 8, 10, 15, 20};

  TH2* ret=0;
  if(axis==1){
    ret=new TH2D(uniq(), TString::Format(";%s", hists1D[0]->GetXaxis()->GetTitle()),
                 hists1D[0]->GetXaxis()->GetNbins(), hists1D[0]->GetXaxis()->GetXbins()->GetArray(),
                 nyBins, yBins);
  }
  else{
    ret=new TH2D(uniq(), TString::Format(";;%s", hists1D[0]->GetXaxis()->GetTitle()),
                 nxBins, xBins,
                 hists1D[0]->GetXaxis()->GetNbins(), hists1D[0]->GetXaxis()->GetXbins()->GetArray());
  }

  ret->SetLineColor(hists1D[0]->GetLineColor());
  ret->SetLineStyle(hists1D[0]->GetLineStyle());

  for(unsigned int iHist=0; iHist<hists1D.size(); ++iHist){
    for(int j=0; j<hists1D[0]->GetXaxis()->GetNbins()+1; ++j){
      int ixBin=axis==1 ? j       : iHist+1;
      int iyBin=axis==1 ? iHist+1 : j;
      double content=hists1D[iHist]->GetBinContent(j);
      ret->SetBinContent(ixBin, iyBin, content);
    }
  }

  return ret;
}

// =====================================================================
vector<std::pair<TH2*, const char*> > getSystHistsAndOpts(MnvH2D* data, bool pt)
{
  MnvPlotter plotter;
  plotter.ApplyStyle(kCCQENuStyle);
  
  // For each bin in the other variable, make a vector of the
  // systematic histograms
  const int nBins=pt ? 12 : 13;
  vector<vector<TH1*> > histsPT;
  histsPT.resize(nBins);

  // Get MnvPlotter to plot all the histograms, and slurp them into histsPT
  for(int i=0; i<nBins; ++i){
    // First plot the histograms in the dummy canvas...
    TCanvas c;
    MnvH1D* proj=pt ? data->ProjectionY(uniq(), i+1, i+1) : data->ProjectionX(uniq(), i+1, i+1);
    plotter.DrawErrorSummary(proj, "N", true, true, -1, false, "", true);
    std::vector<TH1*> padHists=getPadHists(&c);
    histsPT[i]=padHists;
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
  vector<vector<TH1*> > histsPT_transpose;
  int nSyst=histsPT[0].size();
  cout << "There are " << nSyst << " systematics" << endl;
  histsPT_transpose.resize(nSyst);

  for(int iSyst=0; iSyst<nSyst; ++iSyst){
    for(unsigned int iBin=0; iBin<histsPT.size(); ++iBin){  
      histsPT_transpose[iSyst].push_back(histsPT[iBin][iSyst]);
    }
  }

  vector<std::pair<TH2*, const char*> > histsPT2D;
  // TODO: Figure out why the last systematic is crashing
  for(int iSyst=0; iSyst<histsPT_transpose.size(); ++iSyst){
    TH2* h2d=concatenateHists(histsPT_transpose[iSyst], pt ? 2 : 1);
    // We want to draw all of these histograms as graphs, and exclude
    // the zero bins, to get rid of ROOT artifacts. The "graph0" draw
    // option does that (and I made it safe to pass all graphs)
    histsPT2D.push_back(std::make_pair(h2d, "graph0 l"));
  }

  return histsPT2D;
}

// =====================================================================
void makePlots(bool pt, string location)
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

  TFile f(Form("%s_CV/CrossSection_per_nucleon_iterations_4_SpecialSampleIncluded.root_pzmu.root",location.c_str()));
  MnvH2D* data=(MnvH2D*)f.Get("h_pzmu_ptmu_data_nobck_unfold_effcor_cross_section");

  MnvH2D* zerobin = (MnvH2D*)data->Clone("zero");
  zerobin->Reset();
  zerobin->ClearAllErrorBands();
  for(int i=0;i<zerobin->GetNbinsX();i++){
    for(int j=0;j<zerobin->GetNbinsY();j++){
      if(i==5 && j==12){
	zerobin->SetBinContent(i+1,j+1,0.0);
	zerobin->SetBinError(i+1,j+1,0.0);
      }
      else{
	zerobin->SetBinContent(i+1,j+1,1.0);
	zerobin->SetBinError(i+1,j+1,0.0);
      }
      
    }
  }
  zerobin->AddMissingErrorBandsAndFillWithCV(*data);
  zerobin->SaveAs("test.root");
  data->Multiply(data,zerobin);
  vector<std::pair<TH2*, const char*> > histsPT2D=getSystHistsAndOpts(data,  pt);


  vector<string> Sysnames;
  Sysnames.push_back("Total Uncertainty");
  Sysnames.push_back("Statisical");
  Sysnames.push_back("Low Recoil Fits");
  Sysnames.push_back("FSI Model");
  Sysnames.push_back("Flux");
  Sysnames.push_back("Muon Reconstruction");
  Sysnames.push_back("Others");
  Sysnames.push_back("Cross Section Models");

  

  TLegend *leg = new TLegend(0.17, 0.7, 0.31, 0.9);
  leg->SetLineColor(kWhite);
  leg->SetFillColor(kWhite);
  for(int i=0;i<histsPT2D.size();i++)leg->AddEntry(histsPT2D[i].first,Sysnames[i].c_str(),"l");
  if(pt){
    GridCanvas* gcPT=plotpT1D(histsPT2D);
    gcPT->SetYTitle("Fractional uncertainty");
    gcPT->SetYLimits(0, 0.49);
    leg->Draw("SAME");
    gcPT->Modified();
    gcPT->Print("errors-pt.eps");
    gcPT->Print("errors-pt.png");
    gcPT->Print("errors-pt.C");
  }
  else{
    GridCanvas* gcPT=plotpz1D(histsPT2D);
    gcPT->SetYTitle("Fractional uncertainty");
    gcPT->SetYLimits(0, 0.49);
    gcPT->Modified();
    gcPT->Print("errors-pz.eps");
    gcPT->Print("errors-pz.png");
    gcPT->Print("errors-pz.C");
  }
}

int main(int argc, char* argv[])
{
  makePlots(true,argv[1]);
  makePlots(false,argv[1]);
  return 0;
}
