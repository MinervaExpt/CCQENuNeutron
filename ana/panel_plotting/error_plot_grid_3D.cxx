//  #include "../util/plot/myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvPlotter.h"
#include "PlotUtils/HyperDimLinearizer.h"//THIS HAS TO CHANGE TO BE INCLUDED IN THE MAKE FILE EVENTUALLY.
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

    
  // hdptbins.push_back(0.0);
  // hdptbins.push_back(0.075);//added ME
  // hdptbins.push_back(0.15);
  // hdptbins.push_back(0.25);//added ME
  // hdptbins.push_back(0.4);
  // hdptbins.push_back(0.7);//added ME
  // hdptbins.push_back(1.0);
  // hdptbins.push_back(1.5);//added ME
  // hdptbins.push_back(2.5);

  // hdpzbins.push_back(1.5);
  // hdpzbins.push_back(3.5);//added ME
  // hdpzbins.push_back(5.0);
  // hdpzbins.push_back(10.0);//added ME
  // hdpzbins.push_back(20.0);

  // hdrecoilbins_short.push_back(0.);
  // hdrecoilbins_short.push_back(1500.);

  // hdrecoilbins.push_back(0.);
  // hdrecoilbins.push_back(60.);
  // hdrecoilbins.push_back(140.);
  // hdrecoilbins.push_back(240.);
  // hdrecoilbins.push_back(480.);
  // hdrecoilbins.push_back(1000.);


  // cout << "concatenateHists with " << hists1D.size() << " hists, axis=" << axis << endl;
  const int nyBins=6;
  const double yBins[nyBins+1]={0,0.15,0.25,0.4,0.7,1.0,2.5};

  const int nxBins=13;
  const double xBins[nxBins+1]={0,40,80,120,160,200,240,280,320,360,400,600,800,1000};

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
vector<std::pair<TH2*, const char*> > getSystHistsAndOpts(MnvH2D* data, bool pt, string group = "")
{
  MnvPlotter plotter;
  plotter.ApplyStyle(kCCQENuStyle);

  vector<string> vertnames = data->GetVertErrorBandNames();
  vector<string> latnames = data->GetLatErrorBandNames();
  // For each bin in the other variable, make a vector of the
  // systematic histograms
  const int nBins=pt ? 13 : 8;
  vector<vector<TH1*> > histsPT;
  histsPT.resize(nBins);

  // Get MnvPlotter to plot all the histograms, and slurp them into histsPT
  for(int i=0; i<nBins; ++i){
    // First plot the histograms in the dummy canvas...
    TCanvas c;
    MnvH1D* proj=pt ? data->ProjectionY(uniq(), i+1, i+1) : data->ProjectionX(uniq(), i+1, i+1);
    if(group=="") plotter.DrawErrorSummary(proj, "N", true, true, -1, false,"", true);
    else{
      if(std::count(vertnames.begin(),vertnames.end(),group)){
	TH1D err = proj->GetVertErrorBand(group)->GetErrorBand(true);
	err.DrawClone();
      }
      else{
	TH1D err = proj->GetLatErrorBand(group)->GetErrorBand(true);
	err.DrawClone();
      }
    }
    
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
int makePlots(bool pt, string location, bool drawGroups, string histoname)
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

  TFile f(Form("%s_CV/CrossSection_per_nucleon_3D_pzptreco_iterations_4_CombinedPlaylists.root_big3d.root",location.c_str()));
  MnvH2D* dataMnv=(MnvH2D*)f.Get(histoname.c_str());
  //  MnvH2D* dataMnv=(MnvH2D*)f.Get("h_pzptrec_data_nobck_unfold_effcor_cross_section");
  //  MnvH2D* dataMnv=(MnvH2D*)f.Get("h_pzptrec_mc_nobck_unfold_effcor");
  //MnvH2D* dataMnv=(MnvH2D*)f.Get("h_pzptrec_mc_nobck_unfold");


  if(!dataMnv){
    cout << "Failed to get the histogram" << endl;
    return 1;
  }
  cout << dataMnv << endl;

  /*
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
  */

  vector<double> recoil3Dbins;
  vector<double> pt3Dbins;
  vector<double> pz3Dbins;

  pt3Dbins.push_back(0.0);
  pt3Dbins.push_back(0.15);
  pt3Dbins.push_back(0.25);//added ME
  pt3Dbins.push_back(0.4);
  pt3Dbins.push_back(0.7);//added ME
  pt3Dbins.push_back(1.0);
  pt3Dbins.push_back(2.5);

  pz3Dbins.push_back(1.5);
  pz3Dbins.push_back(3.5);//added ME
  pz3Dbins.push_back(8.0);
  pz3Dbins.push_back(20.0);

  for(int i=0;i<10;i++)recoil3Dbins.push_back(i*40);
  for(int i=0;i<4;i++)recoil3Dbins.push_back(i*200+400);




  std::vector<std::vector<double> > full3D;
  full3D.push_back(recoil3Dbins);
  full3D.push_back(pt3Dbins);
  full3D.push_back(pz3Dbins);
  
  HyperDimLinearizer *my3d = new HyperDimLinearizer(full3D,0);

  std::cout << "Starting up getting the projections" << std::endl;
  std::vector<MnvH2D*> dataresults = my3d->Get2DMnvHistos(dataMnv,true);
  
  for(int i=1;i<dataresults.size()-1;i++){
    MnvH2D *data = dataresults[i];
    vector<string> names = data->GetErrorBandNames();
    if(!drawGroups){
      vector<std::pair<TH2*, const char*> > histsPT2D=getSystHistsAndOpts(data,  pt, "");   
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
      for(int j=0;j<histsPT2D.size();j++)leg->AddEntry(histsPT2D[j].first,Sysnames[j].c_str(),"l");
      if(pt){
	GridCanvas* gcPT=plotpT1D(histsPT2D);
	gcPT->SetYTitle("Fractional uncertainty");
	gcPT->SetYLimits(0, 0.49);
	leg->Draw("SAME");
	gcPT->Modified();
	gcPT->Print(Form("errors-pt_%d.eps",i));
	gcPT->Print(Form("errors-pt_%d.png",i));
	gcPT->Print(Form("errors-pt_%d.C",i));
      }
      else{
	GridCanvas* gcPT=plotpz1D(histsPT2D,NULL,true);
	gcPT->SetYTitle("Fractional uncertainty");
	gcPT->SetYLimits(0, 0.49);
	gcPT->Modified();
	gcPT->Print(Form("errors-pz_%d.eps",i));
	gcPT->Print(Form("errors-pz_%d.png",i));
	gcPT->Print(Form("errors-pz_%d.C",i));
      }
    }
    else{
      for(int n=0;n<names.size();n++){
	string group = names[n];
	
	vector<std::pair<TH2*, const char*> > histsPT2D=getSystHistsAndOpts(data,  pt, group);   
	TLegend *leg = new TLegend(0.17, 0.7, 0.31, 0.9);
	leg->SetLineColor(kWhite);
	leg->SetFillColor(kWhite);
      	for(int j=0;j<histsPT2D.size();j++)leg->AddEntry(histsPT2D[j].first,group.c_str(),"l");
	if(pt){
	  GridCanvas* gcPT=plotpT1D(histsPT2D);
	  gcPT->SetYTitle("Fractional uncertainty");
	  gcPT->SetYLimits(0, 0.25);
	  leg->Draw("SAME");
	  gcPT->Modified();
	  gcPT->Print(Form("errors-%s-pt_%d.eps",group.c_str(),i));
	  gcPT->Print(Form("errors-%s-pt_%d.png",group.c_str(),i));
	  gcPT->Print(Form("errors-%s-pt_%d.C",group.c_str(),i));
	}
	else{
	  GridCanvas* gcPT=plotpz1D(histsPT2D,NULL,true);
	  gcPT->SetYTitle("Fractional uncertainty");
	  gcPT->SetYLimits(0, 0.25);
	  gcPT->Modified();
	  gcPT->Print(Form("errors-%s-pz_%d.eps",group.c_str(),i));
	  gcPT->Print(Form("errors-%s-pz_%d.png",group.c_str(),i));
	  gcPT->Print(Form("errors-%s-pz_%d.C",group.c_str(),i));
	}
      }
    }
  }
}
int main(int argc, char* argv[])
{
  makePlots(true,argv[1],false, argv[2]);
  makePlots(false,argv[1],false, argv[2]);

  makePlots(true,argv[1],true, argv[2]);
  makePlots(false,argv[1],true, argv[2]);
  return 0;
}
