#ifndef PLOT_H
#define PLOT_H

#include "TClass.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TGraphErrors.h"

#include <iostream>
#include <vector>
#include <math.h>
#include <iostream>

using namespace std;

//======================================================================

// Return a unique string on every call, so we can name root hists
// without it clobbering them
TString uniq()
{
  static int i=0;
  return TString::Format("uniq%d", i++);
}


//======================================================================

// Convert a histogram to a graph so we can draw it without ROOT
// drawing extra lines at the end and other such nonsense
TGraphErrors* histToGraph(TH1* h, double multiplier=1, bool includeZeros=true)
{
  TGraphErrors* grE=new TGraphErrors;
  grE->SetLineColor(h->GetLineColor());
  grE->SetLineStyle(h->GetLineStyle());
  grE->SetLineWidth(h->GetLineWidth());

  for(int i=0; i<h->GetNbinsX(); ++i){
    int bin=i+1;
    double binVal=h->GetBinContent(bin);
    if(includeZeros || binVal!=0){
      int binNo=grE->GetN();
      grE->SetPoint(binNo, h->GetBinCenter(bin), multiplier*binVal);
      grE->SetPointError(binNo, 0, multiplier*h->GetBinError(bin));
    }
  }
  return grE;
}

//======================================================================

TH1* getPadAxisHist(TPad* pad)
{
  TIter next(pad->GetListOfPrimitives());
  TObject *obj;
  while (( obj=next() )) {
    // cout << obj->GetName() << endl;
    if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
      // cout << "getPadAxisHist returning hist with name " << obj->GetName() << endl;
      return (TH1*)obj;
    }
  }
  return NULL;
}

//======================================================================
std::vector<TH1*> getPadHists(TPad* pad)
{
  std::vector<TH1*> ret;
  std::vector<string> names;
  TIter next(pad->GetListOfPrimitives());
  TObject *obj;
  while (( obj=next() )) {
    //    cout << obj->GetName() << endl;
    if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
      string name = obj->GetName();
      //      cout << "getPadAxisHist returning hist with name " << obj->GetName() <<"t"<< obj->GetTitle() << endl;
      bool isFound = false;
      for(unsigned int i=0;i<names.size();i++){
	if(names[i]==name){
	  isFound = true;
	  break;
	}
      }
      if(!isFound){
	ret.push_back((TH1*)obj);
	names.push_back(name);
      }

    }
  }
  return ret;
}

//======================================================================
double getPadMax(TPad* pad)
{
  TIter next(pad->GetListOfPrimitives());
  TObject *obj;
  Double_t runningMax=-9e99;//Hparam.ymax;
  while (( obj=next() )) {
    if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
      TH1* curHist=(TH1*)obj;
      const double thisMax=curHist->GetBinContent(curHist->GetMaximumBin());
      if (thisMax > runningMax) {
        runningMax=thisMax;
      }
    }
  }
  return runningMax;
}

//======================================================================
void reMax(TPad* pad, double ymin=0, double headroom=1.1)
{
  TH1* firstHist=getPadAxisHist(pad);
  if(!firstHist){
    std::cerr << "reMax: No histogram in pad " << pad->GetName() << ". Can't reMax" << std::endl;
    return;
  }
  firstHist->GetYaxis()->SetRangeUser(ymin, headroom*getPadMax(pad));
  pad->Update();
}

//======================================================================
void drawBinRange(TH2* h, int axis, int bin, const char* varName, const char* numFormatStr=".2f", bool left=false)
{
  double varmin=axis==1 ? h->GetXaxis()->GetBinLowEdge(bin) : h->GetYaxis()->GetBinLowEdge(bin);
  double varmax=axis==1 ? h->GetXaxis()->GetBinUpEdge(bin) :  h->GetYaxis()->GetBinUpEdge(bin);
  
  TString formatStr(TString::Format("%%%s < %%s < %%%s", numFormatStr, numFormatStr));

  TLatex* la=0;
  TString text(TString::Format(formatStr.Data(), varmin, varName, varmax));
  if(left){
    la=new TLatex(gPad->GetLeftMargin()+0.02,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(13); // top left
  }
  else{
    la=new TLatex(1-gPad->GetRightMargin()-0.01,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(33); // top right
  }

  la->SetNDC();
  la->SetTextFont(42);
  la->SetTextSize(0.02);
  la->Draw();
}

//======================================================================
void drawVarBinRange(TH2* h, int axis, int minbin, int maxbin, const char* varName, const char* numFormatStr=".2f", bool left=false)
{
  double varmin=axis==1 ? h->GetXaxis()->GetBinLowEdge(minbin) : h->GetYaxis()->GetBinLowEdge(minbin);
  double varmax=axis==1 ? h->GetXaxis()->GetBinUpEdge(maxbin) :  h->GetYaxis()->GetBinUpEdge(maxbin);
  
  TString formatStr(TString::Format("%%%s < %%s < %%%s", numFormatStr, numFormatStr));

  TLatex* la=0;
  TString text(TString::Format(formatStr.Data(), varmin, varName, varmax));
  if(left){
    la=new TLatex(gPad->GetLeftMargin()+0.02,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(13); // top left
  }
  else{
    la=new TLatex(1-gPad->GetRightMargin()-0.01,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(33); // top right
  }

  la->SetNDC();
  la->SetTextFont(42);
  la->SetTextSize(0.03);
  la->Draw();
}

//======================================================================
GridCanvas* plotpT1D(std::vector<std::pair<TH2*, const char*> > histAndOpts,
                     double* multipliers=NULL)
{
  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();


  int grid_x = sqrt(nbins_pz)+1;
  int grid_y = nbins_pz/(grid_x-1);
  cout << "Plotting plotpT1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pz; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionY(uniq(), i+1, i+1);
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }

        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 1, i+1, "var", ".1f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle("Muon transverse momentum (GeV)");
  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotYAxis1D(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle, double* multipliers=NULL, int xmin=0, int xmax = 0)
{
  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();
  int imin = 0, imax = nbins_pz;
  if(xmax>xmin)
  {
    nbins_pz = xmax-xmin+1;
    imin = xmin;
    imax = xmax+1;
  }


  int grid_x = sqrt(nbins_pz)+1;
  int grid_y = nbins_pz/(grid_x-1);
  if(grid_x*grid_y-nbins_pz==grid_x) grid_y-=1;
  cout << nbins_pz - grid_x*grid_y << endl;
  cout << "Plotting plotYAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=imin; i<imax; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1-imin);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionY(uniq(), i+1, i+1);
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }

        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 1, i+1, celltitle.c_str(), ".1f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle(xaxistitle.c_str());
  //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
GridCanvas* plotpT1DAntiNu(std::vector<std::pair<TH2*, const char*> > histAndOpts,
                     double* multipliers=NULL)
{
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), 4, 3, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<11; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionY(uniq(), i+1, i+1);
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }

        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 1, i+1, "p_{||}/GeV", ".1f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.0f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle("Muon transverse momentum (GeV)");
  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/nucleon)");
  gc->SetTitleSize(22.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}

//======================================================================
TH1* rebinpz(TH1* h)
{
  const int nBins=12;
  // Less aggressive version
  // double newBins[nBins+1]={1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 8, 10, 12, 14};
  double newBins[nBins+1]={1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10};

  TH1* ret=new TH1D(uniq(), TString::Format("%s;%s;%s", h->GetTitle(), h->GetXaxis()->GetTitle(), h->GetYaxis()->GetTitle()),
                    nBins, newBins);

  ret->SetLineColor(h->GetLineColor());
  ret->SetLineStyle(h->GetLineStyle());
  ret->SetLineWidth(h->GetLineWidth());
  ret->SetMarkerStyle(h->GetMarkerStyle());
  ret->SetMarkerSize(h->GetMarkerSize());
  
  // Do under/overflow too
  for(int i=0; i<nBins+2; ++i){
    ret->SetBinContent(i, h->GetBinContent(i));
    ret->SetBinError(i, h->GetBinError(i));
  }
  return ret;
}

//======================================================================
GridCanvas* plotpz1D(std::vector<std::pair<TH2*, const char*> > histAndOpts,
                     double* multipliers=NULL, bool is3D=false)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  int grid_x = sqrt(nbins_pt)+1;
  int grid_y = nbins_pt/(grid_x-1)-1;
  cout << "Plotting plotpZ1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<nbins_pt; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* projtmp=h2d->ProjectionX(uniq(), i+1, i+1);
      TH1* proj=is3D?(TH1*)projtmp->Clone("clone"):rebinpz(projtmp);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, "p_{T}/GeV", ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  const int nLabels=3;
  const double positions[nLabels]={5, 8, 10};
  const char* valueStrings[nLabels]={"5", "10", "20"};
  
  if(!is3D){
    gc->SetManualXLabels(nLabels, positions, valueStrings);
    gc->SetXTitle("Muon longitudinal momentum (GeV)");
  }
  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotXAxis1D(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle,double* multipliers=NULL,int ymin=0, int ymax=0)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();
  int imin = 0, imax = nbins_pt;
  if(ymax>ymin) 
  {
    nbins_pt = ymax-ymin+1;
    imin = ymin;
    imax = ymax+1;
  }

  int grid_x = sqrt(nbins_pt)+1;
  int grid_y = nbins_pt/(grid_x-1);
  if(grid_x*grid_y-nbins_pt==grid_x) grid_y-=1;//6*
  cout << "Plotting plotXAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=imin; i<imax; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1-imin);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotXAxis1D_Nue(std::vector<std::pair<TH2*, const char*> > histAndOpts, string xaxistitle, string celltitle,
                     double* multipliers=NULL)
{

  
  int nbins_pz = histAndOpts[0].first->GetNbinsX();
  int nbins_pt = histAndOpts[0].first->GetNbinsY();

  int grid_x = sqrt(nbins_pt)+1;
  int grid_y = nbins_pt/(grid_x-1);
  if(grid_x*grid_y-nbins_pt==grid_x) grid_y-=1;//6*
  cout << "Plotting plotXAxis1D with a grid of " << nbins_pz << "\t" << nbins_pt << "\t" << grid_x << "\t" << grid_y << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), 2, 3, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pt bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<5; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);
      proj->GetXaxis()->SetRangeUser(0,100);
      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, celltitle.c_str(), ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.1f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }


  gc->SetXTitle(xaxistitle.c_str());
    //  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/C^{12})");
  gc->SetTitleSize(20.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotpz1DAntiNu(std::vector<std::pair<TH2*, const char*> > histAndOpts,
                     double* multipliers=NULL)
{
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=new GridCanvas(uniq(), 3, 2, 800, 500);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();

  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them
  for(int i=0; i<6; ++i){
    // Change to the appropriate pad in the canvas
    TPad* pad=(TPad*)gc->cd(i+1);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();

      TH1* proj=h2d->ProjectionX(uniq(), i+1, i+1);
      //      TH1* proj=rebinpz(projtmp);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
        proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
      }

    }

    drawBinRange(histAndOpts[0].first, 2, i+1, "p_{T}/GeV", ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.0f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle("Muon longitudinal momentum (GeV)");
  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/nucleon)");
  gc->SetTitleSize(22.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotvtx1D(std::vector<std::pair<TH2*, const char*> > histAndOpts,
		      double* multipliers=NULL, TH2* denominator=NULL, bool track2=false)
{
  cout << "Doing ratios" << endl;
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=track2? new GridCanvas(uniq(), 3, 2, 4000, 2000): new GridCanvas(uniq(), 3, 3, 4000, 2000);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();
  
  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them


  const int nbins = track2 ? 4:8;  
  int minb2[4]={1,3,6,9};
  int maxb2[4]={2,5,8,13};

  int minb[8]={1,3,4,5,6,7,8,9};
  int maxb[8]={2,3,4,5,6,7,8,13};

  
  for(int i=0; i<nbins; ++i){
    // Change to the appropriate pad in the canvas
    
    int val = i>1?i+2:i+1;
    if(!track2)val=i+1;
    TPad* pad=(TPad*)gc->cd(val);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      opt.ToLower();
      TH1* projtmp = track2? h2d->ProjectionX(uniq(), minb2[i],maxb2[i]):h2d->ProjectionX(uniq(), minb[i],maxb[i]);
      TH1* demproj = track2? denominator->ProjectionX(uniq(), minb2[i],maxb2[i]):denominator->ProjectionX(uniq(), minb[i],maxb[i]);
      projtmp->Divide(demproj);

      TH1* proj=(TH1*)projtmp->Clone();
      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);

      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
	if(opt.Contains("histe")){
	  proj->SetFillColor(kWhite);
	  proj->SetFillStyle(0);
	  proj->DrawClone("histc");
	  proj->SetLineWidth(2);
	  proj->SetFillColor(kRed);
	  proj->SetFillStyle(3003);
	  proj->Draw("E5SAME");
	}
	else proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
	  
      }

    }

    if(track2) drawVarBinRange(histAndOpts[0].first, 2, minb2[i],maxb2[i], "p_{T}/GeV", ".2f", false);
    else  drawVarBinRange(histAndOpts[0].first, 2, minb[i],maxb[i], "p_{T}/GeV", ".2f", false);

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.0f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle("Vertex Energy");
  gc->SetYTitle("Events");
  gc->SetTitleSize(100.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}


//======================================================================
GridCanvas* plotvtxrate1D(std::vector<std::pair<TH2*, const char*> > histAndOpts,
			  double* multipliers=NULL,bool track2=false, bool binbybin =false)
{
  // Make a canvas of 4x3 plots with pixel size 800x500
  GridCanvas* gc=track2? new GridCanvas(uniq(), 3, 2, 4000, 2000): new GridCanvas(uniq(), 3, 3, 4000, 2000);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  gc->ResetPads();
  
  // Loop over the pz bins. For each one we'll get the 1D projections
  // in that bin and draw them


  const int nbins = track2 ? 4:8;  
  int minb2[4]={1,3,6,9};
  int maxb2[4]={2,5,8,13};

  int minb[8]={1,3,4,5,6,7,8,9};
  int maxb[8]={2,3,4,5,6,7,8,13};


  
  for(int i=0; i<nbins; ++i){
    // Change to the appropriate pad in the canvas
    
    int val = i>1?i+2:i+1;
    if(!track2)val=i+1;
    TPad* pad=(TPad*)gc->cd(val);

    for(unsigned int j=0; j<histAndOpts.size(); ++j){
      TH2* h2d=histAndOpts[j].first;
      TString opt(histAndOpts[j].second);
      //      h2d->SaveAs(Form("%s-%s-%d.root",uniq().Data(),h2d->GetName(),track2));
      opt.ToLower();
      TH1* proj= NULL;
      if(!binbybin) proj = track2? h2d->ProjectionX(uniq(), minb2[i],maxb2[i]):h2d->ProjectionX(uniq(), minb[i],maxb[i]);
      else proj = h2d->ProjectionX(uniq(),i+1,i+1);

      proj->SetLineColor(h2d->GetLineColor());
      proj->SetLineStyle(h2d->GetLineStyle());
      proj->SetLineWidth(h2d->GetLineWidth());
      proj->SetMarkerStyle(h2d->GetMarkerStyle());
      proj->SetMarkerSize(h2d->GetMarkerSize());
      proj->GetXaxis()->SetNdivisions(4);
      proj->GetYaxis()->SetNdivisions(4);
      if(multipliers) proj->Scale(multipliers[i]);

      if(opt.Contains("graph")){
        bool includeZeros=!opt.Contains("graph0");
        opt.ReplaceAll("graph0", "");
        opt.ReplaceAll("graph", "");
        TGraphErrors* gr=histToGraph(proj, 1, includeZeros);

        // If the very first histogram to be drawn is given the
        // "graph" option, draw a zeroed-out copy of the histogram
        // into the canvas so the histogram owns the axes, not the
        // graph
        if(j==0){
          TH1* axes=(TH1*)proj->Clone(uniq());
          axes->Reset();
          axes->SetLineColor(kWhite);
          axes->Draw();
        }
        gr->Draw(opt);
      }
      else{
	if(opt.Contains("histe")){
	  proj->SetFillColor(kWhite);
	  proj->SetFillStyle(0);
	  proj->DrawClone("histc");
	  proj->SetLineWidth(2);
	  proj->SetFillColor(kRed);
	  proj->SetFillStyle(3003);
	  proj->Draw("E5SAME");
	}
	else proj->Draw(j==0 ? TString(histAndOpts[j].second) : TString(histAndOpts[j].second)+" same");
	  
      }

    }
    if(!binbybin){
      if(track2) drawVarBinRange(histAndOpts[0].first, 2, minb2[i],maxb2[i], "p_{T}/GeV", ".2f", false);
      else  drawVarBinRange(histAndOpts[0].first, 2, minb[i],maxb[i], "p_{T}/GeV", ".2f", false);
    }
    else{
      drawBinRange(histAndOpts[0].first, 2, i+1, "p_{T}/GeV", ".2f", false);
    }

    // Do the same for the multiplier text, if necessary
    if(multipliers && multipliers[i]!=1){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
                             1-pad->GetTopMargin()-0.08,
                             TString::Format("#times %.0f", multipliers[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    }
  }

  gc->SetXTitle("Vertex Energy");
  gc->SetYTitle("Events");
  gc->SetTitleSize(100.0);
  // These two lines are the magic incantation to synchronize everything, put all the pads in the right place, etc
  gc->ResetPads();
  gc->Draw();

  return gc;
}





#endif
