// GridCanvas.cxx -
//
// Philip Rodrigues
// Thursday, February  7 2013
//

#include "GridCanvas.h"
#include "TLatex.h"
#include "TPad.h"
#include "TH1.h"
#include "TAxis.h"
#include "TList.h"
#include "TGraph.h"
#include "TStyle.h"

#include "plot.h"

#include <iostream>

GridCanvas::GridCanvas(const char* name, int nPadsX, int nPadsY, int ww, int wh)
  : TCanvas(name, "title", ww, wh), fNPadsX(nPadsX), fNPadsY(nPadsY),
    fInterpadSpace(0.001),
    fXTitleLatex(new TLatex), fYTitleLatex(new TLatex),
    fTitleAlignment(kAlignCenter),
    fXTitle(""), fYTitle(""),
    fXTitleDrawn(false), fYTitleDrawn(false),
    fTitleFont(-1), fTitleSize(-1),
    fManualXLabels(false)
{
  fPads.resize(fNPadsX*fNPadsY);
  fPads2D.resize(fNPadsX);

  //const double padWidth=(1-fLeftMargin-fRightMargin)/fNPadsX;
  //const double padHeight=(1-fTopMargin-fBottomMargin)/fNPadsY;

  for(int i=0; i<fNPadsX; ++i){
    fPads2D[i].resize(fNPadsY);
    for(int j=0; j<fNPadsY; ++j){
      int counter=fNPadsX*(fNPadsY-1-j)+i;

      TPad *pad = new TPad(TString::Format("pad%d", counter), "foo",
                           0, 0, 1, 1);
      pad->SetNumber(counter+1);
      fPads[counter]=pad;
      fPads2D[i][j]=pad;

      TCanvas::cd();
      pad->Draw();
    }
  }
  ResetPads();
}

GridCanvas::~GridCanvas()
{
  delete fXTitleLatex;
  delete fYTitleLatex;
  for(unsigned int i=0; i<fPads.size(); ++i) delete fPads[i];
}

void GridCanvas::ResetPads()
{
  const double gridWidth=(1-fLeftMargin-fRightMargin);
  const double gridHeight=(1-fTopMargin-fBottomMargin);

  const double frameWidth=gridWidth/fNPadsX;
  const double frameHeight=gridHeight/fNPadsY;
  const double aspectRatio=frameWidth/frameHeight;

  for(int i=0; i<fNPadsX; ++i){
    fPads2D[i].resize(fNPadsY);
    for(int j=0; j<fNPadsY; ++j){
      int counter=fNPadsX*(fNPadsY-1-j)+i;

      // const double thisPadHeight = j==0 ? padHeight+fBottomMargin  : padHeight;
      // const double thisPadWidth  = i==0 ? padWidth+fLeftMargin : padWidth;

      TPad *pad = fPads[counter];

      pad->SetFillStyle(4000);

      const double left=fLeftMargin + i*frameWidth + fInterpadSpace/aspectRatio;
      const double right=fRightMargin + (fNPadsX-i-1)*frameWidth + fInterpadSpace/aspectRatio;
      const double bottom=fBottomMargin + j*frameHeight + fInterpadSpace;
      const double top=fTopMargin + (fNPadsY-j-1)*frameHeight + fInterpadSpace;

      pad->SetLeftMargin(left);
      pad->SetRightMargin(right);
      pad->SetBottomMargin(bottom);
      pad->SetTopMargin(top);

      //printf("(%d, %d) c=%d l=%.2f r=%.2f b=%.2f t=%.2f\n", i, j, counter, left, right, bottom, top);
      // pad->SetBottomMargin(j==0 ? (thisPadHeight-padHeight)/thisPadHeight : fInterpadSpace);
      // pad->SetTopMargin(fInterpadSpace);
      // pad->SetLeftMargin(i==0 ? (thisPadWidth-padWidth)/thisPadWidth : fInterpadSpace);
      // pad->SetRightMargin(fInterpadSpace);
    }
  }
}

TH1* GridCanvas::GetPadHist(TPad* pad)
{
  TList* prims=pad->GetListOfPrimitives();
  for(int i=0; i<prims->GetEntries(); ++i){
    TObject* obj=prims->At(i);
    // std::cout << "Object name=" << obj->GetName() << " classname=" << obj->ClassName() << std::endl;
    TH1* h=dynamic_cast<TH1*>(obj);
    if(h) return h;
    TGraph* gr=dynamic_cast<TGraph*>(obj);
    if(gr) return gr->GetHistogram();
  }
  return 0;
}

void GridCanvas::Paint(Option_t* option)
{
  bool anyPadModified=false;
  for(unsigned int i=0; i<fPads.size(); ++i) anyPadModified=anyPadModified || fPads[i]->IsModified();
  if(anyPadModified || IsModified()){
    SetHistTexts();
  }
  TCanvas::Paint(option);
}

void GridCanvas::SetHistTexts()
{
  for(int i=0; i<fNPadsX; ++i){
    fPads2D[i].resize(fNPadsY);
    for(int j=0; j<fNPadsY; ++j){
      int counter=fNPadsX*(fNPadsY-1-j)+i;

      // const double thisPadHeight = j==0 ? padHeight+fBottomMargin  : padHeight;
      // const double thisPadWidth  = i==0 ? padWidth+fLeftMargin : padWidth;

      TPad *pad = fPads[counter];

      TH1* hist=GetPadHist(pad);
      if(!hist) continue;

      hist->GetXaxis()->SetTitleSize(0);
      hist->GetYaxis()->SetTitleSize(0);

      if(i!=0) hist->GetYaxis()->SetLabelSize(0);
      if(j!=0 || fManualXLabels) hist->GetXaxis()->SetLabelSize(0);
    }
  }
  if(fXTitle=="" && GetPadHist(fPads[0])) SetXTitle(GetPadHist(fPads[0])->GetXaxis()->GetTitle());
  if(fYTitle=="" && GetPadHist(fPads[0])) SetYTitle(GetPadHist(fPads[0])->GetYaxis()->GetTitle());

  DrawTitles();
}

void GridCanvas::SetManualXLabels(int nLabels, const double* positions, const char** valueStrings,
                                  double yoffset)
{
  fManualXLabels=true;
  for(int i=0; i<fNPadsX; ++i){
    cout << "PAD " << i << endl;
    fPads2D[i].resize(fNPadsY);
    int j=0;
    int counter=fNPadsX*(fNPadsY-1-j)+i;
    
    TPad *pad = fPads[counter];
    pad->cd();
    pad->Update();

    TH1* hist=GetPadHist(pad);
    if(hist==NULL) continue;
    double x1=pad->GetUxmin();
    double x2=pad->GetUxmax();
    double y1=pad->GetUymin();
    double y2=pad->GetUymax();

    double lmarg=pad->GetLeftMargin();
    double rmarg=pad->GetRightMargin();
    double tmarg=pad->GetTopMargin();
    double bmarg=pad->GetBottomMargin();

    // cout << "y2=" << y2 << " y1=" << y1 << " ypos=" << (y1-yoffset*(y2-y1)) << endl;

    for(int i=0; i<nLabels; ++i){
      // We have to place the TLatex in NDC so it doesn't move when
      // the user changes the y axis limits. This'll still break if
      // the user changes the x limits, but you can't have
      // everything...
      // 
      // Presumably the real fix is to paint() the latex in
      // GridCanvas::Paint() instead of doing this, but laziness

      double y=y1-yoffset*(y2-y1);

      double xndc= lmarg + (positions[i]-x1)*(1-rmarg-lmarg)/(x2-x1);
      double yndc= bmarg + (y-y1)*(tmarg-bmarg)/(y2-y1);

      TLatex* la=new TLatex(xndc, yndc, /* positions[i], y1-yoffset*(y2-y1), */
                            valueStrings[i]);
      la->SetNDC();
      la->SetTextAlign(23);
      
      la->SetTextFont(hist->GetXaxis()->GetLabelFont());
      // Can't do this because we set the label size to zero so ROOT's
      // labels aren't shown
      // la->SetTextSize(axis->GetLabelSize());
      la->SetTextSize(gStyle->GetLabelSize());
      la->Draw();
    }
    cout << "DONE" << endl;
  }
}

void GridCanvas::SetXTitle(const char* xtitle)
{
  fXTitle=xtitle;
}

void GridCanvas::SetYTitle(const char* ytitle)
{
  fYTitle=ytitle;
}

TLatex* GridCanvas::GetXTitle()
{
  return fXTitleLatex;
}

TLatex* GridCanvas::GetYTitle()
{
  return fYTitleLatex;
}

void GridCanvas::SetXLimits(double xmin, double xmax)
{
  for(unsigned int i=0; i<fPads.size(); ++i){
    TH1* h=GetPadHist(fPads[i]);
    if(!h) continue;
    h->GetXaxis()->SetRangeUser(xmin, xmax);
  }
}

void GridCanvas::SetYLimits(double ymin, double ymax)
{
  for(unsigned int i=0; i<fPads.size(); ++i){
    TH1* h=GetPadHist(fPads[i]);
    if(!h) continue;
    h->GetYaxis()->SetRangeUser(ymin, ymax);
  }
}

void GridCanvas::DrawTitles()
{
  fXTitleLatex->SetTitle(fXTitle);
  fYTitleLatex->SetTitle(fYTitle);

  if(GetPadHist(fPads[0])){
    // TH1* h=GetPadHist(fPads[0]);
    // fXTitleLatex->SetTextFont(h->GetXaxis()->GetTitleFont());
    // fYTitleLatex->SetTextFont(h->GetYaxis()->GetTitleFont());

    // fXTitleLatex->SetTextSize(h->GetXaxis()->GetTitleSize());
    // fYTitleLatex->SetTextSize(h->GetYaxis()->GetTitleSize());

    fXTitleLatex->SetTextFont(fTitleFont==-1 ? 43 : fTitleFont );
    fYTitleLatex->SetTextFont(fTitleFont==-1 ? 43 : fTitleFont);

    fXTitleLatex->SetTextSize(fTitleSize==-1 ? 30 : fTitleSize);
    fYTitleLatex->SetTextSize(fTitleSize==-1 ? 30 : fTitleSize);
  }

  fYTitleLatex->SetTextAngle(90);

  double xposx, xposy, yposx, yposy;
  switch(fTitleAlignment){
  case kAlignRight:
    xposx=1-fRightMargin;
    yposx=0.02;

    xposy=0.02;
    yposy=1-fTopMargin;

    fXTitleLatex->SetTextAlign(31);
    fYTitleLatex->SetTextAlign(33);
    break;
  case kAlignCenter:
    xposx=fLeftMargin+0.5*(1-fLeftMargin-fRightMargin);
    yposx=0.02;

    xposy=0.01;
    yposy=fBottomMargin+0.5*(1-fLeftMargin-fRightMargin);

    fXTitleLatex->SetTextAlign(21);
    fYTitleLatex->SetTextAlign(23);
    break;
  default:
    std::cerr << "Unknown alignment type in GridCanvas::DrawTitles()" << std::endl;
    return;
  }
  fXTitleLatex->SetX(xposx);   fXTitleLatex->SetY(yposx);
  fYTitleLatex->SetX(xposy);   fYTitleLatex->SetY(yposy);

  fXTitleLatex->SetNDC();
  fYTitleLatex->SetNDC();

  if(!fXTitleDrawn){
    std::cout << "Drawing x title" << std::endl;
    TCanvas::cd();
    fXTitleLatex->Draw();
    fXTitleDrawn=true;
  }

  if(!fYTitleDrawn){
    std::cout << "Drawing y title" << std::endl;
    TCanvas::cd();
    fYTitleLatex->Draw();
    fYTitleDrawn=true;
  }
}

void GridCanvas::SetLogx(Int_t value)
{
  for(unsigned int i=0; i<fPads.size(); ++i) fPads[i]->SetLogx(value);
}

void GridCanvas::SetLogy(Int_t value)
{
  for(unsigned int i=0; i<fPads.size(); ++i) fPads[i]->SetLogy(value);
}

void GridCanvas::SetLogz(Int_t value)
{
  for(unsigned int i=0; i<fPads.size(); ++i) fPads[i]->SetLogz(value);
}

void GridCanvas::SetGridx(Int_t value)
{
  for(unsigned int i=0; i<fPads.size(); ++i) fPads[i]->SetGridx(value);
}

void GridCanvas::SetGridy(Int_t value)
{
  for(unsigned int i=0; i<fPads.size(); ++i) fPads[i]->SetGridy(value);
}

void GridCanvas::SetTicksy(const char* option)
{
  for(unsigned int i=0; i<fPads.size(); ++i){
    TH1* h=GetPadHist(fPads[i]);
    h->GetYaxis()->SetTicks(option);
    if(strcmp(option, "+")==0){
      // TODO: Work out what the hell this number means, why it's so
      // stupid, and how to make it leave the numbers unchanged
      h->GetYaxis()->SetLabelOffset(-0.015);
    }
  }
}

void GridCanvas::SetTicksx(const char* option)
{
  for(unsigned int i=0; i<fPads.size(); ++i){
    TH1* h=GetPadHist(fPads[i]);
    h->GetXaxis()->SetTicks(option);
    if(strcmp(option, "+")==0){
      // TODO: Work out what the hell this number means, why it's so
      // stupid, and how to make it leave the numbers unchanged
      h->GetXaxis()->SetLabelOffset(-0.015);
    }
  }
}

void GridCanvas::Remax(double ymin)
{
  double allMax=-9e99;
  for(unsigned int i=0; i<fPads.size(); ++i) allMax=std::max(allMax, getPadMax(fPads[i]));
  SetYLimits(ymin, 1.1*allMax);
}

void GridCanvas::SetLeftMargin(Float_t margin)
{
  TCanvas::SetLeftMargin(margin);
  ResetPads();
}

void GridCanvas::SetRightMargin(Float_t margin)
{
  TCanvas::SetRightMargin(margin);
  ResetPads();
}

void GridCanvas::SetTopMargin(Float_t margin)
{
  TCanvas::SetTopMargin(margin);
  ResetPads();
}

void GridCanvas::SetBottomMargin(Float_t margin)
{
  TCanvas::SetBottomMargin(margin);
  ResetPads();
}


ClassImp(GridCanvas)
