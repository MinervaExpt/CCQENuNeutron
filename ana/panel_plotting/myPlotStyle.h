#ifndef MYPLOTSTYLE_H
#define MYPLOTSTYLE_H

#include "TStyle.h"
#include "TColor.h"

#include "TList.h"
#include "TPad.h"
#include "TH1.h"
#include "TCanvas.h"

#include <iostream>

//======================================================================
void setRedPalette()
{
  const int NRGBs = 9;
  const int NCont = 99;
  static bool initialized=false;
  static int colors[99];

  if(!initialized){
    // White -> red
    Double_t stops[NRGBs] = { 0.00, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000};
    Double_t red[NRGBs]   = { 1.00, 1.00, 0.99, 0.99, 0.98, 0.94, 0.80, 0.65, 0.40 };
    Double_t green[NRGBs] = { 0.96, 0.88, 0.73, 0.57, 0.42, 0.23, 0.09, 0.06, 0.00 };
    Double_t blue[NRGBs]  = { 0.94, 0.82, 0.63, 0.45, 0.29, 0.17, 0.11, 0.08, 0.05 };
    int colmin=TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    for(int i=0; i<NCont; ++i) colors[i]=colmin+i;

    initialized=true;
  }
  gStyle->SetNumberContours(NCont);
  gStyle->SetPalette(NCont, colors);

}

//======================================================================
void setCorrelationPalette(double whiteFrac=0.5)
{
  // A colour palette that goes blue->white->red, useful for
  // correlation matrices
  const int NRGBs = 3;
  const int NCont = 99;
  static bool initialized=false;
  static int colors[99];

  if(!initialized){
    gStyle->SetNumberContours(NCont);
    Double_t stops[NRGBs] = { 0.00, whiteFrac, 1.00};
    Double_t red[NRGBs]   = { 0.00, 1.00,      1.00};
    Double_t green[NRGBs] = { 0.00, 1.00,      0.00};
    Double_t blue[NRGBs]  = { 1.00, 1.00,      0.00};
    int colmin=TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    for(int i=0; i<NCont; ++i) colors[i]=colmin+i;

    initialized=true;
  }
  gStyle->SetNumberContours(NCont);
  gStyle->SetPalette(NCont, colors);
}

//======================================================================
void setRainbowToWhitePalette()
{
  // Matt Strait's colour palette that fades out to white at zero
  const int NRGBs = 7, NCont = 999;
  gStyle->SetNumberContours(NCont);
  Double_t stops[NRGBs] = { 0.00, 0.05, 0.23, 0.45, 0.60, 0.85, 1.00 };
  Double_t red[NRGBs]   = { 1.00, 0.00, 0.00, 0.00, 1.00, 1.00, 0.33 };
  Double_t green[NRGBs] = { 1.00, 1.00, 0.30, 0.40, 1.00, 0.00, 0.00 };
  Double_t blue[NRGBs]  = { 1.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
}

//======================================================================
void setBlackbodyPalette()
{
  // Available as gStyle->SetPalette(56) in sufficiently recent versions of ROOT, but not ours
  const int nRGBs = 5;
  const int NCont = 99;
  static bool initialized=false;
  static int colors[99];

  if(!initialized){
    gStyle->SetNumberContours(NCont);
    Double_t stops[nRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00};
    Double_t red[nRGBs]   = { 1.00, 1.00, 1.00, 0.50, 0.00};
    Double_t green[nRGBs] = { 1.00, 1.00, 0.55, 0.00, 0.00};
    Double_t blue[nRGBs]  = { 1.00, 0.00, 0.00, 0.00, 0.00};
    int colmin=TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, NCont);
    for(int i=0; i<NCont; ++i) colors[i]=colmin+i;

    initialized=true;
  }
  gStyle->SetNumberContours(NCont);
  gStyle->SetPalette(NCont, colors);
}

//======================================================================
void myPlotStyle()
{
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.15);

  // use large fonts

  gStyle->SetLabelFont(42);
  gStyle->SetTitleFont(42);

  gStyle->SetTextSize(0.08);

  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.06,"z");
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleX(0.25);
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleOffset(1.2, "Y");

  // use bold lines and markers
  gStyle->SetMarkerStyle(1);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars and y error bar caps
  //gStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  gStyle->SetTickLength(0.01, "Y");
  gStyle->SetTickLength(0.02, "X");

  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetStripDecimals(false);
}

#endif
