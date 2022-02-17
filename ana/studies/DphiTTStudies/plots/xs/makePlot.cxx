#include <iostream>

#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLatex.h"
//using namespace std;
//using namespace ROOT;
#include <vector>
#include "TGraph.h"

void format( TGraph *g1, TGraph *g2, int color )
{
  g1->SetLineStyle(2);
  g1->SetLineWidth(2);

  g2->SetLineStyle(1);
  g2->SetLineWidth(2);

  g1->SetLineColor( color );
  g2->SetLineColor( color+2 );

  g1->SetMarkerColor( color );
  g2->SetMarkerColor( color+2 );

  g1->GetXaxis()->SetRangeUser(0,100);
  g2->GetXaxis()->SetRangeUser(0,100);
}

vector<double> modify( TGraph* ela, TGraph *inela, TGraph *total, double min, double max, double inelaFrac, double x )
{
  double vtotal = total->Eval(x);
  double vinela = inela->Eval(x)*inelaFrac;
  double vela = 0;
  if(vinela > vtotal) vinela = vtotal;
  else vela = vtotal - vinela;
  vector<double> ret(3,0);
  ret[0] = vela;
  ret[1] = vinela;
  ret[2] = vtotal;
  return ret;

}

vector<TGraph*> modify( TGraph* ela, TGraph *inela, TGraph *total, double min, double max, double inelaFrac )
{
  vector<TGraph*>ret;
  ret.push_back( new TGraph(*ela) );
  ret.push_back( new TGraph(*inela )); 
  ret.push_back( new TGraph(*total) );
  int nPoints = ela->GetN();
  for( int i = 0 ; i< nPoints; i++ )
  {
    double x,y;
    ela->GetPoint(i, x, y);
    if(x<min) continue;
    if(x>max) break;
    vector<double> eit = modify( ela, inela, total, min, max, inelaFrac, x );
    ret[0]->SetPoint(i, x, eit[0]);
    ret[1]->SetPoint(i, x, eit[1]);
    ret[2]->SetPoint(i, x, eit[2]);
  }
  return ret;
}


void makePlots(double min, double max, double rw)
{
  TFile *fnew = new TFile("new_geant4_xsec.root");
  TFile *fold = new TFile("cross_section_neutron_carbon_hd.root");

  TGraph *new_elastic = (TGraph*) fnew->Get("elastic_xsec");
  TGraph *new_inelastic = (TGraph*) fnew->Get("inelastic_xsec");
  TGraph *new_total = (TGraph*) fnew->Get("total_xsec");

  vector<TGraph*> modified = modify( new_elastic, new_inelastic, new_total, min,max, rw );

  format( new_elastic, modified[0], kBlue );
  format( new_inelastic, modified[1],kRed );
  format( new_total, modified[2], kBlack );

  TCanvas *c = new TCanvas("c1", "c1");
  new_total->Draw("");
  new_elastic->Draw("same");
  modified[0]->Draw("same");
  new_inelastic->Draw("same");
  modified[1]->Draw("same");


  TLegend *l = new TLegend(0.5,0.6,0.85,0.91);
  l->AddEntry( new_total, "Total Neutron XS" );
  l->AddEntry( new_elastic, "Elastic XS" );
  l->AddEntry( new_inelastic, "Inelastic XS" );
  l->AddEntry( modified[0], Form( "Modified Elastic XS" ));
  l->AddEntry( modified[1], Form( "Modified Inelastic XS"));
  l->Draw();

  TLatex *latex = new TLatex();
  latex->SetTextSize( latex->GetTextSize()*0.75 );
  latex->SetTextAlign(22);
  latex->DrawLatexNDC(0.5,0.95, Form("Inelastic Modified between (%.0f, %.0f) by %.1f",min,max,rw) );
  c->Print(Form("plot_%.0f_%.0f_%.2f.pdf",min,max,rw) );

  fnew->Close();
  fold->Close();
}

void makeOldNewPlot()
{
  TFile *fnew = new TFile("new_geant4_xsec.root");
  TFile *fold = new TFile("cross_section_neutron_carbon_hd.root");

  TGraph *new_elastic = (TGraph*) fnew->Get("elastic_xsec");
  TGraph *new_inelastic = (TGraph*) fnew->Get("inelastic_xsec");
  TGraph *new_total = (TGraph*) fnew->Get("total_xsec");

  TGraph *old_elastic = (TGraph*) fold->Get("elXS");
  TGraph *old_inelastic = (TGraph*) fold->Get("inelXS");
  TGraph *old_total = (TGraph*) fold->Get("totalXS");
  TGraphErrors *abfalterer = (TGraphErrors*) fold->Get("abfalterer");
  abfalterer->SetMarkerColor(kBlue-10);
  abfalterer->SetLineColor(kBlue-10);




  format( old_elastic, new_elastic, kBlue );
  format( old_inelastic, new_inelastic,kRed );
  format( old_total, new_total, kBlack );

  TCanvas *c = new TCanvas("c1", "c1");
  new_total->Draw("");
  old_total->Draw("same");
  new_elastic->Draw("same");
  old_elastic->Draw("same");
  new_inelastic->Draw("same");
  old_inelastic->Draw("same");
  abfalterer->Draw("same");


  TLegend *l = new TLegend(0.5,0.6,0.85,0.91);
  l->AddEntry( new_total, "Total XS - 4.10.3.p03 FTFP_BERT" );
  l->AddEntry( new_elastic, "Elastic XS - 4.10.3.p03 FTFP_BERT" );
  l->AddEntry( new_inelastic, "Inelastic XS - 4.10.3.p03 FTFP_BERT" );
  l->AddEntry( old_total, "Total XS" );
  l->AddEntry( old_elastic, Form( "Elastic XS" ));
  l->AddEntry( old_inelastic, Form( "Inelastic XS"));
  l->AddEntry( abfalterer, Form( "Abfalterer"));
  l->Draw();

  TLatex *latex = new TLatex();
  latex->SetTextSize( latex->GetTextSize()*0.75 );
  latex->SetTextAlign(22);
  latex->DrawLatexNDC(0.5,0.95, "Comparing old and new XS" );
  c->Print(Form("plot_comparison.pdf") );


}

void makeCarbonHydrogen()
{
  TFile *fnew = new TFile("new_geant4_xsec.root");
  TFile *fold = new TFile("cross_section_neutron_hydrogen_hd.root");

  TGraph *new_elastic = (TGraph*) fnew->Get("elastic_xsec");
  TGraph *new_inelastic = (TGraph*) fnew->Get("inelastic_xsec");
  TGraph *new_total = (TGraph*) fnew->Get("total_xsec");

  TGraph *old_elastic = (TGraph*) fold->Get("elXS");
  TGraph *old_inelastic = (TGraph*) fold->Get("inelXS");
  TGraph *old_total = (TGraph*) fold->Get("totalXS");
  TGraphErrors *abfalterer = (TGraphErrors*) fold->Get("abfalterer");
  abfalterer->SetMarkerColor(kBlue-10);
  abfalterer->SetLineColor(kBlue-10);




  format( old_elastic, new_elastic, kBlue );
  format( old_inelastic, new_inelastic,kRed );
  format( old_total, new_total, kBlack );

  TCanvas *c = new TCanvas("c1", "c1");
  new_total->Draw("");
  old_total->Draw("same");
  new_elastic->Draw("same");
  old_elastic->Draw("same");
  new_inelastic->Draw("same");
  old_inelastic->Draw("same");
  abfalterer->Draw("same");


  TLegend *l = new TLegend(0.5,0.6,0.85,0.91);
  l->AddEntry( new_total, "Total C XS - 4.10.3.p03 FTFP_BERT" );
  l->AddEntry( new_elastic, "Elastic C XS - 4.10.3.p03 FTFP_BERT" );
  l->AddEntry( new_inelastic, "Inelastic C XS - 4.10.3.p03 FTFP_BERT" );
  l->AddEntry( old_total, "Total H XS" );
  l->AddEntry( old_elastic, Form( "Elastic H XS" ));
  l->AddEntry( old_inelastic, Form( "Inelastic H XS"));
  l->AddEntry( abfalterer, Form( "Abfalterer Hydrogen"));
  l->Draw();

  TLatex *latex = new TLatex();
  latex->SetTextSize( latex->GetTextSize()*0.75 );
  latex->SetTextAlign(22);
  latex->DrawLatexNDC(0.5,0.95, "Comparing old and new XS" );
  c->Print(Form("plot_comparison_CH.pdf") );


}

void makePlot()
{
  const int nRW = 4, nRanges = 4;
  double reweights[nRW] = {0.1,0.3,1.5,2};
  double mins[nRanges] = {0,0,40,70};
  double maxs[nRanges] = {20,40,70,100};
  for( int i = 0; i< nRanges; i++ )
  {
    for( int j = 0; j<nRW; j++ ) makePlots(mins[i], maxs[i], reweights[j] );
  }
  makeOldNewPlot();
  makeCarbonHydrogen();
}
