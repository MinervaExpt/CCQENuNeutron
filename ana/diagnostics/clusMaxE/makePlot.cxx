#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"


void FormatHisto( TH1D* h, int color )
{
  h->SetLineWidth(2);
  h->SetLineColor(color);
}

void FindMaxPurityCut( TH1D* htotal, TH1D* hpart, TH1D*&purity, TH1D*& efficiency )
{
  purity = (TH1D*) hpart->Clone("purity");
  purity->Reset();
  purity->SetTitle("purity");

  efficiency= (TH1D*) hpart->Clone("efficiency");
  efficiency->Reset();
  efficiency->SetTitle("efficiency");

  int max_bin = purity->GetNbinsX()+1;
  double total_total = htotal->Integral(1,max_bin);
  double part_total = hpart->Integral(1,max_bin);
  for( int i = 0; i< max_bin; i++ )
  {
    int bin = i+1;

    double part_sum = hpart->Integral(bin, max_bin );
    double total_sum = htotal->Integral(bin, max_bin );
    double pur = part_sum/total_sum;
    double eff = part_sum/part_total;
    purity->SetBinContent(bin, pur );
    efficiency->SetBinContent( bin, eff );
  }
}
int makePlot()
{
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  TChain* CCQENu= new TChain("CCQENu");
  CCQENu->Add("/minerva/data/users/tejinc/CCQENu/test/Merged/Merged*root");


  TH1D* hproton = new TH1D("hproton","hproton",25,0,100 );
  TH1D* hpi0 = (TH1D*) hproton->Clone("hpi0");
  TH1D* hoth = (TH1D*) hproton->Clone("hoth");
  //CCQENu->Draw("CCQENu_BlobClusterMaxE/(pow((CCQENu_BlobBegX-CCQENu_BlobEndX)^2+(CCQENu_BlobBegY-CCQENu_BlobEndY)^2+(CCQENu_BlobBegZ-CCQENu_BlobEndZ)^2,.5)/abs(CCQENu_BlobBegZ-CCQENu_BlobEndZ))>>hproton","multiplicity==1 && recoil_energy_nonmuon_nonvtx100mm < 500 && CCQENu_nuHelicity == 2 
// Max$(CCQENu_BlobTotalE)==CCQENu_BlobTotalE && \
  //
  CCQENu->Draw("CCQENu_BlobClusterMaxE>>hproton",
      "multiplicity==1 && CCQENu_nuHelicity == 2 \
      && phys_vertex_is_fiducial && CCQENu_BlobIs3D==1 && \
      CCQENu_BlobMCPID==2212 && CCQENu_RecoPattern==3 && CCQENu_BlobTotalE>0 && \
      Max$(CCQENu_BlobTotalE)==CCQENu_BlobTotalE && \
      ( (recoil_energy_nonmuon_nonvtx100mm<150 ) && \
         recoil_energy_nonmuon_nonvtx100mm/1000.<0.05+0.3*CCQENu_Q2/pow(10,6) && \
         recoil_energy_nonmuon_nonvtx100mm < 500 )" );


  CCQENu->Draw("CCQENu_BlobClusterMaxE>>hoth",
      "multiplicity==1 && CCQENu_nuHelicity == 2 \
      && phys_vertex_is_fiducial && CCQENu_BlobIs3D==1 && \
      CCQENu_BlobMCPID!=111 && CCQENu_BlobMCPID!=2212 && CCQENu_RecoPattern==3 && CCQENu_BlobTotalE>0 && \
      Max$(CCQENu_BlobTotalE)==CCQENu_BlobTotalE && \
      ( (recoil_energy_nonmuon_nonvtx100mm<150 ) && \
         recoil_energy_nonmuon_nonvtx100mm/1000.<0.05+0.3*CCQENu_Q2/pow(10,6) && \
         recoil_energy_nonmuon_nonvtx100mm < 500 )" );

  CCQENu->Draw("CCQENu_BlobClusterMaxE>>hpi0",
      "multiplicity==1  && CCQENu_nuHelicity == 2 \
      && phys_vertex_is_fiducial && CCQENu_BlobIs3D==1 && \
      CCQENu_BlobMCPID==111 && CCQENu_RecoPattern==3 && CCQENu_BlobTotalE>0 && \
      Max$(CCQENu_BlobTotalE)==CCQENu_BlobTotalE && \
      ( (recoil_energy_nonmuon_nonvtx100mm<150 ) && \
         recoil_energy_nonmuon_nonvtx100mm/1000.<0.05+0.3*CCQENu_Q2/pow(10,6) && \
         recoil_energy_nonmuon_nonvtx100mm < 500 )" );


  cout<<"hproton:"<<hproton->Integral()<<endl;
  cout<<"hpi0:"<<hpi0->Integral()<<endl;
  cout<<"hoth:"<<hoth->Integral()<<endl;


  TH1D* htotal = (TH1D*) hproton->Clone("htotal");
  htotal->Add(hpi0);
  htotal->Add(hoth);

  FormatHisto(htotal, kBlack );
  FormatHisto(hproton, kRed );
  FormatHisto(hpi0, kBlue );
  FormatHisto(hoth, kGray );


  TCanvas *c1 = new TCanvas("c1","c1" );
  c1->SetGrid();
  c1->cd();
  htotal->SetTitle("Leading Cluster Energy");
  htotal->GetXaxis()->SetTitle("Energy (MeV)" );

  htotal->Draw("hist");
  hpi0->Draw("histsame");
  hproton->Draw("histsame");
  hoth->Draw("histsame");
  TLegend*leg = new TLegend(.7,.7,.9,.9);
  leg->AddEntry(htotal, "all blobs");
  leg->AddEntry(hproton, "proton blobs");
  leg->AddEntry(hpi0, "#pi^{0} blobs");
  leg->AddEntry(hoth, "other blobs");
  leg->Draw();
  c1->Print("Blob_ClusMaxE.pdf");




  TH1D* purity,*efficiency;
  FindMaxPurityCut( htotal, hproton, purity, efficiency );
  cout<<"found"<<endl;
  efficiency->Draw();

  TH1D* mult = (TH1D*)purity->Clone("mult");
  mult->Multiply(efficiency);
  mult->GetXaxis()->SetTitle( "Energy (MeV)");
  mult->SetTitle( "Efficiency/Purity against Leading Cluster Energy");

  FormatHisto( purity, kRed );
  FormatHisto( efficiency, kBlue );
  FormatHisto( mult, kBlack );
  mult->GetYaxis()->SetRangeUser(0,1);

  mult->Draw("hist");
  purity->Draw("histsame");
  efficiency->Draw("histsame");
  TLegend*leg2 = new TLegend(.7,.7,.9,.9);
  leg2->AddEntry( purity, "proton purity" );
  leg2->AddEntry( efficiency, "proton efficiency" );
  leg2->AddEntry( mult, "product = purity #times eff");
  leg2->Draw();
  c1->SetGrid();
  c1->Print("Blob_effpur.pdf");

  return 0;

}

