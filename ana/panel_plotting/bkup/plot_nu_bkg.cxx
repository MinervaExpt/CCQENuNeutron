//#include "myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
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
#include "localColor.h"
#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"

using namespace PlotUtils;

void makePlots(bool doMultipliers, string location,string sample)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);


  std::vector<double> hdptbins;
  std::vector<double> hdpzbins;
  std::vector<double> hdrecoilbins_short;//to test 3D code via a dummy 2D 
  std::vector<double> hdrecoilbins;
  
  hdptbins.push_back(0.0);
  hdptbins.push_back(0.15);
  hdptbins.push_back(0.4);
  hdptbins.push_back(1.0);
  hdptbins.push_back(2.5);

  hdpzbins.push_back(1.5);
  hdpzbins.push_back(5.0);
  hdpzbins.push_back(20.0);

  hdrecoilbins_short.push_back(0.);
  hdrecoilbins_short.push_back(1500.);

  hdrecoilbins.push_back(0.);
  hdrecoilbins.push_back(60.);
  hdrecoilbins.push_back(140.);
  hdrecoilbins.push_back(240.);
  hdrecoilbins.push_back(480.);
  hdrecoilbins.push_back(1000.);
  
  std::vector<std::vector<double> > full3D;
  full3D.push_back(hdrecoilbins);
  full3D.push_back(hdptbins);
  full3D.push_back(hdpzbins);

  
  HyperDimLinearizer *my3d = new HyperDimLinearizer(full3D,0);



  TFile f1(Form("%s_CV/SideBand3D_SpecialSampleIncluded.root",location.c_str()));//scale
  TFile f2(Form("%s_CV/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-%s_SpecialSampleIncluded.root",location.c_str(),sample.c_str()));//sample 1 trk
  TFile f3(Form("%s_CV/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-%s_SpecialSampleIncluded.root",location.c_str(),sample.c_str()));//sample 2 trk


  string scaleform = "";
  if(sample=="MichelSideBand") scaleform = "_michel";
  if(sample=="BlobSideBand") scaleform = "_blobs";
  if(sample=="MicBlobSideBand") scaleform = "_micblob";

  MnvH2D* scale_1track = (MnvH2D*)f1.Get(Form("h_weights_1track_pzptrecbins_qelikenot%s",scaleform.c_str()));
  MnvH2D* scale_2track = (MnvH2D*)f1.Get(Form("h_weights_1track_pzptrecbins_qelikenot%s",scaleform.c_str()));

  MnvH2D* dataMnv_1trk=(MnvH2D*)f2.Get("h_pzptrec_data");
  MnvH2D* dataMnv_2trk=(MnvH2D*)f3.Get("h_pzptrec_data");
  
  MnvH2D* mcsigMnv_1trk=(MnvH2D*)f2.Get("h_pzptrec_qelike");
  MnvH2D* mcsigMnv_2trk=(MnvH2D*)f3.Get("h_pzptrec_qelike");

  MnvH2D* mcbkgMnv_1trk=(MnvH2D*)f2.Get("h_pzptrec_qelikenot");
  MnvH2D* mcbkgMnv_2trk=(MnvH2D*)f3.Get("h_pzptrec_qelikenot");
  
  //unscaled
  MnvH2D* mcMnv_1trk=(MnvH2D*)f2.Get("h_pzptrec_mc");
  MnvH2D* mcMnv_2trk=(MnvH2D*)f3.Get("h_pzptrec_mc");



  //apply scale
  mcbkgMnv_1trk->Multiply(scale_1track,mcbkgMnv_1trk);
  mcbkgMnv_2trk->Multiply(scale_2track,mcbkgMnv_2trk);

  //add samples
  mcsigMnv_1trk->Add(mcbkgMnv_1trk);
  mcsigMnv_2trk->Add(mcbkgMnv_2trk);
  

  std::cout << "Starting up getting the projections" << std::endl;
  std::vector<TH2D*> dataresults_1trk = my3d->Get2DHistos(dataMnv_1trk,true);
  std::vector<TH2D*> mcresults_1trk = my3d->Get2DHistos(mcMnv_1trk,true);
  std::vector<TH2D*> scaledmcresults_1trk = my3d->Get2DHistos(mcsigMnv_1trk,true);

  std::vector<TH2D*> dataresults_2trk = my3d->Get2DHistos(dataMnv_2trk,true);
  std::vector<TH2D*> mcresults_2trk = my3d->Get2DHistos(mcMnv_2trk,true);
  std::vector<TH2D*> scaledmcresults_2trk = my3d->Get2DHistos(mcsigMnv_2trk,true);
  
  



  
  for(unsigned int i=1;i<hdpzbins.size();i++){
    double pzwidth = hdpzbins[i]-hdpzbins[i-1];
    dataresults_1trk[i]->Scale(1/pzwidth,"width");
    mcresults_1trk[i]->Scale(1/pzwidth,"width");
    scaledmcresults_1trk[i]->Scale(1/pzwidth,"width");

    dataresults_2trk[i]->Scale(1/pzwidth,"width");
    mcresults_2trk[i]->Scale(1/pzwidth,"width");
    scaledmcresults_2trk[i]->Scale(1/pzwidth,"width");
  }
  
  vector<int> mycolors = getColors(2);
  for(unsigned int i=1;i<hdpzbins.size();i++){
    cout << "DOING BIN " << i << endl;
    mcresults_1trk[i]->SetLineColor(kRed);
    mcresults_1trk[i]->SetLineWidth(2);
    mcresults_2trk[i]->SetLineColor(kRed);
    mcresults_2trk[i]->SetLineWidth(2);

    scaledmcresults_1trk[i]->SetLineColor(kBlue);
    scaledmcresults_1trk[i]->SetLineWidth(2);
    scaledmcresults_2trk[i]->SetLineColor(kBlue);
    scaledmcresults_2trk[i]->SetLineWidth(2);

    dataresults_1trk[i]->SetMarkerStyle(kFullCircle);
    dataresults_1trk[i]->SetMarkerSize(0.7);
    dataresults_1trk[i]->SetLineColor(kBlack);
    dataresults_1trk[i]->SetLineWidth(2);
    dataresults_2trk[i]->SetMarkerStyle(kFullCircle);
    dataresults_2trk[i]->SetMarkerSize(0.7);
    dataresults_2trk[i]->SetLineColor(kBlack);
    dataresults_2trk[i]->SetLineWidth(2);




    // Make a list of the histograms we want to draw, along with the
    // draw options we want to use for them. You can add "graph" to the
    // draw options if you want the histogram to be converted to a graph
    // and then drawn. In that case the draw options are interpreted as
    // options to TGraphErrors::Draw().
    //
    // I don't know what happens if you put a "graph" first in the list,
    // so don't do that. Make sure the first item doesn't have "graph"
    // in its options
    std::vector<std::pair<TH2*, const char*> > histAndOpts_1trk;
    histAndOpts_1trk.push_back(std::make_pair(dataresults_1trk[i], "histpe"));
    histAndOpts_1trk.push_back(std::make_pair(mcresults_1trk[i],       "histl"));
    histAndOpts_1trk.push_back(std::make_pair(scaledmcresults_1trk[i],       "graph le1"));


    std::vector<std::pair<TH2*, const char*> > histAndOpts_2trk;
    histAndOpts_2trk.push_back(std::make_pair(dataresults_2trk[i], "histpe"));
    histAndOpts_2trk.push_back(std::make_pair(mcresults_2trk[i],       "histl"));
    histAndOpts_2trk.push_back(std::make_pair(scaledmcresults_2trk[i],       "graph le1"));


    cout << "Mults" << endl;


    // ----------------------------------------------------------------------------------
    //
    // First make pt in bins of pz
      
    // Values to multiply each bin by to get them on a similar range
    //These are for the XY!! 
    double multipliers1[5];
    double multipliers2[5];
    /*
    if(sample=="MichelSideBand"){
      multipliers1 = {2,1,1,3,20};
      multipliers2 = {50,10,20,50,200};
    }

    if(sample=="BlobSideBand"){
      multipliers1 = {500,20,4,2,3};
      multipliers2 = {10000,500,75,25,30};
    }

   if(sample=="MicBlobSideBand"){
      multipliers1 = {10000,200,75,10,20};
      multipliers2 = {10000,1000,1000,200,200};
    }
    */


    // Example of adding a legend. The co-ordinate system is NDC on the
    // entire canvas, ie (0,0) in the bottom left corner of the canvas
    // (not the individual pad), and (1,1) in the top right
    TLegend* leg=new TLegend(0.7, 0.2, 1, 0.5);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(dataresults_1trk[i], "MINERvA data", "lpe");
    leg->AddEntry(mcresults_1trk[i], "MnvGENIE", "l");
    leg->AddEntry(scaledmcresults_1trk[i], "MnvGENIE w/ sc. bkg.", "l");
    

    TLatex mytex;
    mytex.SetTextSize(0.05);
    string mystring =     Form("%.2f < P_{||} [GeV] < %.2f",hdpzbins[i-1]/1000,hdpzbins[i]/1000);
    mytex.DrawLatex(0.35,0.96,mystring.c_str());


    GridCanvas* gc1=NULL;
    if(i==1)gc1=plotpT1D(histAndOpts_1trk, doMultipliers ? multipliers1 : NULL);
    if(i==2)gc1=plotpT1D(histAndOpts_1trk, doMultipliers ? multipliers2 : NULL);
    // Set the y range manually. Can also use gc1->Remax() to guess automatically
    gc1->SetYLimits(0, 10);
    //gc1->SetYLimits(0, 1e-39);
    gc1->SetYTitle("d^{3}N/dp_{T}dp_{||}d_{E_{vis}} (cm^{2}/GeV^{2}MeV)");
    gc1->Modified();
  
    leg->Draw("SAME");
    gc1->Print(doMultipliers ? Form("nu-3d-bkg_%s-comps-pt-multiplier_bin_%d-1track.eps",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-1track.eps",sample.c_str(),i));
    gc1->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-1track.png",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-1track.png",sample.c_str(),i));
    gc1->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-1track.C",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-1track.C",sample.c_str(),i));




    //Now 2 track

    /*
    if(sample=="MichelSideBand"){
      multipliers1 = {400,7,2,2,5};
      multipliers2 = {500,100,40,40,100};
    }

    if(sample=="BlobSideBand"){
      multipliers1 = {50000,200,80,5,3};
      multipliers2 = {10000,50000,1500,75,20};
    }

   if(sample=="MicBlobSideBand"){
      multipliers1 = {10000,200,150,30,10};
      multipliers2 = {10000,1000,2000,400,100};
    }
    */

    GridCanvas* gc2=NULL;
    if(i==1)gc2=plotpT1D(histAndOpts_2trk, doMultipliers ? multipliers1 : NULL);
    if(i==2)gc2=plotpT1D(histAndOpts_2trk, doMultipliers ? multipliers2 : NULL);

    // Set the y range manually. Can also use gc2->Remax() to guess automatically
    gc2->SetYLimits(0, 10);
    //gc2->SetYLimits(0, 1e-39);
    gc2->SetYTitle("d^{3}N/dp_{T}dp_{||}d_{E_{vis}} (cm^{2}/GeV^{2}MeV)");
    gc2->Modified();
  
    leg->Draw("SAME");
    gc2->Print(doMultipliers ? Form("nu-3d-bkg_%s-comps-pt-multiplier_bin_%d-2track.eps",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-2track.eps",sample.c_str(),i));
    gc2->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-2track.png",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-2track.png",sample.c_str(),i));
    gc2->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-2track.C",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-2track.C",sample.c_str(),i));
  }
}

int main(int argc, char* argv[])
{
  makePlots(true,argv[1],"MichelSideBand");
  makePlots(true,argv[1],"BlobSideBand");
  makePlots(true,argv[1],"MicBlobSideBand");


  return 0;
}
