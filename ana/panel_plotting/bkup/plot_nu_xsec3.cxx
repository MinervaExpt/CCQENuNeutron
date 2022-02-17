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

void makePlots(bool doMultipliers,bool doGenies, string location)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  
  TFile f1(Form("%s_CV/CrossSection_per_nucleon_3D_pzptreco_iterations_4_CombinedPlaylists.root_big3d.root",location.c_str()));//final
  MnvH2D* dataMnv=(MnvH2D*)f1.Get("h_pzptrec_data_nobck_unfold_effcor_cross_section");
  MnvH2D* mcMnv=(MnvH2D*)f1.Get("h_pzptrec_mc_nobck_unfold_effcor_cross_section");
  /*
  dataMnv->PopVertErrorBand("GENIE_MaRES");
  dataMnv->PopVertErrorBand("GENIE_MvRES");
  dataMnv->PopVertErrorBand("GENIE_NormCCRES");
  dataMnv->PopVertErrorBand("GENIE_NormDISCC");
  dataMnv->PopVertErrorBand("Reweight_Neutron");
  dataMnv->PopVertErrorBand("Reweight_Pion");
  dataMnv->PopVertErrorBand("Reweight_Proton");

  dataMnv->PopLatErrorBand("Birks_Response_Proton");
  */
  //Model Component
  MnvH2D* mcMnv_qelike_qe = (MnvH2D*)f1.Get("h_pzptrec_cross_section_qelike_qe");//Get from N track
  MnvH2D* mcMnv_qelike_res = (MnvH2D*)f1.Get("h_pzptrec_cross_section_qelike_res");//Get from N track
  MnvH2D* mcMnv_qelike_dis = (MnvH2D*)f1.Get("h_pzptrec_cross_section_qelike_dis");//Get from N track
  MnvH2D* mcMnv_qelike_2p2h = (MnvH2D*)f1.Get("h_pzptrec_cross_section_qelike_2p2h");//Get from N track


  //  MnvH2D* nomGenieMnv = (MnvH2D*)f2.Get("h_pzptrec_mc_nobck_unfold_effcor_cross_section");
  //  MnvH2D* bestGenieMnv = (MnvH2D*)f3.Get("h_pzptrec_mc_nobck_unfold_effcor_cross_section");
 
  vector<double> recoil3Dbins;
  vector<double> pt3Dbins;
  vector<double> pz3Dbins;

  pt3Dbins.push_back(0.0);
  //  pt3Dbins.push_back(0.075);//added ME
  pt3Dbins.push_back(0.15);
  pt3Dbins.push_back(0.25);//added ME
  pt3Dbins.push_back(0.4);
  pt3Dbins.push_back(0.7);//added ME
  pt3Dbins.push_back(1.0);
  //  pt3Dbins.push_back(1.5);//added ME
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
  std::vector<TH2D*> dataresults = my3d->Get2DHistos(dataMnv,true);
  std::vector<TH2D*> dataresults_statonly = my3d->Get2DHistos(dataMnv,false);
  std::vector<TH2D*> mcresults = my3d->Get2DHistos(mcMnv,false);
  std::cout << "Starting up getting the projections of the subcomponents" << std::endl;
  std::vector<TH2D*> mcresults_qelike_2p2h = my3d->Get2DHistos(mcMnv_qelike_2p2h,false);
  std::vector<TH2D*> mcresults_qelike_qe = my3d->Get2DHistos(mcMnv_qelike_qe,false);
  std::vector<TH2D*> mcresults_qelike_res = my3d->Get2DHistos(mcMnv_qelike_res,false);
  std::vector<TH2D*> mcresults_qelike_dis = my3d->Get2DHistos(mcMnv_qelike_dis,false);

  cout << "I found " << mcresults.size() << " histograms to work with " << endl;
  for(unsigned int i=1;i<pz3Dbins.size();i++){
    double pzwidth = pz3Dbins[i]-pz3Dbins[i-1];
    dataresults[i]->Scale(1e42/pzwidth,"width");
    dataresults_statonly[i]->Scale(1e42/pzwidth,"width");
    mcresults[i]->Scale(1e42/pzwidth,"width");
    mcresults_qelike_2p2h[i]->Scale(1e42/pzwidth,"width");
    mcresults_qelike_qe[i]->Scale(1e42/pzwidth,"width");
    mcresults_qelike_res[i]->Scale(1e42/pzwidth,"width");
    mcresults_qelike_dis[i]->Scale(1e42/pzwidth,"width");
  }

  vector<int> mycolors = getColors(2);
  for(unsigned int i=1;i<pz3Dbins.size();i++){//skip underflow and overflow pz bins
    cout << "DOING BIN " << i << endl;
    // Get the data histogram with stat error and with total error
    // separately so we can plot them both for inner and outer ticks
    //  TH2* nomGenie=new TH2D(nomGenieMnv->GetCVHistoWithStatError());
    //  TH2* bestGenie = new TH2D(bestGenieMnv->GetCVHistoWithStatError());
    /*
      TH2* mc_qelike_qe = new TH2D(mcMnv_qelike_qe->GetCVHistoWithStatError());
      TH2* mc_qelike_res = new TH2D(mcMnv_qelike_res->GetCVHistoWithStatError());
      TH2* mc_qelike_dis = new TH2D(mcMnv_qelike_dis->GetCVHistoWithStatError());
      TH2* mc_qelike_2p2h = new TH2D(mcMnv_qelike_2p2h->GetCVHistoWithStatError());
      TH2* mc_qelike_2p2h_no_lowrec = new TH2D(mcMnv_qelike_2p2h_no_lowrec->GetCVHistoWithStatError());
    */
    // These line and marker styles will be propagated to the 1D plots
    

      mcresults[i]->SetLineColor(kRed);
      mcresults[i]->SetLineWidth(2);
  
      //      nomGenie->SetLineColor(kBlue);
      //      nomGenie->SetLineWidth(2);

      //      bestGenie->SetLineColor(mycolors[10]);
      //      bestGenie->SetLineWidth(2);

      mcresults_qelike_qe[i]->SetLineColor(mycolors[3]);
      mcresults_qelike_res[i]->SetLineColor(mycolors[4]);
      mcresults_qelike_dis[i]->SetLineColor(mycolors[5]);
      mcresults_qelike_2p2h[i]->SetLineColor(mycolors[6]);
      // mcresults_qelike_2p2h_no_lowrec->SetLineColor(mycolors[6]);
      // mcresults_qelike_2p2h_no_lowrec->SetLineStyle(3);

      // These line and marker styles will be propagated to the 1D plots
      dataresults[i]->SetMarkerStyle(kFullCircle);
      dataresults[i]->SetMarkerSize(0.7);
      dataresults[i]->SetLineColor(kBlack);
      dataresults[i]->SetLineWidth(2);

      dataresults_statonly[i]->SetMarkerStyle(1);
      dataresults_statonly[i]->SetMarkerSize(100);
      dataresults_statonly[i]->SetLineColor(kViolet);
      dataresults_statonly[i]->SetLineWidth(10);

      // dataStat->SetMarkerStyle(1);
      // dataStat->SetLineColor(kBlack);
      // dataStat->SetLineWidth(2);
    

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
    histAndOpts.push_back(std::make_pair(dataresults[i], "histpe"));
    histAndOpts.push_back(std::make_pair(mcresults[i],       "histl"));
    histAndOpts.push_back(std::make_pair(mcresults_qelike_qe[i],       "histl"));
    histAndOpts.push_back(std::make_pair(mcresults_qelike_res[i],       "histl"));
    histAndOpts.push_back(std::make_pair(mcresults_qelike_dis[i],       "histl"));
    histAndOpts.push_back(std::make_pair(mcresults_qelike_2p2h[i],       "histl"));
    histAndOpts.push_back(std::make_pair(dataresults_statonly[i], "histpe"));

    histAndOpts.push_back(std::make_pair(dataresults[i],     "graph0 ep"));

    cout << "Mults" << endl;


    // ----------------------------------------------------------------------------------
    //
    // First make pt in bins of pz

    // Values to multiply each bin by to get them on a similar range
    //These are for the XY!! 
    double multipliers1[]={1,1,1,1.5,2,2,2,2,2,4,4,8,15};
    double multipliers2[]={0.5,0.8,0.8,1,1,1,1,1,1,2,2,4,10};
    double multipliers3[]={10,10,10,10,10,10,25,25,40,50,75,150,200};
    double multipliers4[]={25,25,40,80,200};
    double multipliers5[]={4,100,1,1,1};
    double multipliers6[]={50,20,1,1,1};
    double multipliers7[]={10,4,1,1,1};

    GridCanvas* gc=NULL;
    if(i==1)gc=plotYAxis1D(histAndOpts, "Muon Transverse Momentum (GeV)", "#Sigma T_{p}", doMultipliers ? multipliers1 : NULL);
    if(i==2)gc=plotYAxis1D(histAndOpts, "Muon Transverse Momentum (GeV)", "#Sigma T_{p}", doMultipliers ? multipliers2 : NULL);
    if(i==3)gc=plotYAxis1D(histAndOpts, "Muon Transverse Momentum (GeV)", "#Sigma T_{p}", doMultipliers ? multipliers3 : NULL);
    if(i==4)gc=plotYAxis1D(histAndOpts, "Muon Transverse Momentum (GeV)", "#Sigma T_{p}", doMultipliers ? multipliers4 : NULL);
    if(i==5)gc=plotYAxis1D(histAndOpts, "Muon Transverse Momentum (GeV)", "#Sigma T_{p}", doMultipliers ? multipliers5 : NULL);
    if(i==6)gc=plotYAxis1D(histAndOpts, "Muon Transverse Momentum (GeV)", "#Sigma T_{p}", doMultipliers ? multipliers6 : NULL);
    if(i==7)gc=plotYAxis1D(histAndOpts, "Muon Transverse Momentum (GeV)", "#Sigma T_{p}", doMultipliers ? multipliers7 : NULL);
    // Set the y range manually. Can also use gc->Remax() to guess automatically
    gc->SetYLimits(0.01, 4.59);
    //gc->SetYLimits(0, 1e-39);
    gc->SetYTitle("d^{3}#sigma/dp_{T}dp_{||}d_{E_{vis}} (x10^{-42} cm^{2}/GeV^{3}/c^{3}/Nucleon)");
    gc->SetLogy(true);
    gc->Modified();
    // Example of adding a legend. The co-ordinate system is NDC on the
    // entire canvas, ie (0,0) in the bottom left corner of the canvas
    // (not the individual pad), and (1,1) in the top right
    TLegend* leg=new TLegend(0.8, 0.1, 1, 0.3);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(dataresults[i], "MINERvA data", "lpe");
    leg->AddEntry(mcresults[i], "MnvGENIE", "l");
    leg->AddEntry(mcresults_qelike_qe[i], "QELike-QE", "l");
    leg->AddEntry(mcresults_qelike_res[i], "QELike-RES", "l");
    leg->AddEntry(mcresults_qelike_dis[i], "QELike-DIS", "l");
    leg->AddEntry(mcresults_qelike_2p2h[i], "QELike-2p2h", "l");
    

    TLatex mytex;
    mytex.SetTextSize(0.05);
    string mystring =     Form("%.2f < P_{||} [GeV] < %.2f",pz3Dbins[i-1],pz3Dbins[i]);
    mytex.DrawLatex(0.35,0.96,mystring.c_str());
    leg->Draw("SAME");
    gc->Print(doMultipliers ? Form("nu-3d-xsec-comps-pt-multiplier_bin_%d.eps",i) : Form("nu-3d-xsec-comps-pt_bin_%d.eps",i));
    gc->Print(doMultipliers ? Form("nu-3d-xsec-comps-pt-multiplier_bin_%d.png",i) : Form("nu-3d-xsec-comps-pt_bin_%d.png",i));
    gc->Print(doMultipliers ? Form("nu-3d-xsec-comps-pt-multiplier_bin_%d.C",i) : Form("nu-3d-xsec-comps-pt_bin_%d.C",i));
    

    // // ----------------------------------------------------------------------------------
    // //
    // // Now make pz in bins of pt. It's all the same

    // // Values to multiply each bin by to get them on a similar range
    double multiplierspz1[]={20,4,1.5,1,2,4,50,1000000};
    double multiplierspz2[]={10,2,1,0.5,1,2,20,2000};
    double multiplierspz3[]={100,30,15,15,20,40,350,10000};
    double multiplierspz4[]={100,100,20,40,50,100,400,10000};
    double multiplierspz5[]={100,20,5,200,1};
    double multiplierspz6[]={40,10,7.5,100,7500};
    double multiplierspz7[]={40,10,7.5,100,7500};
    
    TLegend* leg2=new TLegend(0.8, 0.1, 1, 0.4);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.03);
    leg2->AddEntry(dataresults[i], "MINERvA data", "lpe");
    leg2->AddEntry(mcresults[i], "MnvGENIE", "l");
    leg2->AddEntry(mcresults_qelike_qe[i], "QELike-QE", "l");
    leg2->AddEntry(mcresults_qelike_res[i], "QELike-RES", "l");
    leg2->AddEntry(mcresults_qelike_dis[i], "QELike-DIS", "l");
    leg2->AddEntry(mcresults_qelike_2p2h[i], "QELike-2p2h", "l");

    // // plotXAxis1D fiddles the x axis values to squash up the tail so it
    // // doesn't take up all the horizontal space.
    GridCanvas* gc2= NULL;
    if(i==1)gc2=plotXAxis1D(histAndOpts, "#Sigma T_{p} (MeV)", "P_{t,muon}", doMultipliers ? multiplierspz1 : NULL);
    if(i==2)gc2=plotXAxis1D(histAndOpts, "#Sigma T_{p} (MeV)", "P_{t,muon}", doMultipliers ? multiplierspz2 : NULL);
    if(i==3)gc2=plotXAxis1D(histAndOpts, "#Sigma T_{p} (MeV)", "P_{t,muon}", doMultipliers ? multiplierspz3 : NULL);
    if(i==4)gc2=plotXAxis1D(histAndOpts, "#Sigma T_{p} (MeV)", "P_{t,muon}", doMultipliers ? multiplierspz4 : NULL);
    if(i==5)gc2=plotXAxis1D(histAndOpts, "#Sigma T_{p} (MeV)", "P_{t,muon}", doMultipliers ? multiplierspz5 : NULL);
    if(i==6)gc2=plotXAxis1D(histAndOpts, "#Sigma T_{p} (MeV)", "P_{t,muon}", doMultipliers ? multiplierspz6 : NULL);
    if(i==7)gc2=plotXAxis1D(histAndOpts, "#Sigma T_{p} (MeV)", "P_{t,muon}", doMultipliers ? multiplierspz7 : NULL);
    gc2->SetYLimits(0.01, 5.49);
    gc2->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
    gc2->SetLogy(true);
    gc2->Modified();
    mytex.DrawLatex(0.35,0.96,mystring.c_str());
    leg2->Draw("SAME");
    if(doGenies){
      gc2->Print(doMultipliers ? Form("nu-3d-xsec-genies-pz-multiplier_bin_%d.eps",i) : Form("nu-3d-xsec-genies-pz_bin_%d.eps",i));
      gc2->Print(doMultipliers ? Form("nu-3d-xsec-genies-pz-multiplier_bin_%d.png",i) : Form("nu-3d-xsec-genies-pz_bin_%d.png",i));
      gc2->Print(doMultipliers ? Form("nu-3d-xsec-genies-pz-multiplier_bin_%d.C",i) : Form("nu-3d-xsec-genies-pz_bin_%d.C",i));
    }
    else{
      gc2->Print(doMultipliers ? Form("nu-3d-xsec-comps-pz-multiplier_bin_%d.eps",i) : Form("nu-3d-xsec-comps-pz_bin_%d.eps",i));
      gc2->Print(doMultipliers ? Form("nu-3d-xsec-comps-pz-multiplier_bin_%d.png",i) : Form("nu-3d-xsec-comps-pz_bin_%d.png",i));
      gc2->Print(doMultipliers ? Form("nu-3d-xsec-comps-pz-multiplier_bin_%d.C",i) : Form("nu-3d-xsec-comps-pz_bin_%d.C",i));
    }

  }
}

int main(int argc, char* argv[])
{
  //  makePlots(true,true);
  makePlots(true,false,argv[1]);
  //  makePlots(false,true);
  makePlots(false,false,argv[1]);

  return 0;
}
