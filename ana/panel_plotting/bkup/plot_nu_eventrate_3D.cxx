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

void makePlots(bool doMultipliers, bool doTracks, bool doBkgSub, string location)
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
  hdptbins.push_back(0.075);//added ME
  hdptbins.push_back(0.15);
  hdptbins.push_back(0.25);//added ME
  hdptbins.push_back(0.4);
  hdptbins.push_back(0.7);//added ME
  hdptbins.push_back(1.0);
  hdptbins.push_back(1.5);//added ME
  hdptbins.push_back(2.5);

  hdpzbins.push_back(1.5);
  hdpzbins.push_back(3.5);//added ME
  hdpzbins.push_back(5.0);
  hdpzbins.push_back(10.0);//added ME
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




  //three files 1-track, 2+track, N-track
  //CV
  TFile f1(Form("%s_CV/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal_CombinedPlaylists.root",location.c_str()));//1track
  TFile f2(Form("%s_CV/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_CombinedPlaylists.root",location.c_str()));//2track
  TFile f3(Form("%s_CV/MuonEventSelection_MakeFlux-1_Multiplicity-0_Sample-Signal_CombinedPlaylists.root",location.c_str()));//Ntrack
  //CV without low recoil tune
  TFile f4(Form("%s_pion_rpa/MuonEventSelection_MakeFlux-1_Multiplicity-0_Sample-Signal_CombinedPlaylists.root",location.c_str()));//NTrack
  //Constraint File
  TFile f5(Form("%s_CV/SideBand3D_CombinedPlaylists.root",location.c_str()));

  //need pzmuptmu
  MnvH2D* mcMnv=(MnvH2D*)f3.Get("h_pzptrec_mc");//Get from N track
  MnvH2D* dataMnv = (MnvH2D*)f3.Get("h_pzptrec_data");//Get from N track
  MnvH2D* SingleTrackMnv = (MnvH2D*)f1.Get("h_pzptrec_data");//Get from 1 track
  MnvH2D* MultiTrackMnv = (MnvH2D*)f2.Get("h_pzptrec_data");//Get from 2+ track
  MnvH2D* SingleTrackMCMnv = (MnvH2D*)f1.Get("h_pzptrec_mc");//Get from 1 track
  MnvH2D* MultiTrackMCMnv = (MnvH2D*)f2.Get("h_pzptrec_mc");//Get from 2+ track

  MnvH2D* mcMnv_qelike_qe = (MnvH2D*)f3.Get("h_pzptrec_qelike_qe");//Get from N track
  MnvH2D* mcMnv_qelike_res = (MnvH2D*)f3.Get("h_pzptrec_qelike_res");//Get from N track
  MnvH2D* mcMnv_qelike_dis = (MnvH2D*)f3.Get("h_pzptrec_qelike_dis");//Get from N track
  MnvH2D* mcMnv_qelike_2p2h = (MnvH2D*)f3.Get("h_pzptrec_qelike_oth");//Get from N track
  MnvH2D* mcMnv_qelike_2p2h_no_lowrec = (MnvH2D*)f4.Get("h_pzptrec_qelike_oth");//Get from N track
  //bkgs
  MnvH2D* mcMnv_1track_bkg = (MnvH2D*)f1.Get("h_pzptrec_qelikenot");//Get from 1 track
  MnvH2D* mcMnv_2track_bkg = (MnvH2D*)f2.Get("h_pzptrec_qelikenot");//Get from 2+ track
  //the constraint
  MnvH2D* bkgConstraint_1track = (MnvH2D*)f5.Get("h_weights_1track_pzptrecbins_qelikenot");//get from 1 track
  MnvH2D* bkgConstraint_2track = (MnvH2D*)f5.Get("h_weights_2track_pzptrecbins_qelikenot");//get from 2 track
  //apply constraint
  mcMnv_1track_bkg->Multiply(mcMnv_1track_bkg,bkgConstraint_1track);
  mcMnv_2track_bkg->Multiply(mcMnv_2track_bkg,bkgConstraint_2track);
  
  MnvH2D* mcMnv_bkg = (MnvH2D*)mcMnv_1track_bkg->Clone("tmp");
  mcMnv_bkg->Add(mcMnv_2track_bkg);
  
  dataMnv->GetXaxis()->SetTitle("p_{||} (GeV)");

  if(doBkgSub){
    dataMnv->Add(mcMnv_bkg,-1);
    mcMnv->Add(mcMnv_bkg,-1);
    SingleTrackMnv->Add(mcMnv_1track_bkg,-1);
    MultiTrackMnv->Add(mcMnv_2track_bkg,-1);
    SingleTrackMCMnv->Add(mcMnv_1track_bkg,-1);
    MultiTrackMCMnv->Add(mcMnv_2track_bkg,-1);
  }
  

  std::vector<TH2D*> dataStat = my3d->Get2DHistos(dataMnv,false);
  std::vector<TH2D*> data = my3d->Get2DHistos(dataMnv,true);
  std::vector<TH2D*> mc = my3d->Get2DHistos(mcMnv,false);
  std::vector<TH2D*> mc_qelike_qe = my3d->Get2DHistos(mcMnv_qelike_qe,false);
  std::vector<TH2D*> mc_qelike_res = my3d->Get2DHistos(mcMnv_qelike_res,false);
  std::vector<TH2D*> mc_qelike_dis = my3d->Get2DHistos(mcMnv_qelike_dis,false);
  std::vector<TH2D*> mc_qelike_2p2h = my3d->Get2DHistos(mcMnv_qelike_2p2h,false);
  std::vector<TH2D*> mc_qelike_2p2h_no_lowrec = my3d->Get2DHistos(mcMnv_qelike_2p2h_no_lowrec,false);

  std::vector<TH2D*> mc_bkg = my3d->Get2DHistos(mcMnv_bkg,false);
  std::vector<TH2D*> data_1track = my3d->Get2DHistos(SingleTrackMnv,false);
  std::vector<TH2D*> data_2track = my3d->Get2DHistos(MultiTrackMnv,false);
  std::vector<TH2D*> mc_1track = my3d->Get2DHistos(SingleTrackMCMnv,false);
  std::vector<TH2D*> mc_2track = my3d->Get2DHistos(MultiTrackMCMnv,false);

  // // Get the data histogram with stat error and with total error
  // // separately so we can plot them both for inner and outer ticks
  // TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  // TH2* data=new TH2D(dataMnv->GetCVHistoWithError());
  // TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());
  // TH2* mc_qelike_qe = new TH2D(mcMnv_qelike_qe->GetCVHistoWithStatError());
  // TH2* mc_qelike_res = new TH2D(mcMnv_qelike_res->GetCVHistoWithStatError());
  // TH2* mc_qelike_dis = new TH2D(mcMnv_qelike_dis->GetCVHistoWithStatError());
  // TH2* mc_qelike_2p2h = new TH2D(mcMnv_qelike_2p2h->GetCVHistoWithStatError());
  // //  TH2* mc_qelike_2p2h_no_lowrec = new TH2D(mcMnv_qelike_2p2h_no_lowrec->GetCVHistoWithStatError());
  // TH2* mc_bkg = new TH2D(mcMnv_bkg->GetCVHistoWithStatError());
  // TH2* data_1track = new TH2D(SingleTrackMnv->GetCVHistoWithStatError());
  // TH2* data_2track = new TH2D(MultiTrackMnv->GetCVHistoWithStatError());
  // TH2* mc_1track = new TH2D(SingleTrackMCMnv->GetCVHistoWithStatError());
  // TH2* mc_2track = new TH2D(MultiTrackMCMnv->GetCVHistoWithStatError());
  for(int i=1;i<hdpzbins.size();i++){

    double pzwidth = hdpzbins[i]-hdpzbins[i-1];

    dataStat[i]->Scale(1e-2/pzwidth,"width");
    data[i]->Scale(1e-2/pzwidth,"width");
    mc[i]->Scale(1e-2/pzwidth,"width");
    mc_qelike_qe[i]->Scale(1e-2/pzwidth,"width");
    mc_qelike_res[i]->Scale(1e-2/pzwidth,"width");
    mc_qelike_dis[i]->Scale(1e-2/pzwidth,"width");
    mc_qelike_2p2h[i]->Scale(1e-2/pzwidth,"width");
    mc_bkg[i]->Scale(1e-2/pzwidth,"width");
    data_1track[i]->Scale(1e-2/pzwidth,"width");
    data_2track[i]->Scale(1e-2/pzwidth,"width");
    mc_1track[i]->Scale(1e-2/pzwidth,"width");
    mc_2track[i]->Scale(1e-2/pzwidth,"width");
    mc_qelike_2p2h_no_lowrec[i]->Scale(1e-2/pzwidth,"width");


    // These line and marker styles will be propagated to the 1D plots
    vector<int> mycolors = getColors(2);
    mc[i]->SetLineColor(kRed);
    mc[i]->SetLineWidth(2);
    
    mc_1track[i]->SetLineColor(kViolet);
    mc_1track[i]->SetLineStyle(2);
    mc_2track[i]->SetLineColor(kViolet);
    mc_2track[i]->SetLineStyle(1);
    //need to add signal and bkg colors
    mc_qelike_qe[i]->SetLineColor(mycolors[3]);
    mc_qelike_res[i]->SetLineColor(mycolors[4]);
    mc_qelike_dis[i]->SetLineColor(mycolors[5]);
    mc_qelike_2p2h[i]->SetLineColor(mycolors[6]);
    mc_qelike_2p2h_no_lowrec[i]->SetLineColor(mycolors[6]);
    mc_qelike_2p2h_no_lowrec[i]->SetLineStyle(2);
    mc_bkg[i]->SetLineColor(mycolors[10]);
    //need to add 1track and 2 track
    
    // These line and marker styles will be propagated to the 1D plots
    data[i]->SetMarkerStyle(kFullCircle);
    data[i]->SetMarkerSize(0.7);
    data[i]->SetLineColor(kBlack);
    data[i]->SetLineWidth(2);
    
    dataStat[i]->SetMarkerStyle(1);
    dataStat[i]->SetLineColor(kBlack);
    dataStat[i]->SetLineWidth(2);
    
    data_1track[i]->SetLineColor(kBlack);
    data_1track[i]->SetLineStyle(2);
    data_2track[i]->SetLineColor(kBlack);
    data_2track[i]->SetLineStyle(1);
    
    
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
    histAndOpts.push_back(std::make_pair(mc[i],       "hist"));
    
    //Do by track breakdown
    if(doTracks){
      histAndOpts.push_back(std::make_pair(mc_1track[i],       "histl"));
      histAndOpts.push_back(std::make_pair(mc_2track[i],       "histl"));
      histAndOpts.push_back(std::make_pair(data_1track[i],       "histl"));
      histAndOpts.push_back(std::make_pair(data_2track[i],       "histl"));
    }
    
    //do by physics process
    else{
      histAndOpts.push_back(std::make_pair(mc_qelike_qe[i],       "histl"));
      histAndOpts.push_back(std::make_pair(mc_qelike_res[i],       "histl"));
      histAndOpts.push_back(std::make_pair(mc_qelike_dis[i],       "histl"));
      histAndOpts.push_back(std::make_pair(mc_qelike_2p2h[i],       "histl"));
      histAndOpts.push_back(std::make_pair(mc_qelike_2p2h_no_lowrec[i],       "histl"));
      histAndOpts.push_back(std::make_pair(mc_bkg[i],       "histl"));
    }
    
    histAndOpts.push_back(std::make_pair(dataStat[i], "graph e"));
    histAndOpts.push_back(std::make_pair(data[i],     "graph ep"));
    
    
    
    // ----------------------------------------------------------------------------------
    //
    // First make pt in bins of pz
    
    // Values to multiply each bin by to get them on a similar range
    double multipliers[]={1,1,1,2,5,1,1,1,1,1,1,1,1,1,1};
    double multipliers2[]={0.5,0.5,0.5,1,2,1,1,1,1,1,1,1,1,1,1};
    double multipliers3[]={0.5,1,1,2,7,1,1,1,1,1,1,1,1,1};
    double multipliers4[]={10,10,10,50,100,1,1,1,1,1,1,1,1,1,1};
    GridCanvas* gc=NULL;
    if(i==1) gc = plotpT1D(histAndOpts, doMultipliers ? multipliers : NULL);
    if(i==2) gc = plotpT1D(histAndOpts, doMultipliers ? multipliers2 : NULL);
    if(i==3) gc = plotpT1D(histAndOpts, doMultipliers ? multipliers3 : NULL);
    if(i==4) gc = plotpT1D(histAndOpts, doMultipliers ? multipliers4 : NULL);
    if(i==5) gc = plotpT1D(histAndOpts, doMultipliers ? multipliers2 : NULL);
    // Set the y range manually. Can also use gc->Remax() to guess automatically
    gc->SetYLimits(0, 5);
    gc->SetYTitle("Event Rate(x10^{-2}) per GeV^{2}");
    gc->Modified();
    // Example of adding a legend. The co-ordinate system is NDC on the
    // entire canvas, ie (0,0) in the bottom left corner of the canvas
    // (not the individual pad), and (1,1) in the top right
    TLegend* leg=new TLegend(0.23, 0.7, 0.39, 0.9);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(data[i], "MINERvA data", "lpe");
    leg->AddEntry(mc[i], "MINERvA Tune", "l");
    if(doTracks){
      leg->AddEntry(mc_1track[i], "MC 1-track","l");
      leg->AddEntry(mc_2track[i], "MC Multi-track","l");
      leg->AddEntry(data_1track[i], "data 1-track","l");
      leg->AddEntry(data_2track[i], "data Multi-track","l");
    }
    else{
      leg->AddEntry(mc_qelike_qe[i],"QE","l");
      leg->AddEntry(mc_qelike_res[i],"Resonant","l");
      leg->AddEntry(mc_qelike_dis[i],"DIS","l");
      leg->AddEntry(mc_qelike_2p2h[i],"2p2h","l");
      leg->AddEntry(mc_qelike_2p2h_no_lowrec[i],"2p2h without fit","l");
      leg->AddEntry(mc_bkg[i],"Background","l");
      
    }
    
    leg->Draw("SAME");
    if(doTracks){
      if(doBkgSub){
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-bkgsub-trackmult-pt-multiplier-bin-%d.eps",i) : Form("nu-2d-evtrate-bkgsub-trackmult-pt-bin-%d.eps",i));
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-bkgsub-trackmult-pt-multiplier-bin-%d.png",i) : Form("nu-2d-evtrate-bkgsub-trackmult-pt-bin-%d.png",i));
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-bkgsub-trackmult-pt-multiplier-bin-%d.C",i) : Form("nu-2d-evtrate-bkgsub-trackmult-pt-bin-%d.C",i));
      }
      else{
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-trackmult-pt-multiplier-bin-%d.eps",i) : Form("nu-2d-evtrate-trackmult-pt-bin-%d.eps",i));
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-trackmult-pt-multiplier-bin-%d.png",i) : Form("nu-2d-evtrate-trackmult-pt-bin-%d.png",i));
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-trackmult-pt-multiplier-bin-%d.C",i) : Form("nu-2d-evtrate-trackmult-pt-bin-%d.C",i));
      }
    }
    else{
      gc->Print(doMultipliers ? Form("nu-2d-evtrate-model-pt-multiplier-bin-%d.eps",i) : Form("nu-2d-evtrate-model-pt-bin-%d.eps",i));
      gc->Print(doMultipliers ? Form("nu-2d-evtrate-model-pt-multiplier-bin-%d.png",i) : Form("nu-2d-evtrate-model-pt-bin-%d.png",i));
      gc->Print(doMultipliers ? Form("nu-2d-evtrate-model-pt-multiplier-bin-%d.C",i) : Form("nu-2d-evtrate-model-pt-bin-%d.C",i));
    }
    
    // ----------------------------------------------------------------------------------
    //
    // Now make pz in bins of pt. It's all the same
    
    // Values to multiply each bin by to get them on a similar range
    double  multipliers5[]={10,2,1,1,1,5,75,1,1,1,1,1,1,1};
    double  multipliers6[]={5,1,0.5,0.5,0.5,1,5,500,1,1,1,1,1,1};
    double  multipliers7[]={10,2,1,0.5,1,2,5,100,1,1,1,1,1,1};
    double  multipliers8[]={100,50,25,25,25,50,100,500,1,1,1,1,1,1};
			    
    
    // plotpz1D fiddles the x axis values to squash up the tail so it
    // doesn't take up all the horizontal space.
    GridCanvas* gc2=NULL;
    if(i==1) gc2 = plotpz1D(histAndOpts, doMultipliers ? multipliers5 : NULL,true);
    if(i==2) gc2 = plotpz1D(histAndOpts, doMultipliers ? multipliers6 : NULL,true);
    if(i==3) gc2 = plotpz1D(histAndOpts, doMultipliers ? multipliers7 : NULL,true);
    if(i==4) gc2 = plotpz1D(histAndOpts, doMultipliers ? multipliers8 : NULL,true);
    if(i==5) gc2 = plotpz1D(histAndOpts, doMultipliers ? multipliers4 : NULL,true);

    gc2->SetYLimits(0, 5);
    gc2->SetYTitle("Event Rate(x10^{-2}) per GeV^{2}");
    gc2->Modified();
    if(doTracks){
      if(doBkgSub){
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-bkgsub-trackmult-pz-multiplier-bin-%d.eps",i) : Form("nu-2d-evtrate-bkgsub-trackmult-pz-bin-%d.eps",i));
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-bkgsub-trackmult-pz-multiplier-bin-%d.png",i) : Form("nu-2d-evtrate-bkgsub-trackmult-pz-bin-%d.png",i));
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-bkgsub-trackmult-pz-multiplier-bin-%d.C",i) : Form("nu-2d-evtrate-bkgsub-trackmult-pz-bin-%d.C",i));
      }
      else{
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-trackmult-pz-multiplier-bin-%d.eps",i) : Form("nu-2d-evtrate-trackmult-pz-bin-%d.eps",i));
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-trackmult-pz-multiplier-bin-%d.png",i) : Form("nu-2d-evtrate-trackmult-pz-bin-%d.png",i));
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-trackmult-pz-multiplier-bin-%d.C",i) : Form("nu-2d-evtrate-trackmult-pz-bin-%d.C",i));
      }
    }
    else{
      gc2->Print(doMultipliers ? Form("nu-2d-evtrate-model-pz-multiplier-bin-%d.eps",i) : Form("nu-2d-evtrate-model-pz-bin-%d.eps",i));
      gc2->Print(doMultipliers ? Form("nu-2d-evtrate-model-pz-multiplier-bin-%d.png",i) : Form("nu-2d-evtrate-model-pz-bin-%d.png",i));
      gc2->Print(doMultipliers ? Form("nu-2d-evtrate-model-pz-multiplier-bin-%d.C",i) : Form("nu-2d-evtrate-model-pz-bin-%d.C",i));
    }
  }
}

int main(int argc, char* argv[])
{

  string location = argv[1];
  //multipliers
  makePlots(true,true,false,location);
  makePlots(true,true,true,location);
  makePlots(true,false,false,location);
  //standard
  makePlots(false,true,false,location);
  makePlots(false,true,true,location);
  makePlots(false,false,false,location);

  return 0;
}
