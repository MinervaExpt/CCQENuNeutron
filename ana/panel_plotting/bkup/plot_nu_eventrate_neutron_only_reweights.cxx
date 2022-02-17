//#include "myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvH1D.h"
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
#include "TROOT.h"
#include "localColor.h"
#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"


#include "include/CCQENuPlotUtils.h"

#include "include/CCQENuUtils.h"
#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/HyperDimLinearizer.h"

#include "include/CommonBins.h"
#include "include/GeneralFunc.h"



using namespace PlotUtils;
//forward declaration
//gROOT->SetBatch();

int nPlaylists = 7;
void drawPlots(string savename, TH2* data, TH2* dataStat, TH2* mc, TH2* mc_qelike_qe_h, TH2* mc_qelike_qe_oth, TH2* mc_qelike_res, TH2* mc_qelike_dis, TH2* mc_qelike_2p2h, string xtitle, string ytitle, string celltitle, string axis="x", double scale =1, double ymin=0, double ymax=2.4, double * multipliers=NULL, bool ratio = false, bool setlogx = false , bool setlogy = false);

void draw1DPlots(CCQENuPlotUtils* utils, string savename, MnvH2D** hists, double pot_data,  double pot_mc, int bmin=0, int bmax=-1, string axis="x", bool bkSub = false);

void makePlots(bool doMultipliers, bool doBkgSub, string fname_signal, string fname_bkg_weight, string var, string sideband="Signal", bool drawRCut=false )
{
  bool doCenter = (var == "dthetaPq2qe" && drawRCut == true );
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  //three files 2-track, 2-track bkg_fitted
  //CV
  TFile* f1 = new TFile( fname_signal.c_str(), "READ");// original signal histogram
  TFile* f2 = new TFile( fname_bkg_weight.c_str(), "READ");// weighted histogram


  TVector2 *pot = (TVector2*)f1->Get("pot;1");
  for( int i = 2; i<=nPlaylists; i++ ) (*pot)+= *((TVector2*)f1->Get(Form( "pot;%d", i) ));


  double pot_data = pot->X();
  double pot_mc = pot->Y();
  double pot_norm = pot_data/pot_mc;

  bool fluxHistoExists = true;
  CCQENuPlotUtils *putils = new CCQENuPlotUtils( fluxHistoExists );
  CCQENuUtils *utils = new CCQENuUtils( false, fluxHistoExists );
  axis_binning xbins = (var=="dthetaRq2qe")?dthetaReactbins:dthetaPerpbins;
  axis_binning ybins = muonPtbins;
  axis_binning zbins = Q2bins;
  HyperDimLinearizer* hdl = GetHDL(xbins , ybins, zbins);

  // Get CV histos
  MnvH2D *h_signal[nHistos];
  if (doCenter )
  {
    utils->bookHistos( f1, h_signal, Form( "h_%s_ptmu_center", var.c_str() ) );
  }
  else 
  {
    utils->bookHistos( f1, h_signal, Form( "h_%s_ptmu", var.c_str() ) );
  }
  putils->scaleMCHistos( h_signal, pot_norm );

  // Get Background weights
  MnvH2D *h_signal_bkg_weights = (MnvH2D*) f2->Get( Form( "h_weights_%sptbins_qelikenot_%s", var.c_str(), sideband.c_str() ));

  // Get Tuned Backgrounds:
  MnvH2D* h_signal_bkg = (MnvH2D*) h_signal[kQELikeNot]->Clone( Form("h_%s_bkg", var.c_str() ) );
  h_signal_bkg->Multiply( h_signal_bkg, h_signal_bkg_weights );

  // Get Tuned Data Signal
  MnvH2D* h_signal_mc_nobkg = (MnvH2D*) h_signal[kMC]->Clone( Form("h_%s_mc_nobkg", var.c_str() ));
  if (doBkgSub ) h_signal_mc_nobkg->Add( h_signal[kQELikeNot], -1 );

  MnvH2D* h_signal_data_nobkg = (MnvH2D*) h_signal[kData]->Clone( Form("h_%s_data_nobkg", var.c_str() ));
  h_signal_data_nobkg->AddMissingErrorBandsAndFillWithCV( *h_signal_bkg );
  if (doBkgSub ) h_signal_data_nobkg->Add( h_signal_bkg, -1 );

  vector<MnvH2D*> vec_data = hdl->Get2DMnvHistos( h_signal_data_nobkg, true );
  vector<MnvH2D*> vec_mc = hdl->Get2DMnvHistos( h_signal_mc_nobkg, true );
  vector<MnvH2D*> vec_mc_qelike = hdl->Get2DMnvHistos( h_signal[kQELike], true);
  vector<MnvH2D*> vec_mc_qelike_qe = hdl->Get2DMnvHistos( h_signal[kQELike_QE], true);
  vector<MnvH2D*> vec_mc_qelike_qe_h = hdl->Get2DMnvHistos( h_signal[kQELike_QE_H], true);
  vector<MnvH2D*> vec_mc_qelike_qe_oth = hdl->Get2DMnvHistos( h_signal[kQELike_QE_OTH], true);
  vector<MnvH2D*> vec_mc_qelike_res = hdl->Get2DMnvHistos( h_signal[kQELike_RES], true);
  vector<MnvH2D*> vec_mc_qelike_dis = hdl->Get2DMnvHistos( h_signal[kQELike_DIS], true);
  vector<MnvH2D*> vec_mc_qelike_2p2h = hdl->Get2DMnvHistos( h_signal[kQELike_2p2h], true);
  vector<MnvH2D*> vec_mc_qelikenot = hdl->Get2DMnvHistos( h_signal[kQELikeNot], true);
  //MnvH2D* h_data =            ProjectionXZ( vec_data, var+"_data" );
  //MnvH2D* h_mc =              ProjectionXZ( vec_mc, var+"_mc");
  //MnvH2D* h_mc_qelike =       ProjectionXZ( vec_mc_qelike, var+"_qelike");
  //MnvH2D* h_mc_qelike_qe =    ProjectionXZ( vec_mc_qelike_qe, var+"_qelike_qe");
  //MnvH2D* h_mc_qelike_res =   ProjectionXZ( vec_mc_qelike_res, var+"_qelike_res");
  //MnvH2D* h_mc_qelike_dis =   ProjectionXZ( vec_mc_qelike_dis, var+"_qelike_dis");
  //MnvH2D* h_mc_qelike_2p2h =  ProjectionXZ( vec_mc_qelike_2p2h, var+"_qelike_2p2h");
  //MnvH2D* h_mc_qelikenot =  ProjectionXZ( vec_mc_qelikenot, var+"_qelikenot");
  //
  cout<<"ProjectionHDL"<<endl;
  string dims = "xz";
  MnvH2D* h_data =            ProjectionHDL( vec_data, var+"_data",hdl, xbins, ybins,zbins,dims );
  MnvH2D* h_mc =              ProjectionHDL( vec_mc, var+"_mc", hdl, xbins, ybins,zbins,dims );
  MnvH2D* h_mc_qelike =       ProjectionHDL( vec_mc_qelike, var+"_qelike", hdl, xbins, ybins,zbins,dims );
  MnvH2D* h_mc_qelike_qe =    ProjectionHDL( vec_mc_qelike_qe_h, var+"_qelike_qe", hdl, xbins, ybins,zbins,dims );
  MnvH2D* h_mc_qelike_qe_h =    ProjectionHDL( vec_mc_qelike_qe_h, var+"_qelike_qe_h", hdl, xbins, ybins,zbins,dims );
  MnvH2D* h_mc_qelike_qe_oth =    ProjectionHDL( vec_mc_qelike_qe_oth, var+"_qelike_qe_oth", hdl, xbins, ybins,zbins,dims );
  MnvH2D* h_mc_qelike_res =   ProjectionHDL( vec_mc_qelike_res, var+"_qelike_res", hdl, xbins, ybins,zbins,dims );
  MnvH2D* h_mc_qelike_dis =   ProjectionHDL( vec_mc_qelike_dis, var+"_qelike_dis", hdl, xbins, ybins,zbins,dims );
  MnvH2D* h_mc_qelike_2p2h =  ProjectionHDL( vec_mc_qelike_2p2h, var+"_qelike_2p2h", hdl, xbins, ybins,zbins,dims );
  MnvH2D* h_mc_qelikenot =    ProjectionHDL( vec_mc_qelikenot, var+"_qelikenot", hdl, xbins, ybins,zbins,dims );




  cout<<"h_signal_bkg has "<<h_signal_bkg->GetVertErrorBandNames().size()<<" VertErrBand"<<endl;
  cout<<"h_signal_data_nobkg has "<<h_signal_data_nobkg->GetVertErrorBandNames().size()<<" VertErrBand"<<endl;
  cout<<"h_signal_mc_nobkg has "<<h_signal_mc_nobkg->GetVertErrorBandNames().size()<<" VertErrBand"<<endl;

  cout<<"h_data has "<<vec_data[0]->GetVertErrorBandNames().size()<<" VertErrBand"<<endl;
  cout<<"h_mc has "<<vec_mc[0]->GetVertErrorBandNames().size()<<" VertErrBand"<<endl;

  cout<<"h_mc_qelike_qe_h has "<<h_mc_qelike_qe_h->GetVertErrorBandNames().size()<<endl;
  cout<<"vec_mc_qelike_qe has "<<vec_mc_qelike_qe[0]->GetVertErrorBandNames().size()<<endl;
  cout<<"h_signal[kQELike_QE] has "<<h_signal[kQELike_QE]->GetVertErrorBandNames().size()<<endl;
  cout<<"h_signal[kQELike_RES] has "<<h_signal[kQELike_RES]->GetVertErrorBandNames().size()<<endl;

  MnvH2D* h_mc_qelike_notqe = (MnvH2D*) h_mc->Clone( (var+"_qelike_notqe").c_str() );
  h_mc_qelike_notqe->Add( h_mc_qelike, -1 );

  // get ratio in dthetaP, scaling all mc distribution
  MnvH1D* h_ratio_data_mc = (MnvH1D*) h_data->ProjectionX( (var+"_ratio_data_mc").c_str() ); 
  MnvH2D* h_ratio_data_mc_2D = (MnvH2D*) h_data->Clone( (var+"_ratio_data_mc_2D").c_str() );
  MnvH2D* h_mc_weighted = (MnvH2D*) h_mc->Clone( (var+"_mc_weighted").c_str() );


  h_ratio_data_mc->Divide( h_ratio_data_mc, h_mc->ProjectionX() );
  ExpandHisto( h_ratio_data_mc, h_ratio_data_mc_2D );
  h_mc_weighted->Multiply( h_mc_weighted, h_ratio_data_mc_2D );


  // get ratio in dthetaP, scaling only qelike QE distribution
  MnvH1D* h_ratio_data_mc_qelike_qe = (MnvH1D*) h_data->ProjectionX( (var+"_ratio_data_mc_qelike_qe").c_str() );
  MnvH2D* h_ratio_data_mc_qelike_qe_2D = (MnvH2D*) h_data->Clone( (var+"_ratio_data_mc_qelike_qe_2D").c_str() );
  MnvH2D* h_mc_weighted_qelike_qe= (MnvH2D*) h_mc_qelike_qe->Clone( (var+"_mc_weighted_qelike_qe").c_str() );
  MnvH2D* h_mc_weighted_qelike_all = (MnvH2D*) h_mc_qelike_notqe->Clone( (var+"_mc_weighted_qelike_all").c_str() );

  h_ratio_data_mc_qelike_qe->Add( h_mc_qelike_notqe->ProjectionX(), -1 ); // get only QelikeQE
  h_ratio_data_mc_qelike_qe->Divide( h_ratio_data_mc_qelike_qe,  h_mc_qelike->ProjectionX() ); // Divide data-noqe by mc_qelike
  ExpandHisto( h_ratio_data_mc_qelike_qe, h_ratio_data_mc_qelike_qe_2D );//expanding weights into 2D histos
 
  h_mc_weighted_qelike_qe->Multiply( h_mc_weighted_qelike_qe, h_ratio_data_mc_qelike_qe_2D );

  h_mc_weighted_qelike_all->Add( h_mc_weighted_qelike_qe );

  double multipliers[]={ 80, 32, 8, 4, 
                         1 , 1 , 1, 1,
                         1, 1, 2, 4,
                         40};



  //Saving unweighted data, MC into TH2D
  TH2* dataStat = new TH2D( h_data->GetCVHistoWithStatError() );
  TH2* data =     new TH2D( h_data->GetCVHistoWithError() );
  TH2* mc  =      new TH2D( h_mc->GetCVHistoWithStatError() );
  TH2* mc_qelike_qe_h= new TH2D( h_mc_qelike_qe_h->GetCVHistoWithStatError() );
  TH2* mc_qelike_qe_oth= new TH2D( h_mc_qelike_qe_oth->GetCVHistoWithStatError() );
  TH2* mc_qelike_res= new TH2D( h_mc_qelike_res->GetCVHistoWithStatError() );
  TH2* mc_qelike_dis= new TH2D( h_mc_qelike_dis->GetCVHistoWithStatError() );
  TH2* mc_qelike_2p2h= new TH2D( h_mc_qelike_2p2h->GetCVHistoWithStatError() );

  TH2* mc_equalweighted = new TH2D( h_mc_weighted->GetCVHistoWithStatError() );
  TH2* mc_weighted_qelike_qe = new TH2D( h_mc_weighted_qelike_qe->GetCVHistoWithStatError() );
  TH2* mc_weighted_qelike_all = new TH2D( h_mc_weighted_qelike_all->GetCVHistoWithStatError() );

  
  // ----------------------------------------------------------------------------------
  // First make dthetaP in bins of q2qe
  // ----------------------------------------------------------------------------------

  // scaling histos by width
  string savename = "plots/nu-2d-evtrate-model-q2qe-"+var+"-"+sideband;
  string ext = ".pdf";
  string doBkg=(doBkgSub)? "_bkgsub": "";
  string xtitle = "#delta#theta_{P}", ytitle = "Event Rate (#times 10^{-3}) per Degree";
  if( var == "dthetaRq2qe" ) xtitle = "#delta#theta_{R}";
  string celltitle = "Q^{2}_{QE}";
  double scale = 1e-3;
  double ymin = 0, ymax = 2.4;

  drawPlots(savename+doBkg+ext, data,  dataStat,  mc,  mc_qelike_qe_h, mc_qelike_qe_oth,  mc_qelike_res,  mc_qelike_dis,  mc_qelike_2p2h, xtitle, ytitle, celltitle, "x", scale, ymin,ymax, multipliers, false );

  savename =  "plots/nu-2d-evtrate-ratio-model-q2qe-"+var+"-"+sideband;
  ytitle = "Ratio to MnvTune";
  drawPlots(savename+doBkg+ext, data,  dataStat,  mc,  mc_qelike_qe_h, mc_qelike_qe_oth,  mc_qelike_res,  mc_qelike_dis,  mc_qelike_2p2h, xtitle, ytitle, celltitle, "x", scale, -.4,2.4, NULL, true );
  //gc->Print("plots/nu-2d-evtrate-model-q2qe-.png");
  //

  //draw Q2 Distributions
  TH1* data_qe = (TH1*) data->ProjectionY("h_data_qe");
  TH1* dataStat_qe = (TH1*) dataStat->ProjectionY("h_dataStat_qe");
  TH1* mc_qe = (TH1*) mc->ProjectionY("h_mc_qe");
  TH1* mc_qe_equalweighted = (TH1*) mc_equalweighted->ProjectionY("h_mc_equalweighted");

  data_qe->Scale(1,"width"); dataStat_qe->Scale(1,"width");
  mc_qe->Scale(1,"width"); mc_qe_equalweighted->Scale(1,"width");

  mc_qe->SetLineColor(kBlue);
  mc_qe_equalweighted->SetLineColor(kRed);
  TCanvas* c = new TCanvas("c", "c", 800,800 );

  mc_qe->GetXaxis()->SetTitle("Q^{2}_{QE} (GeV^{2})");
  mc_qe->Draw("hist");
  mc_qe_equalweighted->Draw("hist same");
  data_qe->Draw("pe same");
  dataStat_qe->Draw("pe same");
  TLegend* leg=new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->AddEntry( data_qe, "Data" );
  leg->AddEntry( mc_qe, "MC" );
  leg->AddEntry( mc_qe_equalweighted, "MC #delta#theta_{P} - weighted" );
  leg->Draw();
  c->Print(Form("plots/q2qe_%s.pdf", sideband.c_str() ));

  
  //Draw dthetaP weights
  h_ratio_data_mc_qelike_qe->Draw("pe");
  c->Print("plots/dthetaP_weights.pdf");


  //Now let's reweight in Q2qe bins

  double multipliers2[]={ 40, 40, 16, 8, 4,4, 2, 1,
                          1, 2, 4, 4, 8, 16,40,40};



  // scaling histos by width
  savename = "plots/nu-2d-evtrate-model-"+sideband+"-"+var;
  if(doCenter) savename+="-center";
  xtitle = "Q^{2}_{QE} (GeV^{2})"; ytitle = "Event Rate (#times 10^{-3}) per GeV^{2}";
  celltitle = "#delta#theta_{P}";
  scale = 1e-3;
  ymin = 0; ymax = 2.8;

  drawPlots(savename+doBkg+ext, data,  dataStat,  mc,  mc_qelike_qe_h, mc_qelike_qe_oth,  mc_qelike_res,  mc_qelike_dis,  mc_qelike_2p2h, xtitle, ytitle, celltitle, "y", scale, ymin,ymax, multipliers2, false,true );

  savename =  "plots/nu-2d-evtrate-ratio-model-"+sideband+"-"+var; 
  if(doCenter) savename+="-center";
  ytitle = "Ratio to MnvTune";
  drawPlots(savename+doBkg+ext, data,  dataStat,  mc,  mc_qelike_qe_h, mc_qelike_qe_oth,  mc_qelike_res,  mc_qelike_dis,  mc_qelike_2p2h, xtitle, ytitle, celltitle, "y", scale, -.4,2.4, NULL, true,true );
  //gc->Print("plots/nu-2d-evtrate-model-q2qe-.png");
  //

  savename =  "plots/nu-2d-evtrate-model-equalweighted-"+sideband+"-"+var;
  if(doCenter) savename+="-center";
  ytitle = "Event Rate (#times 10^{-3}) per GeV^{2}";
  drawPlots(savename+doBkg+ext, data,  dataStat,  mc_equalweighted,  mc_qelike_qe_h, mc_qelike_qe_oth,  mc_qelike_res,  mc_qelike_dis,  mc_qelike_2p2h, xtitle, ytitle, celltitle, "y", scale, ymin,ymax, multipliers2, false ,true);

  savename =  "plots/nu-2d-evtrate-ratio-model-equalweighted-"+sideband+"-"+var;
  if(doCenter) savename+="-center";
  ytitle =  "Ratio to MnvTune (EqualWeight)";
  drawPlots(savename+doBkg+ext, data,  dataStat,  mc_equalweighted,  mc_qelike_qe_h, mc_qelike_qe_oth,  mc_qelike_res,  mc_qelike_dis,  mc_qelike_2p2h, xtitle, ytitle, celltitle, "y", scale, -.4,2.4, NULL, true ,true);


  //savename =  "plots/nu-2d-evtrate-model-q2qe-qeweighted.pdf";
  //ytitle =  "Event Rate (#times 10^{3}) per Degree";
  //drawPlots(savename, data,  dataStat,  mc_weighted_qelike_all,  mc_weighted_qelike_qe,  mc_qelike_res,  mc_qelike_dis,  mc_qelike_2p2h, xtitle, ytitle, celltitle, "x", scale, ymin,ymax, multipliers, false );

  //savename =  "plots/nu-2d-evtrate-ratio-model-q2qe-qeweighted.pdf";
  //ytitle =  "Ratio to MnvTune (EqualWeight)";
  //drawPlots(savename, data,  dataStat,  mc_weighted_qelike_all,  mc_weighted_qelike_qe,  mc_qelike_res,  mc_qelike_dis,  mc_qelike_2p2h, xtitle, ytitle, celltitle, "x", scale, 0,2, multipliers, true );
}

void makePlot1D(string fname_signal, string fname_bkg_weight, string sideband="Signal", string var = "dthetaP", string varname = "#delta#theta_{P} (degree)", bool doCenter=false)
{
  cout<<"makePlot1D"<<endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  //three files 2-track, 2-track bkg_fitted
  //CV
  TFile* f1 = new TFile( fname_signal.c_str(), "READ");// original signal histogram
  TFile* f2 = new TFile( fname_bkg_weight.c_str(), "READ");// weighted histogram
  

  MnvH2D *h_qelikenot_weights = (MnvH2D*) f2->Get( Form("h_weights_%sptbins_qelikenot_%s", var.c_str(),sideband.c_str() ));
  cout<<"Read h_qelikenot_weights? "<<(h_qelikenot_weights!=0)<<endl;

  TVector2 *pot = (TVector2*)f1->Get("pot;1");
  //for( int i = 2; i<8; i++ ) (*pot)+= *((TVector2*)f1->Get(Form( "pot;%d", i) ));
  double pot_data = pot->X();
  double pot_mc = pot->Y();
  double pot_norm = pot_data/pot_mc;

  bool fluxHistoExists = true;
  CCQENuPlotUtils *putils = new CCQENuPlotUtils( fluxHistoExists );
  CCQENuUtils *utils = new CCQENuUtils( false, fluxHistoExists );

  MnvH2D *h_var_ptmu[nHistos];
  if(doCenter) utils->bookHistos( f1, h_var_ptmu, Form("h_%s_ptmu_center", var.c_str()) );
  else utils->bookHistos( f1, h_var_ptmu, Form("h_%s_ptmu", var.c_str()) );
  for( int i = 0; i< nHistos; i++ )
  {
    h_var_ptmu[i]->GetXaxis()->SetTitle( varname.c_str());
    h_var_ptmu[i]->GetYaxis()->SetTitle("#it{p}_{T}^{#mu} (GeV/#it{c})");
  }
  cout<<"bookHistos done"<<endl;
  putils->scaleMCHistos( h_var_ptmu, pot_norm );

  string ext = ".pdf";
  string savename = "plots/nu-1d-evtrate-"+sideband+"-ptmu";
  draw1DPlots(putils, savename+ext, h_var_ptmu, pot_data,  pot_mc,0,-1, "y");
  savename = "plots/nu-1d-evtrate-"+sideband+"-"+var;
  draw1DPlots(putils, savename+ext, h_var_ptmu, pot_data,  pot_mc,0,-1, "x");

  h_var_ptmu[kQELikeNot]->Multiply( h_var_ptmu[kQELikeNot], h_qelikenot_weights );
  h_var_ptmu[kQELikeNot_SingleChargedPion]->Multiply( h_var_ptmu[kQELikeNot_SingleChargedPion], (MnvH2D*) f2->Get(Form("hs_weights_%sptbins_bgType_SingleChargedPion", var.c_str() ) ) );
  h_var_ptmu[kQELikeNot_SingleNeutralPion]->Multiply( h_var_ptmu[kQELikeNot_SingleNeutralPion], (MnvH2D*) f2->Get(Form("hs_weights_%sptbins_bgType_SingleNeutralPion", var.c_str() ) ) );
  h_var_ptmu[kQELikeNot_MultiPion]->Multiply( h_var_ptmu[kQELikeNot_MultiPion], (MnvH2D*) f2->Get(Form("hs_weights_%sptbins_bgType_MultiPions", var.c_str() ) ) );

  savename = "plots/nu-1d-evtrate-weighted-"+sideband+"-ptmu";
  draw1DPlots(putils, savename+ext, h_var_ptmu, pot_data,  pot_mc,0,-1, "y");
  savename = "plots/nu-1d-evtrate-weighted-"+sideband+"-"+var;
  draw1DPlots(putils, savename+ext, h_var_ptmu, pot_data,  pot_mc,0,-1, "x");

  savename = "plots/nu-1d-evtrate-weighted-nobck-"+sideband+"-ptmu";
  draw1DPlots(putils, savename+ext, h_var_ptmu, pot_data,  pot_mc,0,-1, "y", true);
  savename = "plots/nu-1d-evtrate-weighted-nobck-"+sideband+"-"+var;
  draw1DPlots(putils, savename+ext, h_var_ptmu, pot_data,  pot_mc,0,-1, "x", true);




  MnvPlotter *plotter = new MnvPlotter;
  savename = "plots/nu-1d-evtrate-ptmu-qelikenot-weight.pdf";
  TCanvas* c = new TCanvas("c","c", 1200,800 );
  MnvH1D* h_qelikenot_weights1D = (MnvH1D*) h_qelikenot_weights->ProjectionY("h_qelikenot_weights1D",1,1);
  h_qelikenot_weights1D->GetXaxis()->SetTitle("p_{T}^{#mu} (GeV/#it{c})" );
  h_qelikenot_weights1D->GetYaxis()->SetTitle( "MC Weights" );

  plotter->DrawMCWithErrorBand( h_qelikenot_weights1D );
  c->Print(savename.c_str() );
}

void drawPlots(string savename, TH2* data, TH2* dataStat, TH2* mc, TH2* mc_qelike_qe_h, TH2* mc_qelike_qe_oth, TH2* mc_qelike_res, TH2* mc_qelike_dis, TH2* mc_qelike_2p2h, string xtitle, string ytitle, string celltitle, string axis, double scale , double ymin, double ymax, double * multipliers, bool ratio , bool setlogx  , bool setlogy )
{
  //Plotting unweighed data/MC ratio
  vector<int> mycolors = getColors(2);
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);

  //need to add signal and bkg colors
  mc_qelike_qe_h->SetLineColor( mycolors[2] );
  mc_qelike_qe_oth->SetLineColor( mycolors[3] );
  mc_qelike_res->SetLineColor( mycolors[4] );
  mc_qelike_dis->SetLineColor( mycolors[5] );
  mc_qelike_2p2h->SetLineColor( mycolors[6] );

  data->SetMarkerStyle( 1 );
  data->SetMarkerSize(0.2);
  data->SetLineColor(kBlack);
  data->SetLineWidth(1);

  dataStat->SetMarkerStyle(1);
  dataStat->SetMarkerSize(0.2);
  dataStat->SetLineColor(kBlack);
  dataStat->SetLineWidth(1);

  //set options
  std::vector<std::pair<TH2*, const char*> > histAndOpts;


  histAndOpts.push_back(std::make_pair( (TH2*) dataStat->Clone("dataStat"), "e"));
  histAndOpts.push_back(std::make_pair( (TH2*) data->Clone("data"),     "ep"));
  histAndOpts.push_back(std::make_pair( (TH2*) mc->Clone("mc"),       "hist"));
  histAndOpts.push_back(std::make_pair( (TH2*) mc_qelike_qe_h->Clone("mc_qqeh"),       "histl"));
  histAndOpts.push_back(std::make_pair( (TH2*) mc_qelike_qe_oth->Clone("mc_qqeoth"),   "histl"));
  histAndOpts.push_back(std::make_pair( (TH2*) mc_qelike_res->Clone("mc_qres"),       "histl"));
  histAndOpts.push_back(std::make_pair( (TH2*) mc_qelike_dis->Clone("mc_qdis"),       "histl"));
  histAndOpts.push_back(std::make_pair( (TH2*) mc_qelike_2p2h->Clone("mc_q2p2h"),       "histl"));


  //for( int i = 0; i< histAndOpts.size();i++) histAndOpts[i].first->Scale(1,"width");

  TH2* mcRatio;
  if(ratio) mcRatio = (TH2*) mc->Clone("mcRatio");
 
  for( unsigned int i = 0 ; i< histAndOpts.size(); i++ ) 
  {
    if(ratio) histAndOpts[i].first->Divide(histAndOpts[i].first ,mcRatio);
    //else histAndOpts[i].first->Scale(scale, "width");
  }
  GridCanvas* gc = (axis == "x")? plotXAxis1D( histAndOpts, xtitle, celltitle, multipliers ): plotYAxis1D( histAndOpts, xtitle, celltitle, multipliers ) ;
  gc->SetLogx(setlogx);
  gc->SetLogy(setlogy);
  gc->SetGridx(true);
  gc->SetYLimits(ymin, ymax);
  gc->SetYTitle(ytitle.c_str());
  gc->Modified();


  TLegend* leg=new TLegend(0.7, 0.1, 0.9, 0.3);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(data, "MINERvA data", "lpe");
  leg->AddEntry(mc, "MINERvA Tune", "l");
  leg->AddEntry(mc_qelike_qe_h,"QE-H","l");
  leg->AddEntry(mc_qelike_qe_oth,"QE-Oth","l");
  leg->AddEntry(mc_qelike_res,"Resonant","l");
  leg->AddEntry(mc_qelike_dis,"DIS","l");
  leg->AddEntry(mc_qelike_2p2h,"2p2h","l");
  //leg->AddEntry(mc_qelike_2p2h_no_lowrec,"2p2h without fit","l");
  leg->Draw("SAME");

  gc->Print(savename.c_str());
}

void draw1DPlots(CCQENuPlotUtils* utils, string savename, MnvH2D** hists, double pot_data,  double pot_mc, int bmin, int bmax, string axis, bool bkSub)
{
  MnvH1D* h[nHistos];
  for( int i = 0; i< nHistos; i++ )
  {
    if( axis == "x") h[i]= (MnvH1D*) hists[i]->ProjectionX( Form( "%s_projX", hists[i]->GetName() ), bmin, bmax );
    else h[i]= (MnvH1D*) hists[i]->ProjectionY( Form( "%s_projY", hists[i]->GetName() ), bmin, bmax );
  }
  //mainly stacked plots
  //

  if(bkSub)
  {
    h[kData]->Add( h[kQELikeNot], -1 );
    h[kQELikeNot]->Reset();
    h[kQELikeNot_NoPions]->Reset();
    h[kQELikeNot_SingleNeutralPion]->Reset();
    h[kQELikeNot_SingleChargedPion]->Reset();
    h[kQELikeNot_MultiPion]->Reset();
  }
  string xlabel = ( axis == "x")? "#delta#theta_{P} (Degree)": "#it{p}_{T}^{#mu} (GeV/#it{c})";
  string ylabel = ( axis == "x")? "Event Rate (Degree)": "Event Rate (GeV/#it{c})";
  TCanvas* c = new TCanvas("c","c", 1200,800 );
  c->cd();
  cout<<"Start DrawStacked"<<endl;
  utils->drawStacked( h, "QELike_split_PionInFS", false, pot_data, pot_mc, true, -1, xlabel,ylabel );
  cout<<"Start Print"<<endl;
  c->Print(savename.c_str());

  //utils->Print( c, savename );
}


int main(int argc, char* argv[])
{

  string f_orig = argv[1];
  string f_bkg_weighted = argv[2];
  string sideband = argv[3];
  //multipliers
  string var = "dthetaPq2qe";
  //makePlots(true,false,f_orig, f_bkg_weighted, var, sideband);
  //makePlots(true,true,f_orig, f_bkg_weighted, var, sideband);

  //makePlots(true,false,f_orig, f_bkg_weighted, var, sideband, true);
  //makePlots(true,true,f_orig, f_bkg_weighted, var, sideband, true);

  var = "dthetaRq2qe";
  //makePlots(true,false,f_orig, f_bkg_weighted, var, sideband);
  //makePlots(true,true,f_orig, f_bkg_weighted, var, sideband);

  ////dthetaP and ptmu
  makePlot1D( f_orig, f_bkg_weighted, sideband,"dthetaP", "#delta#theta_{P} (degree)");
  makePlot1D( f_orig, f_bkg_weighted, sideband,"dthetaR", "#delta#theta_{R} (degree)");
  return 0;
}
