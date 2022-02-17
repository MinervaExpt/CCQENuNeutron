#include "include/CCQENuPlotUtils.h"
#include "include/CCQENuUtils.h"
#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/TargetUtils.h"

using namespace CCQENU_ANA;

void ScaleBackground(MnvH2D** h, MnvH2D* h_bg_weights, MnvH2D* &h_data_nobck, MnvH2D* &h_mc_nobck, MnvH2D*& h_tuned_bkg);//This is the new sideband method


int SideBandFitHists( string output_filename, string filename_1track_signal,string filename_2track_signal, string filename_1track_blobs,string filename_2track_blobs, string filename_1track_michel,string filename_2track_michel, string filename_bg_weights, bool makeFluxConstraintHisto, bool applyFluxConstraint)
{
  // read file to get histograms
  TFile *f_1track_signal = new TFile( filename_1track_signal.c_str(), "READ" );
  if (f_1track_signal->IsZombie() || f_1track_signal->GetListOfKeys()->IsEmpty()){
    Error("SideBandFitPlots","Could not get histogram ROOT file or it was empty.");
    return 1;
  }
  TFile *f_2track_signal = new TFile( filename_2track_signal.c_str(), "READ" );
  if (f_2track_signal->IsZombie() || f_2track_signal->GetListOfKeys()->IsEmpty()){
    Error("SideBandFitPlots","Could not get histogram ROOT file or it was empty.");
    return 1;
  }


  TFile *f_1track_blobs = new TFile( filename_1track_blobs.c_str(), "READ" );
  if (f_1track_blobs->IsZombie() || f_1track_blobs->GetListOfKeys()->IsEmpty()){
    Error("SideBandFitPlots","Could not get histogram ROOT file or it was empty.");
    return 1;
  }
  TFile *f_2track_blobs = new TFile( filename_2track_blobs.c_str(), "READ" );
  if (f_2track_blobs->IsZombie() || f_2track_blobs->GetListOfKeys()->IsEmpty()){
    Error("SideBandFitPlots","Could not get histogram ROOT file or it was empty.");
    return 1;
  }


  TFile *f_1track_michel = new TFile( filename_1track_michel.c_str(), "READ" );
  if (f_1track_michel->IsZombie() || f_1track_michel->GetListOfKeys()->IsEmpty()){
    Error("SideBandFitPlots","Could not get histogram ROOT file or it was empty.");
    return 1;
  }
  TFile *f_2track_michel = new TFile( filename_2track_michel.c_str(), "READ" );
  if (f_2track_michel->IsZombie() || f_2track_michel->GetListOfKeys()->IsEmpty()){
    Error("SideBandFitPlots","Could not get histogram ROOT file or it was empty.");
    return 1;
  }

  TFile *f_bg_weights = new TFile( filename_bg_weights.c_str(), "READ" );
  if (f_bg_weights->IsZombie() || f_bg_weights->GetListOfKeys()->IsEmpty()){
    Error("SideBandFitPlots","Could not get histogram ROOT file or it was empty.");
    return 1;
  }
  CCQENuUtils *utils = new CCQENuUtils( false, makeFluxConstraintHisto );
  utils->setPlaylist("minerva13");//Dominant LE playlist. Need to figure out if there is some correction needed here...
  utils->setFluxReweighterPlaylist();
  CCQENuPlotUtils *plotUtils = new CCQENuPlotUtils();

  //-------------------------------------------
  // Grab event counts before cuts for norm
  //-------------------------------------------
  TVector2 *evt = (TVector2*)f_1track_signal->Get("n_events");
  double data_events = evt->X();
  double mc_events   = evt->Y();
  
  cout<< "Number of Data Events (1-track) = " << data_events << endl;
  cout<< "Number of MC Events (1-track) = " <<   mc_events << endl;
  
  evt = (TVector2*)f_2track_signal->Get("n_events");
  data_events = evt->X();
  mc_events   = evt->Y();
  
  cout<< "Number of Data Events (2-track) = " << data_events << endl;
  cout<< "Number of MC Events (2-track) = " <<   mc_events << endl;

  //--------------------------------------------
  // Load histos to plot 
  //--------------------------------------------
  //2-D
  //signal
  MnvH2D *h_pzmu_ptmu_1track_signal[nHistos], *h_pzmu_ptmu_2track_signal[nHistos];
  MnvH2D *h_enu_ptmu_1track_signal[nHistos], *h_enu_ptmu_2track_signal[nHistos];
  MnvH2D *h_q2_ptmu_1track_signal[nHistos], *h_q2_ptmu_2track_signal[nHistos];
  MnvH2D *h_recoil_ptmu_1track_signal[nHistos], *h_recoil_ptmu_2track_signal[nHistos];
  MnvH2D *h_recoil_inc_ptmu_1track_signal[nHistos], *h_recoil_inc_ptmu_2track_signal[nHistos];

  //Blobs
  MnvH2D *h_pzmu_ptmu_1track_blobs[nHistos], *h_pzmu_ptmu_2track_blobs[nHistos];
  MnvH2D *h_enu_ptmu_1track_blobs[nHistos], *h_enu_ptmu_2track_blobs[nHistos];
  MnvH2D *h_q2_ptmu_1track_blobs[nHistos], *h_q2_ptmu_2track_blobs[nHistos];
  MnvH2D *h_recoil_ptmu_1track_blobs[nHistos], *h_recoil_ptmu_2track_blobs[nHistos];
  MnvH2D *h_recoil_inc_ptmu_1track_blobs[nHistos], *h_recoil_inc_ptmu_2track_blobs[nHistos];

  //Michel
  MnvH2D *h_pzmu_ptmu_1track_michel[nHistos], *h_pzmu_ptmu_2track_michel[nHistos];
  MnvH2D *h_enu_ptmu_1track_michel[nHistos], *h_enu_ptmu_2track_michel[nHistos];
  MnvH2D *h_q2_ptmu_1track_michel[nHistos], *h_q2_ptmu_2track_michel[nHistos];
  MnvH2D *h_recoil_ptmu_1track_michel[nHistos], *h_recoil_ptmu_2track_michel[nHistos];
  MnvH2D *h_recoil_inc_ptmu_1track_michel[nHistos], *h_recoil_inc_ptmu_2track_michel[nHistos];


  //signal
  plotUtils->bookHistos( f_1track_signal, h_pzmu_ptmu_1track_signal, "h_pzmu_ptmu");
  plotUtils->bookHistos( f_2track_signal, h_pzmu_ptmu_2track_signal, "h_pzmu_ptmu");
  plotUtils->bookHistos( f_1track_signal, h_q2_ptmu_1track_signal, "h_q2_ptmu");
  plotUtils->bookHistos( f_2track_signal, h_q2_ptmu_2track_signal, "h_q2_ptmu");
  plotUtils->bookHistos( f_1track_signal, h_enu_ptmu_1track_signal, "h_enu_ptmu");
  plotUtils->bookHistos( f_2track_signal, h_enu_ptmu_2track_signal, "h_enu_ptmu");
  // plotUtils->bookHistos( f_1track_signal, h_recoil_ptmu_1track_signal,"h_recoil_ptmu");
  // plotUtils->bookHistos( f_2track_signal, h_recoil_ptmu_2track_signal,"h_recoil_ptmu");
  // plotUtils->bookHistos( f_1track_signal, h_recoil_inc_ptmu_1track_signal,"h_recoil_inc_ptmu");
  // plotUtils->bookHistos( f_2track_signal, h_recoil_inc_ptmu_2track_signal,"h_recoil_inc_ptmu");
  //blobs
  plotUtils->bookHistos( f_1track_blobs, h_pzmu_ptmu_1track_blobs, "h_pzmu_ptmu");
  plotUtils->bookHistos( f_2track_blobs, h_pzmu_ptmu_2track_blobs, "h_pzmu_ptmu");
  plotUtils->bookHistos( f_1track_blobs, h_q2_ptmu_1track_blobs, "h_q2_ptmu");
  plotUtils->bookHistos( f_2track_blobs, h_q2_ptmu_2track_blobs, "h_q2_ptmu");
  plotUtils->bookHistos( f_1track_blobs, h_enu_ptmu_1track_blobs, "h_enu_ptmu");
  plotUtils->bookHistos( f_2track_blobs, h_enu_ptmu_2track_blobs, "h_enu_ptmu");
  // plotUtils->bookHistos( f_1track_blobs, h_recoil_ptmu_1track_blobs,"h_recoil_ptmu");
  // plotUtils->bookHistos( f_2track_blobs, h_recoil_ptmu_2track_blobs,"h_recoil_ptmu");
  // plotUtils->bookHistos( f_1track_blobs, h_recoil_inc_ptmu_1track_blobs,"h_recoil_inc_ptmu");
  // plotUtils->bookHistos( f_2track_blobs, h_recoil_inc_ptmu_2track_blobs,"h_recoil_inc_ptmu");
  //michel
  plotUtils->bookHistos( f_1track_michel, h_pzmu_ptmu_1track_michel, "h_pzmu_ptmu");
  plotUtils->bookHistos( f_2track_michel, h_pzmu_ptmu_2track_michel, "h_pzmu_ptmu");
  plotUtils->bookHistos( f_1track_michel, h_q2_ptmu_1track_michel, "h_q2_ptmu");
  plotUtils->bookHistos( f_2track_michel, h_q2_ptmu_2track_michel, "h_q2_ptmu");
  plotUtils->bookHistos( f_1track_michel, h_enu_ptmu_1track_michel, "h_enu_ptmu");
  plotUtils->bookHistos( f_2track_michel, h_enu_ptmu_2track_michel, "h_enu_ptmu");
  // plotUtils->bookHistos( f_1track_michel, h_recoil_ptmu_1track_michel,"h_recoil_ptmu");
  // plotUtils->bookHistos( f_2track_michel, h_recoil_ptmu_2track_michel,"h_recoil_ptmu");
  // plotUtils->bookHistos( f_1track_michel, h_recoil_inc_ptmu_1track_michel,"h_recoil_inc_ptmu");
  // plotUtils->bookHistos( f_2track_michel, h_recoil_inc_ptmu_2track_michel,"h_recoil_inc_ptmu");



  //2D BG Weights
  //pz
  //Signal
  MnvH2D *h_pzmu_ptmu_bg_weights_1track = (MnvH2D*) f_bg_weights->Get("h_weights_1track_pzptbins_qelikenot");
  MnvH2D *h_pzmu_ptmu_bg_weights_2track = (MnvH2D*) f_bg_weights->Get("h_weights_2track_pzptbins_qelikenot");

  MnvH2D *h_q2_ptmu_bg_weights_1track = (MnvH2D*) f_bg_weights->Get("h_weights_1track_q2ptbins_qelikenot");
  MnvH2D *h_q2_ptmu_bg_weights_2track = (MnvH2D*) f_bg_weights->Get("h_weights_2track_q2ptbins_qelikenot");

  MnvH2D *h_enu_ptmu_bg_weights_1track = (MnvH2D*) f_bg_weights->Get("h_weights_1track_enuptbins_qelikenot");
  MnvH2D *h_enu_ptmu_bg_weights_2track = (MnvH2D*) f_bg_weights->Get("h_weights_2track_enuptbins_qelikenot");

  // MnvH2D *h_recoil_ptmu_bg_weights_1track = (MnvH2D*) f_bg_weights->Get("h_weights_1track_recoilptbins_qelikenot");
  // MnvH2D *h_recoil_ptmu_bg_weights_2track = (MnvH2D*) f_bg_weights->Get("h_weights_2track_recoilptbins_qelikenot");

  // MnvH2D *h_recoil_inc_ptmu_bg_weights_1track = (MnvH2D*) f_bg_weights->Get("h_weights_1track_recoilincptbins_qelikenot");
  // MnvH2D *h_recoil_inc_ptmu_bg_weights_2track = (MnvH2D*) f_bg_weights->Get("h_weights_2track_recoilincptbins_qelikenot");

  //blobs
  MnvH2D *h_pzmu_ptmu_bg_weights_1track_blobs = (MnvH2D*) f_bg_weights->Get("h_weights_1track_pzptbins_qelikenot_blobs");
  MnvH2D *h_pzmu_ptmu_bg_weights_2track_blobs = (MnvH2D*) f_bg_weights->Get("h_weights_2track_pzptbins_qelikenot_blobs");

  MnvH2D *h_q2_ptmu_bg_weights_1track_blobs = (MnvH2D*) f_bg_weights->Get("h_weights_1track_q2ptbins_qelikenot_blobs");
  MnvH2D *h_q2_ptmu_bg_weights_2track_blobs = (MnvH2D*) f_bg_weights->Get("h_weights_2track_q2ptbins_qelikenot_blobs");

  MnvH2D *h_enu_ptmu_bg_weights_1track_blobs = (MnvH2D*) f_bg_weights->Get("h_weights_1track_enuptbins_qelikenot_blobs");
  MnvH2D *h_enu_ptmu_bg_weights_2track_blobs = (MnvH2D*) f_bg_weights->Get("h_weights_2track_enuptbins_qelikenot_blobs");

  // MnvH2D *h_recoil_ptmu_bg_weights_1track_blobs = (MnvH2D*) f_bg_weights->Get("h_weights_1track_recoilptbins_qelikenot_blobs");
  // MnvH2D *h_recoil_ptmu_bg_weights_2track_blobs = (MnvH2D*) f_bg_weights->Get("h_weights_2track_recoilptbins_qelikenot_blobs");

  // MnvH2D *h_recoil_inc_ptmu_bg_weights_1track_blobs = (MnvH2D*) f_bg_weights->Get("h_weights_1track_recoilincptbins_qelikenot_blobs");
  // MnvH2D *h_recoil_inc_ptmu_bg_weights_2track_blobs = (MnvH2D*) f_bg_weights->Get("h_weights_2track_recoilincptbins_qelikenot_blobs");


  //michel
  MnvH2D *h_pzmu_ptmu_bg_weights_1track_michel = (MnvH2D*) f_bg_weights->Get("h_weights_1track_pzptbins_qelikenot_michel");
  MnvH2D *h_pzmu_ptmu_bg_weights_2track_michel = (MnvH2D*) f_bg_weights->Get("h_weights_2track_pzptbins_qelikenot_michel");

  MnvH2D *h_q2_ptmu_bg_weights_1track_michel = (MnvH2D*) f_bg_weights->Get("h_weights_1track_q2ptbins_qelikenot_michel");
  MnvH2D *h_q2_ptmu_bg_weights_2track_michel = (MnvH2D*) f_bg_weights->Get("h_weights_2track_q2ptbins_qelikenot_michel");

  MnvH2D *h_enu_ptmu_bg_weights_1track_michel = (MnvH2D*) f_bg_weights->Get("h_weights_1track_enuptbins_qelikenot_michel");
  MnvH2D *h_enu_ptmu_bg_weights_2track_michel = (MnvH2D*) f_bg_weights->Get("h_weights_2track_enuptbins_qelikenot_michel");

  // MnvH2D *h_recoil_ptmu_bg_weights_1track_michel = (MnvH2D*) f_bg_weights->Get("h_weights_1track_recoilptbins_qelikenot_michel");
  // MnvH2D *h_recoil_ptmu_bg_weights_2track_michel = (MnvH2D*) f_bg_weights->Get("h_weights_2track_recoilptbins_qelikenot_michel");

  // MnvH2D *h_recoil_inc_ptmu_bg_weights_1track_michel = (MnvH2D*) f_bg_weights->Get("h_weights_1track_recoilincptbins_qelikenot_michel");
  // MnvH2D *h_recoil_inc_ptmu_bg_weights_2track_michel = (MnvH2D*) f_bg_weights->Get("h_weights_2track_recoilincptbins_qelikenot_michel");

  
  //--------------------------------------------
  // Get POT Normalization factor:
  // This is necessary because Background Subtraction
  // uses POT Normalized MC
  //--------------------------------------------
  
  double pot_data = plotUtils->getPOTData( f_1track_signal );
  double pot_mc = plotUtils->getPOTMC( f_1track_signal );
  double pot_scale = plotUtils->getPOTNormFactor( f_1track_signal );
  
  if(pot_scale==0) pot_scale = 1.0;
  cout<< "---------------------------------------"  << endl;
  cout<< "POT INFORMATION:"  << endl;
  cout<< "POT Data = " << pot_data << endl;
  cout<< "POT MC   = " << pot_mc << endl;
  cout<< "POT Scale factor = " << pot_scale << endl;
  cout<< "---------------------------------------"  << endl;

  //2-D
  //pz
  plotUtils->scaleMCHistos( h_pzmu_ptmu_1track_signal, pot_scale);
  plotUtils->scaleMCHistos( h_pzmu_ptmu_2track_signal, pot_scale);
  plotUtils->scaleMCHistos( h_q2_ptmu_1track_signal, pot_scale);
  plotUtils->scaleMCHistos( h_q2_ptmu_2track_signal, pot_scale);
  plotUtils->scaleMCHistos( h_enu_ptmu_1track_signal, pot_scale);
  plotUtils->scaleMCHistos( h_enu_ptmu_2track_signal, pot_scale);
  // plotUtils->scaleMCHistos( h_recoil_ptmu_1track_signal, pot_scale);
  // plotUtils->scaleMCHistos( h_recoil_ptmu_2track_signal, pot_scale);
  // plotUtils->scaleMCHistos( h_recoil_inc_ptmu_1track_signal, pot_scale);
  // plotUtils->scaleMCHistos( h_recoil_inc_ptmu_2track_signal, pot_scale);

  plotUtils->scaleMCHistos( h_pzmu_ptmu_1track_blobs, pot_scale);
  plotUtils->scaleMCHistos( h_pzmu_ptmu_2track_blobs, pot_scale);
  plotUtils->scaleMCHistos( h_q2_ptmu_1track_blobs, pot_scale);
  plotUtils->scaleMCHistos( h_q2_ptmu_2track_blobs, pot_scale);
  plotUtils->scaleMCHistos( h_enu_ptmu_1track_blobs, pot_scale);
  plotUtils->scaleMCHistos( h_enu_ptmu_2track_blobs, pot_scale);
  // plotUtils->scaleMCHistos( h_recoil_ptmu_1track_blobs, pot_scale);
  // plotUtils->scaleMCHistos( h_recoil_ptmu_2track_blobs, pot_scale);
  // plotUtils->scaleMCHistos( h_recoil_inc_ptmu_1track_blobs, pot_scale);
  // plotUtils->scaleMCHistos( h_recoil_inc_ptmu_2track_blobs, pot_scale);

  plotUtils->scaleMCHistos( h_pzmu_ptmu_1track_michel, pot_scale);
  plotUtils->scaleMCHistos( h_pzmu_ptmu_2track_michel, pot_scale);
  plotUtils->scaleMCHistos( h_q2_ptmu_1track_michel, pot_scale);
  plotUtils->scaleMCHistos( h_q2_ptmu_2track_michel, pot_scale);
  plotUtils->scaleMCHistos( h_enu_ptmu_1track_michel, pot_scale);
  plotUtils->scaleMCHistos( h_enu_ptmu_2track_michel, pot_scale);
  // plotUtils->scaleMCHistos( h_recoil_ptmu_1track_michel, pot_scale);
  // plotUtils->scaleMCHistos( h_recoil_ptmu_2track_michel, pot_scale);
  // plotUtils->scaleMCHistos( h_recoil_inc_ptmu_1track_michel, pot_scale);
  // plotUtils->scaleMCHistos( h_recoil_inc_ptmu_2track_michel, pot_scale);

  //--------------------------------------------
  // Background Subtraction
  //--------------------------------------------
  cout<< "---------------------------------------"  << endl;
  cout<< "Background Subtraction"  << endl;
  cout<< "---------------------------------------"  << endl;

  MnvH2D *h_data_pzmu_ptmu_weighted_bck_1track_signal, *h_data_pzmu_ptmu_predicted_bck_1track_signal, *h_mc_pzmu_ptmu_weighted_bck_1track_signal, *h_mc_pzmu_ptmu_tuned_bck_1track_signal;
  MnvH2D *h_data_pzmu_ptmu_weighted_bck_2track_signal, *h_data_pzmu_ptmu_predicted_bck_2track_signal, *h_mc_pzmu_ptmu_weighted_bck_2track_signal, *h_mc_pzmu_ptmu_tuned_bck_2track_signal;

  MnvH2D *h_data_pzmu_ptmu_weighted_bck_1track_blobs, *h_data_pzmu_ptmu_predicted_bck_1track_blobs, *h_mc_pzmu_ptmu_weighted_bck_1track_blobs, *h_mc_pzmu_ptmu_tuned_bck_1track_blobs;
  MnvH2D *h_data_pzmu_ptmu_weighted_bck_2track_blobs, *h_data_pzmu_ptmu_predicted_bck_2track_blobs, *h_mc_pzmu_ptmu_weighted_bck_2track_blobs, *h_mc_pzmu_ptmu_tuned_bck_2track_blobs;

  MnvH2D *h_data_pzmu_ptmu_weighted_bck_1track_michel, *h_data_pzmu_ptmu_predicted_bck_1track_michel, *h_mc_pzmu_ptmu_weighted_bck_1track_michel, *h_mc_pzmu_ptmu_tuned_bck_1track_michel;
  MnvH2D *h_data_pzmu_ptmu_weighted_bck_2track_michel, *h_data_pzmu_ptmu_predicted_bck_2track_michel, *h_mc_pzmu_ptmu_weighted_bck_2track_michel, *h_mc_pzmu_ptmu_tuned_bck_2track_michel;

  MnvH2D *bkg_tuned_pzmu_ptmu_1track_signal,*bkg_tuned_pzmu_ptmu_1track_blobs,*bkg_tuned_pzmu_ptmu_1track_michel;
  MnvH2D *bkg_tuned_pzmu_ptmu_2track_signal,*bkg_tuned_pzmu_ptmu_2track_blobs,*bkg_tuned_pzmu_ptmu_2track_michel;
  
  cout<<"Background Scale of 1 Track Sub-Sample"<<endl;
  ScaleBackground( h_pzmu_ptmu_1track_signal, h_pzmu_ptmu_bg_weights_1track, h_data_pzmu_ptmu_weighted_bck_1track_signal, h_mc_pzmu_ptmu_weighted_bck_1track_signal,bkg_tuned_pzmu_ptmu_1track_signal);
  ScaleBackground( h_pzmu_ptmu_1track_blobs, h_pzmu_ptmu_bg_weights_1track_blobs, h_data_pzmu_ptmu_weighted_bck_1track_blobs, h_mc_pzmu_ptmu_weighted_bck_1track_blobs,bkg_tuned_pzmu_ptmu_1track_blobs);
  ScaleBackground( h_pzmu_ptmu_1track_michel, h_pzmu_ptmu_bg_weights_1track_michel, h_data_pzmu_ptmu_weighted_bck_1track_michel, h_mc_pzmu_ptmu_weighted_bck_1track_michel,bkg_tuned_pzmu_ptmu_1track_michel);
  cout<<"Background Scale of >=2 Track Sub-Sample"<<endl;
  ScaleBackground( h_pzmu_ptmu_2track_signal, h_pzmu_ptmu_bg_weights_2track, h_data_pzmu_ptmu_weighted_bck_2track_signal, h_mc_pzmu_ptmu_weighted_bck_2track_signal,bkg_tuned_pzmu_ptmu_2track_signal);
  ScaleBackground( h_pzmu_ptmu_2track_blobs, h_pzmu_ptmu_bg_weights_2track_blobs, h_data_pzmu_ptmu_weighted_bck_2track_blobs, h_mc_pzmu_ptmu_weighted_bck_2track_blobs,bkg_tuned_pzmu_ptmu_2track_blobs);
  ScaleBackground( h_pzmu_ptmu_2track_michel, h_pzmu_ptmu_bg_weights_2track_michel, h_data_pzmu_ptmu_weighted_bck_2track_michel, h_mc_pzmu_ptmu_weighted_bck_2track_michel,bkg_tuned_pzmu_ptmu_2track_michel);

  /*
  MnvH2D *h_data_recoil_ptmu_weighted_bck_1track_signal, *h_data_recoil_ptmu_predicted_bck_1track_signal, *h_mc_recoil_ptmu_weighted_bck_1track_signal, *h_mc_recoil_ptmu_tuned_bck_1track_signal;
  MnvH2D *h_data_recoil_ptmu_weighted_bck_2track_signal, *h_data_recoil_ptmu_predicted_bck_2track_signal, *h_mc_recoil_ptmu_weighted_bck_2track_signal, *h_mc_recoil_ptmu_tuned_bck_2track_signal;

  MnvH2D *h_data_recoil_ptmu_weighted_bck_1track_blobs, *h_data_recoil_ptmu_predicted_bck_1track_blobs, *h_mc_recoil_ptmu_weighted_bck_1track_blobs, *h_mc_recoil_ptmu_tuned_bck_1track_blobs;
  MnvH2D *h_data_recoil_ptmu_weighted_bck_2track_blobs, *h_data_recoil_ptmu_predicted_bck_2track_blobs, *h_mc_recoil_ptmu_weighted_bck_2track_blobs, *h_mc_recoil_ptmu_tuned_bck_2track_blobs;

  MnvH2D *h_data_recoil_ptmu_weighted_bck_1track_michel, *h_data_recoil_ptmu_predicted_bck_1track_michel, *h_mc_recoil_ptmu_weighted_bck_1track_michel, *h_mc_recoil_ptmu_tuned_bck_1track_michel;
  MnvH2D *h_data_recoil_ptmu_weighted_bck_2track_michel, *h_data_recoil_ptmu_predicted_bck_2track_michel, *h_mc_recoil_ptmu_weighted_bck_2track_michel, *h_mc_recoil_ptmu_tuned_bck_2track_michel;

  MnvH2D *bkg_tuned_recoil_ptmu_1track_signal,*bkg_tuned_recoil_ptmu_1track_blobs,*bkg_tuned_recoil_ptmu_1track_michel;
  MnvH2D *bkg_tuned_recoil_ptmu_2track_signal,*bkg_tuned_recoil_ptmu_2track_blobs,*bkg_tuned_recoil_ptmu_2track_michel;
  
  cout<<"Background Scale of 1 Track Sub-Sample"<<endl;
  ScaleBackground( h_recoil_ptmu_1track_signal, h_recoil_ptmu_bg_weights_1track, h_data_recoil_ptmu_weighted_bck_1track_signal, h_mc_recoil_ptmu_weighted_bck_1track_signal,bkg_tuned_recoil_ptmu_1track_signal);
  ScaleBackground( h_recoil_ptmu_1track_blobs, h_recoil_ptmu_bg_weights_1track_blobs, h_data_recoil_ptmu_weighted_bck_1track_blobs, h_mc_recoil_ptmu_weighted_bck_1track_blobs,bkg_tuned_recoil_ptmu_1track_blobs);
  ScaleBackground( h_recoil_ptmu_1track_michel, h_recoil_ptmu_bg_weights_1track_michel, h_data_recoil_ptmu_weighted_bck_1track_michel, h_mc_recoil_ptmu_weighted_bck_1track_michel,bkg_tuned_recoil_ptmu_1track_michel);
  cout<<"Background Scale of >=2 Track Sub-Sample"<<endl;
  ScaleBackground( h_recoil_ptmu_2track_signal, h_recoil_ptmu_bg_weights_2track, h_data_recoil_ptmu_weighted_bck_2track_signal, h_mc_recoil_ptmu_weighted_bck_2track_signal,bkg_tuned_recoil_ptmu_2track_signal);
  ScaleBackground( h_recoil_ptmu_2track_blobs, h_recoil_ptmu_bg_weights_2track_blobs, h_data_recoil_ptmu_weighted_bck_2track_blobs, h_mc_recoil_ptmu_weighted_bck_2track_blobs,bkg_tuned_recoil_ptmu_2track_blobs);
  ScaleBackground( h_recoil_ptmu_2track_michel, h_recoil_ptmu_bg_weights_2track_michel, h_data_recoil_ptmu_weighted_bck_2track_michel, h_mc_recoil_ptmu_weighted_bck_2track_michel,bkg_tuned_recoil_ptmu_2track_michel);


  MnvH2D *h_data_recoil_inc_ptmu_weighted_bck_1track_signal, *h_data_recoil_inc_ptmu_predicted_bck_1track_signal, *h_mc_recoil_inc_ptmu_weighted_bck_1track_signal, *h_mc_recoil_inc_ptmu_tuned_bck_1track_signal;
  MnvH2D *h_data_recoil_inc_ptmu_weighted_bck_2track_signal, *h_data_recoil_inc_ptmu_predicted_bck_2track_signal, *h_mc_recoil_inc_ptmu_weighted_bck_2track_signal, *h_mc_recoil_inc_ptmu_tuned_bck_2track_signal;

  MnvH2D *h_data_recoil_inc_ptmu_weighted_bck_1track_blobs, *h_data_recoil_inc_ptmu_predicted_bck_1track_blobs, *h_mc_recoil_inc_ptmu_weighted_bck_1track_blobs, *h_mc_recoil_inc_ptmu_tuned_bck_1track_blobs;
  MnvH2D *h_data_recoil_inc_ptmu_weighted_bck_2track_blobs, *h_data_recoil_inc_ptmu_predicted_bck_2track_blobs, *h_mc_recoil_inc_ptmu_weighted_bck_2track_blobs, *h_mc_recoil_inc_ptmu_tuned_bck_2track_blobs;

  MnvH2D *h_data_recoil_inc_ptmu_weighted_bck_1track_michel, *h_data_recoil_inc_ptmu_predicted_bck_1track_michel, *h_mc_recoil_inc_ptmu_weighted_bck_1track_michel, *h_mc_recoil_inc_ptmu_tuned_bck_1track_michel;
  MnvH2D *h_data_recoil_inc_ptmu_weighted_bck_2track_michel, *h_data_recoil_inc_ptmu_predicted_bck_2track_michel, *h_mc_recoil_inc_ptmu_weighted_bck_2track_michel, *h_mc_recoil_inc_ptmu_tuned_bck_2track_michel;

  MnvH2D *bkg_tuned_recoil_inc_ptmu_1track_signal,*bkg_tuned_recoil_inc_ptmu_1track_blobs,*bkg_tuned_recoil_inc_ptmu_1track_michel;
  MnvH2D *bkg_tuned_recoil_inc_ptmu_2track_signal,*bkg_tuned_recoil_inc_ptmu_2track_blobs,*bkg_tuned_recoil_inc_ptmu_2track_michel;
  
  cout<<"Background Scale of 1 Track Sub-Sample"<<endl;
  ScaleBackground( h_recoil_inc_ptmu_1track_signal, h_recoil_inc_ptmu_bg_weights_1track, h_data_recoil_inc_ptmu_weighted_bck_1track_signal, h_mc_recoil_inc_ptmu_weighted_bck_1track_signal,bkg_tuned_recoil_inc_ptmu_1track_signal);
  ScaleBackground( h_recoil_inc_ptmu_1track_blobs, h_recoil_inc_ptmu_bg_weights_1track_blobs, h_data_recoil_inc_ptmu_weighted_bck_1track_blobs, h_mc_recoil_inc_ptmu_weighted_bck_1track_blobs,bkg_tuned_recoil_inc_ptmu_1track_blobs);
  ScaleBackground( h_recoil_inc_ptmu_1track_michel, h_recoil_inc_ptmu_bg_weights_1track_michel, h_data_recoil_inc_ptmu_weighted_bck_1track_michel, h_mc_recoil_inc_ptmu_weighted_bck_1track_michel,bkg_tuned_recoil_inc_ptmu_1track_michel);
  cout<<"Background Scale of >=2 Track Sub-Sample"<<endl;
  ScaleBackground( h_recoil_inc_ptmu_2track_signal, h_recoil_inc_ptmu_bg_weights_2track, h_data_recoil_inc_ptmu_weighted_bck_2track_signal, h_mc_recoil_inc_ptmu_weighted_bck_2track_signal,bkg_tuned_recoil_inc_ptmu_2track_signal);
  ScaleBackground( h_recoil_inc_ptmu_2track_blobs, h_recoil_inc_ptmu_bg_weights_2track_blobs, h_data_recoil_inc_ptmu_weighted_bck_2track_blobs, h_mc_recoil_inc_ptmu_weighted_bck_2track_blobs,bkg_tuned_recoil_inc_ptmu_2track_blobs);
  ScaleBackground( h_recoil_inc_ptmu_2track_michel, h_recoil_inc_ptmu_bg_weights_2track_michel, h_data_recoil_inc_ptmu_weighted_bck_2track_michel, h_mc_recoil_inc_ptmu_weighted_bck_2track_michel,bkg_tuned_recoil_inc_ptmu_2track_michel);
  
  */
  MnvH2D *h_data_q2_ptmu_weighted_bck_1track_signal, *h_data_q2_ptmu_predicted_bck_1track_signal, *h_mc_q2_ptmu_weighted_bck_1track_signal, *h_mc_q2_ptmu_tuned_bck_1track_signal;
  MnvH2D *h_data_q2_ptmu_weighted_bck_2track_signal, *h_data_q2_ptmu_predicted_bck_2track_signal, *h_mc_q2_ptmu_weighted_bck_2track_signal, *h_mc_q2_ptmu_tuned_bck_2track_signal;

  MnvH2D *h_data_q2_ptmu_weighted_bck_1track_blobs, *h_data_q2_ptmu_predicted_bck_1track_blobs, *h_mc_q2_ptmu_weighted_bck_1track_blobs, *h_mc_q2_ptmu_tuned_bck_1track_blobs;
  MnvH2D *h_data_q2_ptmu_weighted_bck_2track_blobs, *h_data_q2_ptmu_predicted_bck_2track_blobs, *h_mc_q2_ptmu_weighted_bck_2track_blobs, *h_mc_q2_ptmu_tuned_bck_2track_blobs;

  MnvH2D *h_data_q2_ptmu_weighted_bck_1track_michel, *h_data_q2_ptmu_predicted_bck_1track_michel, *h_mc_q2_ptmu_weighted_bck_1track_michel, *h_mc_q2_ptmu_tuned_bck_1track_michel;
  MnvH2D *h_data_q2_ptmu_weighted_bck_2track_michel, *h_data_q2_ptmu_predicted_bck_2track_michel, *h_mc_q2_ptmu_weighted_bck_2track_michel, *h_mc_q2_ptmu_tuned_bck_2track_michel;

  MnvH2D *bkg_tuned_q2_ptmu_1track_signal,*bkg_tuned_q2_ptmu_1track_blobs,*bkg_tuned_q2_ptmu_1track_michel;
  MnvH2D *bkg_tuned_q2_ptmu_2track_signal,*bkg_tuned_q2_ptmu_2track_blobs,*bkg_tuned_q2_ptmu_2track_michel;
  
  cout<<"Background Scale of 1 Track Sub-Sample"<<endl;
  ScaleBackground( h_q2_ptmu_1track_signal, h_q2_ptmu_bg_weights_1track, h_data_q2_ptmu_weighted_bck_1track_signal, h_mc_q2_ptmu_weighted_bck_1track_signal,bkg_tuned_q2_ptmu_1track_signal);
  ScaleBackground( h_q2_ptmu_1track_blobs, h_q2_ptmu_bg_weights_1track_blobs, h_data_q2_ptmu_weighted_bck_1track_blobs, h_mc_q2_ptmu_weighted_bck_1track_blobs,bkg_tuned_q2_ptmu_1track_blobs);
  ScaleBackground( h_q2_ptmu_1track_michel, h_q2_ptmu_bg_weights_1track_michel, h_data_q2_ptmu_weighted_bck_1track_michel, h_mc_q2_ptmu_weighted_bck_1track_michel,bkg_tuned_q2_ptmu_1track_michel);
  cout<<"Background Scale of >=2 Track Sub-Sample"<<endl;
  ScaleBackground( h_q2_ptmu_2track_signal, h_q2_ptmu_bg_weights_2track, h_data_q2_ptmu_weighted_bck_2track_signal, h_mc_q2_ptmu_weighted_bck_2track_signal,bkg_tuned_q2_ptmu_2track_signal);
  ScaleBackground( h_q2_ptmu_2track_blobs, h_q2_ptmu_bg_weights_2track_blobs, h_data_q2_ptmu_weighted_bck_2track_blobs, h_mc_q2_ptmu_weighted_bck_2track_blobs,bkg_tuned_q2_ptmu_2track_blobs);
  ScaleBackground( h_q2_ptmu_2track_michel, h_q2_ptmu_bg_weights_2track_michel, h_data_q2_ptmu_weighted_bck_2track_michel, h_mc_q2_ptmu_weighted_bck_2track_michel,bkg_tuned_q2_ptmu_2track_michel);

  MnvH2D *h_data_enu_ptmu_weighted_bck_1track_signal, *h_data_enu_ptmu_predicted_bck_1track_signal, *h_mc_enu_ptmu_weighted_bck_1track_signal, *h_mc_enu_ptmu_tuned_bck_1track_signal;
  MnvH2D *h_data_enu_ptmu_weighted_bck_2track_signal, *h_data_enu_ptmu_predicted_bck_2track_signal, *h_mc_enu_ptmu_weighted_bck_2track_signal, *h_mc_enu_ptmu_tuned_bck_2track_signal;

  MnvH2D *h_data_enu_ptmu_weighted_bck_1track_blobs, *h_data_enu_ptmu_predicted_bck_1track_blobs, *h_mc_enu_ptmu_weighted_bck_1track_blobs, *h_mc_enu_ptmu_tuned_bck_1track_blobs;
  MnvH2D *h_data_enu_ptmu_weighted_bck_2track_blobs, *h_data_enu_ptmu_predicted_bck_2track_blobs, *h_mc_enu_ptmu_weighted_bck_2track_blobs, *h_mc_enu_ptmu_tuned_bck_2track_blobs;

  MnvH2D *h_data_enu_ptmu_weighted_bck_1track_michel, *h_data_enu_ptmu_predicted_bck_1track_michel, *h_mc_enu_ptmu_weighted_bck_1track_michel, *h_mc_enu_ptmu_tuned_bck_1track_michel;
  MnvH2D *h_data_enu_ptmu_weighted_bck_2track_michel, *h_data_enu_ptmu_predicted_bck_2track_michel, *h_mc_enu_ptmu_weighted_bck_2track_michel, *h_mc_enu_ptmu_tuned_bck_2track_michel;

  MnvH2D *bkg_tuned_enu_ptmu_1track_signal,*bkg_tuned_enu_ptmu_1track_blobs,*bkg_tuned_enu_ptmu_1track_michel;
  MnvH2D *bkg_tuned_enu_ptmu_2track_signal,*bkg_tuned_enu_ptmu_2track_blobs,*bkg_tuned_enu_ptmu_2track_michel;
  
  cout<<"Background Scale of 1 Track Sub-Sample"<<endl;
  ScaleBackground( h_enu_ptmu_1track_signal, h_enu_ptmu_bg_weights_1track, h_data_enu_ptmu_weighted_bck_1track_signal, h_mc_enu_ptmu_weighted_bck_1track_signal,bkg_tuned_enu_ptmu_1track_signal);
  ScaleBackground( h_enu_ptmu_1track_blobs, h_enu_ptmu_bg_weights_1track_blobs, h_data_enu_ptmu_weighted_bck_1track_blobs, h_mc_enu_ptmu_weighted_bck_1track_blobs,bkg_tuned_enu_ptmu_1track_blobs);
  ScaleBackground( h_enu_ptmu_1track_michel, h_enu_ptmu_bg_weights_1track_michel, h_data_enu_ptmu_weighted_bck_1track_michel, h_mc_enu_ptmu_weighted_bck_1track_michel,bkg_tuned_enu_ptmu_1track_michel);
  cout<<"Background Scale of >=2 Track Sub-Sample"<<endl;
  ScaleBackground( h_enu_ptmu_2track_signal, h_enu_ptmu_bg_weights_2track, h_data_enu_ptmu_weighted_bck_2track_signal, h_mc_enu_ptmu_weighted_bck_2track_signal,bkg_tuned_enu_ptmu_2track_signal);
  ScaleBackground( h_enu_ptmu_2track_blobs, h_enu_ptmu_bg_weights_2track_blobs, h_data_enu_ptmu_weighted_bck_2track_blobs, h_mc_enu_ptmu_weighted_bck_2track_blobs,bkg_tuned_enu_ptmu_2track_blobs);
  ScaleBackground( h_enu_ptmu_2track_michel, h_enu_ptmu_bg_weights_2track_michel, h_data_enu_ptmu_weighted_bck_2track_michel, h_mc_enu_ptmu_weighted_bck_2track_michel,bkg_tuned_enu_ptmu_2track_michel);
  


  //-----------------------------------------------
  // Merge both background subtracted sub-samples
  //-----------------------------------------------
  cout<<"Merging 1 and >=2 Tracks Sub-Samples..."<<endl;
  MnvH2D *h_pzmu_ptmu_data_1track_weighted_bck_signal = (MnvH2D*)h_data_pzmu_ptmu_weighted_bck_1track_signal->Clone("h_pzmu_ptmu_data_1track_weighted_bck_signal");
  MnvH2D *h_pzmu_ptmu_data_2track_weighted_bck_signal = (MnvH2D*)h_data_pzmu_ptmu_weighted_bck_2track_signal->Clone("h_pzmu_ptmu_data_2track_weighted_bck_signal");
  MnvH2D *h_pzmu_ptmu_mc_1track_weighted_bck_signal = (MnvH2D*)h_mc_pzmu_ptmu_weighted_bck_1track_signal->Clone("h_pzmu_ptmu_mc_1track_weighted_bck_signal");
  MnvH2D *h_pzmu_ptmu_mc_2track_weighted_bck_signal = (MnvH2D*)h_mc_pzmu_ptmu_weighted_bck_2track_signal->Clone("h_pzmu_ptmu_mc_2track_weighted_bck_signal");

  MnvH2D *h_pzmu_ptmu_data_1track_weighted_bck_blobs = (MnvH2D*)h_data_pzmu_ptmu_weighted_bck_1track_blobs->Clone("h_pzmu_ptmu_data_1track_weighted_bck_blobs");
  MnvH2D *h_pzmu_ptmu_data_2track_weighted_bck_blobs = (MnvH2D*)h_data_pzmu_ptmu_weighted_bck_2track_blobs->Clone("h_pzmu_ptmu_data_2track_weighted_bck_blobs");
  MnvH2D *h_pzmu_ptmu_mc_1track_weighted_bck_blobs = (MnvH2D*)h_mc_pzmu_ptmu_weighted_bck_1track_blobs->Clone("h_pzmu_ptmu_mc_1track_weighted_bck_blobs");
  MnvH2D *h_pzmu_ptmu_mc_2track_weighted_bck_blobs = (MnvH2D*)h_mc_pzmu_ptmu_weighted_bck_2track_blobs->Clone("h_pzmu_ptmu_mc_2track_weighted_bck_blobs");

  MnvH2D *h_pzmu_ptmu_data_1track_weighted_bck_michel = (MnvH2D*)h_data_pzmu_ptmu_weighted_bck_1track_michel->Clone("h_pzmu_ptmu_data_1track_weighted_bck_michel");
  MnvH2D *h_pzmu_ptmu_data_2track_weighted_bck_michel = (MnvH2D*)h_data_pzmu_ptmu_weighted_bck_2track_michel->Clone("h_pzmu_ptmu_data_2track_weighted_bck_michel");
  MnvH2D *h_pzmu_ptmu_mc_1track_weighted_bck_michel = (MnvH2D*)h_mc_pzmu_ptmu_weighted_bck_1track_michel->Clone("h_pzmu_ptmu_mc_1track_weighted_bck_michel");
  MnvH2D *h_pzmu_ptmu_mc_2track_weighted_bck_michel = (MnvH2D*)h_mc_pzmu_ptmu_weighted_bck_2track_michel->Clone("h_pzmu_ptmu_mc_2track_weighted_bck_michel");

  /*
  MnvH2D *h_recoil_ptmu_data_1track_weighted_bck_signal = (MnvH2D*)h_data_recoil_ptmu_weighted_bck_1track_signal->Clone("h_recoil_ptmu_data_1track_weighted_bck_signal");
  MnvH2D *h_recoil_ptmu_data_2track_weighted_bck_signal = (MnvH2D*)h_data_recoil_ptmu_weighted_bck_2track_signal->Clone("h_recoil_ptmu_data_2track_weighted_bck_signal");
  MnvH2D *h_recoil_ptmu_mc_1track_weighted_bck_signal = (MnvH2D*)h_mc_recoil_ptmu_weighted_bck_1track_signal->Clone("h_recoil_ptmu_mc_1track_weighted_bck_signal");
  MnvH2D *h_recoil_ptmu_mc_2track_weighted_bck_signal = (MnvH2D*)h_mc_recoil_ptmu_weighted_bck_2track_signal->Clone("h_recoil_ptmu_mc_2track_weighted_bck_signal");

  MnvH2D *h_recoil_ptmu_data_1track_weighted_bck_blobs = (MnvH2D*)h_data_recoil_ptmu_weighted_bck_1track_blobs->Clone("h_recoil_ptmu_data_1track_weighted_bck_blobs");
  MnvH2D *h_recoil_ptmu_data_2track_weighted_bck_blobs = (MnvH2D*)h_data_recoil_ptmu_weighted_bck_2track_blobs->Clone("h_recoil_ptmu_data_2track_weighted_bck_blobs");
  MnvH2D *h_recoil_ptmu_mc_1track_weighted_bck_blobs = (MnvH2D*)h_mc_recoil_ptmu_weighted_bck_1track_blobs->Clone("h_recoil_ptmu_mc_1track_weighted_bck_blobs");
  MnvH2D *h_recoil_ptmu_mc_2track_weighted_bck_blobs = (MnvH2D*)h_mc_recoil_ptmu_weighted_bck_2track_blobs->Clone("h_recoil_ptmu_mc_2track_weighted_bck_blobs");

  MnvH2D *h_recoil_ptmu_data_1track_weighted_bck_michel = (MnvH2D*)h_data_recoil_ptmu_weighted_bck_1track_michel->Clone("h_recoil_ptmu_data_1track_weighted_bck_michel");
  MnvH2D *h_recoil_ptmu_data_2track_weighted_bck_michel = (MnvH2D*)h_data_recoil_ptmu_weighted_bck_2track_michel->Clone("h_recoil_ptmu_data_2track_weighted_bck_michel");
  MnvH2D *h_recoil_ptmu_mc_1track_weighted_bck_michel = (MnvH2D*)h_mc_recoil_ptmu_weighted_bck_1track_michel->Clone("h_recoil_ptmu_mc_1track_weighted_bck_michel");
  MnvH2D *h_recoil_ptmu_mc_2track_weighted_bck_michel = (MnvH2D*)h_mc_recoil_ptmu_weighted_bck_2track_michel->Clone("h_recoil_ptmu_mc_2track_weighted_bck_michel");


  MnvH2D *h_recoil_inc_ptmu_data_1track_weighted_bck_signal = (MnvH2D*)h_data_recoil_inc_ptmu_weighted_bck_1track_signal->Clone("h_recoil_inc_ptmu_data_1track_weighted_bck_signal");
  MnvH2D *h_recoil_inc_ptmu_data_2track_weighted_bck_signal = (MnvH2D*)h_data_recoil_inc_ptmu_weighted_bck_2track_signal->Clone("h_recoil_inc_ptmu_data_2track_weighted_bck_signal");
  MnvH2D *h_recoil_inc_ptmu_mc_1track_weighted_bck_signal = (MnvH2D*)h_mc_recoil_inc_ptmu_weighted_bck_1track_signal->Clone("h_recoil_inc_ptmu_mc_1track_weighted_bck_signal");
  MnvH2D *h_recoil_inc_ptmu_mc_2track_weighted_bck_signal = (MnvH2D*)h_mc_recoil_inc_ptmu_weighted_bck_2track_signal->Clone("h_recoil_inc_ptmu_mc_2track_weighted_bck_signal");

  MnvH2D *h_recoil_inc_ptmu_data_1track_weighted_bck_blobs = (MnvH2D*)h_data_recoil_inc_ptmu_weighted_bck_1track_blobs->Clone("h_recoil_inc_ptmu_data_1track_weighted_bck_blobs");
  MnvH2D *h_recoil_inc_ptmu_data_2track_weighted_bck_blobs = (MnvH2D*)h_data_recoil_inc_ptmu_weighted_bck_2track_blobs->Clone("h_recoil_inc_ptmu_data_2track_weighted_bck_blobs");
  MnvH2D *h_recoil_inc_ptmu_mc_1track_weighted_bck_blobs = (MnvH2D*)h_mc_recoil_inc_ptmu_weighted_bck_1track_blobs->Clone("h_recoil_inc_ptmu_mc_1track_weighted_bck_blobs");
  MnvH2D *h_recoil_inc_ptmu_mc_2track_weighted_bck_blobs = (MnvH2D*)h_mc_recoil_inc_ptmu_weighted_bck_2track_blobs->Clone("h_recoil_inc_ptmu_mc_2track_weighted_bck_blobs");

  MnvH2D *h_recoil_inc_ptmu_data_1track_weighted_bck_michel = (MnvH2D*)h_data_recoil_inc_ptmu_weighted_bck_1track_michel->Clone("h_recoil_inc_ptmu_data_1track_weighted_bck_michel");
  MnvH2D *h_recoil_inc_ptmu_data_2track_weighted_bck_michel = (MnvH2D*)h_data_recoil_inc_ptmu_weighted_bck_2track_michel->Clone("h_recoil_inc_ptmu_data_2track_weighted_bck_michel");
  MnvH2D *h_recoil_inc_ptmu_mc_1track_weighted_bck_michel = (MnvH2D*)h_mc_recoil_inc_ptmu_weighted_bck_1track_michel->Clone("h_recoil_inc_ptmu_mc_1track_weighted_bck_michel");
  MnvH2D *h_recoil_inc_ptmu_mc_2track_weighted_bck_michel = (MnvH2D*)h_mc_recoil_inc_ptmu_weighted_bck_2track_michel->Clone("h_recoil_inc_ptmu_mc_2track_weighted_bck_michel");

  */
  MnvH2D *h_q2_ptmu_data_1track_weighted_bck_signal = (MnvH2D*)h_data_q2_ptmu_weighted_bck_1track_signal->Clone("h_q2_ptmu_data_1track_weighted_bck_signal");
  MnvH2D *h_q2_ptmu_data_2track_weighted_bck_signal = (MnvH2D*)h_data_q2_ptmu_weighted_bck_2track_signal->Clone("h_q2_ptmu_data_2track_weighted_bck_signal");
  MnvH2D *h_q2_ptmu_mc_1track_weighted_bck_signal = (MnvH2D*)h_mc_q2_ptmu_weighted_bck_1track_signal->Clone("h_q2_ptmu_mc_1track_weighted_bck_signal");
  MnvH2D *h_q2_ptmu_mc_2track_weighted_bck_signal = (MnvH2D*)h_mc_q2_ptmu_weighted_bck_2track_signal->Clone("h_q2_ptmu_mc_2track_weighted_bck_signal");

  MnvH2D *h_q2_ptmu_data_1track_weighted_bck_blobs = (MnvH2D*)h_data_q2_ptmu_weighted_bck_1track_blobs->Clone("h_q2_ptmu_data_1track_weighted_bck_blobs");
  MnvH2D *h_q2_ptmu_data_2track_weighted_bck_blobs = (MnvH2D*)h_data_q2_ptmu_weighted_bck_2track_blobs->Clone("h_q2_ptmu_data_2track_weighted_bck_blobs");
  MnvH2D *h_q2_ptmu_mc_1track_weighted_bck_blobs = (MnvH2D*)h_mc_q2_ptmu_weighted_bck_1track_blobs->Clone("h_q2_ptmu_mc_1track_weighted_bck_blobs");
  MnvH2D *h_q2_ptmu_mc_2track_weighted_bck_blobs = (MnvH2D*)h_mc_q2_ptmu_weighted_bck_2track_blobs->Clone("h_q2_ptmu_mc_2track_weighted_bck_blobs");

  MnvH2D *h_q2_ptmu_data_1track_weighted_bck_michel = (MnvH2D*)h_data_q2_ptmu_weighted_bck_1track_michel->Clone("h_q2_ptmu_data_1track_weighted_bck_michel");
  MnvH2D *h_q2_ptmu_data_2track_weighted_bck_michel = (MnvH2D*)h_data_q2_ptmu_weighted_bck_2track_michel->Clone("h_q2_ptmu_data_2track_weighted_bck_michel");
  MnvH2D *h_q2_ptmu_mc_1track_weighted_bck_michel = (MnvH2D*)h_mc_q2_ptmu_weighted_bck_1track_michel->Clone("h_q2_ptmu_mc_1track_weighted_bck_michel");
  MnvH2D *h_q2_ptmu_mc_2track_weighted_bck_michel = (MnvH2D*)h_mc_q2_ptmu_weighted_bck_2track_michel->Clone("h_q2_ptmu_mc_2track_weighted_bck_michel");

  
  MnvH2D *h_enu_ptmu_data_1track_weighted_bck_signal = (MnvH2D*)h_data_enu_ptmu_weighted_bck_1track_signal->Clone("h_enu_ptmu_data_1track_weighted_bck_signal");
  MnvH2D *h_enu_ptmu_data_2track_weighted_bck_signal = (MnvH2D*)h_data_enu_ptmu_weighted_bck_2track_signal->Clone("h_enu_ptmu_data_2track_weighted_bck_signal");
  MnvH2D *h_enu_ptmu_mc_1track_weighted_bck_signal = (MnvH2D*)h_mc_enu_ptmu_weighted_bck_1track_signal->Clone("h_enu_ptmu_mc_1track_weighted_bck_signal");
  MnvH2D *h_enu_ptmu_mc_2track_weighted_bck_signal = (MnvH2D*)h_mc_enu_ptmu_weighted_bck_2track_signal->Clone("h_enu_ptmu_mc_2track_weighted_bck_signal");

  MnvH2D *h_enu_ptmu_data_1track_weighted_bck_blobs = (MnvH2D*)h_data_enu_ptmu_weighted_bck_1track_blobs->Clone("h_enu_ptmu_data_1track_weighted_bck_blobs");
  MnvH2D *h_enu_ptmu_data_2track_weighted_bck_blobs = (MnvH2D*)h_data_enu_ptmu_weighted_bck_2track_blobs->Clone("h_enu_ptmu_data_2track_weighted_bck_blobs");
  MnvH2D *h_enu_ptmu_mc_1track_weighted_bck_blobs = (MnvH2D*)h_mc_enu_ptmu_weighted_bck_1track_blobs->Clone("h_enu_ptmu_mc_1track_weighted_bck_blobs");
  MnvH2D *h_enu_ptmu_mc_2track_weighted_bck_blobs = (MnvH2D*)h_mc_enu_ptmu_weighted_bck_2track_blobs->Clone("h_enu_ptmu_mc_2track_weighted_bck_blobs");

  MnvH2D *h_enu_ptmu_data_1track_weighted_bck_michel = (MnvH2D*)h_data_enu_ptmu_weighted_bck_1track_michel->Clone("h_enu_ptmu_data_1track_weighted_bck_michel");
  MnvH2D *h_enu_ptmu_data_2track_weighted_bck_michel = (MnvH2D*)h_data_enu_ptmu_weighted_bck_2track_michel->Clone("h_enu_ptmu_data_2track_weighted_bck_michel");
  MnvH2D *h_enu_ptmu_mc_1track_weighted_bck_michel = (MnvH2D*)h_mc_enu_ptmu_weighted_bck_1track_michel->Clone("h_enu_ptmu_mc_1track_weighted_bck_michel");
  MnvH2D *h_enu_ptmu_mc_2track_weighted_bck_michel = (MnvH2D*)h_mc_enu_ptmu_weighted_bck_2track_michel->Clone("h_enu_ptmu_mc_2track_weighted_bck_michel");


  //==================================================================
  // Create ROOT file to store histograms
  // Write to file all the created histograms 
  //==================================================================
  TFile *f_output = new TFile( output_filename.c_str(), "RECREATE");

  cout<< "--------------------------------------------"  << endl;
  cout<< "Writing histograms to file"  << endl;
  cout<< "--------------------------------------------"  << endl;
  //pz
  h_pzmu_ptmu_data_1track_weighted_bck_signal->Write();
  h_pzmu_ptmu_data_2track_weighted_bck_signal->Write();
  h_pzmu_ptmu_data_1track_weighted_bck_blobs->Write();
  h_pzmu_ptmu_data_2track_weighted_bck_blobs->Write();
  h_pzmu_ptmu_data_1track_weighted_bck_michel->Write();
  h_pzmu_ptmu_data_2track_weighted_bck_michel->Write();


  h_pzmu_ptmu_mc_1track_weighted_bck_signal->Write();
  h_pzmu_ptmu_mc_2track_weighted_bck_signal->Write();
  h_pzmu_ptmu_mc_1track_weighted_bck_blobs->Write();
  h_pzmu_ptmu_mc_2track_weighted_bck_blobs->Write();
  h_pzmu_ptmu_mc_1track_weighted_bck_michel->Write();
  h_pzmu_ptmu_mc_2track_weighted_bck_michel->Write();


  bkg_tuned_pzmu_ptmu_1track_signal->Write("h_tuned_bkg_pzmu_ptmu_1track_signal");
  bkg_tuned_pzmu_ptmu_1track_blobs->Write("h_tuned_bkg_pzmu_ptmu_1track_blobs");
  bkg_tuned_pzmu_ptmu_1track_michel->Write("h_tuned_bkg_pzmu_ptmu_1track_michel");

  bkg_tuned_pzmu_ptmu_2track_signal->Write("h_tuned_bkg_pzmu_ptmu_2track_signal");
  bkg_tuned_pzmu_ptmu_2track_blobs->Write("h_tuned_bkg_pzmu_ptmu_2track_blobs");
  bkg_tuned_pzmu_ptmu_2track_michel->Write("h_tuned_bkg_pzmu_ptmu_2track_michel");
  /*
  h_recoil_ptmu_data_1track_weighted_bck_signal->Write();
  h_recoil_ptmu_data_2track_weighted_bck_signal->Write();
  h_recoil_ptmu_data_1track_weighted_bck_blobs->Write();
  h_recoil_ptmu_data_2track_weighted_bck_blobs->Write();
  h_recoil_ptmu_data_1track_weighted_bck_michel->Write();
  h_recoil_ptmu_data_2track_weighted_bck_michel->Write();


  h_recoil_ptmu_mc_1track_weighted_bck_signal->Write();
  h_recoil_ptmu_mc_2track_weighted_bck_signal->Write();
  h_recoil_ptmu_mc_1track_weighted_bck_blobs->Write();
  h_recoil_ptmu_mc_2track_weighted_bck_blobs->Write();
  h_recoil_ptmu_mc_1track_weighted_bck_michel->Write();
  h_recoil_ptmu_mc_2track_weighted_bck_michel->Write();


  bkg_tuned_recoil_ptmu_1track_signal->Write("h_tuned_bkg_recoil_ptmu_1track_signal");
  bkg_tuned_recoil_ptmu_1track_blobs->Write("h_tuned_bkg_recoil_ptmu_1track_blobs");
  bkg_tuned_recoil_ptmu_1track_michel->Write("h_tuned_bkg_recoil_ptmu_1track_michel");

  bkg_tuned_recoil_ptmu_2track_signal->Write("h_tuned_bkg_recoil_ptmu_2track_signal");
  bkg_tuned_recoil_ptmu_2track_blobs->Write("h_tuned_bkg_recoil_ptmu_2track_blobs");
  bkg_tuned_recoil_ptmu_2track_michel->Write("h_tuned_bkg_recoil_ptmu_2track_michel");

  h_recoil_inc_ptmu_data_1track_weighted_bck_signal->Write();
  h_recoil_inc_ptmu_data_2track_weighted_bck_signal->Write();
  h_recoil_inc_ptmu_data_1track_weighted_bck_blobs->Write();
  h_recoil_inc_ptmu_data_2track_weighted_bck_blobs->Write();
  h_recoil_inc_ptmu_data_1track_weighted_bck_michel->Write();
  h_recoil_inc_ptmu_data_2track_weighted_bck_michel->Write();


  h_recoil_inc_ptmu_mc_1track_weighted_bck_signal->Write();
  h_recoil_inc_ptmu_mc_2track_weighted_bck_signal->Write();
  h_recoil_inc_ptmu_mc_1track_weighted_bck_blobs->Write();
  h_recoil_inc_ptmu_mc_2track_weighted_bck_blobs->Write();
  h_recoil_inc_ptmu_mc_1track_weighted_bck_michel->Write();
  h_recoil_inc_ptmu_mc_2track_weighted_bck_michel->Write();


  bkg_tuned_recoil_inc_ptmu_1track_signal->Write("h_tuned_bkg_recoil_inc_ptmu_1track_signal");
  bkg_tuned_recoil_inc_ptmu_1track_blobs->Write("h_tuned_bkg_recoil_inc_ptmu_1track_blobs");
  bkg_tuned_recoil_inc_ptmu_1track_michel->Write("h_tuned_bkg_recoil_inc_ptmu_1track_michel");

  bkg_tuned_recoil_inc_ptmu_2track_signal->Write("h_tuned_bkg_recoil_inc_ptmu_2track_signal");
  bkg_tuned_recoil_inc_ptmu_2track_blobs->Write("h_tuned_bkg_recoil_inc_ptmu_2track_blobs");
  bkg_tuned_recoil_inc_ptmu_2track_michel->Write("h_tuned_bkg_recoil_inc_ptmu_2track_michel");
  */
  h_q2_ptmu_data_1track_weighted_bck_signal->Write();
  h_q2_ptmu_data_2track_weighted_bck_signal->Write();
  h_q2_ptmu_data_1track_weighted_bck_blobs->Write();
  h_q2_ptmu_data_2track_weighted_bck_blobs->Write();
  h_q2_ptmu_data_1track_weighted_bck_michel->Write();
  h_q2_ptmu_data_2track_weighted_bck_michel->Write();


  h_q2_ptmu_mc_1track_weighted_bck_signal->Write();
  h_q2_ptmu_mc_2track_weighted_bck_signal->Write();
  h_q2_ptmu_mc_1track_weighted_bck_blobs->Write();
  h_q2_ptmu_mc_2track_weighted_bck_blobs->Write();
  h_q2_ptmu_mc_1track_weighted_bck_michel->Write();
  h_q2_ptmu_mc_2track_weighted_bck_michel->Write();


  bkg_tuned_q2_ptmu_1track_signal->Write("h_tuned_bkg_q2_ptmu_1track_signal");
  bkg_tuned_q2_ptmu_1track_blobs->Write("h_tuned_bkg_q2_ptmu_1track_blobs");
  bkg_tuned_q2_ptmu_1track_michel->Write("h_tuned_bkg_q2_ptmu_1track_michel");

  bkg_tuned_q2_ptmu_2track_signal->Write("h_tuned_bkg_q2_ptmu_2track_signal");
  bkg_tuned_q2_ptmu_2track_blobs->Write("h_tuned_bkg_q2_ptmu_2track_blobs");
  bkg_tuned_q2_ptmu_2track_michel->Write("h_tuned_bkg_q2_ptmu_2track_michel");


  h_enu_ptmu_data_1track_weighted_bck_signal->Write();
  h_enu_ptmu_data_2track_weighted_bck_signal->Write();
  h_enu_ptmu_data_1track_weighted_bck_blobs->Write();
  h_enu_ptmu_data_2track_weighted_bck_blobs->Write();
  h_enu_ptmu_data_1track_weighted_bck_michel->Write();
  h_enu_ptmu_data_2track_weighted_bck_michel->Write();


  h_enu_ptmu_mc_1track_weighted_bck_signal->Write();
  h_enu_ptmu_mc_2track_weighted_bck_signal->Write();
  h_enu_ptmu_mc_1track_weighted_bck_blobs->Write();
  h_enu_ptmu_mc_2track_weighted_bck_blobs->Write();
  h_enu_ptmu_mc_1track_weighted_bck_michel->Write();
  h_enu_ptmu_mc_2track_weighted_bck_michel->Write();


  bkg_tuned_enu_ptmu_1track_signal->Write("h_tuned_bkg_enu_ptmu_1track_signal");
  bkg_tuned_enu_ptmu_1track_blobs->Write("h_tuned_bkg_enu_ptmu_1track_blobs");
  bkg_tuned_enu_ptmu_1track_michel->Write("h_tuned_bkg_enu_ptmu_1track_michel");

  bkg_tuned_enu_ptmu_2track_signal->Write("h_tuned_bkg_enu_ptmu_2track_signal");
  bkg_tuned_enu_ptmu_2track_blobs->Write("h_tuned_bkg_enu_ptmu_2track_blobs");
  bkg_tuned_enu_ptmu_2track_michel->Write("h_tuned_bkg_enu_ptmu_2track_michel");



  for(int i=0;i<nHistos;i++){
    h_pzmu_ptmu_1track_signal[i]->Write(Form("h_pzmu_ptmu_1track_signal_%s",names[i].c_str()));
    h_pzmu_ptmu_2track_signal[i]->Write(Form("h_pzmu_ptmu_2track_signal_%s",names[i].c_str()));
    h_pzmu_ptmu_1track_blobs[i]->Write(Form("h_pzmu_ptmu_1track_blobs_%s",names[i].c_str()));
    h_pzmu_ptmu_2track_blobs[i]->Write(Form("h_pzmu_ptmu_2track_blobs_%s",names[i].c_str()));
    h_pzmu_ptmu_1track_michel[i]->Write(Form("h_pzmu_ptmu_1track_michel_%s",names[i].c_str()));
    h_pzmu_ptmu_2track_michel[i]->Write(Form("h_pzmu_ptmu_2track_michel_%s",names[i].c_str()));

    // h_recoil_ptmu_1track_signal[i]->Write(Form("h_recoil_ptmu_1track_signal_%s",names[i].c_str()));
    // h_recoil_ptmu_2track_signal[i]->Write(Form("h_recoil_ptmu_2track_signal_%s",names[i].c_str()));
    // h_recoil_ptmu_1track_blobs[i]->Write(Form("h_recoil_ptmu_1track_blobs_%s",names[i].c_str()));
    // h_recoil_ptmu_2track_blobs[i]->Write(Form("h_recoil_ptmu_2track_blobs_%s",names[i].c_str()));
    // h_recoil_ptmu_1track_michel[i]->Write(Form("h_recoil_ptmu_1track_michel_%s",names[i].c_str()));
    // h_recoil_ptmu_2track_michel[i]->Write(Form("h_recoil_ptmu_2track_michel_%s",names[i].c_str()));

    // h_recoil_inc_ptmu_1track_signal[i]->Write(Form("h_recoil_inc_ptmu_1track_signal_%s",names[i].c_str()));
    // h_recoil_inc_ptmu_2track_signal[i]->Write(Form("h_recoil_inc_ptmu_2track_signal_%s",names[i].c_str()));
    // h_recoil_inc_ptmu_1track_blobs[i]->Write(Form("h_recoil_inc_ptmu_1track_blobs_%s",names[i].c_str()));
    // h_recoil_inc_ptmu_2track_blobs[i]->Write(Form("h_recoil_inc_ptmu_2track_blobs_%s",names[i].c_str()));
    // h_recoil_inc_ptmu_1track_michel[i]->Write(Form("h_recoil_inc_ptmu_1track_michel_%s",names[i].c_str()));
    // h_recoil_inc_ptmu_2track_michel[i]->Write(Form("h_recoil_inc_ptmu_2track_michel_%s",names[i].c_str()));


    h_q2_ptmu_1track_signal[i]->Write(Form("h_q2_ptmu_1track_signal_%s",names[i].c_str()));
    h_q2_ptmu_2track_signal[i]->Write(Form("h_q2_ptmu_2track_signal_%s",names[i].c_str()));
    h_q2_ptmu_1track_blobs[i]->Write(Form("h_q2_ptmu_1track_blobs_%s",names[i].c_str()));
    h_q2_ptmu_2track_blobs[i]->Write(Form("h_q2_ptmu_2track_blobs_%s",names[i].c_str()));
    h_q2_ptmu_1track_michel[i]->Write(Form("h_q2_ptmu_1track_michel_%s",names[i].c_str()));
    h_q2_ptmu_2track_michel[i]->Write(Form("h_q2_ptmu_2track_michel_%s",names[i].c_str()));


    h_enu_ptmu_1track_signal[i]->Write(Form("h_enu_ptmu_1track_signal_%s",names[i].c_str()));
    h_enu_ptmu_2track_signal[i]->Write(Form("h_enu_ptmu_2track_signal_%s",names[i].c_str()));
    h_enu_ptmu_1track_blobs[i]->Write(Form("h_enu_ptmu_1track_blobs_%s",names[i].c_str()));
    h_enu_ptmu_2track_blobs[i]->Write(Form("h_enu_ptmu_2track_blobs_%s",names[i].c_str()));
    h_enu_ptmu_1track_michel[i]->Write(Form("h_enu_ptmu_1track_michel_%s",names[i].c_str()));
    h_enu_ptmu_2track_michel[i]->Write(Form("h_enu_ptmu_2track_michel_%s",names[i].c_str()));

  }


  f_bg_weights->Close();
  f_2track_signal->Close();
  f_1track_signal->Close();
  f_2track_blobs->Close();
  f_1track_blobs->Close();
  f_2track_michel->Close();
  f_1track_michel->Close();
  f_output->Close();

  delete f_output; 
  delete f_1track_signal;
  delete f_2track_signal;

  delete f_1track_blobs;
  delete f_2track_blobs;

  delete f_1track_michel;
  delete f_2track_michel;
  delete f_bg_weights; 

  delete plotUtils; 
  delete utils;

  cout<< "All done..."  << endl;
  return 0;
};




void ScaleBackground(MnvH2D** h, MnvH2D* h_bg_weights, MnvH2D* &h_data_nobck, MnvH2D* &h_mc_nobck, MnvH2D* &h_tuned_bkg){
  h_tuned_bkg = (MnvH2D*)h[kQELikeNot]->Clone( Form( "%s_tunedbck", h[kQELikeNot]->GetName()) );
  h_mc_nobck   = (MnvH2D*)h[kMC]->Clone( Form( "%s_nobck_mc", h[kMC]->GetName()) );
  h_data_nobck = (MnvH2D*)h[kData]->Clone( Form( "%s_nobck_data", h[kData]->GetName()) );
  h_data_nobck->AddMissingErrorBandsAndFillWithCV( *h_mc_nobck );//<---- make sure errors continue on! All data universes = 1.0 wgt.
  h_tuned_bkg->Multiply(h_tuned_bkg,h_bg_weights);//This is the data driven background!
  h_data_nobck->Add(h_tuned_bkg,-1.0);
  h_mc_nobck->Add( h[kQELikeNot], -1.0);

}


int main( int argc, char *argv[])
{
  ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if (argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t./SideBandFitHists  Output_file  Input_file_1TRACK  Input_file_2TRACK  Input_file_BKG  Input_file_MIGRATION  Input_file_EFFICIENCY  Input_file_FLUX  Num_iterations  Make_Flux_Constraint_Histo  Apply_Flux_Constraint \n\n"<<
      "\t-output_file\t =\t Cross Section output file\n"<<
      "\t-input_file_1track\t =\t Muon Event Selection(Signal) for 1 track sub-sample\n"<<
      "\t-input_file_2track\t =\t Muon Event Selection(Signal) for >=2 track sub-sample\n"<<
      "\t-input_file_1track\t =\t Muon Event Selection(Blobs) for 1 track sub-sample\n"<<
      "\t-input_file_2track\t =\t Muon Event Selection(Blobs) for >=2 track sub-sample\n"<<
      "\t-input_file_1track\t =\t Muon Event Selection(Michel) for 1 track sub-sample\n"<<
      "\t-input_file_2track\t =\t Muon Event Selection(Michel) for >=2 track sub-sample\n"<<
      "\t-input_file_bkg\t =\t Background weights filename\n"<<
      "\t-Make_Flux_Constraint_Histo\t =\t If TRUE Enter 1; If FALSE Enter 0\n" << 
      "\t-Apply_Flux_Constraint\t =\t If TRUE Enter 1; If FALSE Enter 0\n" << 
      "\t********************************************************************************************** \n"<<
      "\t Please see : MuonSelectionHists.cxx, BackgroundWeights.cxx, MigrationMatrixHists.cxx, EffPurityHists.cxx and Ana/Flux/python/compute_flux.py for getting the necessary input files"<< std::endl; 
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    return 0; 
  }

  //! Default parameters
  std::vector<std::string> par;
  par.push_back("SideBandFitPlots");
  par.push_back( Form("%s/ana/rootfiles/SideBandFitPlots.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/MuonSelectionHists_Signal_1track.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/MuonSelectionHists_Signal_2track.root",getenv("CCQENUROOT") ) );

  par.push_back( Form("%s/ana/rootfiles/MuonSelectionHists_Blobs_1track.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/MuonSelectionHists_Blobs_2track.root",getenv("CCQENUROOT") ) );

  par.push_back( Form("%s/ana/rootfiles/MuonSelectionHists_Michel_1track.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/MuonSelectionHists_Michel_2track.root",getenv("CCQENUROOT") ) );

  par.push_back( Form("%s/ana/rootfiles/BackgroundWeights.root",getenv("CCQENUROOT") ) );
  par.push_back("0"); 
  par.push_back("0"); 

  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;

  bool fluxConstraintHistoNeeded = ( *(par.end()-2) == "1" ) ? true : false; 

  bool applyFluxConstraint = ( par.back() == "1" ) ? true : false; 

  return SideBandFitHists( par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], fluxConstraintHistoNeeded, applyFluxConstraint );

}
