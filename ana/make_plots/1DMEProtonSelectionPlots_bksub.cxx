#include "include/CCQENuPlotUtils.h"
#include "TParameter.h"

using namespace CCQENU_ANA;

void drawQ2Lines(CCQENuPlotUtils *utils, int color = kBlack, Option_t* opts = "");
void drawEnuLines(CCQENuPlotUtils *utils, int color = kBlue, Option_t* opts = "");

void plottingOneInBinsOfAnother( MnvH2D *de2DHisto[], TString plotVariable, TString norm, bool areaNormalizeEntireHisto, MnvPlotter *plotter, CCQENuPlotUtils *utils, string style, bool includeData, TString sample, TString multiplicityLabel, double pot_data, double pot_mc, bool draw_with_sys_band ); 


int ProtonSelectionPlots( string filename, string filename_orig, string style = "QE")
{
  // read file to get histograms
  cout<<"Begin ProtonSelectionPlots"<<endl;
  bool area_norm  = false;
  TFile *f = new TFile( filename.c_str(), "READ" );
  TFile *f_orig = new TFile( filename_orig.c_str(), "READ" );
  if (f->IsZombie() || f->GetListOfKeys()->IsEmpty()){
    Error("VertexEnergyPlots","Could not get histogram ROOT file or it was empty.");
    return 1;
  }

  //-------------------------------------------
  // Get event counts before cuts for norm
  //-------------------------------------------
  TVector2 *evt = (TVector2*)f_orig->Get("n_events");
  double data_events = evt->X();
  double mc_events = evt->Y();
  
  cout<< " Number of Data Events = " << data_events << endl;
  cout<< " Number of MC Events = " <<   mc_events << endl;

  //-----------------------------------------------------------
  // Will the input ROOT file be used for flux constraints  
  //-----------------------------------------------------------
  bool fluxHistoExists = false; 
  TParameter<bool> *fluxParam = (TParameter<bool>*)f_orig->Get("fluxConstraintHisto"); 
  if( fluxParam==NULL ) { 
    cout << "Input ROOT file doesn't have the necessary histos for the flux constraint !" << endl;
  } else { 
    if ( 0==fluxParam->GetVal() ) { 
      cout << "Input ROOT file might have the necessary flux histos but does not have its flux constrained !" << endl;  
    }
    fluxHistoExists = true; 
  }

  cout << "fluxHistoExists " << fluxHistoExists << endl; 

  CCQENuPlotUtils *utils = new CCQENuPlotUtils( fluxHistoExists );
  //CCQENuPlotUtils *utils = new CCQENuPlotUtils( );
  MnvPlotter *plotter = new MnvPlotter;

  TCanvas* c;

  TString norm = area_norm ? "area" : "pot";

  double pot_data = utils->getPOTData( f_orig );
  double pot_mc = utils->getPOTMC( f_orig );
  double pot_scale = utils->getPOTNormFactor( f_orig );
  if(pot_scale==0) pot_scale = 1.0;
  cout<< "POT DATA: " << pot_data << " POT MC: " << pot_mc << " POT Scale Factor = " << pot_scale << endl; 
 
  //--------------------------------------------
  // Load histos to plot 
  //--------------------------------------------
  //---------------------------------------------
  // Neutron Angular Histos
  //---------------------------------------------
  MnvH2D* h_dthetaR_ptmu_2track_signal[nHistos],  *h_dthetaR_ptmu_2track_blobs[nHistos],  *h_dthetaR_ptmu_2track_michel[nHistos], *h_dthetaR_ptmu_2track_micblob[nHistos]; 
  MnvH2D* h_tuned_bkg_dthetaR_ptmu_2track_signal, *h_tuned_bkg_dthetaR_ptmu_2track_blobs, *h_tuned_bkg_dthetaR_ptmu_2track_michel,*h_tuned_bkg_dthetaR_ptmu_2track_micblob;
  MnvH2D* h_data_nobkg_dthetaR_ptmu_2track_signal,*h_data_bkg_dthetaR_ptmu_2track_blobs,  *h_data_bkg_dthetaR_ptmu_2track_michel, *h_data_bkg_dthetaR_ptmu_2track_micblob; 
  MnvH2D* h_mc_nobkg_dthetaR_ptmu_2track_signal,  *h_mc_bkg_dthetaR_ptmu_2track_blobs,    *h_mc_bkg_dthetaR_ptmu_2track_michel,   *h_mc_bkg_dthetaR_ptmu_2track_micblob;   

  MnvH2D* h_dthetaP_ptmu_2track_signal[nHistos],  *h_dthetaP_ptmu_2track_blobs[nHistos],  *h_dthetaP_ptmu_2track_michel[nHistos], *h_dthetaP_ptmu_2track_micblob[nHistos]; 
  MnvH2D* h_tuned_bkg_dthetaP_ptmu_2track_signal, *h_tuned_bkg_dthetaP_ptmu_2track_blobs, *h_tuned_bkg_dthetaP_ptmu_2track_michel,*h_tuned_bkg_dthetaP_ptmu_2track_micblob;
  MnvH2D* h_data_nobkg_dthetaP_ptmu_2track_signal,*h_data_bkg_dthetaP_ptmu_2track_blobs,  *h_data_bkg_dthetaP_ptmu_2track_michel, *h_data_bkg_dthetaP_ptmu_2track_micblob; 
  MnvH2D* h_mc_nobkg_dthetaP_ptmu_2track_signal,  *h_mc_bkg_dthetaP_ptmu_2track_blobs,    *h_mc_bkg_dthetaP_ptmu_2track_michel,   *h_mc_bkg_dthetaP_ptmu_2track_micblob;   


  MnvH1D* h_dthetaP_signal[nHistos], *h_dthetaR_signal[nHistos];
  MnvH1D* h_dthetaP_blobs[nHistos], *h_dthetaR_blobs[nHistos];
  MnvH1D* h_dthetaP_michel[nHistos], *h_dthetaR_michel[nHistos];
  MnvH1D* h_dthetaP_micblob[nHistos], *h_dthetaR_micblob[nHistos];

  utils->bookHistos( f, h_dthetaR_ptmu_2track_signal,"h_dthetaR_ptmu_2track_signal");
  utils->bookHistos( f, h_dthetaR_ptmu_2track_blobs ,"h_dthetaR_ptmu_2track_blobs");
  utils->bookHistos( f, h_dthetaR_ptmu_2track_michel,"h_dthetaR_ptmu_2track_michel");
  utils->bookHistos( f, h_dthetaR_ptmu_2track_micblob ,"h_dthetaR_ptmu_2track_micblob");

  utils->bookHistos( f, h_dthetaP_ptmu_2track_signal,"h_dthetaP_ptmu_2track_signal");
  utils->bookHistos( f, h_dthetaP_ptmu_2track_blobs ,"h_dthetaP_ptmu_2track_blobs");
  utils->bookHistos( f, h_dthetaP_ptmu_2track_michel,"h_dthetaP_ptmu_2track_michel");
  utils->bookHistos( f, h_dthetaP_ptmu_2track_micblob ,"h_dthetaP_ptmu_2track_micblob");

  
  h_dthetaR_ptmu_2track_signal[kQELikeNot]->Reset();  
  h_dthetaR_ptmu_2track_blobs[kQELikeNot]->Reset();
  h_dthetaR_ptmu_2track_michel[kQELikeNot]->Reset();
  h_dthetaR_ptmu_2track_micblob[kQELikeNot]->Reset();

  h_dthetaP_ptmu_2track_signal[kQELikeNot]->Reset();
  h_dthetaP_ptmu_2track_blobs[kQELikeNot]->Reset();
  h_dthetaP_ptmu_2track_michel[kQELikeNot]->Reset();
  h_dthetaP_ptmu_2track_micblob[kQELikeNot]->Reset();

  //utils->scaleMCHistos( h_dthetaP_ptmu_2track_signal, pot_scale);
  //utils->scaleMCHistos( h_dthetaP_ptmu_2track_blobs, pot_scale);
  //utils->scaleMCHistos( h_dthetaP_ptmu_2track_michel, pot_scale);

  //utils->scaleMCHistos( h_dthetaR_ptmu_2track_signal, pot_scale);
  //utils->scaleMCHistos( h_dthetaR_ptmu_2track_blobs, pot_scale);
  //utils->scaleMCHistos( h_dthetaR_ptmu_2track_michel, pot_scale);

  h_tuned_bkg_dthetaR_ptmu_2track_signal  =   (MnvH2D*) f->Get("h_tuned_bkg_dthetaR_ptmu_2track_signal");  
  h_tuned_bkg_dthetaR_ptmu_2track_blobs   =   (MnvH2D*) f->Get("h_tuned_bkg_dthetaR_ptmu_2track_blobs");
  h_tuned_bkg_dthetaR_ptmu_2track_michel  =   (MnvH2D*) f->Get("h_tuned_bkg_dthetaR_ptmu_2track_michel");
  h_tuned_bkg_dthetaR_ptmu_2track_micblob   =   (MnvH2D*) f->Get("h_tuned_bkg_dthetaR_ptmu_2track_micblob");
                                          
  h_tuned_bkg_dthetaP_ptmu_2track_signal  =   (MnvH2D*) f->Get("h_tuned_bkg_dthetaP_ptmu_2track_signal");  
  h_tuned_bkg_dthetaP_ptmu_2track_blobs   =   (MnvH2D*) f->Get("h_tuned_bkg_dthetaP_ptmu_2track_blobs");
  h_tuned_bkg_dthetaP_ptmu_2track_michel  =   (MnvH2D*) f->Get("h_tuned_bkg_dthetaP_ptmu_2track_michel");
  h_tuned_bkg_dthetaP_ptmu_2track_micblob   =   (MnvH2D*) f->Get("h_tuned_bkg_dthetaP_ptmu_2track_micblob");

  //================================================================
  h_dthetaR_ptmu_2track_signal[kData] = (MnvH2D*) f->Get("h_dthetaR_ptmu_data_2track_weighted_nobck_signal");
  h_dthetaR_ptmu_2track_blobs[kData]= (MnvH2D*) f->Get("h_dthetaR_ptmu_data_2track_weighted_nobck_blobs");
  h_dthetaR_ptmu_2track_michel[kData]= (MnvH2D*) f->Get("h_dthetaR_ptmu_data_2track_weighted_nobck_michel");
  h_dthetaR_ptmu_2track_micblob[kData]= (MnvH2D*) f->Get("h_dthetaR_ptmu_data_2track_weighted_nobck_micblob");

  h_dthetaR_ptmu_2track_signal[kMC] = (MnvH2D*) f->Get("h_dthetaR_ptmu_mc_2track_weighted_nobck_signal");
  h_dthetaR_ptmu_2track_blobs[kMC]= (MnvH2D*) f->Get("h_dthetaR_ptmu_mc_2track_weighted_nobck_blobs");
  h_dthetaR_ptmu_2track_michel[kMC]= (MnvH2D*) f->Get("h_dthetaR_ptmu_mc_2track_weighted_nobck_michel");
  h_dthetaR_ptmu_2track_micblob[kMC]= (MnvH2D*) f->Get("h_dthetaR_ptmu_mc_2track_weighted_nobck_micblob");


  h_dthetaP_ptmu_2track_signal[kData] = (MnvH2D*) f->Get("h_dthetaP_ptmu_data_2track_weighted_nobck_signal");
  h_dthetaP_ptmu_2track_blobs[kData]= (MnvH2D*) f->Get("h_dthetaP_ptmu_data_2track_weighted_nobck_blobs");
  h_dthetaP_ptmu_2track_michel[kData]= (MnvH2D*) f->Get("h_dthetaP_ptmu_data_2track_weighted_nobck_michel");
  h_dthetaP_ptmu_2track_micblob[kData]= (MnvH2D*) f->Get("h_dthetaP_ptmu_data_2track_weighted_nobck_micblob");

  h_dthetaP_ptmu_2track_signal[kMC] = (MnvH2D*) f->Get("h_dthetaP_ptmu_mc_2track_weighted_nobck_signal");
  h_dthetaP_ptmu_2track_blobs[kMC]= (MnvH2D*) f->Get("h_dthetaP_ptmu_mc_2track_weighted_nobck_blobs");
  h_dthetaP_ptmu_2track_michel[kMC]= (MnvH2D*) f->Get("h_dthetaP_ptmu_mc_2track_weighted_nobck_michel");
  h_dthetaP_ptmu_2track_micblob[kMC]= (MnvH2D*) f->Get("h_dthetaP_ptmu_mc_2track_weighted_nobck_micblob");





  for(int i=0;i<nHistos;i++){
    h_dthetaR_signal[i] =  h_dthetaR_ptmu_2track_signal[i]->ProjectionX(Form("%s_projX",h_dthetaR_ptmu_2track_signal[i]->GetName()),0,-1);
    h_dthetaR_blobs[i] =  h_dthetaR_ptmu_2track_blobs[i]->ProjectionX(Form("%s_projX",h_dthetaR_ptmu_2track_blobs[i]->GetName()),0,-1);
    h_dthetaR_michel[i] =  h_dthetaR_ptmu_2track_michel[i]->ProjectionX(Form("%s_projX",h_dthetaR_ptmu_2track_michel[i]->GetName()),0,-1);
    h_dthetaR_micblob[i] =  h_dthetaR_ptmu_2track_micblob[i]->ProjectionX(Form("%s_projX",h_dthetaR_ptmu_2track_micblob[i]->GetName()),0,-1);

    h_dthetaP_signal[i] =  h_dthetaP_ptmu_2track_signal[i]->ProjectionX(Form("%s_projX",h_dthetaP_ptmu_2track_signal[i]->GetName()),0,-1);
    h_dthetaP_blobs[i] =  h_dthetaP_ptmu_2track_blobs[i]->ProjectionX(Form("%s_projX",h_dthetaP_ptmu_2track_blobs[i]->GetName()),0,-1);
    h_dthetaP_michel[i] =  h_dthetaP_ptmu_2track_michel[i]->ProjectionX(Form("%s_projX",h_dthetaP_ptmu_2track_michel[i]->GetName()),0,-1);
    h_dthetaP_micblob[i] =  h_dthetaP_ptmu_2track_micblob[i]->ProjectionX(Form("%s_projX",h_dthetaP_ptmu_2track_micblob[i]->GetName()),0,-1);

    /*
    int bmin = h_neutron_dperp[i]->FindBin(-.02);
    int bmax = h_neutron_dperp[i]->FindBin(.02);
    h_neutron_dreact_dperpcut[i] =  h_neutron_dperp_dreactplane[i]->ProjectionX(Form("%s_projX",h_neutron_dreactplane_qsq[i]->GetName()),bmin,bmax);
    */
  }



  //--------------------------------------------
  // Plot Muon Kinematics...
  //--------------------------------------------
  cout << "Plotting Proton Angulars..." << endl;
  c = new TCanvas("cNeutA","Neutron Angulars"); 
  bool includeData = true;

  //Settings   
  TGaxis::SetMaxDigits(2);
  //2-D
  plotter->SetRedHeatPalette();

  //-----------------------------------------------------------------------------
  //Make plots of the RAW event counts before POT and binwidth normalizations 
  //-----------------------------------------------------------------------------
   
   //--------------------------------------------
   // Get Area or POT Normalization
   //--------------------------------------------

   //1-D
   //double tmu_scale                          = (area_norm)? utils->getAreaNormFactor(h_tmu)   : pot_scale;
   //utils->scaleMCHistos( h_tmu, tmu_scale);
   //utils->scaleMCHistos( h_dthetaP_signal, pot_scale);
   //utils->scaleMCHistos( h_dthetaP_blobs, pot_scale);
   //utils->scaleMCHistos( h_dthetaP_michel, pot_scale);

   //utils->scaleMCHistos( h_dthetaR_signal, pot_scale);
   //utils->scaleMCHistos( h_dthetaR_blobs, pot_scale);
   //utils->scaleMCHistos( h_dthetaR_michel, pot_scale);

   cout<<"Finished scaling histos"<<endl;
  //==============================================================  
  //1-D Data with MC broken down into it's different templates
  //==============================================================
  cout<<h_dthetaP_signal[kData]->GetName()<<endl;
  utils->drawStacked( h_dthetaP_signal, style, area_norm, pot_data, pot_mc, includeData );
  utils->Print( c, Form("ProtonSelection%s/dthetaP_signal_%s_%s", "Signal", "2", "pot") );
  utils->drawStacked( h_dthetaP_blobs, style, area_norm, pot_data, pot_mc, includeData );
  utils->Print( c, Form("ProtonSelection%s/dthetaP_signal_%s_%s", "BlobSideBand", "2", "pot") );
  utils->drawStacked( h_dthetaP_michel, style, area_norm, pot_data, pot_mc, includeData );
  utils->Print( c, Form("ProtonSelection%s/dthetaP_signal_%s_%s", "MichelSideBand", "2", "pot") );
  utils->drawStacked( h_dthetaP_micblob, style, area_norm, pot_data, pot_mc, includeData );
  utils->Print( c, Form("ProtonSelection%s/dthetaP_signal_%s_%s", "MicBlobSideBand", "2", "pot") );

  utils->drawStacked( h_dthetaR_signal, style, area_norm, pot_data, pot_mc, includeData );
  utils->Print( c, Form("ProtonSelection%s/dthetaR_signal_%s_%s", "Signal", "2", "pot") );
  utils->drawStacked( h_dthetaR_blobs, style, area_norm, pot_data, pot_mc, includeData );
  utils->Print( c, Form("ProtonSelection%s/dthetaR_signal_%s_%s", "BlobSideBand", "2", "pot") );
  utils->drawStacked( h_dthetaR_michel, style, area_norm, pot_data, pot_mc, includeData );
  utils->Print( c, Form("ProtonSelection%s/dthetaR_signal_%s_%s", "MichelSideBand", "2", "pot") );
  utils->drawStacked( h_dthetaR_micblob, style, area_norm, pot_data, pot_mc, includeData );
  utils->Print( c, Form("ProtonSelection%s/dthetaR_signal_%s_%s", "MicBlobSideBand", "2", "pot") );


   cout<<"Drawn 1D"<<endl;


  //c->SetLogy(1);
  //utils->drawStacked( h_multiplicity, style, area_norm, pot_data, pot_mc, includeData );
  //utils->Print( c, Form("ProtonSelection%s/multiplicity_%s_%s_logy", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  //c->SetLogy(0);


/*
  const int first_bin = 1;
  const int last_bin  = h_q22[kMC]->GetNbinsX()-2;

  h_q22[kMC]->GetXaxis()->SetRange( first_bin, last_bin );
  h_q22[kData]->GetXaxis()->SetRange( first_bin, last_bin );
  h_q22_old[kMC]->GetXaxis()->SetRange( first_bin, last_bin );
  h_q22_old[kData]->GetXaxis()->SetRange( first_bin, last_bin );

  utils->drawStacked( h_q22, style, area_norm, pot_data, pot_mc, includeData );
  utils->Print( c, Form("ProtonSelection%s/q22_bin_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  utils->drawStacked( h_q22_old, style, area_norm, pot_data, pot_mc, includeData );
  utils->Print( c, Form("ProtonSelection%s/q22_old_bin_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
*/
  //====================================================================
  // For drawing any systematics or error summaries take projections 
  // of 2D histos (1D histos don't have error bands on them any more) 
  //====================================================================
  //MnvH1D *h_pzmu_mc_temp = (MnvH1D*)h_pzmu_ptmu[kMC]->ProjectionX();
  //MnvH1D *h_ptmu_mc_temp = (MnvH1D*)h_pzmu_ptmu[kMC]->ProjectionY();

  //MnvH1D *h_pzmu_data_temp = (MnvH1D*)h_pzmu_ptmu[kData]->ProjectionX();
  //MnvH1D *h_ptmu_data_temp = (MnvH1D*)h_pzmu_ptmu[kData]->ProjectionY();


  //MnvH1D *h_tmu_mc_temp = (MnvH1D*)h_tmu_theta[kMC]->ProjectionX();
  //MnvH1D *h_theta_mc_temp = (MnvH1D*)h_tmu_theta[kMC]->ProjectionY();

  //MnvH1D *h_tmu_data_temp = (MnvH1D*)h_tmu_theta[kData]->ProjectionX();
  //MnvH1D *h_theta_data_temp = (MnvH1D*)h_tmu_theta[kData]->ProjectionY();


  
  //==================================================
  // 1-D Plot ratios 
  //==================================================
  //plotter->DrawDataMCRatio(h_ptmu_data_temp, h_ptmu_mc_temp, 1.0, true, true, -1.5, 1.5, "Data / MC", area_norm);  
  //utils->Print( c, Form("ProtonSelection%s/muon_pt_data_mc_ratio_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) ); 

  plotter->DrawDataMCRatio(h_dthetaP_signal[kData], h_dthetaP_signal[kMC], 1.0, true, true, -1.5, 1.5, "Data / MC", area_norm);  
  utils->Print( c, Form("ProtonSelection%s/h_dthetaP_ratio_%s_%s", "Signal", "2", "pot" )); 
  plotter->DrawDataMCRatio(h_dthetaP_blobs[kData], h_dthetaP_blobs[kMC], 1.0, true, true, -1.5, 1.5, "Data / MC", area_norm);  
  utils->Print( c, Form("ProtonSelection%s/h_dthetaP_ratio_%s_%s", "BlobSideBand", "2", "pot" )); 
  plotter->DrawDataMCRatio(h_dthetaP_michel[kData], h_dthetaP_michel[kMC], 1.0, true, true, -1.5, 1.5, "Data / MC", area_norm);  
  utils->Print( c, Form("ProtonSelection%s/h_dthetaP_ratio_%s_%s", "MichelSideBand", "2", "pot" )); 
  plotter->DrawDataMCRatio(h_dthetaP_micblob[kData], h_dthetaP_micblob[kMC], 1.0, true, true, -1.5, 1.5, "Data / MC", area_norm);  
  utils->Print( c, Form("ProtonSelection%s/h_dthetaP_ratio_%s_%s", "MicBlobSideBand", "2", "pot" )); 

  plotter->DrawDataMCRatio(h_dthetaR_signal[kData], h_dthetaR_signal[kMC], 1.0, true, true, -1.5, 1.5, "Data / MC", area_norm);  
  utils->Print( c, Form("ProtonSelection%s/h_dthetaR_ratio_%s_%s", "Signal", "2", "pot" )); 
  plotter->DrawDataMCRatio(h_dthetaR_blobs[kData], h_dthetaR_blobs[kMC], 1.0, true, true, -1.5, 1.5, "Data / MC", area_norm);  
  utils->Print( c, Form("ProtonSelection%s/h_dthetaR_ratio_%s_%s", "BlobSideBand", "2", "pot" )); 
  plotter->DrawDataMCRatio(h_dthetaR_michel[kData], h_dthetaR_michel[kMC], 1.0, true, true, -1.5, 1.5, "Data / MC", area_norm);  
  utils->Print( c, Form("ProtonSelection%s/h_dthetaR_ratio_%s_%s", "MichelSideBand", "2", "pot" )); 
  plotter->DrawDataMCRatio(h_dthetaR_micblob[kData], h_dthetaR_micblob[kMC], 1.0, true, true, -1.5, 1.5, "Data / MC", area_norm);  
  utils->Print( c, Form("ProtonSelection%s/h_dthetaR_ratio_%s_%s", "MicBlobSideBand", "2", "pot" )); 




  //==================================================
  // 1-D of Data with MC with stats + systematics
  // Plots of T_mu, Theta_mu, pT_mu and pZ_mu 
  //==================================================
  plotter->DrawDataMCWithErrorBand(h_dthetaR_signal[kData], h_dthetaR_signal[kMC], 1.0, "TR", false, NULL, NULL, area_norm, true);
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/h_dthetaR_signal_sys_%s_%s",  "Signal", "2", "pot" ) );
  plotter->DrawDataMCWithErrorBand(h_dthetaR_blobs[kData], h_dthetaR_blobs[kMC], 1.0, "TR", false, NULL, NULL, area_norm, true);
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/h_dthetaR_blobs_sys_%s_%s",  "BlobSideBand", "2", "pot" ) );
  plotter->DrawDataMCWithErrorBand(h_dthetaR_michel[kData], h_dthetaR_michel[kMC], 1.0, "TR", false, NULL, NULL, area_norm, true);
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/h_dthetaR_signal_sys_%s_%s",  "MichelSideBand", "2", "pot" ) );
  plotter->DrawDataMCWithErrorBand(h_dthetaR_micblob[kData], h_dthetaR_micblob[kMC], 1.0, "TR", false, NULL, NULL, area_norm, true);
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/h_dthetaR_micblob_sys_%s_%s",  "MicBlobSideBand", "2", "pot" ) );

  plotter->DrawDataMCWithErrorBand(h_dthetaP_signal[kData], h_dthetaP_signal[kMC], 1.0, "TR", false, NULL, NULL, area_norm, true);
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/h_dthetaP_signal_sys_%s_%s",  "Signal", "2", "pot" ) );
  plotter->DrawDataMCWithErrorBand(h_dthetaP_blobs[kData], h_dthetaP_blobs[kMC], 1.0, "TR", false, NULL, NULL, area_norm, true);
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/h_dthetaP_blobs_sys_%s_%s",  "BlobSideBand", "2", "pot" ) );
  plotter->DrawDataMCWithErrorBand(h_dthetaP_michel[kData], h_dthetaP_michel[kMC], 1.0, "TR", false, NULL, NULL, area_norm, true);
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/h_dthetaP_signal_sys_%s_%s",  "MichelSideBand", "2", "pot" ) );
  plotter->DrawDataMCWithErrorBand(h_dthetaP_micblob[kData], h_dthetaP_micblob[kMC], 1.0, "TR", false, NULL, NULL, area_norm, true);
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/h_dthetaP_micblob_sys_%s_%s",  "MicBlobSideBand", "2", "pot" ) );



  //======================================================
  // Plot Fractional Systematic Error Summary 
  // Plots of pT_mu and pZ_mu (errors on MC right now) 
  //======================================================
  plotter->error_summary_group_map = utils->getSystematicGroupMap();
  plotter->error_color_map         = utils->getSystematicGroupMapColors();
  TGaxis::SetMaxDigits(5);
  //--------------------------------------------
  // DATA-MC Error Summary for Muon pT and pZ 
  //--------------------------------------------
  unsigned int main_hist = kData;
  plotter->DrawErrorSummary(h_dthetaR_signal[main_hist], "TR", true, true, 0.00001, area_norm );
  plotter->AddPlotLabel( Form("%s #bullet%s", "dthetaR" , (area_norm)? "Shape Errors" : "Absolute Errors"), 0.33, 0.94, 0.025 );
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/h_dthetaR_signal_with_sys_errorsummary_%s_%s", "Signal", "2", "pot" ) );

  plotter->DrawErrorSummary(h_dthetaP_signal[main_hist], "TR", true, true, 0.00001, area_norm );
  plotter->AddPlotLabel( Form("%s #bullet%s", "dthetaP" , (area_norm)? "Shape Errors" : "Absolute Errors"), 0.33, 0.94, 0.025 );
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/h_dthetaP_signal_with_sys_errorsummary_%s_%s",   "Signal", "2", "pot" ) );


  plotter->DrawErrorSummary(h_dthetaR_blobs[main_hist], "TR", true, true, 0.00001, area_norm );
  plotter->AddPlotLabel( Form("%s #bullet%s", "dthetaR" , (area_norm)? "Shape Errors" : "Absolute Errors"), 0.33, 0.94, 0.025 );
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/h_dthetaR_blobs_with_sys_errorsummary_%s_%s", "BlobSideBand", "2", "pot" ) );

  plotter->DrawErrorSummary(h_dthetaP_blobs[main_hist], "TR", true, true, 0.00001, area_norm );
  plotter->AddPlotLabel( Form("%s #bullet%s", "dthetaP" , (area_norm)? "Shape Errors" : "Absolute Errors"), 0.33, 0.94, 0.025 );
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/h_dthetaP_blobs_with_sys_errorsummary_%s_%s",   "BlobSideBand", "2", "pot" ) );
  

  plotter->DrawErrorSummary(h_dthetaR_michel[main_hist], "TR", true, true, 0.00001, area_norm );
  plotter->AddPlotLabel( Form("%s #bullet%s", "dthetaR" , (area_norm)? "Shape Errors" : "Absolute Errors"), 0.33, 0.94, 0.025 );
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/h_dthetaR_michel_with_sys_errorsummary_%s_%s", "MichelSideBand", "2", "pot" ) );

  plotter->DrawErrorSummary(h_dthetaP_michel[main_hist], "TR", true, true, 0.00001, area_norm );
  plotter->AddPlotLabel( Form("%s #bullet%s", "dthetaP" , (area_norm)? "Shape Errors" : "Absolute Errors"), 0.33, 0.94, 0.025 );
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/h_dthetaP_michel_with_sys_errorsummary_%s_%s",   "MichelSideBand", "2", "pot" ) );


  plotter->DrawErrorSummary(h_dthetaR_micblob[main_hist], "TR", true, true, 0.00001, area_norm );
  plotter->AddPlotLabel( Form("%s #bullet%s", "dthetaR" , (area_norm)? "Shape Errors" : "Absolute Errors"), 0.33, 0.94, 0.025 );
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/h_dthetaR_micblob_with_sys_errorsummary_%s_%s", "MicBlobSideBand", "2", "pot" ) );

  plotter->DrawErrorSummary(h_dthetaP_micblob[main_hist], "TR", true, true, 0.00001, area_norm );
  plotter->AddPlotLabel( Form("%s #bullet%s", "dthetaP" , (area_norm)? "Shape Errors" : "Absolute Errors"), 0.33, 0.94, 0.025 );
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/h_dthetaP_micblob_with_sys_errorsummary_%s_%s",   "MicBlobSideBand", "2", "pot" ) );





   //------------------------------------------
  // Plot individual systematic groups now 
  //------------------------------------------
  //double old_headroom = plotter->headroom;
  //int old_legend_n_columns = plotter->legend_n_columns;
  plotter->headroom = 1.0;
  plotter->legend_n_columns = 2;
  plotter->axis_maximum = 0.4; 
  plotter->axis_minimum = 0.00001; 

  //--------------------------------------------------------------------------------------------------
  // Drawing error summaries in the following groups "as fractional uncertainties" 
  // This was developed for testing the effects (or not) of the flux constraint implementations  
  // Nowadays it is used for viewing the components of an error summary band 
  //--------------------------------------------------------------------------------------------------
  //plotter->axis_maximum_group = 0.2; 
  //plotter->DrawErrorSummary(h_pzmu_mc_temp, "TL", false, true, 0.00001, area_norm, "Recoil Reconstruction" );
  //utils->Print( c, Form("ProtonSelection%s/muon_pz_with_frac_sys_errorsummary_RecoilRecons_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  

  plotter->axis_maximum_group = 0.2; 
  std::vector<std::string> vertErrNames = h_dthetaR_signal[main_hist]->GetVertErrorBandNames();
  std::vector<std::string> latErrNames = h_dthetaR_signal[main_hist]->GetLatErrorBandNames();
  std::vector<std::string>::iterator errName;
  std::cout<<"Vertical Error Bands: "<<std::endl;
  for( errName = vertErrNames.begin(); errName != vertErrNames.end();++errName ) std::cout<<*errName<<std::endl;

  std::cout<<"Lateral Error Bands: "<<std::endl;
  for( errName = latErrNames.begin(); errName != latErrNames.end();++errName ) std::cout<<*errName<<std::endl;




  for( MnvPlotter::ErrorSummaryGroupMap::iterator it = plotter->error_summary_group_map.begin(); it!=plotter->error_summary_group_map.end();++it )
  {
    std::string groupname = it->first;
    cout<<"Plotting Systematic Group: "<<groupname<<endl;
    for(errName = it->second.begin(); errName != it->second.end(); ++ errName ) std::cout<<*errName<<std::endl;

    plotter->DrawErrorSummary(h_dthetaR_signal[main_hist], "TL", false, true, 0.00001, area_norm, groupname );
    utils->Print( c, Form("ProtonSelection%s/h_dthetaR_signal_with_sys_errorsummary_frac_%s_%s_%s", "Signal", groupname.c_str(), "2", "pot" ) );
    plotter->DrawErrorSummary(h_dthetaR_blobs[main_hist], "TL", false, true, 0.00001, area_norm, groupname );
    utils->Print( c, Form("ProtonSelection%s/h_dthetaR_blobs_with_sys_errorsummary_frac_%s_%s_%s", "BlobSideBand", groupname.c_str(), "2", "pot" ) );
    plotter->DrawErrorSummary(h_dthetaR_michel[main_hist], "TL", false, true, 0.00001, area_norm, groupname );
    utils->Print( c, Form("ProtonSelection%s/h_dthetaR_michel_with_sys_errorsummary_frac_%s_%s_%s", "MichelSideBand", groupname.c_str(), "2", "pot" ) );
    plotter->DrawErrorSummary(h_dthetaR_micblob[main_hist], "TL", false, true, 0.00001, area_norm, groupname );
    utils->Print( c, Form("ProtonSelection%s/h_dthetaR_micblob_with_sys_errorsummary_frac_%s_%s_%s", "MicBlobSideBand", groupname.c_str(), "2", "pot" ) );

    plotter->DrawErrorSummary(h_dthetaP_signal[main_hist], "TL", false, true, 0.00001, area_norm, groupname );
    utils->Print( c, Form("ProtonSelection%s/h_dthetaP_signal_with_sys_errorsummary_frac_%s_%s_%s", "Signal", groupname.c_str(), "2", "pot" ) );
    plotter->DrawErrorSummary(h_dthetaP_blobs[main_hist], "TL", false, true, 0.00001, area_norm, groupname );
    utils->Print( c, Form("ProtonSelection%s/h_dthetaP_blobs_with_sys_errorsummary_frac_%s_%s_%s", "BlobSideBand", groupname.c_str(), "2", "pot" ) );
    plotter->DrawErrorSummary(h_dthetaP_michel[main_hist], "TL", false, true, 0.00001, area_norm, groupname );
    utils->Print( c, Form("ProtonSelection%s/h_dthetaP_michel_with_sys_errorsummary_frac_%s_%s_%s", "MichelSideBand", groupname.c_str(), "2", "pot" ) );
    plotter->DrawErrorSummary(h_dthetaP_micblob[main_hist], "TL", false, true, 0.00001, area_norm, groupname );
    utils->Print( c, Form("ProtonSelection%s/h_dthetaP_micblob_with_sys_errorsummary_frac_%s_%s_%s", "MicBlobSideBand", groupname.c_str(), "2", "pot" ) );



  }

  bool hasMuonEnergySideBand = h_dthetaR_signal[main_hist]->HasLatErrorBand("Muon_Energy");




  //------------------------------------------------------------------------------------------------
  // Drawing error summaries in the following groups "as absolute uncertainties" 
  // This was developed for testing the effects (or not) of the flux constraint implementations
  // You might need to tune the Y-axis to see all the absolute uncertainties  
  //------------------------------------------------------------------------------------------------
  //=====================================================================================
  // Plotting one variable in bins of the other (presently for pT_mu and pZ_mu) 
  // Also make plots for kData and kMC histos with the sys error band on the MC 
  //=====================================================================================
  //plottingOneInBinsOfAnother(h_pzmu_ptmu, "pZ", norm, true, plotter, utils, style, includeData, sample, multiplicity_label, pot_data, pot_mc, true);
  
  //plottingOneInBinsOfAnother(h_pzmu_ptmu, "pT", norm, true, plotter, utils, style, includeData, sample, multiplicity_label, pot_data, pot_mc, true); 

  //====================================
  //full projection split
  //====================================
  /*
  mnvh1d *recoil_proj[nhistos];
  mnvh1d *recoil_inc_proj[nhistos];
  for(int i=0;i<nhistos;i++){
    recoil_proj[i] =  h_recoil_ptmu[i]->projectionx(form("%s_projx",h_recoil_ptmu[i]->getname()),0,-1);
    recoil_inc_proj[i] =  h_recoil_inc_ptmu[i]->projectionx(form("%s_projx",h_recoil_inc_ptmu[i]->getname()),0,-1);
  }
  utils->drawstacked( recoil_proj, style, area_norm, pot_data, pot_mc, includedata );
  utils->print( c, form("neutronselection%s/recoil_projection_%s_%s", sample.data(), multiplicity_label.data(), norm.data()) );

  utils->drawstacked( recoil_inc_proj, style, area_norm, pot_data, pot_mc, includedata );
  utils->print( c, form("neutronselection%s/recoil_inc_projection_%s_%s", sample.data(), multiplicity_label.data(), norm.data()) );
  */
  //=====================================================================================
  // Plotting one variable in bins of the other (presently for pT_mu and pZ_mu) 
  // Make stacked plots of the MC split into its various styles 
  //=====================================================================================
  //plottingOneInBinsOfAnother(h_pzmu_ptmu, "pZ", norm, true, plotter, utils, style, includeData, sample, multiplicity_label, pot_data, pot_mc, false);

  //plottingOneInBinsOfAnother(h_pzmu_ptmu, "pT", norm, true, plotter, utils, style, includeData, sample, multiplicity_label, pot_data, pot_mc, false);

  //Normalize by binwidth
  f->Close();
  delete f;

  delete utils;
  delete plotter; 

  return 0;
};

void drawEnuLines(CCQENuPlotUtils *utils, int color, Option_t* opts){

  TF1 *fenu;
  fenu = utils->getEnuLine(2.0, 0.0, 30.0, kBlue+1);
  fenu->Draw(opts);
  fenu = utils->getEnuLine(4.0, 0.0, 30.0, kBlue+1);
  fenu->Draw(opts);
  fenu = utils->getEnuLine(5.99, 0.0, 30.0, kBlue+1);
  fenu->Draw(opts);
  fenu = utils->getEnuLine(8.0, 0.0, 30.0, kBlue+1);
  fenu->Draw(opts);
  fenu = utils->getEnuLine(10.0, 0.0, 30.0, kBlue+1);
  fenu->Draw(opts);
  fenu = utils->getEnuLine(15.0, 0.0, 30.0, kBlue+1);
  fenu->Draw(opts);
  fenu = utils->getEnuLine(20.0, 0.0, 30.0, kBlue+1);
  fenu->Draw(opts);
  fenu = utils->getEnuLine(25.0, 0.0, 30.0, kBlue+1);
  fenu->Draw(opts);
  fenu = utils->getEnuLine(30.0, 0.0, 30.0, kBlue+1);
  fenu->Draw(opts);
}

void drawQ2Lines(CCQENuPlotUtils *utils, int color, Option_t* opts){
  TF1 *fq2;
  fq2 = utils->getQ2Line(0.1, 0.0, 30.0, kBlack);
  fq2->Draw(opts);
  //fq2 = utils->getQ2Line(0.2, 0.0, 10.0, kBlack);
  //fq2->Draw(opts);
  fq2 = utils->getQ2Line(0.4, 0.0, 30.0, kBlack);
  fq2->Draw(opts);
  fq2 = utils->getQ2Line(0.8, 0.0, 30.0, kBlack);
  fq2->Draw(opts);
  fq2 = utils->getQ2Line(1.2, 0.0, 30.0, kBlack);
  fq2->Draw(opts);
  fq2 = utils->getQ2Line(2.0, 0.0, 30.0, kBlack);
  fq2->Draw(opts);
  fq2 = utils->getQ2Line(4.0, 0.0, 30.0, kBlack);
  fq2->Draw(opts);
  fq2 = utils->getQ2Line(6.0, 0.0, 30.0, kBlack);
  fq2->Draw(opts);
}

void plottingOneInBinsOfAnother( MnvH2D *de2DHisto[], TString plotVariable, TString norm, bool areaNormalizeEntireHisto, MnvPlotter *plotter, CCQENuPlotUtils *utils, string style, bool includeData, TString sample, TString multiplicityLabel, double pot_data, double pot_mc, bool draw_with_sys_band ) {  
  
  TCanvas* canv;

  bool area_norm = (norm.EqualTo("area")) ? true: false; 

  bool plotXInBinsOfY = false; 
  if( plotVariable == "pz" || plotVariable == "PZ" || plotVariable == "pZ" || plotVariable == "Pz" ) { 
    std::cout << "You will be plotting " << plotVariable << " in bins of PT" << std::endl; 
    plotXInBinsOfY = true; 
  } else if( plotVariable == "pt" || plotVariable == "PT" || plotVariable == "pT" || plotVariable == "Pt" ) {
    std::cout << "You will be plotting " << plotVariable << " in bins of PZ" << std::endl; 
  } else if( plotVariable == "dperp" ) {
    std::cout << "You will be plotting " << plotVariable << " in bins of Qsq" << std::endl; 
  } else {  
    std::cout << "Could not find the variable you want to plot " << std::endl; 
  }
  
  double areaNorm_scale = 0;  
  if( area_norm && areaNormalizeEntireHisto ) { 
    areaNorm_scale = de2DHisto[kData]->Integral() / de2DHisto[kMC]->Integral(); 
    de2DHisto[kMC]->Scale( areaNorm_scale ); 
    de2DHisto[kQELikeNot]->Scale( areaNorm_scale ); 
  }
  
  int lastBinNeeded = 0; 
  if ( plotXInBinsOfY ) { 
    lastBinNeeded = de2DHisto[kData]->GetNbinsY()+1;  
  } else { 
    lastBinNeeded = de2DHisto[kData]->GetNbinsX()+1;  
  }
  
  TString binVariable = plotXInBinsOfY ? "PT" : "PZ"; 

  //--------------------------------------------------------------
  // Taking projections of pZ for each bin of pT or vice versa 
  //--------------------------------------------------------------  
  for( int i=1; i<lastBinNeeded; ++i ) {
    
    canv  = new TCanvas("canv","Muon Kinematics in Bins"); 
    MnvH1D *h_projections[nHistos]; 
    
    cout << "Value of " << binVariable << " bin " << i << endl; 
    if( plotXInBinsOfY ) { 
      for( unsigned int j=kData; j<nHistos; ++j ) { 
	h_projections[j] = (MnvH1D*)de2DHisto[j]->ProjectionX(Form("h_%s_In_%sBin_%i_%s", plotVariable.Data(), binVariable.Data(), i, names[j].c_str()), i, i ); 
      }
      h_projections[kMC]->GetYaxis()->SetTitle(Form("Events / %.1f GeV", h_projections[kMC]->GetXaxis()->GetBinWidth(1))); 
    } else { 
      for( unsigned int j=kData; j<nHistos; ++j ) { 
	h_projections[j] = (MnvH1D*)de2DHisto[j]->ProjectionY(Form("h_%s_%sBin_%i_%s", plotVariable.Data(), binVariable.Data(), i, names[j].c_str()), i, i ); 
      }
      h_projections[kMC]->GetYaxis()->SetTitle(Form("Events / %.1f GeV", h_projections[kMC]->GetXaxis()->GetBinWidth(1)));  
    }
    
    //------------------------------------------------------------------
    // This is area normalizing individually inside each pT/pZ bin 
    // after taking the projections of MnvH2Ds. 
    // Don't know if this is correct or not, but it might be worth 
    // it for comparison purposes.    
    //------------------------------------------------------------------
    if( area_norm && !areaNormalizeEntireHisto ) { 
      areaNorm_scale = h_projections[kData]->Integral() / h_projections[kMC]->Integral(); 
      h_projections[kMC]->Scale( areaNorm_scale ); 
      h_projections[kQELikeNot]->Scale( areaNorm_scale ); 
    }
    
    double lowerBinEdge, upperBinEdge; 
    if( plotXInBinsOfY ) { 
      lowerBinEdge = de2DHisto[kData]->GetYaxis()->GetBinLowEdge(i); 
      upperBinEdge = de2DHisto[kData]->GetYaxis()->GetBinUpEdge(i); 
    } else { 
      lowerBinEdge = de2DHisto[kData]->GetXaxis()->GetBinLowEdge(i); 
      upperBinEdge = de2DHisto[kData]->GetXaxis()->GetBinUpEdge(i); 
    }
    
    //---------------------------------------------------------
    // Draw with systematics error band or the MC split up 
    //---------------------------------------------------------
    if( draw_with_sys_band ) { 
      cout << "Drawing Data-MC with systematic band " << endl; 
      plotter->DrawDataMCWithErrorBand( h_projections[kData], h_projections[kMC], 1.0, "TR", false, (sample.EqualTo("Signal"))? h_projections[kQELikeNot]: NULL, NULL, area_norm ); 
      plotter->AddPlotLabel( Form("%.1f #leq %s < %.1f", lowerBinEdge, binVariable.Data(), upperBinEdge), 0.5, 0.93 );
      //(area_norm) ? plotter->WriteNorm("Area Normalized", 0.3, 0.88, 0.045, pot_data) : plotter->AddPOTNormBox(pot_data,pot_mc, 0.3, 0.88); 
      utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
      utils->Print( canv, Form("ProtonSelection%s/muon_%s_In_%sBin_%i_with_sys_%s_%s", sample.Data(), plotVariable.Data(), binVariable.Data(), i, multiplicityLabel.Data(), norm.Data()) ); 
    } else {
      cout << "Drawing Data-MC with the MC split up " << endl; 
      utils->drawStacked( h_projections, style, area_norm, pot_data, pot_mc, includeData );
      plotter->AddPlotLabel( Form("%.1f #leq %s < %.1f", lowerBinEdge, binVariable.Data(), upperBinEdge), 0.5, 0.93 );
      utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
      utils->Print( canv, Form("ProtonSelection%s/muon_%s_In_%sBin_%i_splitMC_%s_%s", sample.Data(), plotVariable.Data(), binVariable.Data(), i, multiplicityLabel.Data(), norm.Data()) );  
    }
    
    //Do you want to plot ratio too (e.g. AcceptanceStudiesPlots) ? 

    if(canv) canv->Close();  
  } //End of for loop over bins 
}



int main( int argc, char *argv[])
{
  ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);
  if (argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
        "\t-./ProtonSelectionPlots path/to/filename Sample Area_Normalize\n\n"<<
        "\t-path/to/filename\t =\t Path to the input root file\n"<<
        "\t-path/to/filename\t =\t Path to the no background sub input root file\n"<<
        "\t-Area_Normalize\t\t =\t Area normalize Sample\n"<<
        "\t-MC_Plot_Style\t\t =\t QE, QE_PionInFS, QELike, QELike_split, QELike_PionInFS, QELike_split_PionInFS"<<std::endl;
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    return 0; 
  }

  //! Default parameters
  std::vector<std::string> par;
  par.push_back("ProtonSelectionPlots");
  par.push_back( Form("%s/ana/rootfiles/ProtonSelectionHists.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/ProtonSelectionHists_nobgSub.root",getenv("CCQENUROOT") ) );
  par.push_back("QE");

  //! Set user parameters
  for( int i=0; i<argc; ++i){
    std::cout<<"Input parameter: "<< argv[i]<<std::endl;
    par.at(i) = argv[i];
  }

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;
  
  return ProtonSelectionPlots(par[1], par[2],  par[3]);
}
