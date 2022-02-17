#include "include/CCQENuPlotUtils.h"
#include "TParameter.h"

using namespace CCQENU_ANA;

void drawQ2Lines(CCQENuPlotUtils *utils, int color = kBlack, Option_t* opts = "");
void drawEnuLines(CCQENuPlotUtils *utils, int color = kBlue, Option_t* opts = "");

void plottingOneInBinsOfAnother( MnvH2D *de2DHisto[], TString plotVariable, TString norm, bool areaNormalizeEntireHisto, MnvPlotter *plotter, CCQENuPlotUtils *utils, string style, bool includeData, TString sample, TString multiplicityLabel, double pot_data, double pot_mc, bool draw_with_sys_band ); 


int ProtonSelectionPlots( string filename, bool area_norm = false, string style = "QE")
{
  // read file to get histograms
  TFile *f = new TFile( filename.c_str(), "READ" );
  if (f->IsZombie() || f->GetListOfKeys()->IsEmpty()){
    Error("VertexEnergyPlots","Could not get histogram ROOT file or it was empty.");
    return 1;
  }

  //-------------------------------------------
  // Get event counts before cuts for norm
  //-------------------------------------------
  TVector2 *evt = (TVector2*)f->Get("n_events");
  double data_events = evt->X();
  double mc_events = evt->Y();
  
  cout<< " Number of Data Events = " << data_events << endl;
  cout<< " Number of MC Events = " <<   mc_events << endl;

  //-----------------------------------------------------------
  // Will the input ROOT file be used for flux constraints  
  //-----------------------------------------------------------
  bool fluxHistoExists = false; 
  TParameter<bool> *fluxParam = (TParameter<bool>*)f->Get("appliedFluxConstraint"); 
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

  //--------------------------------------------
  // Load histos to plot 
  //--------------------------------------------
  //---------------------------------------------
  // Neutron Angular Histos
  //---------------------------------------------
  MnvH2D *h_neutron_dperp_qsq[nHistos], *h_neutron_phit_qsq[nHistos], *h_neutron_dreactplane_qsq[nHistos], *h_neutron_dperp_dreactplane[nHistos], *h_neutron_dreactplane_vtxE[nHistos];
  MnvH2D *h_neutron_dreactplane_protonKE[nHistos], *h_neutron_dreactplane_protonAngle[nHistos], *h_neutron_dreactplane_expectedAngle[nHistos];
  MnvH1D *h_neutron_dperp[nHistos], *h_neutron_phit[nHistos], *h_neutron_dreactplane[nHistos];
  MnvH1D* h_neutron_dreact_dperpcut[nHistos], *h_neutron_vtxE[nHistos];
  MnvH1D* h_neutron_protonKE[nHistos], *h_neutron_protonAngle[nHistos], *h_neutron_expectedAngle[nHistos];

  utils->bookHistos( f, h_neutron_dreactplane_vtxE, "h_neutron_dreactplane_vtxE");
  utils->bookHistos( f, h_neutron_dperp_qsq, "h_neutron_dperp_qsq");
  utils->bookHistos( f, h_neutron_dperp_dreactplane, "h_neutron_dperp_dreactplane");
  utils->bookHistos( f, h_neutron_dreactplane_qsq, "h_neutron_dreactplane_qsq");
  utils->bookHistos( f, h_neutron_phit_qsq, "h_neutron_phit_qsq");
  utils->bookHistos( f, h_neutron_dreactplane_protonAngle, "h_neutron_dreactplane_protonAngle");
  utils->bookHistos( f, h_neutron_dreactplane_protonKE, "h_neutron_dreactplane_protonKE");
  utils->bookHistos( f, h_neutron_dreactplane_expectedAngle, "h_neutron_dreactplane_expectedAngle");
  for(int i=0;i<nHistos;i++){
    h_neutron_dperp[i] =  h_neutron_dperp_qsq[i]->ProjectionX(Form("%s_projX",h_neutron_dperp_qsq[i]->GetName()),0,-1);
    h_neutron_dreactplane[i] =  h_neutron_dreactplane_qsq[i]->ProjectionX(Form("%s_projX",h_neutron_dreactplane_qsq[i]->GetName()),0,-1);
    h_neutron_phit[i] =  h_neutron_phit_qsq[i]->ProjectionX(Form("%s_projX",h_neutron_phit_qsq[i]->GetName()),0,-1);
    h_neutron_vtxE[i] =  h_neutron_dreactplane_vtxE[i]->ProjectionY( Form("%s_projY",h_neutron_dreactplane_vtxE[i]->GetName()),0,-1);
    h_neutron_protonKE[i] = h_neutron_dreactplane_protonKE[i]->ProjectionY( Form("%s_projY",h_neutron_dreactplane_protonKE[i]->GetName()),0,-1);
    h_neutron_protonAngle[i] = h_neutron_dreactplane_protonAngle[i]->ProjectionY( Form("%s_projY",h_neutron_dreactplane_protonAngle[i]->GetName()),0,-1);
    h_neutron_expectedAngle[i] = h_neutron_dreactplane_expectedAngle[i]->ProjectionY( Form("%s_projY",h_neutron_dreactplane_expectedAngle[i]->GetName()),0,-1);
    int bmin = h_neutron_dperp[i]->FindBin(-.02);
    int bmax = h_neutron_dperp[i]->FindBin(.02);
    h_neutron_dreact_dperpcut[i] =  h_neutron_dperp_dreactplane[i]->ProjectionX(Form("%s_projX",h_neutron_dreactplane_qsq[i]->GetName()),bmin,bmax);
    cout<< h_neutron_dreactplane_qsq[i]->GetName()<<endl;
    cout<< h_neutron_dreact_dperpcut[i]->Integral();
  }

  //--------------------------------------------
  // Plot Muon Kinematics...
  //--------------------------------------------
  cout << "Plotting Neutron Angulars..." << endl;
  c = new TCanvas("cNeutA","Neutron Angulars"); 
  bool includeData = true;

  TMap *sampleMap = (TMap*)f->Get("sample");
  TObjString *multiplicity_obj = (TObjString*)sampleMap->GetValue("multiplicity");
  TString multiplicity_label = multiplicity_obj->GetString();

  TObjString *sample_obj = (TObjString*)sampleMap->GetValue("sample");
  TString sample = sample_obj->GetString();

  cout << "-------------------------------------------------------"<<endl;
  cout << "\tMultiplicity = " << multiplicity_label.Data() << endl;
  cout << "\tSample = "       << sample.Data() << endl;
  cout << "\tMC plot style = "       << style.c_str() << endl;
  cout << "-------------------------------------------------------"<<endl;
  
  TString toplabel;
  if( sample.EqualTo("Signal") )
    toplabel = "#nu_{#mu} Tracker #rightarrow CCQE-like";
  else if( sample.EqualTo("SideBand") )
    toplabel = "#nu_{#mu} Tracker #bullet Recoil Sideband";
  else if( sample.EqualTo("MichelSideBand") )
    toplabel = "#nu_{#mu} Tracker #bullet Michel Sample";

  //Settings   
  TGaxis::SetMaxDigits(2);
  //2-D
  plotter->SetRedHeatPalette();

  //-----------------------------------------------------------------------------
  //Make plots of the RAW event counts before POT and binwidth normalizations 
  //-----------------------------------------------------------------------------
  MnvH2D *h_neutron_dperp_qsq_mc_noBWNorm = (MnvH2D*) h_neutron_dperp_qsq[kMC]->Clone("h_neutron_dperp_qsq_mc_noBWNorm"); 
  h_neutron_dperp_qsq_mc_noBWNorm->GetZaxis()->SetTitle("Events");
  c->SetLogz(1);
  h_neutron_dperp_qsq_mc_noBWNorm->Draw("colz");
  plotter->AddHistoTitle("MC");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dperp_qsq_mc_rawEvts_colz_logz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  

  MnvH2D *h_neutron_dperp_qsq_mc_noBWNorm_qeh_ratio = (MnvH2D*) h_neutron_dperp_qsq[kQE_H]->Clone("h_neutron_dperp_qsq_mc_noBWNorm_qeh_ratio"); 
  MnvH2D *h_neutron_dperp_qsq_mc_noBWNorm_qeh = (MnvH2D*) h_neutron_dperp_qsq[kQE_H]->Clone("h_neutron_dperp_qsq_mc_noBWNorm_qeh"); 
  c->SetLogz(0);
  h_neutron_dperp_qsq_mc_noBWNorm_qeh_ratio->GetZaxis()->SetTitle("Ratio");
  h_neutron_dperp_qsq_mc_noBWNorm_qeh_ratio->Divide( h_neutron_dperp_qsq_mc_noBWNorm_qeh, h_neutron_dperp_qsq_mc_noBWNorm );
  h_neutron_dperp_qsq_mc_noBWNorm_qeh_ratio->Draw("colz");
  plotter->AddHistoTitle("MC");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dperp_qsq_mc_rawEvts_qeh_ratio_colz_logz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
 


  MnvH2D *h_neutron_dperp_qsq_data_noBWNorm = (MnvH2D*) h_neutron_dperp_qsq[kData]->Clone("h_neutron_dperp_qsq_data_noBWNorm"); 
  h_neutron_dperp_qsq_data_noBWNorm->GetZaxis()->SetTitle("Events");
  c->SetLogz(1);
  h_neutron_dperp_qsq_data_noBWNorm->Draw("colz");
  plotter->AddHistoTitle("DATA");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dperp_qsq_data_rawEvts_colz_logz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );

  MnvH2D *h_neutron_dperp_qsq_data_fracStatError = (MnvH2D*) h_neutron_dperp_qsq[kData]->Clone("h_neutron_dperp_qsq_data_fracStatError");
  h_neutron_dperp_qsq_data_fracStatError->Clear(); 
  int dperp_maxbin = h_neutron_dperp_qsq[kData]->GetNbinsX();
  int qsq_maxbin = h_neutron_dperp_qsq[kData]->GetNbinsY();
  for( int perp=1; perp<=dperp_maxbin; perp++ ) { 
    for( int qsq=1; qsq<=qsq_maxbin; qsq++ ) { 
      double fracStatError = 1.0; 
      if( h_neutron_dperp_qsq[kData]->GetBinContent(perp, qsq) > 0.0 ) fracStatError = 1.0/ sqrt( h_neutron_dperp_qsq[kData]->GetBinError(perp,qsq) ); 
      //cout << "Bin Content " << h_pzmu_ptmu[kData]->GetBinContent(pZ, pT) << " Bin Error " << h_pzmu_ptmu[kData]->GetBinError(pZ, pT) << endl;  
      h_neutron_dperp_qsq_data_fracStatError->SetBinContent(perp,qsq, fracStatError);
      h_neutron_dperp_qsq_data_fracStatError->SetBinError(perp,qsq, 1.0); 
    }
  }
  c->SetLogz(0); 
  h_neutron_dperp_qsq_data_fracStatError->Draw("colz");
  plotter->AddHistoTitle("DATA");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dperp_qsq_data_fracStatError_colz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );

  //======================

  MnvH2D *h_neutron_dperp_dreactplane_mc_noBWNorm = (MnvH2D*) h_neutron_dperp_dreactplane[kMC]->Clone("h_neutron_dperp_dreactplane_mc_noBWNorm"); 
  h_neutron_dperp_dreactplane_mc_noBWNorm->GetZaxis()->SetTitle("Events");
  c->SetLogz(1);
  h_neutron_dperp_dreactplane_mc_noBWNorm->Draw("colz");
  plotter->AddHistoTitle("MC");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dperp_dreactplane_mc_rawEvts_colz_logz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  
  MnvH2D *h_neutron_dperp_dreactplane_data_noBWNorm = (MnvH2D*) h_neutron_dperp_dreactplane[kData]->Clone("h_neutron_dperp_dreactplane_data_noBWNorm"); 
  h_neutron_dperp_dreactplane_data_noBWNorm->GetZaxis()->SetTitle("Events");
  c->SetLogz(1);
  h_neutron_dperp_dreactplane_data_noBWNorm->Draw("colz");
  plotter->AddHistoTitle("DATA");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dperp_dreactplane_data_rawEvts_colz_logz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );

  MnvH2D *h_neutron_dperp_dreactplane_data_fracStatError = (MnvH2D*) h_neutron_dperp_dreactplane[kData]->Clone("h_neutron_dperp_dreactplane_data_fracStatError");
  h_neutron_dperp_dreactplane_data_fracStatError->Clear(); 
  int dreactplane_maxbin = h_neutron_dperp_dreactplane[kData]->GetNbinsY();
  for( int perp=1; perp<=dperp_maxbin; perp++ ) { 
    for( int dreactplane=1; dreactplane<=dreactplane_maxbin; dreactplane++ ) { 
      double fracStatError = 1.0; 
      if( h_neutron_dperp_dreactplane[kData]->GetBinContent(perp, dreactplane) > 0.0 ) fracStatError = 1.0/ sqrt( h_neutron_dperp_dreactplane[kData]->GetBinError(perp,dreactplane) ); 
      //cout << "Bin Content " << h_pzmu_ptmu[kData]->GetBinContent(pZ, pT) << " Bin Error " << h_pzmu_ptmu[kData]->GetBinError(pZ, pT) << endl;  
      h_neutron_dperp_dreactplane_data_fracStatError->SetBinContent(perp,dreactplane, fracStatError);
      h_neutron_dperp_dreactplane_data_fracStatError->SetBinError(perp,dreactplane, 1.0); 
    }
  }
  c->SetLogz(0); 
  h_neutron_dperp_dreactplane_data_fracStatError->Draw("colz");
  plotter->AddHistoTitle("DATA");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dperp_dreactplane_data_fracStatError_colz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );

  //===================
  MnvH2D *h_neutron_dreactplane_qsq_mc_noBWNorm = (MnvH2D*) h_neutron_dreactplane_qsq[kMC]->Clone("h_neutron_dreactplane_qsq_mc_noBWNorm"); 
  h_neutron_dreactplane_qsq_mc_noBWNorm->GetZaxis()->SetTitle("Events");
  c->SetLogz(1);
  h_neutron_dreactplane_qsq_mc_noBWNorm->Draw("colz");
  plotter->AddHistoTitle("MC");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dreactplane_qsq_mc_rawEvts_colz_logz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  
  MnvH2D *h_neutron_dreactplane_qsq_data_noBWNorm = (MnvH2D*) h_neutron_dreactplane_qsq[kData]->Clone("h_neutron_dreactplane_qsq_data_noBWNorm"); 
  h_neutron_dreactplane_qsq_data_noBWNorm->GetZaxis()->SetTitle("Events");
  c->SetLogz(1);
  h_neutron_dreactplane_qsq_data_noBWNorm->Draw("colz");
  plotter->AddHistoTitle("DATA");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dreactplane_qsq_data_rawEvts_colz_logz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );

  MnvH2D *h_neutron_dreactplane_qsq_data_fracStatError = (MnvH2D*) h_neutron_dreactplane_qsq[kData]->Clone("h_neutron_dreactplane_qsq_data_fracStatError");
  h_neutron_dreactplane_qsq_data_fracStatError->Clear(); 
  for( int perp=1; perp<=dreactplane_maxbin; perp++ ) { 
    for( int qsq=1; qsq<=qsq_maxbin; qsq++ ) { 
      double fracStatError = 1.0; 
      if( h_neutron_dreactplane_qsq[kData]->GetBinContent(perp, qsq) > 0.0 ) fracStatError = 1.0/ sqrt( h_neutron_dreactplane_qsq[kData]->GetBinError(perp,qsq) ); 
      //cout << "Bin Content " << h_pzmu_ptmu[kData]->GetBinContent(pZ, pT) << " Bin Error " << h_pzmu_ptmu[kData]->GetBinError(pZ, pT) << endl;  
      h_neutron_dreactplane_qsq_data_fracStatError->SetBinContent(perp,qsq, fracStatError);
      h_neutron_dreactplane_qsq_data_fracStatError->SetBinError(perp,qsq, 1.0); 
    }
  }
  c->SetLogz(0); 
  h_neutron_dreactplane_qsq_data_fracStatError->Draw("colz");
  plotter->AddHistoTitle("DATA");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dreactplane_qsq_data_fracStatError_colz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );

  //===================
  MnvH2D *h_neutron_dreactplane_protonKE_mc_noBWNorm = (MnvH2D*) h_neutron_dreactplane_protonKE[kMC]->Clone("h_neutron_dreactplane_protonKE_mc_noBWNorm"); 
  h_neutron_dreactplane_qsq_mc_noBWNorm->GetZaxis()->SetTitle("Events");
  c->SetLogz(1);
  h_neutron_dreactplane_protonKE_mc_noBWNorm->Draw("colz");
  plotter->AddHistoTitle("MC");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dreactplane_protonKE_mc_rawEvts_colz_logz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  
  MnvH2D *h_neutron_dreactplane_protonKE_data_noBWNorm = (MnvH2D*) h_neutron_dreactplane_protonKE[kData]->Clone("h_neutron_dreactplane_protonKE_data_noBWNorm"); 
  h_neutron_dreactplane_protonKE_data_noBWNorm->GetZaxis()->SetTitle("Events");
  c->SetLogz(1);
  h_neutron_dreactplane_protonKE_data_noBWNorm->Draw("colz");
  plotter->AddHistoTitle("DATA");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dreactplane_protonKE_data_rawEvts_colz_logz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );

  MnvH2D *h_neutron_dreactplane_protonKE_data_fracStatError = (MnvH2D*) h_neutron_dreactplane_protonKE[kData]->Clone("h_neutron_dreactplane_protonKE_data_fracStatError");
  h_neutron_dreactplane_protonKE_data_fracStatError->Clear(); 
  int protonKE_maxbin = h_neutron_dreactplane_protonKE[kData]->GetNbinsY();
  for( int perp=1; perp<=dreactplane_maxbin; perp++ ) { 
    for( int protonKE=1; protonKE<=protonKE_maxbin; protonKE++ ) { 
      double fracStatError = 1.0; 
      if( h_neutron_dreactplane_protonKE[kData]->GetBinContent(perp, protonKE) > 0.0 ) fracStatError = 1.0/ sqrt( h_neutron_dreactplane_protonKE[kData]->GetBinError(perp,protonKE) ); 
      //cout << "Bin Content " << h_pzmu_ptmu[kData]->GetBinContent(pZ, pT) << " Bin Error " << h_pzmu_ptmu[kData]->GetBinError(pZ, pT) << endl;  
      h_neutron_dreactplane_protonKE_data_fracStatError->SetBinContent(perp,protonKE, fracStatError);
      h_neutron_dreactplane_protonKE_data_fracStatError->SetBinError(perp,protonKE, 1.0); 
    }
  }
  c->SetLogz(0); 
  h_neutron_dreactplane_protonKE_data_fracStatError->Draw("colz");
  plotter->AddHistoTitle("DATA");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dreactplane_protonKE_data_fracStatError_colz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );

  //===================
  MnvH2D *h_neutron_dreactplane_protonAngle_mc_noBWNorm = (MnvH2D*) h_neutron_dreactplane_protonAngle[kMC]->Clone("h_neutron_dreactplane_protonAngle_mc_noBWNorm"); 
  h_neutron_dreactplane_qsq_mc_noBWNorm->GetZaxis()->SetTitle("Events");
  c->SetLogz(1);
  h_neutron_dreactplane_protonAngle_mc_noBWNorm->Draw("colz");
  plotter->AddHistoTitle("MC");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dreactplane_protonAngle_mc_rawEvts_colz_logz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  
  MnvH2D *h_neutron_dreactplane_protonAngle_data_noBWNorm = (MnvH2D*) h_neutron_dreactplane_protonAngle[kData]->Clone("h_neutron_dreactplane_protonAngle_data_noBWNorm"); 
  h_neutron_dreactplane_protonAngle_data_noBWNorm->GetZaxis()->SetTitle("Events");
  c->SetLogz(1);
  h_neutron_dreactplane_protonAngle_data_noBWNorm->Draw("colz");
  plotter->AddHistoTitle("DATA");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dreactplane_protonAngle_data_rawEvts_colz_logz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );

  MnvH2D *h_neutron_dreactplane_protonAngle_data_fracStatError = (MnvH2D*) h_neutron_dreactplane_protonAngle[kData]->Clone("h_neutron_dreactplane_protonAngle_data_fracStatError");
  h_neutron_dreactplane_protonAngle_data_fracStatError->Clear(); 
  int protonAngle_maxbin = h_neutron_dreactplane_protonAngle[kData]->GetNbinsY();
  for( int perp=1; perp<=dreactplane_maxbin; perp++ ) { 
    for( int protonAngle=1; protonAngle<=protonAngle_maxbin; protonAngle++ ) { 
      double fracStatError = 1.0; 
      if( h_neutron_dreactplane_protonAngle[kData]->GetBinContent(perp, protonAngle) > 0.0 ) fracStatError = 1.0/ sqrt( h_neutron_dreactplane_protonAngle[kData]->GetBinError(perp,protonAngle) ); 
      //cout << "Bin Content " << h_pzmu_ptmu[kData]->GetBinContent(pZ, pT) << " Bin Error " << h_pzmu_ptmu[kData]->GetBinError(pZ, pT) << endl;  
      h_neutron_dreactplane_protonAngle_data_fracStatError->SetBinContent(perp,protonAngle, fracStatError);
      h_neutron_dreactplane_protonAngle_data_fracStatError->SetBinError(perp,protonAngle, 1.0); 
    }
  }
  c->SetLogz(0); 
  h_neutron_dreactplane_protonAngle_data_fracStatError->Draw("colz");
  plotter->AddHistoTitle("DATA");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dreactplane_protonAngle_data_fracStatError_colz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );

  //===================
  MnvH2D *h_neutron_dreactplane_expectedAngle_mc_noBWNorm = (MnvH2D*) h_neutron_dreactplane_expectedAngle[kMC]->Clone("h_neutron_dreactplane_expectedAngle_mc_noBWNorm"); 
  h_neutron_dreactplane_qsq_mc_noBWNorm->GetZaxis()->SetTitle("Events");
  c->SetLogz(1);
  h_neutron_dreactplane_expectedAngle_mc_noBWNorm->Draw("colz");
  plotter->AddHistoTitle("MC");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dreactplane_expectedAngle_mc_rawEvts_colz_logz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  
  MnvH2D *h_neutron_dreactplane_expectedAngle_data_noBWNorm = (MnvH2D*) h_neutron_dreactplane_expectedAngle[kData]->Clone("h_neutron_dreactplane_expectedAngle_data_noBWNorm"); 
  h_neutron_dreactplane_expectedAngle_data_noBWNorm->GetZaxis()->SetTitle("Events");
  c->SetLogz(1);
  h_neutron_dreactplane_expectedAngle_data_noBWNorm->Draw("colz");
  plotter->AddHistoTitle("DATA");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dreactplane_expectedAngle_data_rawEvts_colz_logz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );

  MnvH2D *h_neutron_dreactplane_expectedAngle_data_fracStatError = (MnvH2D*) h_neutron_dreactplane_expectedAngle[kData]->Clone("h_neutron_dreactplane_expectedAngle_data_fracStatError");
  h_neutron_dreactplane_expectedAngle_data_fracStatError->Clear(); 
  int expectedAngle_maxbin = h_neutron_dreactplane_expectedAngle[kData]->GetNbinsY();
  for( int perp=1; perp<=dreactplane_maxbin; perp++ ) { 
    for( int expectedAngle=1; expectedAngle<=expectedAngle_maxbin; expectedAngle++ ) { 
      double fracStatError = 1.0; 
      if( h_neutron_dreactplane_expectedAngle[kData]->GetBinContent(perp, expectedAngle) > 0.0 ) fracStatError = 1.0/ sqrt( h_neutron_dreactplane_expectedAngle[kData]->GetBinError(perp,expectedAngle) ); 
      //cout << "Bin Content " << h_pzmu_ptmu[kData]->GetBinContent(pZ, pT) << " Bin Error " << h_pzmu_ptmu[kData]->GetBinError(pZ, pT) << endl;  
      h_neutron_dreactplane_expectedAngle_data_fracStatError->SetBinContent(perp,expectedAngle, fracStatError);
      h_neutron_dreactplane_expectedAngle_data_fracStatError->SetBinError(perp,expectedAngle, 1.0); 
    }
  }
  c->SetLogz(0); 
  h_neutron_dreactplane_expectedAngle_data_fracStatError->Draw("colz");
  plotter->AddHistoTitle("DATA");
  plotter->AddPlotLabel( toplabel.Data(), 0.285, 0.93, 0.033 );
  plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.93, 0.033, kBlue, 32, 11 );
  utils->Print(c, Form("ProtonSelection%s/muon_neutron_dreactplane_expectedAngle_data_fracStatError_colz_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );



 
 
 
   //--------------------------------------------
   // Get Area or POT Normalization
   //--------------------------------------------
   double pot_data = utils->getPOTData( f );
   double pot_mc = utils->getPOTMC( f );
   double pot_scale = utils->getPOTNormFactor( f );
   if(pot_scale==0) pot_scale = 1.0;
   cout<< "POT DATA: " << pot_data << " POT MC: " << pot_mc << " POT Scale Factor = " << pot_scale << endl; 
 
   //1-D
   //double tmu_scale                          = (area_norm)? utils->getAreaNormFactor(h_tmu)   : pot_scale;
   //utils->scaleMCHistos( h_tmu, tmu_scale);
   double dperp_scale                          = (area_norm)? utils->getAreaNormFactor(h_neutron_dperp ): pot_scale;
   utils->scaleMCHistos( h_neutron_dperp, dperp_scale);
   double phit_scale                           = (area_norm)? utils->getAreaNormFactor(h_neutron_phit ): pot_scale;
   utils->scaleMCHistos( h_neutron_phit, phit_scale);
   double dreactplane_scale                    = (area_norm)? utils->getAreaNormFactor(h_neutron_dreactplane ): pot_scale;
   utils->scaleMCHistos( h_neutron_dreactplane, dreactplane_scale);
   utils->scaleMCHistos( h_neutron_dreact_dperpcut, dreactplane_scale);

   double vtxE_scale                          =(area_norm)? utils->getAreaNormFactor(h_neutron_vtxE ): pot_scale;
   utils->scaleMCHistos( h_neutron_vtxE, vtxE_scale);

   utils->scaleMCHistos( h_neutron_protonKE, vtxE_scale );
   utils->scaleMCHistos( h_neutron_protonAngle, vtxE_scale );
   utils->scaleMCHistos( h_neutron_expectedAngle, vtxE_scale );
   cout<<"dreact scale:"<<dreactplane_scale<<endl;
 
   //2-D
   //double tmu_theta_scale       = (area_norm)? utils->getAreaNormFactor(h_tmu_theta)      : pot_scale;
   //double tmu_costheta_scale    = (area_norm)? utils->getAreaNormFactor(h_tmu_costheta)   : pot_scale;
   //double pzmu_ptmu_scale       = (area_norm)? utils->getAreaNormFactor(h_pzmu_ptmu)      : pot_scale;
   //double pmu_ptmu_scale        = (area_norm)? utils->getAreaNormFactor(h_pmu_ptmu)       : pot_scale;
   //utils->scaleMCHistos( h_tmu_theta, tmu_theta_scale);
   //utils->scaleMCHistos( h_tmu_costheta, tmu_costheta_scale);
   //utils->scaleMCHistos( h_pzmu_ptmu, pzmu_ptmu_scale);
   //utils->scaleMCHistos( h_pmu_ptmu, pmu_ptmu_scale);
   double dperp_qsq_scale                          = (area_norm)? utils->getAreaNormFactor(h_neutron_dperp_qsq ): pot_scale;
   //utils->scaleMCHistos( h_neutron_dperp, dperp_scale);
   double phit_qsq_scale                           = (area_norm)? utils->getAreaNormFactor(h_neutron_phit_qsq ): pot_scale;
   //utils->scaleMCHistos( h_neutron_phit_qsq, phit_qsq_scale);
   double dreactplane_qsq_scale                    = (area_norm)? utils->getAreaNormFactor(h_neutron_dreactplane_qsq ): pot_scale;
   //utils->scaleMCHistos( h_neutron_dreactplane_qsq, dreactplane_qsq_scale);


  //Exlude last theta bin (25 to 40 degrees bin)
  /*
  h_theta[kMC]->GetXaxis()->SetRange( 1, theta_maxbin );
  h_tmu_theta[kMC]->GetYaxis()->SetRange( 1, theta_maxbin );
  h_theta[kData]->GetXaxis()->SetRange( 1, theta_maxbin );
  h_tmu_theta[kData]->GetYaxis()->SetRange( 1, theta_maxbin );
  */

  //Exclude last bin in pz/pt 
  /*
  h_pzmu[kMC]->GetXaxis()->SetRange(1, pz_maxbin );
  h_pzmu[kData]->GetXaxis()->SetRange(1, pz_maxbin );
  h_ptmu[kMC]->GetXaxis()->SetRange(1, pt_maxbin );
  h_ptmu[kData]->GetXaxis()->SetRange(1, pt_maxbin );
  h_pzmu_ptmu[kMC]->GetXaxis()->SetRange(1,pz_maxbin );
  h_pzmu_ptmu[kMC]->GetYaxis()->SetRange(1,pt_maxbin );
  h_pzmu_ptmu[kData]->GetXaxis()->SetRange(1,pz_maxbin );
  h_pzmu_ptmu[kData]->GetYaxis()->SetRange(1,pt_maxbin );
  */


  /* Temporary settings 
  h_pzmu[kMC]->GetYaxis()->SetTitle("Events / 1.0 GeV");
  h_pzmu[kData]->GetYaxis()->SetTitle("Events / 1.0 GeV");
  h_ptmu[kMC]->GetYaxis()->SetTitle("Events / 0.125 GeV");
  h_ptmu[kData]->GetYaxis()->SetTitle("Events / 0.125 GeV");
  */

  //==============================================================  
  //1-D Data with MC broken down into it's different templates
  //==============================================================
  utils->drawStacked( h_neutron_dperp, style, area_norm, pot_data, pot_mc, includeData );
  utils->Print( c, Form("ProtonSelection%s/neutron_dperp_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  utils->drawStacked( h_neutron_phit, style, area_norm, pot_data, pot_mc, includeData );
  utils->Print( c, Form("ProtonSelection%s/neutron_phit_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  utils->drawStacked( h_neutron_dreactplane, style, area_norm, pot_data, pot_mc, includeData );
  utils->Print( c, Form("ProtonSelection%s/neutron_dreactplane_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  utils->drawStacked( h_neutron_vtxE, style, area_norm, pot_data, pot_mc, includeData );
  utils->Print( c, Form("ProtonSelection%s/neutron_vtxE_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  utils->drawStacked( h_neutron_dreact_dperpcut, style, area_norm, pot_data, pot_mc, includeData );
  utils->Print( c, Form("ProtonSelection%s/neutron_dreact_dperpcut_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );

  utils->drawStacked( h_neutron_protonKE, style, area_norm, pot_data, pot_mc, includeData );
  utils->Print( c, Form("ProtonSelection%s/neutron_protonKE_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  utils->drawStacked( h_neutron_protonAngle, style, area_norm, pot_data, pot_mc, includeData,-1, "Angle", "#", 0,180 );
  utils->Print( c, Form("ProtonSelection%s/neutron_protonAngle_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  utils->drawStacked( h_neutron_expectedAngle, style, area_norm, pot_data, pot_mc, includeData,-1, "Angle", "#", 0,180 );
  utils->Print( c, Form("ProtonSelection%s/neutron_expectedAngle_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
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

  plotter->DrawDataMCRatio(h_neutron_dperp[kData], h_neutron_dperp[kMC], 1.0, true, true, -1.5, 1.5, "Data / MC", area_norm);  
  utils->Print( c, Form("ProtonSelection%s/neutron_dperp_mc_ratio_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) ); 
  plotter->DrawDataMCRatio(h_neutron_phit[kData], h_neutron_phit[kMC], 1.0, true, true, -1.5, 1.5, "Data / MC", area_norm);  
  utils->Print( c, Form("ProtonSelection%s/neutron_phit_mc_ratio_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) ); 
  plotter->DrawDataMCRatio(h_neutron_dreactplane[kData], h_neutron_dreactplane[kMC], 1.0, true, true, -1.5, 1.5, "Data / MC", area_norm);  
  utils->Print( c, Form("ProtonSelection%s/neutron_dreactplane_mc_ratio_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) ); 



  //==================================================
  // 1-D of Data with MC with stats + systematics
  // Plots of T_mu, Theta_mu, pT_mu and pZ_mu 
  //==================================================
  plotter->DrawDataMCWithErrorBand(h_neutron_dperp[kData], h_neutron_dperp[kMC], 1.0, "TR", false, (sample.EqualTo("Signal"))? h_neutron_dperp[kQELikeNot]: NULL, NULL, area_norm);
  //(area_norm)? plotter->WriteNorm("Area Normalized", 0.3, 0.88, 0.045, pot_data) : plotter->AddPOTNormBox( pot_data, pot_mc, 0.3, 0.88);
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/neutron_dperp_with_sys_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );

  //======================================================
  // Plot Fractional Systematic Error Summary 
  // Plots of pT_mu and pZ_mu (errors on MC right now) 
  //======================================================
  plotter->error_summary_group_map = utils->getSystematicGroupMap();
  plotter->error_color_map         = utils->getSystematicGroupMapColors();
  TGaxis::SetMaxDigits(5);
/*  
  //--------------------------------------------
  // DATA-MC Error Summary for Muon pT and pZ 
  //--------------------------------------------
  plotter->DrawDataMCErrorSummary(h_pzmu_mc_temp, h_pzmu_data_temp, "TR", false, true, 0.00001, area_norm );
  plotter->AddPlotLabel( Form("%s #bullet%s", toplabel.Data(), (area_norm)? "Shape Errors" : "Absolute Errors"), 0.33, 0.94, 0.025 );
  //plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.94, 0.03, kBlue, 32, 11 );
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/muon_pz_with_sys_errorsummary_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  
  plotter->DrawDataMCErrorSummary(h_ptmu_mc_temp, h_ptmu_data_temp, "TR", false, true, 0.00001, area_norm );
  plotter->AddPlotLabel( Form("%s #bullet%s", toplabel.Data(), (area_norm)? "Shape Errors" : "Absolute Errors"), 0.33, 0.94, 0.025 );
  //plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.94, 0.03, kBlue, 32, 11 );
  utils->writeNorm( area_norm, pot_data, pot_mc, true ); 
  utils->Print( c, Form("ProtonSelection%s/muon_pt_with_sys_errorsummary_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  
  plotter->DrawDataMCErrorSummary(h_q22[kMC], h_q22[kData], "TR", false, true, 0.00001, area_norm );
  plotter->AddPlotLabel( Form("%s #bullet%s", toplabel.Data(), (area_norm)? "Shape Errors" : "Absolute Errors"), 0.33, 0.94, 0.025 );
  //plotter->AddPlotLabel( "MINER#nuA Preliminary", 0.6, 0.94, 0.03, kBlue, 32, 11 );
  utils->writeNorm( area_norm, pot_data, pot_mc, true );
  utils->Print( c, Form("ProtonSelection%s/q22_with_sys_errorsummary_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
*/

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
  plotter->DrawErrorSummary(h_neutron_dperp[kMC], "TL", false, true, 0.00001, area_norm, "Muon Reconstruction");
  utils->Print( c, Form("ProtonSelection%s/neutron_dperp_with_frac_sys_errorsummary_MuonRecons_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
  
  plotter->axis_maximum_group = 0.2; 
  plotter->DrawErrorSummary(h_neutron_phit[kMC], "TL", false, true, 0.00001, area_norm, "Muon Reconstruction");
  utils->Print( c, Form("ProtonSelection%s/neutron_phit_with_frac_sys_errorsummary_MuonRecons_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );

  plotter->axis_maximum_group = 0.2; 
  plotter->DrawErrorSummary(h_neutron_dreactplane[kMC], "TL", false, true, 0.00001, area_norm, "Muon Reconstruction");
  utils->Print( c, Form("ProtonSelection%s/neutron_dreactplane_with_frac_sys_errorsummary_MuonRecons_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );

  //------------------------------------------------------------------------------------------------
  // Drawing error summaries in the following groups "as absolute uncertainties" 
  // This was developed for testing the effects (or not) of the flux constraint implementations
  // You might need to tune the Y-axis to see all the absolute uncertainties  
  //------------------------------------------------------------------------------------------------
  plotter->axis_maximum_group = 2000.0; /*510 for NoFluxConstraint */
  plotter->DrawErrorSummary(h_neutron_dperp[kMC], "TL", false, true, 0.00001, area_norm, "Flux", false );
  utils->Print( c, Form("ProtonSelection%s/neutron_dperp_with_abs_sys_errorsummary_Flux_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
 
  plotter->axis_maximum_group = 2000.0; /*510 for NoFluxConstraint */
  plotter->DrawErrorSummary(h_neutron_phit[kMC], "TL", false, true, 0.00001, area_norm, "Flux", false );
  utils->Print( c, Form("ProtonSelection%s/neutron_phit_with_abs_sys_errorsummary_Flux_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );
 
  plotter->axis_maximum_group = 2000.0; /*510 for NoFluxConstraint */
  plotter->DrawErrorSummary(h_neutron_dreactplane[kMC], "TL", false, true, 0.00001, area_norm, "Flux", false );
  utils->Print( c, Form("ProtonSelection%s/neutron_dreactplane_with_abs_sys_errorsummary_Flux_%s_%s", sample.Data(), multiplicity_label.Data(), norm.Data()) );

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
        "\t-Area_Normalize\t\t =\t Area normalize Sample\n"<<
        "\t-MC_Plot_Style\t\t =\t QE, QE_PionInFS, QELike, QELike_split, QELike_PionInFS, QELike_split_PionInFS"<<std::endl;
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    return 0; 
  }

  //! Default parameters
  std::vector<std::string> par;
  par.push_back("ProtonSelectionPlots");
  par.push_back( Form("%s/ana/rootfiles/ProtonSelectionHists.root",getenv("CCQENUROOT") ) );
  par.push_back("0");
  par.push_back("QE");

  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;
  
  return ProtonSelectionPlots(par[1], atoi(par[2].c_str()),  par[3]);
}
