#include "include/GlobalParameters.h"
#include "include/CCQENuUtils.h"
#include "TParameter.h" 
#include "PlotUtils/HyperDimLinearizer.h"

#include "include/CCQENuPlotUtils.h"

#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/HyperDimLinearizer.h"

using namespace CCQENU_ANA;

MnvH1D* ScaleBackground(MnvH1D** hists, vector<int> &fit_hists, vector<MnvH1D*> &fits, MnvH1D* bck)
{
  vector<int> full_hists({kData,kMC, kQELike, kQELikeNot});
  vector<int> qelike_hists({kQELike_QE_H, kQELike_QE_OTH, kQELike_2p2h, kQELike_RES, kQELike_DIS});
  vector<int> qelikenot_hists({kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions});

  for( auto ptc = fit_hists.begin(); ptc != fit_hists.end(); ++ptc )
  {
    hists[*ptc]->Multiply( hists[*ptc], fits[ptc-fit_hists.begin() ] );
  }

  hists[kMC]->Reset();
  hists[kQELikeNot]->Reset();
  hists[kQELike]->Reset();
  for(auto c : qelike_hists) hists[kQELike]->Add( hists[c] );
  for(auto c : qelikenot_hists) hists[kQELikeNot]->Add( hists[c] );
  hists[kMC]->Add( hists[kQELikeNot] );
  hists[kMC]->Add( hists[kQELike] );

  

  bck = hists[kMC]->Clone("bck");
  bck->Add(hists[kQELike_QE_H],-1);
  MnvH1D* data = hists[kData]->Clone(Form("%s_bcksub", hists[kData]->GetName()));
  data->AddMissingErrorBandsAndFillWithCV(*hists[kMC]);
  data->Add(bck,-1);
  return data;

}

int SideBandFitHistsCycle( string output_filename,string filename_signal,string filename_bg_weights, bool makeFluxConstraintHisto, bool applyFluxConstraint)
{
  TH1::SetDefaultSumw2();
 
  cout<<"Begin SideBandFitHists"<<endl;
  // read file to get histograms
  TFile *f_signal = new TFile( filename_signal.c_str(), "READ" );
  if (f_signal->IsZombie() || f_signal->GetListOfKeys()->IsEmpty()){
    Error("SideBandFitPlots","Could not get Signal histogram ROOT file or it was empty.");
    return 1;
  }

  TFile *f_bg_weights = new TFile( filename_bg_weights.c_str(), "READ");
  if (f_bg_weights->IsZombie() || f_bg_weights->GetListOfKeys()->IsEmpty()){
    Error("SideBandFitPlots","Could not get Background histogram ROOT file or it was empty.");
    return 1;
  }

  cout<<"Define CCQENuUtils"<<endl;
  CCQENuUtils *utils = new CCQENuUtils( false, makeFluxConstraintHisto );
  cout<<"CCQENuUtils setplaylist"<<endl;
  utils->setPlaylist("minervame6A");//Dominant LE playlist. Need to figure out if there is some correction needed here...
  cout<<"CCQENuUtils setfluxreweighter"<<endl;
  utils->setFluxReweighterPlaylist();
  CCQENuPlotUtils *plotUtils = new CCQENuPlotUtils();

  //--------------------------------------------
  // Get POT Normalization factor:
  // This is necessary because Background Subtraction
  // uses POT Normalized MC
  //--------------------------------------------
  double pot_data = plotUtils->getPOTData( f_signal );
  double pot_mc = plotUtils->getPOTMC( f_signal );
  double pot_scale = plotUtils->getPOTNormFactor( f_signal );
  cout<< "POT of Data Events = " << pot_data << endl;
  cout<< "POT of MC Events = " << pot_mc << endl;
   //--------------------------------------------
  // Creating the output file
  //--------------------------------------------
  TFile *f_output = new TFile( output_filename.c_str(), "RECREATE");
  cout<<"Created output file"<<endl;

  //--------------------------------------------
  // Creating holder MnvH1D/2D Objects
  //--------------------------------------------
  // distribution histso
  MnvH1D *h_region_00[nHistos],*h_region_01[nHistos];
  // weights histos
  MnvH1D *h_bg_weights_qelike_qe_oth,*h_bg_weights_qelike_2p2h,*h_bg_weights_qelikenot_scp,*h_bg_weights_qelikenot_snp,*h_bg_weights_qelikenot_mp; 
  // background tuning histos
  MnvH1D *h_data_bck_region_00, *h_data_bck_region_01;
  MnvH1D *h_data_bcksub_region_00, *h_data_bcksub_region_01;

  //--------------------------------------------
  // Looping over variables fitted
  //--------------------------------------------

  plotUtils->bookHistos( f_signal, h_region_00, Form("h_q2qe_region_00") );
  plotUtils->bookHistos( f_signal, h_region_01, Form("h_q2qe_region_01") );
  cout<<"booked angular histograms"<<endl;
  
  h_bg_weights_qelike_qe_oth = (MnvH1D*) f_bg_weights->Get("hs_weights_yvar_bgType_qe_oth");
  h_bg_weights_qelike_2p2h =   (MnvH1D*) f_bg_weights->Get("hs_weights_yvar_bgType_2p2h");
  h_bg_weights_qelikenot_scp = (MnvH1D*) f_bg_weights->Get("hs_weights_yvar_bgType_qelikenot_scp");
  h_bg_weights_qelikenot_snp = (MnvH1D*) f_bg_weights->Get("hs_weights_yvar_bgType_qelikenot_snp");
  h_bg_weights_qelikenot_mp =  (MnvH1D*) f_bg_weights->Get("hs_weights_yvar_bgType_qelikenot_mp");
  cout<<"booked  bck weights"<<endl;
  
  vector<int> fit_hists_id({ kQELike_QE_OTH, kQELike_2p2h, kQELikeNot_SingleNeutralPion, kQELikeNot_SingleChargedPion, kQELikeNot_MultiPion });
  vector<MnvH1D*> fit_hists({ h_bg_weights_qelike_qe_oth,h_bg_weights_qelike_2p2h,h_bg_weights_qelikenot_snp,h_bg_weights_qelikenot_scp,h_bg_weights_qelikenot_mp});

  h_data_bcksub_region_00 = ScaleBackground( h_region_00, fit_hists_id, fit_hists, h_data_bck_region_00 );
  h_data_bcksub_region_01 = ScaleBackground( h_region_01, fit_hists_id, fit_hists, h_data_bck_region_01 );
  cout<<"scaled histograms"<<endl;

  h_data_bcksub_region_00->Write();
  h_data_bcksub_region_01->Write();
  cout<<"saved"<<endl;
  


  f_bg_weights->Close();
  f_signal->Close();
  f_output->Close();

  //delete f_output; 
  //delete f_signal;
  //delete f_bg_weights; 

  //delete plotUtils; 
  //delete utils;

  cout<< "All done..."  << endl;
  return 0;
};





int main( int argc, char *argv[])
{
  ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if (argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t./SideBandFitHists  Output_file Input_file_2TRACK  Input_file_BKG  Input_file_MIGRATION  Input_file_EFFICIENCY  Input_file_FLUX  Num_iterations  Make_Flux_Constraint_Histo  Apply_Flux_Constraint \n\n"<<
      "\t-output_file\t =\t Cross Section output file\n"<<
      "\t-input_file_angles\t =\t Muon Event Selection(Signal) for >=2 track sub-sample\n"<<
      "\t-input_file_angle_bkg\t =\t Background weights filename\n"<<
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
  par.push_back( Form("%s/ana/rootfiles/SideBandFitHists.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/MuonSelectionHists_Angles.root",getenv("CCQENUROOT") ) );
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

  return SideBandFitHistsCycle( par[1], par[2], par[3], fluxConstraintHistoNeeded, applyFluxConstraint );

}
