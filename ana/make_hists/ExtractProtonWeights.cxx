#include "include/CCQENuPlotUtils.h"

#include "include/CCQENuUtils.h"
#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/HyperDimLinearizer.h"

#include "include/CommonBins.h"
#include "include/GeneralFunc.h"

using namespace CCQENU_ANA;


int ExpandHDL( string output_filename,string input_filename)
{
  TFile *input_file = TFile::Open(input_filename.c_str(), "READ");

  if ( input_file->IsZombie() || input_file->GetListOfKeys()->IsEmpty() ){
    Error("ExpandHDL","Could not get histogram ROOT file or it was empty.");
    return 1;
  }

  CCQENuPlotUtils *plotUtils = new CCQENuPlotUtils();

  //MnvH2D* h_dthetaPq2qe_ptmu[nHistos];
  MnvH2D* h_bkg_signal = (MnvH2D*) input_file->Get("h_tuned_bkg_dthetaPq2qe_ptmu_signal");
  MnvH2D* h_mc_weighted_nobck_signal = (MnvH2D*) input_file->Get("h_dthetaPq2qe_ptmu_mc_weighted_nobck_signal");
  MnvH2D* h_data_weighted_nobck_signal = (MnvH2D*) input_file->Get("h_dthetaPq2qe_ptmu_data_weighted_nobck_signal");
  HyperDimLinearizer *hdl = GetHDL( dthetaPerpbins, muonPtbins, Q2bins ); 

  vector<MnvH2D*> bkg= hdl->Get2DMnvHistos(  h_bkg_signal, true);
  vector<MnvH2D*> mc =hdl->Get2DMnvHistos( h_mc_weighted_nobck_signal, true );
  vector<MnvH2D*> data =hdl->Get2DMnvHistos( h_data_weighted_nobck_signal, true );

  MnvH2D* bkg_sum = (MnvH2D*) bkg[0]->Clone("h_dthetaPq2qe_ptmu_bkg_sum");     bkg_sum->SetTitle("h_dthetaPq2qe_ptmu_sum");
  MnvH2D* mc_sum = (MnvH2D*) mc[0]->Clone("h_dthetaPq2qe_ptmu_mc_sum");        mc_sum->SetTitle("h_dthetaPq2qe_ptmu_sum");
  MnvH2D* data_sum = (MnvH2D*) data[0]->Clone("h_dthetaPq2qe_ptmu_data_sum");  data_sum->SetTitle("h_dthetaPq2qe_ptmu_sum");

  TFile* output_file = TFile::Open(output_filename.c_str() , "RECREATE");
  for( unsigned int i = 0; i< bkg.size(); ++i )
  {
    double ptmu = muonPtbins.bin_edges[i];
    bkg[i]->SetTitle(Form("h_dthetaPq2qe_ptmu_bkg_%03d_%.3f",i, ptmu));
    bkg[i]->Write(Form("h_dthetaPq2qe_ptmu_bkg_%03d",i));
    if(i>0) bkg_sum->Add(bkg[i]);
  }
  for( unsigned int i = 0; i< mc.size(); ++i )
  {
    double ptmu = muonPtbins.bin_edges[i];
    mc[i]->SetTitle(Form("h_dthetaPq2qe_ptmu_mc_%03d_%.3f",i, ptmu));
    mc[i]->Write(Form("h_dthetaPq2qe_ptmu_mc_%03d",i));
    if(i>0) mc_sum->Add(mc[i]);
  }
  for( unsigned int i = 0; i< data.size(); ++i )
  {
    double ptmu = muonPtbins.bin_edges[i];
    data[i]->SetTitle(Form("h_dthetaPq2qe_ptmu_data_%03d_%.3f",i, ptmu));
    data[i]->Write(Form("h_dthetaPq2qe_ptmu_data_%03d",i));
    if(i>0) data_sum->Add(data[i]);
  }

  bkg_sum->Write();
  mc_sum->Write();
  data_sum->Write();

  MnvH1D* dataMCRatio = (MnvH1D*) data_sum->ProjectionX("data_mc_ratio_dthetaP");
  dataMCRatio->Divide(dataMCRatio, mc_sum->ProjectionX() );
  MnvH2D* dataMCRatio2D = (MnvH2D*) mc_sum->Clone("data_mc_ratio_dthetaP_q2qe");
  dataMCRatio2D->ClearAllErrorBands();
  dataMCRatio2D->ClearSysErrorMatrices();
  dataMCRatio2D->Reset();
  ExpandHisto( dataMCRatio, dataMCRatio2D, 0 );

  dataMCRatio->Write();
  dataMCRatio2D->Write();
  MnvH2D* mc_sum_rw = (MnvH2D*) mc_sum->Clone("h_dthetaP_q2qe_mc_rw");
  mc_sum_rw->Multiply( mc_sum_rw, dataMCRatio2D );
  mc_sum_rw->Write();




  output_file->Close();
  input_file->Close();


  return 0;
}





int main( int argc, char *argv[])
{
  ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if (argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t./SideBandFitHists  Output_file Input_file_2TRACK  Input_file_BKG  Input_file_MIGRATION  Input_file_EFFICIENCY  Input_file_FLUX  Num_iterations  Make_Flux_Constraint_Histo  Apply_Flux_Constraint \n\n"<<
      "\t-output_file\t =\t Cross Section output file\n"<<
      "\t-input_bg_subtracted\t =\t bkgd subtracted file\n"<<
       "\t********************************************************************************************** \n"<<
      "\t Please see : MuonSelectionHists.cxx, BackgroundWeights.cxx, MigrationMatrixHists.cxx, EffPurityHists.cxx and Ana/Flux/python/compute_flux.py for getting the necessary input files"<< std::endl; 
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    return 0; 
  }

  //! Default parameters
  std::vector<std::string> par;
  par.push_back("SideBandFitPlots");
  par.push_back( Form("%s/ana/rootfiles/minervame1D/ProtonWeights.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/minervame1D/bkgd_sub_hists_ptmu_2D_optimized_bins.root",getenv("CCQENUROOT") ) );

  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }
  return ExpandHDL( par[1], par[2] );

}



