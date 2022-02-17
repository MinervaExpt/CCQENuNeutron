#include "include/CCQENuPlotUtils.h"
#include "include/CCQENuUtils.h"
#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/TargetUtils.h"
#include "include/HDLFunc.h"
#include "include/GeneralFunc.h"
#include "include/CommonBins.h"

using namespace CCQENU_ANA;

const int nRegions = 6;
int getRegion(double dthetaP, double dthetaR)
{
  if( abs(dthetaP) > 55 || dthetaR> 10 || -55 > dthetaR )          return -1;
  if( abs(dthetaP) < 10 && 10 > dthetaR && dthetaR > -10 )           return 0; //hydrogen
  else if ( abs(dthetaP) < 20 && abs(dthetaR)<20 )               return 1; //QE-C 1
  else if ( abs(dthetaP) < 30 && abs(dthetaR)<30 )              return 2; //QE-C 3
  else if ( abs(dthetaP) < 20 && dthetaR<-30 && dthetaR >=-40  ) return 3; //2p2h/RES 4
  else if ( abs(dthetaP) > 20 && dthetaR<-30 && dthetaR >=-40  ) return 4; //2p2h/RES 1
  else if ( abs(dthetaP) < 55 && dthetaR>-55 && dthetaR < -40  ) return 5; //2p2h/RES 2
  else return -1;
}

int GenerateSideBandHists( string signal_files, string background_weights, string signal_sideband, int nPlaylists = 12, bool applyFluxConstraint = true )
{
  TFile *f1 = new TFile( signal_files.c_str(), "READ");
  TFile *f2 = new TFile( background_weights.c_str(), "READ");

  TVector2 *pot = (TVector2*)f1->Get("pot;1");
  for( int i = 2; i<nPlaylists+1; i++ ) (*pot)+= *((TVector2*)f1->Get(Form( "pot;%d", i) ));


  //TVector2 *pot = (TVector2*)f1->Get("pot");
  double pot_data = pot->X();
  double pot_mc = pot->Y();
  double pot_norm = pot_data/pot_mc;
  cout<<"===================== POT ==================="<<endl;
  cout<<"Data: "<<pot_data<<endl;
  cout<<"MC: "<<pot_mc<<endl;
  cout<<"Ratio: "<<pot_norm<<endl;

  CCQENuPlotUtils *putils = new CCQENuPlotUtils( applyFluxConstraint );
  CCQENuUtils *utils = new CCQENuUtils( false, applyFluxConstraint );

  axis_binning xbins = dthetaPerpbins;
  axis_binning ybins = Q2bins;
  axis_binning zbins = dthetaReactbins;

  HyperDimLinearizer* hdl = GetHDL(xbins , ybins, zbins);

  // Get CV histos
  MnvH2D *h_hdl_signal[nHistos];
  cout<<"booking h2d"<<endl;
  utils->bookHistos( f1, h_hdl_signal, Form( "h_dthetaPdthetaR_q2qe") );
  putils->scaleMCHistos( h_hdl_signal, pot_norm );

  MnvH1D *h_q2qe_signal[nHistos];
  utils->bookHistos( f1, h_q2qe_signal, Form( "h_q2qe") );
  putils->scaleMCHistos( h_q2qe_signal, pot_norm );


  MnvH2D *h_bck_weights_qelikenot = (MnvH2D*) f2->Get("h_weights_dthetaPdthetaR_yvarbins_qelikenot_Signal");
  MnvH2D *h_bck_weights_scp = (MnvH2D*) f2->Get("hs_weights_dthetaPdthetaR_yvarbins_bgType_SingleChargedPion");
  MnvH2D *h_bck_weights_snp = (MnvH2D*) f2->Get("hs_weights_dthetaPdthetaR_yvarbins_bgType_SingleNeutralPion");
  MnvH2D *h_bck_weights_mp  = (MnvH2D*) f2->Get("hs_weights_dthetaPdthetaR_yvarbins_bgType_MultiPions");

  cout<<"booked weights"<<endl;
  cout<<"h_bck_weights_qelikenot integral: "<<h_bck_weights_qelikenot->Integral()<<endl;

  h_hdl_signal[kMC]->Add(h_hdl_signal[kQELikeNot],-1 );

  h_hdl_signal[kQELikeNot_SingleChargedPion]->Multiply( h_hdl_signal[kQELikeNot_SingleChargedPion], h_bck_weights_scp );
  h_hdl_signal[kQELikeNot_SingleNeutralPion]->Multiply( h_hdl_signal[kQELikeNot_SingleNeutralPion], h_bck_weights_snp );
  h_hdl_signal[kQELikeNot_MultiPion]->Multiply( h_hdl_signal[kQELikeNot_MultiPion], h_bck_weights_mp );

  h_hdl_signal[kQELikeNot]->Reset();
  h_hdl_signal[kQELikeNot]->Add( h_hdl_signal[kQELikeNot_SingleChargedPion] );
  h_hdl_signal[kQELikeNot]->Add( h_hdl_signal[kQELikeNot_SingleNeutralPion] );
  h_hdl_signal[kQELikeNot]->Add( h_hdl_signal[kQELikeNot_MultiPion] );
  h_hdl_signal[kQELikeNot]->Add( h_hdl_signal[kQELikeNot_NoPions] );

  //h_bck_weights_qelikenot->Reset();
  //h_bck_weights_qelikenot->Add( h_hdl_signal[kQELikeNot_SingleChargedPion] );
  //h_bck_weights_qelikenot->Add( h_hdl_signal[kQELikeNot_SingleNeutralPion] );
  //h_bck_weights_qelikenot->Add( h_hdl_signal[kQELikeNot_MultiPion] );
  //h_bck_weights_qelikenot->Add( h_hdl_signal[kQELikeNot_NoPions] );
  //h_bck_weights_qelikenot->Divide( h_bck_weights_qelikenot, h_hdl_signal[kQELikeNot] );

  //h_hdl_signal[kQELikeNot]->Multiply( h_hdl_signal[kQELikeNot], h_bck_weights_qelikenot );


  h_hdl_signal[kMC]->Add(h_hdl_signal[kQELikeNot] );
  h_hdl_signal[kData]->Reset();
  TH2D* data_qe = (TH2D*) h_hdl_signal[kQELike_QE_OTH]->Clone("data_qe");
  TH2D* data_res = (TH2D*) h_hdl_signal[kQELike_RES]->Clone("data_res");
  TH2D* data_2p2h = (TH2D*) h_hdl_signal[kQELike_2p2h]->Clone("data_2p2h");
  TH2D* data_dis = (TH2D*) h_hdl_signal[kQELike_DIS]->Clone("data_dis");
  TH2D* data_oth = (TH2D*) h_hdl_signal[kQELike_OTH]->Clone("data_oth");
  TH2D* data_bck = (TH2D*) h_hdl_signal[kQELikeNot]->Clone("data_qelikenot");
  h_hdl_signal[kData]->Add( data_qe,1 );
  h_hdl_signal[kData]->Add( data_res,1.2 );
  h_hdl_signal[kData]->Add( data_2p2h,.8 );
  h_hdl_signal[kData]->Add( data_dis,1 );
  h_hdl_signal[kData]->Add( data_oth,1 );
  h_hdl_signal[kData]->Add( data_bck,1 );

  //h_hdl_signal[kData]->AddMissingErrorBandsAndFillWithCV( *h_hdl_signal[kMC] );
  //h_hdl_signal[kData]->Add(h_hdl_signal[kQELikeNot],-1 );

  HDLHistos hdlHists( hdl, xbins, ybins, zbins, h_hdl_signal, true );

  //define region histo
  vector<vector<MnvH1D*>>  h_regions(nRegions);
  vector<vector<MnvH2D*>> slices = hdlHists.GetHistos();
  for( unsigned int i = 0; i< h_regions.size(); i++ )
  {
    for( unsigned int j = 0 ;j < nHistos; j++ )
    {
      string name = Form("h_q2qe_region_%02d_%s", i, names[j].c_str() );
      MnvH1D* h = (MnvH1D*) h_q2qe_signal[j]->Clone( name.c_str() ) ;
      h->Reset();
      h_regions[i].push_back( h );
    }
  }

  for( int x = 1; x <= hdlHists.GetAxisHisto()->GetNbinsX(); x++ )
  {
    double dthetaP = hdlHists.GetXaxis()->GetBinCenter(x);
    if( abs(dthetaP)>55 ) continue;
    for( int z = 1; z <= hdlHists.GetAxisHisto()->GetNbinsZ(); z++ )
    {
      double dthetaR = hdlHists.GetZaxis()->GetBinCenter(z);
      int region = getRegion(dthetaR, dthetaP);
      if (region == -1 ) continue;

      vector<MnvH1D*> projected_histos = hdlHists.ProjectionsY( slices, string("tmp_hists"), x,x,z,z);
      hdlHists.Add( h_regions[region], projected_histos);
    }
  }
  cout<<"h_regions[0]: "<<h_regions[0][0]->Integral()<<endl;

  cout<< "--------------------------------------------"  << endl;
  cout<< "Writing histograms to file"  << endl;
  cout<< "--------------------------------------------"  << endl;

  TFile *f_output = new TFile( signal_sideband.c_str(), "RECREATE" );
  f_output->cd();
  for( int i = 0; i<nRegions; i++ )
  {
    for(unsigned int j = 0; j< nHistos; j++ )
    {
      h_regions[i][j]->Write();
      cout<<"Writing "<<h_regions[i][j]->GetName()<<endl;
    }
  }
  pot->Write("pot");
  f_output->Write();
  f_output->Close();


  return 0;
};



int main( int argc, char *argv[])
{
  ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if (argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t./GenerateSideBandHists  Input_file  Input_background_weights Output_Signal_SideBand nPlaylist applyFluxConstraints\n\n"<<
      "\t********************************************************************************************** \n"<<
      "\t Please see : MuonSelectionHists.cxx, BackgroundWeights.cxx, MigrationMatrixHists.cxx, EffPurityHists.cxx and Ana/Flux/python/compute_flux.py for getting the necessary input files"<< std::endl; 
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    return 0; 
  }

  //! Default parameters
  std::vector<std::string> par;
  par.push_back("SideBandFitPlots");
  par.push_back( Form("%s/ana/rootfiles/MuonSelectionHist.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/SideBandFit.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/SignalSideBand.root",getenv("CCQENUROOT") ) );
  par.push_back( "12" );
  par.push_back( "1" );
  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;


  bool applyFluxConstraint = ( par.back() == "1" ) ? true : false; 

  return GenerateSideBandHists( par[1], par[2], par[3], stoi(par[4]), applyFluxConstraint);

}
