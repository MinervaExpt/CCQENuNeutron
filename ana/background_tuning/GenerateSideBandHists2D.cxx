#include "include/CCQENuPlotUtils.h"
#include "include/CCQENuUtils.h"
#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/TargetUtils.h"
#include "include/HDLFunc.h"
#include "include/GeneralFunc.h"
#include "include/CommonBins.h"

//#include "malloc.h"

using namespace CCQENU_ANA;


bool doRW=false;
int ifunc = 2;
const int nRegions = 6;


int getRegion2D(double dtheta2D, double qsq)
{
  if( abs(dtheta2D) > 40 ) return -1;
  if( abs(dtheta2D) <4 ) return 0;
  if( abs(dtheta2D) <12 ) return 1;
  if( dtheta2D <20 && dtheta2D > 12 ) return 2;
  if( dtheta2D <40 && dtheta2D > 20 ) return 3;
  if( -dtheta2D <20 && -dtheta2D > 12 ) return 4;
  if( -dtheta2D <40 && -dtheta2D > 20 ) return 5;
  return -1;
}


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


TH2D* convert3DToHDL( MnvH2D* h, TH3D* h3, HyperDimLinearizer* hdl, bool resetWeight = false )
{
  cout<<"entering convert3DToHDL"<<endl;
  //h3-> dthetaP vs dthetaR vs Q2qe
  //h -> dthetaP vs Q2qe vs dthetaR
  TH2D* ret = new TH2D( h->GetCVHistoWithStatError() );
  ret->Reset();
  for( int bx = 1; bx < h3->GetNbinsX()+1; bx++)
  {
    double x = h3->GetXaxis()->GetBinCenter( bx );
    for( int by = 1; by < h3->GetNbinsY()+1; by++ )
    {
      double y = h3->GetYaxis()->GetBinCenter( by );
      for( int bz = 1; bz < h3->GetNbinsZ()+1; bz++ )
      {
        double z = h3->GetZaxis()->GetBinCenter( bz );

        double globalx = hdl->GetBin( vector<double>( {x,z,y} ) ).first+0.0001;
        double w = h3->GetBinContent(bx,by,bz);
        int binx = ret->GetXaxis()->FindBin( globalx);
        int biny = ret->GetYaxis()->FindBin( z );
        if(w==0 && resetWeight) w = 1;
        ret->SetBinContent(binx, biny, w );


        //pair<int,int> xy = hdl->GetBin( vector<double>( {x,z,y} ) );
        //double w = h3->GetBinContent(bx,by,bz);
        //ret->SetBinContent( xy.first, xy.second, w );
      }
    }
  }
  return ret;
}

TH2D* reweighQELike_QE_OTH( MnvH2D** hist, TH3D* h3w, HyperDimLinearizer* hdl )
{
  TH2D* hw = convert3DToHDL( hist[kMC], h3w, hdl, true );
  
  MnvH2D* h = hist[kQELike_QE_OTH];
  hist[kMC]->Add( h, -1 );
  hist[kQELike_QE]->Add( h, -1 );

  cout<<"Multiplying Single"<<endl;
  h->MultiplySingle( h, hw );

  hist[kMC]->Add( h );
  hist[kQELike_QE]->Add( h );
  cout<<"done reweighting"<<endl;

  return hw;
}



//void reweighQELike_QE_OTH( MnvH2D** hist, TH3D* h3w, HyperDimLinearizer* hdl )
//{
//  cout<<"entering reweighQELike_QE_OTH"<<endl;
//  MnvH2D* h = hist[kQELike_QE_OTH];
//  TH2D* hw = new TH2D( h->GetCVHistoWithStatError() );
//  cout<<"created hw"<<endl;
//  hw->Reset();
//  for( int bx = 1; bx < h3w->GetNbinsX()+1; bx++)
//  {
//    double x = h3w->GetXaxis()->GetBinCenter( bx );
//    for( int by = 1; by < h3w->GetNbinsY()+1; by++ )
//    {
//      double y = h3w->GetYaxis()->GetBinCenter( by );
//      for( int bz = 1; bz < h3w->GetNbinsZ()+1; bz++ )
//      {
//        double z = h3w->GetZaxis()->GetBinCenter( bz );
//        pair<int,int> xy = hdl->GetBin( vector<double>( {x,z,y} ) );
//        double w = h3w->GetBinContent(bx,by,bz);
//        if(w==0) w = 1;
//        hw->SetBinContent( xy.first, xy.second, w );
//       
//        cout<<Form("%d, %d, %d --- %d,%d--- %f", bx,by,bz,xy.first,xy.second, w )<<endl;;
//      }
//    }
//  }
//
//  //for( int gx = 1; gx < h->GetNbinsX()+1; gx++ )
//  //{
//  //  cout<<"global x: "<< gx <<endl;
//  //  vector<int> bins = hdl->GetValues(gx);
//  //  int xbin = bins[0];
//  //  int ybin = bins[1];
//  //  for( int zbin = 1; zbin != h->GetNbinsY(); zbin++ )
//  //  {
//  //    double w = h3w->GetBinContent( xbin, ybin, zbin );
//  //    if(w == 0) w = 1;
//  //    hw->SetBinContent( gx, zbin, w );
//  //    cout<<Form("%d, %d, %d, %f", xbin, ybin, zbin, w )<<endl;;
//  //  }
//  //}
//  cout<<"Multiplying Single"<<endl;
//  hist[kMC]->Add( h, -1 );
//  h->MultiplySingle( h, hw );
//  hist[kMC]->Add( h );
//  hist[kQELike_QE]->MultiplySingle( hist[kQELike_QE], hw );
//  cout<<"done reweighting"<<endl;
//}

int GenerateSideBandHists( string signal_files, string signal_sideband, int nPlaylists = 12, bool applyFluxConstraint = true )
{
  TFile *f1 = new TFile( signal_files.c_str(), "READ");

  TVector2 *pot = (TVector2*)f1->Get("pot");
  TVector2 *n_events = (TVector2*)f1->Get("n_events");
  //for( int i = 2; i<nPlaylists+1; i++ ) (*pot)+= *((TVector2*)f1->Get(Form( "pot;%d", i) ));


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
  MnvH2D *h_dtheta2D_q2qe_signal[nHistos];
  cout<<"booking theta2d"<<endl;
  utils->bookHistos( f1, h_dtheta2D_q2qe_signal, Form( "h_dtheta2D_q2qe") );
  cout<<"booked"<<endl;

  //putils->scaleMCHistos( h_dtheta2D_q2qe_signal, pot_norm );

  MnvH1D *h_q2qe_signal[nHistos];
  utils->bookHistos( f1, h_q2qe_signal, Form( "h_q2qe") );
  //putils->scaleMCHistos( h_q2qe_signal, pot_norm );

  h_dtheta2D_q2qe_signal[kData]->AddMissingErrorBandsAndFillWithCV( *h_dtheta2D_q2qe_signal[kMC] );


  //define region histo
  vector<vector<MnvH1D*>>  h_regions(nRegions+1);

  cout<<"Begin cloning"<<endl;

  for( unsigned int i = 0; i< h_regions.size(); i++ )
  {
    int iRegion = i;
    if( i == nRegions ) iRegion = 99;
    for( unsigned int j = 0 ;j < nHistos; j++ )
    {
      string name = Form("h_q2qe_region_%02d_%s", iRegion, names[j].c_str() );
      MnvH1D* h = (MnvH1D*) h_q2qe_signal[j]->Clone( name.c_str() ) ;
      h->Reset();
      h_regions[i].push_back( h );
    }
  }
  cout<<"Begin regioning"<<endl;

  for( int x = 1; x <= h_dtheta2D_q2qe_signal[kData]->GetNbinsX(); x++ )
  {
    double dtheta2D = h_dtheta2D_q2qe_signal[kData]->GetXaxis()->GetBinCenter(x);
    int region = getRegion2D( dtheta2D, 0 );
    if(region == -1 ) region = nRegions;
    cout<<"Region "<<region<<" at "<<dtheta2D<<endl;
    //vector<MnvH1D*> projected_histos; //= hdlHists.ProjectionsY( slices, string("tmp_hists"), x,x,z,z);
    //for( int ih = 0; ih < nHistos; ih ++ ) projected_histos.push_back( (MnvH1D* ) h_dtheta2D_q2qe_signal[ih]->ProjectionY( Form("h_dtheta2D_q2qe_signal_%s_%02d", names[ih].c_str(), x ), x,x ) );
    //for( int ih = 0; ih < nHistos; ih ++ ) h_regions[region][ih]->Add( projected_histos[ih] );
    for( int ih = 0; ih < nHistos; ih ++ ) h_regions[region][ih]->Add( (MnvH1D* ) h_dtheta2D_q2qe_signal[ih]->ProjectionY( Form("h_dtheta2D_q2qe_signal_%s_%02d", names[ih].c_str(), x ), x,x ) );
  }

  cout<< "--------------------------------------------"  << endl;
  cout<< "Writing histograms to file"  << endl;
  cout<< "--------------------------------------------"  << endl;

  TFile *f_output = new TFile( signal_sideband.c_str(), "RECREATE" );
  f_output->cd();

  for(unsigned int i = 0; i< nHistos; i++ )
  {
    for( unsigned int r = 0; r<h_regions.size(); r++ )
    {
      h_regions[r][i]->Write();
      cout<<"Writing "<<h_regions[r][i]->GetName()<<endl;
    }
    //h_q2qe_angle_01[i]->Write();
    //h_q2qe_angle_02[i]->Write();
    //h_q2qe_angle_03[i]->Write();
    //h_q2qe_angle_04[i]->Write();
  }
  //Writing small angles
  cout<<"Write POT"<<endl;
  pot->Write("pot");
  cout<<"Write n_events"<<endl;
  n_events->Write("n_events");
  cout<<"Write f_output"<<endl;
  f_output->Write();
  cout<<"Close f_output"<<endl;
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
  //par.push_back( Form("%s/ana/rootfiles/SideBandFit.root",getenv("CCQENUROOT") ) );
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

  return GenerateSideBandHists( par[1], par[2],  stoi(par[3]), applyFluxConstraint);

}
