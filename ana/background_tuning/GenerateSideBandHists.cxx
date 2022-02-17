#include "include/CCQENuPlotUtils.h"
#include "include/CCQENuUtils.h"
#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/TargetUtils.h"
#include "include/HDLFunc.h"
#include "include/GeneralFunc.h"
#include "include/CommonBins.h"

using namespace CCQENU_ANA;

bool doRW=false;
int ifunc = 2;
const int nRegions = 6;
//int getRegion(double dthetaP, double dthetaR)
//{
//  if( abs(dthetaP) > 55 || dthetaR> 10 || -55 > dthetaR )          return -1;
//  if( abs(dthetaP) < 10 && 10 > dthetaR && dthetaR > -10 )           return 0; //hydrogen
//  else if ( abs(dthetaP) < 20 && abs(dthetaR)<20 )               return 1; //QE-C 1
//  else if ( abs(dthetaP) < 30 && abs(dthetaR)<30 )              return 2; //QE-C 3
//  else if ( abs(dthetaP) < 20 && dthetaR<-30 && dthetaR >=-40  ) return 3; //2p2h/RES 4
//  else if ( abs(dthetaP) > 20 && dthetaR<-30 && dthetaR >=-40  ) return 4; //2p2h/RES 1
//  else if ( abs(dthetaP) < 55 && dthetaR>-55 && dthetaR < -40  ) return 5; //2p2h/RES 2
//  else return -1;
//}

int getRegion(double dthetaP, double dthetaR)
{
  if( abs(dthetaP) > 55 || dthetaR> 30 || -55 > dthetaR )          return -1;
  if( abs(dthetaP) < 10 && abs(dthetaR)<10)           return 0; //hydrogen
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
  MnvH2D *h_hdl_signal[nHistos];
  cout<<"booking h2d"<<endl;
  utils->bookHistos( f1, h_hdl_signal, Form( "h_dthetaPdthetaR_q2qe") );
  cout<<"booked"<<endl;

  //putils->scaleMCHistos( h_hdl_signal, pot_norm );
  MnvH1D *h_q2qe_angle_01[nHistos]; utils->bookHistos( f1, h_q2qe_angle_01, "h_q2qe_angle_01" );
  MnvH1D *h_q2qe_angle_02[nHistos]; utils->bookHistos( f1, h_q2qe_angle_02, "h_q2qe_angle_02" );
  MnvH1D *h_q2qe_angle_03[nHistos]; utils->bookHistos( f1, h_q2qe_angle_03, "h_q2qe_angle_03" );
  MnvH1D *h_q2qe_angle_04[nHistos]; utils->bookHistos( f1, h_q2qe_angle_04, "h_q2qe_angle_04" );

  MnvH2D* h_enu_q2qe[nHistos]; utils->bookHistos( f1, h_enu_q2qe, "h_enu_q2qe" );

  MnvH2D *h_q2qe_ptmu_region_00[nHistos]; utils->bookHistos(f1, h_q2qe_ptmu_region_00, "h_q2qe_ptmu_region_00");
  MnvH2D *h_q2qe_ptmu_region_01[nHistos]; utils->bookHistos(f1, h_q2qe_ptmu_region_01, "h_q2qe_ptmu_region_01");
  MnvH2D *h_q2qe_ptmu_region_02[nHistos]; utils->bookHistos(f1, h_q2qe_ptmu_region_02, "h_q2qe_ptmu_region_02");
  MnvH2D *h_q2qe_ptmu_region_03[nHistos]; utils->bookHistos(f1, h_q2qe_ptmu_region_03, "h_q2qe_ptmu_region_03");
  MnvH2D *h_q2qe_ptmu_region_04[nHistos]; utils->bookHistos(f1, h_q2qe_ptmu_region_04, "h_q2qe_ptmu_region_04");
  MnvH2D *h_q2qe_ptmu_region_05[nHistos]; utils->bookHistos(f1, h_q2qe_ptmu_region_05, "h_q2qe_ptmu_region_05");
  MnvH2D *h_q2qe_ptmu_region_99[nHistos]; utils->bookHistos(f1, h_q2qe_ptmu_region_99, "h_q2qe_ptmu_region_99");
  MnvH2D *h_q2qe_ptmu_region_non99[nHistos]; utils->bookHistos(f1, h_q2qe_ptmu_region_non99, "h_q2qe_ptmu_region_non99");


  //MnvH1D *h_q2qe_vtx_region_00[nHistos]; utils->bookHistos( f1, h_q2qe_vtx_region_00, "h_q2qe_vtx_region_00" );
  //MnvH1D *h_q2qe_vtx_region_01[nHistos]; utils->bookHistos( f1, h_q2qe_vtx_region_01, "h_q2qe_vtx_region_01" );
  //MnvH1D *h_q2qe_vtx_region_02[nHistos]; utils->bookHistos( f1, h_q2qe_vtx_region_02, "h_q2qe_vtx_region_02" );
  //MnvH1D *h_q2qe_vtx_region_03[nHistos]; utils->bookHistos( f1, h_q2qe_vtx_region_03, "h_q2qe_vtx_region_03" );
  //MnvH1D *h_q2qe_vtx_region_04[nHistos]; utils->bookHistos( f1, h_q2qe_vtx_region_04, "h_q2qe_vtx_region_04" );
  //MnvH1D *h_q2qe_vtx_region_05[nHistos]; utils->bookHistos( f1, h_q2qe_vtx_region_05, "h_q2qe_vtx_region_05" );
  //MnvH1D *h_q2qe_vtx_region_99[nHistos]; utils->bookHistos( f1, h_q2qe_vtx_region_99, "h_q2qe_vtx_region_99" );

  //MnvH1D *h_q2qe_nonvtx_region_00[nHistos]; utils->bookHistos( f1, h_q2qe_nonvtx_region_00, "h_q2qe_nonvtx_region_00" );
  //MnvH1D *h_q2qe_nonvtx_region_01[nHistos]; utils->bookHistos( f1, h_q2qe_nonvtx_region_01, "h_q2qe_nonvtx_region_01" );
  //MnvH1D *h_q2qe_nonvtx_region_02[nHistos]; utils->bookHistos( f1, h_q2qe_nonvtx_region_02, "h_q2qe_nonvtx_region_02" );
  //MnvH1D *h_q2qe_nonvtx_region_03[nHistos]; utils->bookHistos( f1, h_q2qe_nonvtx_region_03, "h_q2qe_nonvtx_region_03" );
  //MnvH1D *h_q2qe_nonvtx_region_04[nHistos]; utils->bookHistos( f1, h_q2qe_nonvtx_region_04, "h_q2qe_nonvtx_region_04" );
  //MnvH1D *h_q2qe_nonvtx_region_05[nHistos]; utils->bookHistos( f1, h_q2qe_nonvtx_region_05, "h_q2qe_nonvtx_region_05" );
  //MnvH1D *h_q2qe_nonvtx_region_99[nHistos]; utils->bookHistos( f1, h_q2qe_nonvtx_region_99, "h_q2qe_nonvtx_region_99" );


  //if( !h_q2qe_vtx_region_00[0] )
  //{
  //  cout<<"HISTO DOESN'T EXIST!"<<endl;
  //  return 1;
  //}


  MnvH1D *h_q2qe_signal[nHistos], *h_q2qe_non99[nHistos];
  utils->bookHistos( f1, h_q2qe_signal, Form( "h_q2qe") );
  utils->bookHistos( f1, h_q2qe_non99, Form( "h_q2qe") );
  //putils->scaleMCHistos( h_q2qe_signal, pot_norm );

  h_hdl_signal[kData]->AddMissingErrorBandsAndFillWithCV( *h_hdl_signal[kMC] );



  HDLHistos hdlHists( hdl, xbins, ybins, zbins, h_hdl_signal, true );

  //define region histo
  vector<vector<MnvH1D*>>  h_regions(nRegions+1);
  vector<vector<MnvH2D*>> slices = hdlHists.GetHistos();
      //{
      //  MnvH2D* ret = slices[5][kMC];
      //  cout<<"hMC: "<<ret->GetName()<<" "<<ret->Integral()<<endl;
      //  vector<string> vertNames = ret->GetVertErrorBandNames();
      //  for( unsigned int k = 0; k < vertNames.size(); ++k )
      //  {
      //    //cout<<vertNames[k]<<",";
      //    unsigned int nUniverses = ret->GetVertErrorBand( vertNames[k] )->GetNHists();
      //    for (unsigned int iUniv = 0 ;iUniv< nUniverses; iUniv++) {
      //      TH2D* r = ret->GetVertErrorBand( vertNames[k] )->GetHist(iUniv);
      //      cout<<"vert: "<<vertNames[k]<<" : "<<r->Integral()<<endl;
      //    }
      //  }
      //}






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

  for( int x = 1; x <= hdlHists.GetAxisHisto()->GetNbinsX(); x++ )
  {
    double dthetaP = hdlHists.GetXaxis()->GetBinCenter(x);
    //if( abs(dthetaP)>55 ) continue;
    for( int z = 1; z <= hdlHists.GetAxisHisto()->GetNbinsZ(); z++ )
    {
      double dthetaR = hdlHists.GetZaxis()->GetBinCenter(z);
      int region = getRegion(dthetaP, dthetaR);
      //if (region == -1 ) continue;
      if (region == -1) region = nRegions;
      cout<<"Region "<<region<<" at "<<dthetaP<<" "<<dthetaR<<endl;

      vector<MnvH1D*> projected_histos = hdlHists.ProjectionsY( slices, string("tmp_hists"), x,x,z,z);

      hdlHists.Add( h_regions[region], projected_histos);

      for( vector<MnvH1D*>::iterator it = projected_histos.begin(); it!=projected_histos.end();++it)
      {
        (*it)->Reset();
        delete (*it);
      }
      projected_histos.clear();
    }
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

    h_q2qe_non99[i]->Reset();
    h_q2qe_non99[i]->Add( h_regions[0][i] );
    h_q2qe_non99[i]->Add( h_regions[1][i] );
    h_q2qe_non99[i]->Add( h_regions[2][i] );
    h_q2qe_non99[i]->Add( h_regions[3][i] );
    h_q2qe_non99[i]->Add( h_regions[4][i] );
    h_q2qe_non99[i]->Add( h_regions[5][i] );

    h_q2qe_non99[i]->SetName( Form("h_q2qe_non99_%s", names[i].c_str() ) );
    h_q2qe_non99[i]->Write();



    h_q2qe_angle_01[i]->Write();
    h_q2qe_angle_02[i]->Write();
    h_q2qe_angle_03[i]->Write();
    h_q2qe_angle_04[i]->Write();

    //cout<<" Writing "<<h_q2qe_vtx_region_00[i]->GetName()<<endl;
    //h_q2qe_vtx_region_00[i]->Write(); 
    //h_q2qe_vtx_region_01[i]->Write();
    //h_q2qe_vtx_region_02[i]->Write();
    //h_q2qe_vtx_region_03[i]->Write();
    //h_q2qe_vtx_region_04[i]->Write();
    //h_q2qe_vtx_region_05[i]->Write();
    //h_q2qe_vtx_region_99[i]->Write();

    //cout<<" Writing "<<h_q2qe_nonvtx_region_00[i]->GetName()<<endl;
    //h_q2qe_nonvtx_region_00[i]->Write(); 
    //h_q2qe_nonvtx_region_01[i]->Write();
    //h_q2qe_nonvtx_region_02[i]->Write();
    //h_q2qe_nonvtx_region_03[i]->Write();
    //h_q2qe_nonvtx_region_04[i]->Write();
    //h_q2qe_nonvtx_region_05[i]->Write();
    //h_q2qe_nonvtx_region_99[i]->Write();

    h_q2qe_signal[i]->Write();
    h_enu_q2qe[i]->Write();



    h_q2qe_ptmu_region_00[i]->ProjectionY(Form( "h_ptmu_region_00_%s", names[i].c_str() ))->Write();
    h_q2qe_ptmu_region_01[i]->ProjectionY(Form( "h_ptmu_region_01_%s", names[i].c_str() ))->Write(); 
    h_q2qe_ptmu_region_02[i]->ProjectionY(Form( "h_ptmu_region_02_%s", names[i].c_str() ))->Write(); 
    h_q2qe_ptmu_region_03[i]->ProjectionY(Form( "h_ptmu_region_03_%s", names[i].c_str() ))->Write(); 
    h_q2qe_ptmu_region_04[i]->ProjectionY(Form( "h_ptmu_region_04_%s", names[i].c_str() ))->Write(); 
    h_q2qe_ptmu_region_05[i]->ProjectionY(Form( "h_ptmu_region_05_%s", names[i].c_str() ))->Write(); 
    h_q2qe_ptmu_region_99[i]->ProjectionY(Form( "h_ptmu_region_99_%s", names[i].c_str() ))->Write(); 
    h_q2qe_ptmu_region_non99[i]->ProjectionY(Form( "h_ptmu_region_non99_%s", names[i].c_str() ))->Write();



  }
  //Writing small angles
  pot->Write("pot");
  n_events->Write("n_events");
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
