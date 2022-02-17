#include "include/GlobalIncludes.h"
#include "include/CCQENuUtils.h"
#include "include/CCQENuPlotUtils.h"

#include "include/GeneralFunc.h"
#include "include/CommonBins.h"


#include "PlotUtils/HyperDimLinearizer.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TFitter.h"
#include "TMatrixD.h"
using namespace CCQENU_ANA;
using namespace std;


vector<int> categories_to_fit({ kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion} );
vector<string> categories_to_fit_names({"qelike_qe_oth", "qelike_res","qelike_2p2h","qelikenot_scp", "qelikenot_snp"});
vector<string> categories_to_fit_title_names({"QELike && QE OTH", "QELike && RES", "QELike && 2p2h", "QELikeNot && Single Charged Pion", "QELikeNot && Single Neutral Pion"});

vector<string> region_names({"h_q2qe_angle_03","h_q2qe_angle_04","h_q2qe_angle_05", "h_q2qe_angle_06","h_q2qe_angle_07","h_q2qe_angle_08","h_q2qe_angle_09","h_q2qe_angle_10"});

vector<int> used_categories({kData, kMC, kQELike_QE_H});

vector<TH1D> input_histos(nHistos) ;
vector<TVectorD> input_A; //Ax=N-B
vector<TVectorD> input_B;
vector<TVectorD> input_N;

//================Define function names=========
int FitEventRate ( vector<string> filenames ); //signal, signal angles and weights
vector<MnvH2D*> GetFitsHisto( TFile* f, string xvar );
vector<MnvH1D*> GetFitsHisto( TFile* f );
template<class MnvHND> void DoFit( vector<MnvHND*> &h, vector<MnvHND*>&fits );

double normFunction( double *evtrates );
void minuitFunction( int &nDim, double* gout, double &result, double par[], int flg);

MnvH1D*  GetSignalRate( vector<MnvH1D*>&full, vector<vector<MnvH1D*>>&regions );
void SetInput( TH1D &full, vector<vector<TH1D>> &regions );
void PrintInput();
TH1D MinimizeNorm(MnvH1D* h);
void setScaleAndError(TH1D*h, TH1D &scale);




//=================Fitting Functions====================

MnvH1D*  GetSignalRate( vector<MnvH1D*>&full, vector<vector<MnvH1D*>>&regions )
{
  MnvH1D* ret = (MnvH1D*) full[kQELike_QE_H]->Clone("fit_signal");
  ret->ClearAllErrorBands();
  ret->Reset();

  //CV
  TH1D h1d_full = full[kQELike_QE_H]->GetCVHistoWithStatError();
  vector<vector<TH1D>> h1d_regions( regions.size(), vector<TH1D>(nHistos) );
  for( uint i = 0; i< regions.size(); i++)
  {
    for( auto c: used_categories ) h1d_regions[i][c]=regions[i][c]->GetCVHistoWithStatError();
  }
  SetInput( h1d_full, h1d_regions );
  PrintInput();
  TH1D cv = MinimizeNorm( full[kQELike_QE_H] );
  setScaleAndError( (TH1D*) ret, cv );
  cout<<"Scaled CV"<<endl;

  //return ret;
  //Vert Err Band
  cout << "Fitting the Vertical errors" << endl;
  vector<string> vertNames = full[kQELike_QE_H]->GetVertErrorBandNames();
  for( uint k = 0; k< vertNames.size(); k++ )
  {
    cout<<"Working on "<<vertNames[k] <<endl;
    int nUniverses = full[kQELike_QE_H]->GetVertErrorBand( vertNames[k] )->GetNHists();
    ret->AddVertErrorBand( vertNames[k],nUniverses);
    for( int iUniv = 0; iUniv< nUniverses; iUniv++ )
    {
      TH1D* errband = ret->GetVertErrorBand( vertNames[k] )->GetHist(iUniv);
      TH1D h1d_full_errband(*full[kQELike_QE_H]->GetVertErrorBand( vertNames[k] )->GetHist(iUniv) );
      vector<vector<TH1D>> h1d_regions_errband( regions.size(), vector<TH1D>(nHistos) );
      for( uint i = 0; i< regions.size(); i++ )
      {
        for( auto c:used_categories ) h1d_regions_errband[i][c]=TH1D( *regions[i][c]->GetVertErrorBand( vertNames[k] )->GetHist(iUniv) );
      }
      SetInput( h1d_full_errband, h1d_regions_errband );
      TH1D eb = MinimizeNorm( full[kQELike_QE_H] );
      setScaleAndError( errband, eb );
    }
  }

   //Lat Err Band
  cout << "Fitting the Lat errors" << endl;
  vector<string> latNames = full[kQELike_QE_H]->GetLatErrorBandNames();
  for( uint k = 0; k< latNames.size(); k++ )
  {
    cout<<"Working on "<<latNames[k] <<endl;
    int nUniverses = full[kQELike_QE_H]->GetLatErrorBand( latNames[k] )->GetNHists();
    ret->AddLatErrorBand( latNames[k],nUniverses);
    for( int iUniv = 0; iUniv< nUniverses; iUniv++ )
    {
      TH1D* errband = ret->GetLatErrorBand( latNames[k] )->GetHist(iUniv);
      TH1D h1d_full_errband(*full[kQELike_QE_H]->GetLatErrorBand( latNames[k] )->GetHist(iUniv) );
      vector<vector<TH1D>> h1d_regions_errband( regions.size(), vector<TH1D>(nHistos) );
      for( uint i = 0; i< regions.size(); i++ )
      {
        for( auto c:used_categories ) h1d_regions_errband[i][c]=TH1D( *regions[i][c]->GetLatErrorBand( latNames[k] )->GetHist(iUniv) );
      }
      SetInput( h1d_full_errband, h1d_regions_errband );
      TH1D eb = MinimizeNorm( full[kQELike_QE_H] );
      setScaleAndError( errband, eb );
    }
  }
  return ret;
}
void setScaleAndError(TH1D*h, TH1D &scale)
{
  for( int i = 0; i< h->GetNbinsX(); i++ )
  {
    h->SetBinContent(i, scale.GetBinContent(i) );
    h->SetBinError(i, scale.GetBinError(i) );
    //cout<<scale.GetBinContent(i)<<", "<<scale.GetBinError(i)<<endl;
  }
}

TH1D MinimizeNorm(MnvH1D* h)
{
  cout<<"MinimizeNorm"<<endl;
  TH1D ret = h->GetCVHistoWithStatError();
  ret.Reset();

  int nbins = ret.GetNbinsX();
  TFitter* minimizer = new TFitter(nbins);
  double print_val=-1;
  minimizer->ExecuteCommand("SET PRINTOUT",&print_val,1);
  minimizer->SetFCN(minuitFunction);

  //cout<<"Setting Parameters"<<endl;

  for(int i=0;i<nbins;i++)
  {
    double v = h->GetBinContent(i+1);
    minimizer->SetParameter(i,Form("scale_%d",i+1), v,0.00001,v*0,v*10 );
  }  
  minimizer->ExecuteCommand("SIMPLEX",0,0);
  minimizer->ExecuteCommand("MIGRAD",0,0);

  for( int i = 0; i<nbins;i ++ )
  {
    cout<<i+1<<", "<<minimizer->GetParameter(i)<<", "<<minimizer->GetParError(i)<<endl;
    ret.SetBinContent(i+1, minimizer->GetParameter(i));
    ret.SetBinError(  i+1, minimizer->GetParError(i) );
  }
  return ret;
}

void SetInput( TH1D &full, vector<vector<TH1D>> &regions )
{
  cout<<"SetInput"<<endl;
  input_A.clear();
  input_B.clear();
  input_N.clear();
  int nBins = full.GetNbinsX();
  int nRegions = regions.size();
  for( int i = 0; i< nBins; i++ )
  {
    TVectorD A(nRegions);
    TVectorD B(nRegions);
    TVectorD N(nRegions);
    for( int j = 0; j< nRegions; j++ )
    {
      N[j] = regions[j][kData].GetBinContent(i+1);
      B[j] = regions[j][kMC].GetBinContent(i+1) - regions[j][kQELike_QE_H].GetBinContent(i+1);
      double f = full.GetBinContent(i+1);
      A[j] = (f>0)? regions[j][kQELike_QE_H].GetBinContent(i+1)/f : 0;
    }
    input_A.push_back(A);
    input_B.push_back(B);
    input_N.push_back(N);
    //cout<<"SetInput: "<<input_B[i](0)<<endl;
  }
}

void PrintInput()
{
  cout<<"PrintInput"<<endl;
  cout<<"nBins:"<<input_A.size()<<endl;
  for(uint i = 0; i<input_A.size(); i++ )
  {
    TVectorD A = input_A[i];
    TVectorD B = input_B[i];
    TVectorD N = input_N[i];
    cout<<"Bin: "<<i<<endl;
    cout<<"A: ("; for(uint j = 0; j<region_names.size();j++) cout<<A(j)<<","; cout<<")"<<endl;
    cout<<"B: ("; for(uint j = 0; j<region_names.size();j++) cout<<B(j)<<","; cout<<")"<<endl;
    cout<<"N: ("; for(uint j = 0; j<region_names.size();j++) cout<<N(j)<<","; cout<<")"<<endl;
  }
}

//================Minimization=====================

double normFunction( double* input_x )
{
  double norm = 0;
  int nBins = input_A.size();
  //cout<<"normFunction====================="<<endl;
  for( int i = 0; i< nBins; i++ )
  {
    TVectorD Ax = input_A[i];//Nx1
    TVectorD B = input_B[i];
    TVectorD N = input_N[i];
    Ax*= input_x[i];
    //cout<<"bin--"<<i+1<<endl;
    //cout<<Ax(0)<<", "<<B(0)<<", "<<N(0)<<", "<<endl;
    //Ax=N-B
    TVectorD delta= Ax + B - N;
    double dnorm = TMath::Sqrt( delta*delta );
    //cout<<input_x[i]<<", "<<dnorm<<endl;
    norm+=dnorm;
  }
  return norm;
}

void minuitFunction( int &nDim, double* gout, double &result, double par[], int flg)
{
  result = normFunction( par );
}


//================Apply Fit Functions===============
vector<MnvH2D*> GetFitsHisto( TFile* f, string xvar )
{
  vector<MnvH2D*> ret;
  for( unsigned int i = 0; i<categories_to_fit_names.size();i++)
  {
    ret.push_back( (MnvH2D*) f->Get( Form( "hs_weights_%s_yvarbins_bgType_%s", xvar.c_str(), categories_to_fit_names[i].c_str() ) ) );
  }
  return ret;
}
vector<MnvH1D*> GetFitsHisto( TFile* f  )
{
  vector<MnvH1D*> ret;
  for( unsigned int i = 0; i<categories_to_fit_names.size();i++)
  {
    ret.push_back( (MnvH1D*) f->Get( Form( "hs_weights_yvar_bgType_%s", categories_to_fit_names[i].c_str() ) ) );
  }
  return ret;
}


template<class MnvHND> 
void DoFit( vector<MnvHND*> &h, vector<MnvHND*>&fits )
{
  for( uint i = 0; i< categories_to_fit.size(); i++ )
  {
    int cat = categories_to_fit[i];
    cout<<"apply fit: "<<names[cat]<<endl;
    cout<<"weight: "<<fits[i]->GetName()<<endl;
    h[cat]->Multiply( h[cat],fits[i] );
  }

  //h[kQELikeNot_SingleNeutralPion]->Multiply(h[kQELikeNot_SingleNeutralPion], fits[2]);
  //h[kQELikeNot_SingleChargedPion]->Multiply(h[kQELikeNot_SingleChargedPion], fits[2]);
 

  h[kMC]->Reset();
  h[kQELike]->Reset();
  h[kQELikeNot]->Reset();


  h[kQELike_QE]->Reset();
  h[kQELike_QE]->Add(h[kQELike_QE_OTH]);
  h[kQELike_QE]->Add(h[kQELike_QE_H]);
 
  vector<int> qelike({ kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH });
  for( auto cat: qelike ) h[kQELike]->Add( h[cat] );
  vector<int> qelikenot({ kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } );
  for( auto cat: qelikenot ) h[kQELikeNot]->Add( h[cat] );

  h[kMC]->Add( h[kQELike] );
  h[kMC]->Add( h[kQELikeNot] );
}
template void DoFit( vector<MnvH1D*> &h, vector<MnvH1D*>&fits );
template void DoFit( vector<MnvH2D*> &h, vector<MnvH2D*>&fits );



//===============End Helper Functions============
//============================================================================================


int FitEventRate( vector<string> filenames )
{

  //Access Files
  cout << "Loading file: " << filenames[1] << endl;
  TFile* f_signal = new TFile( filenames[1].c_str() , "READ" );
  TFile* f_signal_angles = new TFile( filenames[2].c_str() , "READ" );
  TFile* f_fits = new TFile( filenames[3].c_str(), "READ" );
  

  //Define utils class
  bool fluxHistoExists = true;
  CCQENuPlotUtils *plotutils    = new CCQENuPlotUtils( fluxHistoExists );
  //CCQENuCuts      *cutter       = new CCQENuCuts();
  //MnvPlotter      *plotter      = plotutils->getPlotter();

  //Access histos
  int nRegions = region_names.size();
  vector< vector<MnvH1D*> > hists(nRegions, vector<MnvH1D*>(nHistos) );
  vector<MnvH1D*> h_q2qe_region_00(nHistos);
  plotutils->bookHistos(f_signal_angles, &(h_q2qe_region_00[0]), "h_q2qe_region_00");
  for( int r = 0; r<nRegions; r++ ) plotutils->bookHistos( f_signal, &(hists[r][0]), region_names[r].c_str() );
  //Apply fit
  vector<MnvH1D*> fits = GetFitsHisto( f_fits );
  DoFit(h_q2qe_region_00, fits );
  h_q2qe_region_00[kData]->AddMissingErrorBandsAndFillWithCV( *h_q2qe_region_00[kMC] );
  for(uint i = 0; i< hists.size(); i++ ) 
  {
    DoFit( hists[i], fits );
    hists[i][kData]->AddMissingErrorBandsAndFillWithCV(*hists[i][kMC] );
  }

  cout<<"Done Applying Fit"<<endl;


  //bins( matrix )
  //This is a test on CV only
  MnvH1D* signal = GetSignalRate( h_q2qe_region_00, hists );
  signal->SetName("h_q2qe_region_00_evtrate");

  TFile *output = new TFile( filenames[4].c_str(), "RECREATE" );
  output->cd();
  signal->Write();
  //output->Close();
  return 0;
}

int main( int argc, char *argv[])
{
  ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if (argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t-./SideBandFit  Name_and_path_to_SignalSelectedEvents_file  Name_and_path_to_BlobSideBandSelectedEvents_file Name_and_path_to_MichelSideBandeSelectedEvents_file \n\n"<<
      "\t-Name_and_path_to_background_templates_file\t =\t Name of and path to MuonSelection selecting the signal events 2 track\n"<<
      "\t-Name_and_path_to_Output_file\t =\t Name of and path to the Output ROOT file that will be created \n"<<
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    return 0; 
  }

  //! Default parameters
  std::vector<std::string> par;
  par.push_back("SideBandFit");
  par.push_back( Form("%s/ana/rootfiles/MuonEvent_Signal.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/MuonEvent_Signal_Angle.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/SignalFitWeightsAngle.root",getenv("CCQENUROOT") ) ); 
  par.push_back( Form("%s/ana/rootfiles/MuonEvent_Signal_Fitted.root",getenv("CCQENUROOT") ) );//output file

  const int nArgs = par.size();

  //! Set user parameters
  for( int i=0; i<nArgs; ++i){
    par.at(i) = argv[i];
  }
  for( unsigned int i=0; i<par.size(); ++i)
  std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;
  
  return FitEventRate(par);
}

