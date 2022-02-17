#include "include/CCQENuPlotUtils.h"
#include "include/CCQENuUtils.h"
#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/TargetUtils.h"
#include "include/GeneralFunc.h"


using namespace CCQENU_ANA;


//void SubtractBackground(MnvH2D** h, MnvH2D* h_bg_weights, MnvH2D* &h_data_nobck, MnvH2D* &h_mc_nobck, MnvH2D* &h_mc_tuned_bck, MnvH2D* &h_data_predicted_bck, std::string label);
//void SubtractBackground(MnvH2D** h, MnvH2D* h_bg_weights, MnvH2D* &h_data_nobck, MnvH2D* &h_mc_nobck);
//void SubtractBackground(MnvH1D** hists, MnvH1D* &h_data_nobck, MnvH1D* &h_data_predicted_bck, MnvH1D* &h_mc_nobck, MnvH1D* &h_mc_tuned_bck , std::string label);
template<class MnvHnD> void SubtractBackground(MnvHnD** hists, MnvHnD* &h_data_nobck, MnvHnD* &h_data_predicted_bck, MnvHnD* &h_mc_nobck, MnvHnD* &h_mc_tuned_bck , std::string label);

template<class T>
void CorrectByEfficiency( T *&h_effcor, T* h_unfolded, T** h_effhist ); 
template<class T>
void CorrectByEfficiency( T *&h_effcor, T* h_unfolded, T** h_effhist_num, T** h_effhist_dem ); 

void NormalizeByFluxAndTargets( MnvH1D *&h_normalized, MnvH1D *h_effcor, MnvH1D* h_flux, double pot_scale, bool applyFluxConstraint, CCQENuUtils* utils,MnvH1D *flux_norm= NULL, bool isMC=false );

//void NormalizeByFluxAndTargets( MnvH2D *&h_normalized, MnvH2D *h_effcor, MnvH1D* h_flux, double pot_scale, bool applyFluxConstraint, CCQENuUtils* utils, TMatrixD covUnfold,MnvH2D *flux_norm= NULL );
void NormalizeByFluxAndTargets( MnvH2D *&h_normalized, MnvH2D *h_effcor, MnvH1D* h_flux, double pot_scale, bool applyFluxConstraint, CCQENuUtils* utils,MnvH2D *flux_norm= NULL, bool isMC=false );

MnvH2D* GetFluxNorm( MnvH1D* flux, MnvH2D* hist );

void RebinFluxHist(TH1D* h_flux, TH1D* &h_rebinned_flux);

bool AddNonReweightableGENIESys(MnvH2D *&xs);

void ZeroUnreported(MnvH2D *&xs);
void ZeroUnreported(TMatrixD &matd);

double getBinArea(TH2* h, int globalBin)
{
  int binx, biny, binz;
  h->GetBinXYZ(globalBin, binx, biny, binz);
  return h->GetXaxis()->GetBinWidth(binx)*h->GetYaxis()->GetBinWidth(biny);
}
double getBinArea(TH1* h, int globalBin)
{
  return h->GetXaxis()->GetBinWidth(globalBin);
}

template<class T> 
TMatrixD divideCovByHist(TMatrixD& m, T* h)
{
  TMatrixD ret(m);
  for(int i=0; i<h->fN; ++i){
    for(int j=0; j<h->fN; ++j){
      ret(i,j)=m(i,j)/(h->fArray[i]*h->fArray[j]);
    }
  }
  return ret;
}

template TMatrixD divideCovByHist(TMatrixD& m, TH2D* h);
template TMatrixD divideCovByHist(TMatrixD& m, TH1D* h);

TMatrixD divideCovByHists(TMatrixD m, TH2D* num, TH2D* dem)
{ 

  TH2D *tmp = new TH2D(*num);
  tmp->Divide(dem);
  TMatrixD ret(m);
  cout << "Dividing cov by hists" << num->fN << endl;
  for(int i=0; i<num->fN; ++i){
    for(int j=0; j<num->fN; ++j){
      int x,y,z,usex1,usey1,usex2,usey2;
      num->GetBinXYZ(i,x,y,z);
      usex1=x;
      usey1=y;
      num->GetBinXYZ(j,x,y,z);
      usex2=x;
      usey2=y;
      double eff1 = tmp->GetBinContent(usex1,usey1);
      double eff2 = tmp->GetBinContent(usex2,usey2);
      double val = eff1*eff2;
      if(val!=0)ret(i,j)=m(i,j)/(eff1*eff2);
    }
  }
  return ret;
}

TMatrixD divideCovByBinWidth(TMatrixD& m, TH2D* h)
{
  TMatrixD ret(m);
  for(int i=0; i<h->fN; ++i){
    for(int j=0; j<h->fN; ++j){
      ret(i,j)=m(i,j)/(getBinArea(h, i)*getBinArea(h, j));
    }
  }
  return ret;
}
TMatrixD divideCovByBinWidth(TMatrixD& m, TH1D* h)
{
  TMatrixD ret(m);
  for(int i=0; i<h->fN; ++i){
    for(int j=0; j<h->fN; ++j){
      ret(i,j)=m(i,j)/(getBinArea(h, i)*getBinArea(h, j));
    }
  }
  return ret;
}

//need to clean out result for unreportable bins
void SetBinZero(MnvH2D *&xs, int x, int y){

  xs->SetBinContent(x,y,0);
  vector<string> verterrnames = xs->GetVertErrorBandNames();
  vector<string> laterrnames = xs->GetLatErrorBandNames();
  //Vertical Errors
  for(UInt_t i=0;i<verterrnames.size();i++){
    MnvVertErrorBand2D *tmperr = xs->GetVertErrorBand(verterrnames[i]);
    int nhist = tmperr->GetNHists();
    for(int j=0;j<nhist;j++){
      tmperr->GetHist(j)->SetBinContent(x,y,0);
    }
  }
  //Lateral Errors
  for(UInt_t i=0;i<laterrnames.size();i++){
    MnvLatErrorBand2D *tmperr = xs->GetLatErrorBand(laterrnames[i]);
    int nhist = tmperr->GetNHists();
    for(int j=0;j<nhist;j++){
      tmperr->GetHist(j)->SetBinContent(x,y,0);
    }
  }
}

void ZeroUnreported(MnvH2D *&xs){
  int nBinsX = xs->GetNbinsX()+2;//Include over/under flows
  int nBinsY = xs->GetNbinsY()+2;//Include over/under flows
  for(int i=0;i<nBinsX;i++){
    for(int j=0;j<nBinsY;j++){
  
      if(i==0 || j==0 || i==nBinsX-1 || j==nBinsY-1) SetBinZero(xs,i,j);
      if(i==1 && j>=10) SetBinZero(xs,i,j);
      else if(i==1 && j>=10) SetBinZero(xs,i,j);
      else if(i==2 && j>=11) SetBinZero(xs,i,j);
      else if(i==3 && j>=12) SetBinZero(xs,i,j);
      else if(i==4 && j>=12) SetBinZero(xs,i,j);
      else if(i==5 && j>=13) SetBinZero(xs,i,j);
      else if(i==6 && j>=14) SetBinZero(xs,i,j);
    }
  }
}

void ZeroUnreported(TMatrixD &matd,MnvH2D*xs){
  //Zero the stupid stat covariance matrix we have to cart around
  //mapping of cov matrix to x,y bin
  //Big bins = y axis bins
  //small bins = x axis bins
  int nBinsX = xs->GetNbinsX()+2;//Include over/under flows
  int nBinsY = xs->GetNbinsY()+2;//Include over/under flows
  //20 pt and 30 pz. So bin 35 is pt bin 1 and the 35-30 or bin 5 of pz
  // 35/30 = 1 (with 5 left over) 35-1*30 = 5. 0 to 29 is the first and 30-49 is the second so you subtract 1 off the little bin
  int covMatrixN = nBinsX*nBinsY;
  for(int i=0;i<covMatrixN;i++){
    int bigBin_i = i/nBinsX;
    int littleBin_i = i%nBinsX;
    for(int j=0;j<covMatrixN;j++){
      int bigBin_j = j/nBinsX;
      int littleBin_j = j%nBinsX;
      
      if(littleBin_i==0 || bigBin_i==0 || littleBin_i==nBinsX-1 || bigBin_i==nBinsY-1)matd[i][j]=0;

      if(littleBin_i==1 && bigBin_i>=10) matd[i][j]=0;
      else if(littleBin_i==1 && bigBin_i>=10) matd[i][j]=0;
      else if(littleBin_i==2 && bigBin_i>=11) matd[i][j]=0;
      else if(littleBin_i==3 && bigBin_i>=12) matd[i][j]=0;
      else if(littleBin_i==4 && bigBin_i>=12) matd[i][j]=0;
      else if(littleBin_i==5 && bigBin_i>=13) matd[i][j]=0;
      else if(littleBin_i==6 && bigBin_i>=14) matd[i][j]=0;
      
      if(littleBin_j==0 || bigBin_j==0 || littleBin_j==nBinsX-1 || bigBin_j==nBinsY-1)matd[i][j]=0;

      if(littleBin_j==1 && bigBin_j>=10) matd[i][j]=0;
      else if(littleBin_j==1 && bigBin_j>=10) matd[i][j]=0;
      else if(littleBin_j==2 && bigBin_j>=11) matd[i][j]=0;
      else if(littleBin_j==3 && bigBin_j>=12) matd[i][j]=0;
      else if(littleBin_j==4 && bigBin_j>=12) matd[i][j]=0;
      else if(littleBin_j==5 && bigBin_j>=13) matd[i][j]=0;
      else if(littleBin_j==6 && bigBin_j>=14) matd[i][j]=0;
    }
  }
}

//============Unfolding Functions from TransWarpExtraction ==============
TMatrixD UnfoldDummy( MnvH1D* h_data_nobck, MnvH1D* h_data_unfolded, MnvH1D* h_reco, MnvH1D* h_truth, MnvH2D* h_migration, int num_iter ) 
{
  MinervaUnfold::MnvUnfold unfold;
  cout << "Getting the covariance of the unfolding" << endl;
  TMatrixD unfoldingCovMatrixOrig;
  int correctNbins;
  int matrixRows;  

  TH1D* hUnfoldedDummy  = new TH1D(h_data_unfolded->GetCVHistoWithStatError());
  TH1D* hRecoDummy      = new TH1D(h_reco->GetCVHistoWithStatError());
  TH1D* hTruthDummy     = new TH1D(h_truth->GetCVHistoWithStatError());
  TH1D* hBGSubDataDummy = new TH1D(h_data_nobck->GetCVHistoWithStatError());
  TH2D* hMigrationDummy = new TH2D(h_migration->GetCVHistoWithStatError());
  unfold.UnfoldHisto(hUnfoldedDummy, unfoldingCovMatrixOrig, hMigrationDummy, hRecoDummy, hTruthDummy, hBGSubDataDummy,RooUnfold::kBayes, num_iter);//Stupid RooUnfold.  This is dummy, we don't need iterations

  correctNbins=hUnfoldedDummy->fN;
  matrixRows=unfoldingCovMatrixOrig.GetNrows();

  if(correctNbins!=matrixRows){
    cout << "****************************************************************************" << endl;
    cout << "*  Fixing unfolding matrix size because of RooUnfold bug. From " << matrixRows << " to " << correctNbins << endl;
    cout << "****************************************************************************" << endl;
    // It looks like this, since the extra last two bins don't have any content
    unfoldingCovMatrixOrig.ResizeTo(correctNbins, correctNbins);
  }

  for(int i=0; i<unfoldingCovMatrixOrig.GetNrows(); ++i) unfoldingCovMatrixOrig(i,i)=0;
  delete hUnfoldedDummy;
  delete hMigrationDummy;
  delete hRecoDummy;
  delete hTruthDummy;
  delete hBGSubDataDummy;

  return unfoldingCovMatrixOrig;
}

TMatrixD UnfoldDummy( MnvH2D* h_data_nobck, MnvH2D* h_data_unfolded, MnvH2D* h_reco, MnvH2D* h_truth, MnvH2D* h_migration, int num_iter ) 
{

  MinervaUnfold::MnvUnfold unfold;
  cout << "Getting the covariance of the unfolding" << endl;
  TMatrixD unfoldingCovMatrixOrig;
  int correctNbins;
  int matrixRows;  

  TH2D* hUnfoldedDummy  = new TH2D(h_data_unfolded->GetCVHistoWithStatError());
  TH2D* hRecoDummy      = new TH2D(h_reco->GetCVHistoWithStatError());
  TH2D* hTruthDummy     = new TH2D(h_truth->GetCVHistoWithStatError());
  TH2D* hBGSubDataDummy = new TH2D(h_data_nobck->GetCVHistoWithStatError());
  TH2D* hMigrationDummy = new TH2D(h_migration->GetCVHistoWithStatError());
  unfold.UnfoldHisto2D(hUnfoldedDummy, unfoldingCovMatrixOrig, hMigrationDummy, hRecoDummy, hTruthDummy, hBGSubDataDummy, num_iter);

  correctNbins=hUnfoldedDummy->fN;
  matrixRows=unfoldingCovMatrixOrig.GetNrows();

  if(correctNbins!=matrixRows){
    cout << "****************************************************************************" << endl;
    cout << "*  Fixing unfolding matrix size because of RooUnfold bug. From " << matrixRows << " to " << correctNbins << endl;
    cout << "****************************************************************************" << endl;
    // It looks like this, since the extra last two bins don't have any content
    unfoldingCovMatrixOrig.ResizeTo(correctNbins, correctNbins);
  }

  for(int i=0; i<unfoldingCovMatrixOrig.GetNrows(); ++i) unfoldingCovMatrixOrig(i,i)=0;
  delete hUnfoldedDummy;
  delete hMigrationDummy;
  delete hRecoDummy;
  delete hTruthDummy;
  delete hBGSubDataDummy;

  return unfoldingCovMatrixOrig;
}


bool UnfoldData( MnvH1D* &h_data_unfolded, MnvH1D* h_data_nobck, MnvH1D* h_reco, MnvH1D* h_truth, MnvH2D* h_migration, int num_iter ) 
{
  MinervaUnfold::MnvUnfold unfold;
  bool data_unfolded = false;
  TMatrixD dummyCovMatrix;
  data_unfolded = unfold.UnfoldHisto( h_data_unfolded, dummyCovMatrix, h_migration, h_data_nobck, RooUnfold::kBayes, num_iter, true, true );
  
  if( !data_unfolded )
  {
    std::cout << "Unfolding failed for either data or MC. Please check " << std::endl;
    return false;  
  }

  TMatrixD unfoldingCovMatrixOrig = UnfoldDummy( h_data_nobck, h_data_unfolded, h_reco, h_truth, h_migration, num_iter ); 
  
  h_data_unfolded->PushCovMatrix("unfoldingCov",unfoldingCovMatrixOrig);
 
  return true;
}

bool UnfoldData( MnvH2D* &h_data_unfolded, MnvH2D* h_data_nobck, MnvH2D* h_reco, MnvH2D* h_truth, MnvH2D* h_migration, int num_iter ) 
{
  MinervaUnfold::MnvUnfold unfold;
  bool data_unfolded = false;
  data_unfolded = unfold.UnfoldHisto2D(h_data_unfolded, h_migration, h_reco, h_truth, h_data_nobck, num_iter, true, true );
  
  if( !data_unfolded )
  {
    std::cout << "Unfolding failed for either data or MC. Please check " << std::endl;
    return false;  
  }

  TMatrixD unfoldingCovMatrixOrig = UnfoldDummy( h_data_nobck, h_data_unfolded, h_reco, h_truth, h_migration, num_iter ); 
  h_data_unfolded->PushCovMatrix("unfoldingCov",unfoldingCovMatrixOrig);

  return true;
}

//============Fitting Utils==============
vector<int> fittedCat({ kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion } );
vector<string> fittedCatNames({"qelike_qe_oth", "qelike_res", "qelike_2p2h", "qelikenot_snp", "qelikenot_scp"});
vector<int> qelikeCat({ kQELike_QE_H, kQELike_QE_OTH, kQELike_2p2h, kQELike_RES, kQELike_DIS, kQELike_OTH });
vector<int> qelikenotCat({ kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } );

vector<MnvH1D*> GetFits( TFile* f )
{
  vector<MnvH1D*> ret(fittedCat.size(),NULL);
  for( uint i = 0; i< fittedCat.size(); i++ ) ret[i]=(MnvH1D*) f->Get( Form( "hs_weights_yvar_bgType_%s", fittedCatNames[i].c_str() ) ) ;
  //ret[0] = (MnvH1D*) f->Get("hs_weights_yvar_bgType_qe_oth");
  //ret[1] = (MnvH1D*) f->Get("hs_weights_yvar_bgType_2p2h");
  //ret[2] = (MnvH1D*) f->Get("hs_weights_yvar_bgType_qelikenot_scp");
  //ret[3] = (MnvH1D*) f->Get("hs_weights_yvar_bgType_qelikenot_snp");
  //ret[4] = (MnvH1D*) f->Get("hs_weights_yvar_bgType_qelikenot_mp");
  return ret;
}

vector<MnvH2D*> GetFits( TFile* f, string xvar )
{
  vector<MnvH2D*> ret;
  for( unsigned int i = 0; i<fittedCat.size();i++)
  {
    string hname = Form( "hs_weights_%s_yvarbins_bgType_%s", xvar.c_str(), fittedCatNames[i].c_str() ) ; 
    cout<<"Getting "<<hname<<endl;
    ret.push_back( (MnvH2D*) f->Get( hname.c_str() ) );
  }
  return ret;
}


template<class T>
void ApplyFit( T** hists, vector<T*>&fits )
{
  for( UInt_t i = 0; i< fittedCat.size(); i++ )
  {
    int cat = fittedCat[i];
    hists[cat]->Multiply( hists[cat], fits[i] );
  }
  hists[kMC]->Reset();
  hists[kQELike]->Reset();
  hists[kQELikeNot]->Reset();
  for( auto cat: qelikeCat ) hists[kQELike]->Add( hists[cat] );
  for( auto cat: qelikenotCat ) hists[kQELikeNot]->Add( hists[cat] );
  hists[kMC]->Add( hists[kQELike] );
  hists[kMC]->Add( hists[kQELikeNot] );
}
template void ApplyFit( MnvH1D** hists, vector<MnvH1D*>&fits );
template void ApplyFit( MnvH2D** hists, vector<MnvH2D*>&fits );



int CrossSectionHists( string output_filename, string filename_signal_regions, string filename_bg_weights, string filename_migration_histogram, string filename_eff, int num_iter, bool makeFluxConstraintHisto, bool applyFluxConstraint,int run_type,bool doenu)
{

  cout <<"*********************************" << endl;
  cout << "Starting" << endl;
  cout << output_filename << endl;
  cout <<"*********************************" << endl;

  // read file to get histograms
  TFile *f_signal_regions = new TFile( filename_signal_regions.c_str(), "READ" );
  if (f_signal_regions->IsZombie() || f_signal_regions->GetListOfKeys()->IsEmpty()){
    Error("CrossSectionPlots","Could not get histogram ROOT file or it was empty.");//let's say 20 pt bins 
    return 1;
  }

  TFile *f_bg_weights = new TFile( filename_bg_weights.c_str(), "READ" );
  if (f_bg_weights->IsZombie() || f_bg_weights->GetListOfKeys()->IsEmpty()){
    Error("CrossSectionPlots","Could not get histogram ROOT file or it was empty.");
    return 1;
  }
  TFile *f_migration = new TFile( filename_migration_histogram.c_str(), "READ" );
  if (f_migration->IsZombie() || f_migration->GetListOfKeys()->IsEmpty()){
    Error("CrossSectionPlots","Could not get histogram ROOT file or it was empty.");
    return 1;
  }
  cout << "Loading Effficiency file: " << filename_eff.c_str() << endl;
  TFile *f_efficiency = new TFile( filename_eff.c_str(), "READ" );
  if (f_efficiency->IsZombie() || f_efficiency->GetListOfKeys()->IsEmpty()){
    Error("CrossSectionPlots","Could not get histogram ROOT file or it was empty.");
    return 1;
  }

  CCQENuUtils *utils = new CCQENuUtils( false, makeFluxConstraintHisto );
  utils->setPlaylist("minervame6D");
  utils->setFluxReweighterPlaylist();
  CCQENuPlotUtils *plotUtils = new CCQENuPlotUtils();
  MinervaUnfold::MnvUnfold unfold;

  //-------------------------------------------
  // Grab event counts before cuts for norm
  //-------------------------------------------
  double data_events,mc_events;
  TVector2 *evt = new TVector2();
  evt = (TVector2*)f_signal_regions->Get("n_events");
  data_events = evt->X();
  mc_events   = evt->Y();
  
  cout<< "Number of Data Events (2-track) = " << data_events << endl;
  cout<< "Number of MC Events (2-track) = " <<   mc_events << endl;

  //--------------------------------------------
  // Load histos to plot 
  //--------------------------------------------
  //1-D
  string sname = "h_q2qe_region_00";
  string effpurMig_name = "h_q2qe_angle_10";


  MnvH1D *h_q2qe_signal[nHistos];
  MnvH1D *h_effhist_num[nHistos], *h_effhist_dem[nHistos];
  MnvH2D *h_migration;
  MnvH1D *h_reco, *h_truth;
  //MnvH1D *h_q2qe_angle_01[nHistos];
  //MnvH1D *h_q2qe_angle_02[nHistos];
  //MnvH1D *h_q2qe_angle_03[nHistos];
  //MnvH1D *h_q2qe_angle_04[nHistos];
  cout << "Booking Event Rates" << endl;
  plotUtils->bookHistos( f_signal_regions, h_q2qe_signal, sname.c_str() );

  cout<<" Get Fits "<<endl;
  vector<MnvH1D*> fits = GetFits( f_bg_weights );
  cout << "Apply Fits 1D"<<endl;
  ApplyFit( h_q2qe_signal, fits );

  cout << "Booking Efficiency " << endl;
  plotUtils->bookHistos( f_efficiency, h_effhist_num,effpurMig_name.c_str() );
  plotUtils->bookHistos( f_efficiency, h_effhist_dem,Form("h_truth_q2qe")  );
  for( UInt_t i = 0; i< nHistos; i++ )
  {
    h_effhist_dem[i]->AddMissingErrorBandsAndFillWithCV( *h_q2qe_signal[i] );
    h_effhist_num[i]->AddMissingErrorBandsAndFillWithCV( *h_q2qe_signal[i] );
  }

  cout << "Booking Migration" << endl;
  h_migration = (MnvH2D*) f_migration->Get(Form("%s_migration_qelike_qe_h", (effpurMig_name).c_str() ) );
  h_reco = (MnvH1D*) f_migration->Get(Form("%s_qelike_qe_h", (effpurMig_name+"_reco").c_str() ) );
  h_truth = (MnvH1D*) f_migration->Get(Form("%s_qelike_qe_h", (effpurMig_name+"_truth").c_str() ) );

  MnvH2D* tmp = (MnvH2D*) h_migration->Clone("tmp");
  tmp->Reset();
  tmp->AddMissingErrorBandsAndFillWithCV( *h_q2qe_signal[kMC] );

  cout << h_migration->GetVertErrorBandNames().size() << endl;

  cout << "DONE-q2" << endl;

  //2-D
  string sname2 = "h_enu_q2qe";
  string effpurMig_name2 = "h_enu_q2qe";
  string xvar = "enu";


  MnvH2D *h_enu_q2qe_signal[nHistos];
  MnvH2D *h_enu_effhist_num[nHistos], *h_enu_effhist_dem[nHistos];
  MnvH2D *h_enu_migration;
  MnvH2D *h_enu_reco, *h_enu_truth;
  cout << "Booking Event Rates" << endl;
  plotUtils->bookHistos( f_signal_regions, h_enu_q2qe_signal, sname2.c_str() );

  cout<<" Get Fits "<<endl;
  vector<MnvH2D*> fits_enu = GetFits( f_bg_weights, "muonmomentum" );
  cout<<fits_enu[0]<<endl;
  cout << "Apply Fits 2D"<<endl;
  ApplyFit( h_enu_q2qe_signal, fits_enu );

  cout << "Booking Efficiency " << endl;
  plotUtils->bookHistos( f_efficiency, h_enu_effhist_num,effpurMig_name2.c_str() );
  plotUtils->bookHistos( f_efficiency, h_enu_effhist_dem,Form("h_truth_enu_q2qe")  );
  for( UInt_t i = 0; i< nHistos; i++ )
  {
    h_enu_effhist_dem[i]->AddMissingErrorBandsAndFillWithCV( *h_enu_q2qe_signal[i] );
    h_enu_effhist_num[i]->AddMissingErrorBandsAndFillWithCV( *h_enu_q2qe_signal[i] );
  }

  cout << "Booking Migration" << endl;
  h_enu_migration = (MnvH2D*) f_migration->Get(Form("%s_qelike_qe_h_migration", (effpurMig_name2).c_str() ) );
  h_enu_reco = (MnvH2D*) f_migration->Get(Form("%s_qelike_qe_h_reco", (effpurMig_name2).c_str() ) );
  h_enu_truth = (MnvH2D*) f_migration->Get(Form("%s_qelike_qe_h_truth", (effpurMig_name2).c_str() ) );


  cout << "DONE-enu-q2" << endl;





  //Getting the Flux
  
  //New method for flux
  MnvH1D* h_flux = utils->getReweightedFluxHist();
  h_flux->AddMissingErrorBandsAndFillWithCV( *h_q2qe_signal[kMC] );
  //MnvH1D* h_flux = utils->getUnweightedFluxHist();

  std::vector<std::string> fluxErrorBandNames = h_flux->GetVertErrorBandNames(); 
  for( std::vector<std::string>::iterator itName = fluxErrorBandNames.begin(); itName != fluxErrorBandNames.end(); ++itName ) { 
    cout << "Flux error band name is " << *itName << endl; 
  }

  
  //--------------------------------------------
  // Get POT Normalization factor:
  // This is necessary because Background Subtraction
  // uses POT Normalized MC
  //--------------------------------------------
  
  double pot_data = plotUtils->getPOTData( f_signal_regions );
  double pot_mc = plotUtils->getPOTMC( f_signal_regions );
  double pot_scale = plotUtils->getPOTNormFactor( f_signal_regions );
  
  if(pot_scale==0) pot_scale = 1.0;
  cout<< "---------------------------------------"  << endl;
  cout<< "POT INFORMATION:"  << endl;
  cout<< "POT Data = " << pot_data << endl;
  cout<< "POT MC   = " << pot_mc << endl;
  cout<< "POT Scale factor = " << pot_scale << endl;
  cout<< "---------------------------------------"  << endl;

  //2-D
  cout << "Scale Rate" << endl;
  plotUtils->scaleMCHistos( h_q2qe_signal, pot_scale);
  plotUtils->scaleMCHistos( h_enu_q2qe_signal, pot_scale);
  cout << "Scale migration" << endl;
  h_reco->Scale(pot_scale);
  h_truth->Scale(pot_scale);
  h_migration->Scale(pot_scale);

  h_enu_reco->Scale(pot_scale);
  h_enu_truth->Scale(pot_scale);
  h_enu_migration->Scale(pot_scale);


  //--------------------------------------------
  // Background Subtraction
  //--------------------------------------------
  cout<< "---------------------------------------"  << endl;
  cout<< "Background Subtraction"  << endl;
  cout<< "---------------------------------------"  << endl;

  MnvH1D *h_data_nobck=NULL,*h_data_predicted_bck=NULL,*h_mc_nobck=NULL,*h_mc_tuned_bck=NULL;
  MnvH2D *h_enu_data_nobck=NULL,*h_enu_data_predicted_bck=NULL,*h_enu_mc_nobck=NULL,*h_enu_mc_tuned_bck=NULL;

  cout<<"Background Subtraction of Hydrogen Sub-Sample"<<endl;
  SubtractBackground( h_q2qe_signal, h_data_nobck, h_data_predicted_bck, h_mc_nobck, h_mc_tuned_bck, "hydrogen" );
  SubtractBackground( h_enu_q2qe_signal, h_enu_data_nobck, h_enu_data_predicted_bck, h_enu_mc_nobck, h_enu_mc_tuned_bck, "hydrogen" );
  cout<<h_q2qe_signal[kData]<<", "<<h_data_nobck<<", "<<h_data_predicted_bck<<", "<<h_mc_nobck<<", "<<h_mc_tuned_bck<<endl;
  
  //--------------------------------------------
  // Unfold DATA and MC merged samples
  //--------------------------------------------
  MnvH1D *h_data_unfolded = 0, *h_mc_unfolded = 0;
  TMatrixD unfoldingCovMatrixData;
  TMatrixD unfoldingCovMatrixMC;

  MnvH2D *h_enu_data_unfolded = 0, *h_enu_mc_unfolded = 0;
  cout<< "---------------------------------------"  << endl;
  cout<< "Now Unfolding with " << num_iter << " iterations" << endl;
  cout<< "---------------------------------------"  << endl;
  if (num_iter == 0) {
    //------------------------------------------------------------------
    // Special case code where there are no requested iterations
    // Copies the hist then changes the name to end in _unfold
    //------------------------------------------------------------------
    // DATA
    h_data_unfolded = new MnvH1D(*h_data_nobck);
    h_data_unfolded->SetName(Form( "%s_unfold", h_data_unfolded->GetName() ));
    // MC
    h_mc_unfolded = new MnvH1D(*h_mc_nobck);
    h_mc_unfolded->SetName(Form( "%s_unfold", h_mc_unfolded->GetName() ));

    //cout<<" Pushing unfolding matrix "<<endl;
    //for( int i = 0; i < unfoldingCovMatrixData.GetNrows(); ++i ) unfoldingCovMatrixData(i,i)=0;
    //h_data_unfolded->PushCovMatrix( "unfoldingCov", unfoldingCovMatrixData );
    

  } else {
    //---------------------------------------------------------------------------------
    // Last input variable (= true) ensures that we are unfolding each universe  
    // of the distribution with its corresponding universe in the migration matrix   
    //---------------------------------------------------------------------------------
    //pz
    cout << "Unfolding" << endl;
    cout<< h_data_unfolded<<", "<< unfoldingCovMatrixData.GetName() <<", "<< h_migration <<", "<<h_data_nobck<<endl;
    cout << "data"<<endl;
    bool data_unfolded= UnfoldData( h_data_unfolded, h_data_nobck,h_reco, h_truth, h_migration, num_iter );
    cout << "mc"<<endl;
    bool mc_unfolded= UnfoldData( h_mc_unfolded, h_mc_nobck,h_reco, h_truth, h_migration, num_iter );
    cout<< "Unfolded" <<endl;
    if(!data_unfolded || !mc_unfolded){
      std::cout << "Unfolding failed for either data or MC. Please check " << std::endl;
    }

    cout << "2D data"<<endl;
    data_unfolded= UnfoldData( h_enu_data_unfolded, h_enu_data_nobck,h_enu_reco, h_enu_truth, h_enu_migration, num_iter );
    cout << "2D mc"<<endl;
    mc_unfolded= UnfoldData( h_enu_mc_unfolded, h_enu_mc_nobck,h_enu_reco, h_enu_truth, h_enu_migration, num_iter );
    cout<< "2D Unfolded" <<endl;
    if(!data_unfolded || !mc_unfolded){
      std::cout << "Unfolding failed for either 2D data or MC. Please check " << std::endl;
    }



    //h_data_unfolded = new MnvH1D(*h_data_nobck);
    //h_data_unfolded->SetName(Form( "%s_unfold", h_data_unfolded->GetName() ));
    //h_mc_unfolded = new MnvH1D(*h_mc_nobck);
    //h_mc_unfolded->SetName(Form( "%s_unfold", h_mc_unfolded->GetName() ));

    //bool data_unfolded = unfold.UnfoldHisto(h_data_unfolded, unfoldingCovMatrixData, h_migration, h_data_nobck, RooUnfold::kBayes,num_iter, true,true);
    //bool mc_unfolded = unfold.UnfoldHisto(h_mc_unfolded, unfoldingCovMatrixMC, h_migration, h_mc_nobck, RooUnfold::kBayes,num_iter, true,true);
    
    //NOW we have to get the damn covariance of the unfolding which means doing the unfolding AGAIN!@!
    //Only need this on the data as the MC is used as a CV only.
    //pz
    //cout << "Getting the covariance of the unfolding for PZ" << endl;

    //TH1D* hUnfoldedDummy=new TH1D(h_data_unfolded->GetCVHistoWithStatError());
    //TH2D* hMigrationDummy=new TH2D(h_migration->GetCVHistoWithStatError());
    //TH1D* hRecoDummy=new TH1D(h_reco->GetCVHistoWithStatError());
    //TH1D* hTruthDummy=new TH1D(h_truth->GetCVHistoWithStatError());
    //TH1D* hBGSubDataDummy=new TH1D(h_data_nobck->GetCVHistoWithStatError());

    //unfold.UnfoldHisto(hUnfoldedDummy, unfoldingCovMatrixData, hMigrationDummy, hRecoDummy, hTruthDummy,hBGSubDataDummy, RooUnfold::kBayes,num_iter);
    // There's a bug in RooUnfold that's making it return covariance  // matrices with two extra bins. Kill them here, with a check.  
    // Conveniently, this bug was being hidden by an offsetting bug in  
    // MnvH2D, which is now fixed
    //int correctNbins=hUnfoldedDummy->fN;
    //cout<<"Correct N bins = "<<correctNbins<<endl;
    //int matrixRows=unfoldingCovMatrixData.GetNrows();
    //if(correctNbins!=matrixRows){
    //  cout << "****************************************************************************" << endl;
    //  cout << "*  Fixing unfolding matrix size because of RooUnfold bug. From " << matrixRows << " to " << correctNbins << endl;
    //  cout << "****************************************************************************" << endl;
    //  // It looks like this DTRT, since the extra last two bins don't have any content
    //  unfoldingCovMatrixData.ResizeTo(correctNbins, correctNbins);
    //}


    //delete hUnfoldedDummy;
    //delete hMigrationDummy;
    //delete hRecoDummy;
    //delete hTruthDummy;
    //delete hBGSubDataDummy;
  }


  //--------------------------------------------
  // Pushing unfolding matrix
  //--------------------------------------------
  //cout<<" Pushing unfolding matrix "<<endl;
  //for( int i = 0; i < unfoldingCovMatrixData.GetNrows(); ++i ) unfoldingCovMatrixData(i,i)=0;
  //h_data_unfolded->PushCovMatrix( "unfoldingCov", unfoldingCovMatrixData );
    
    

  //--------------------------------------------
  // Efficiency Correction 
  //--------------------------------------------
  MnvH1D *h_data_effcor = NULL, *h_mc_effcor = NULL;
  cout<< "---------------------------------------"  << endl;
  cout<< "Correcting By Efficiency"  << endl;
  cout<< "---------------------------------------"  << endl;
  //Let's calculate the efficiency on the fly now
  //pz

  cout << "PZ EFF" << endl;
  CorrectByEfficiency(h_data_effcor, h_data_unfolded,h_effhist_num,h_effhist_dem);
  CorrectByEfficiency(h_mc_effcor, h_mc_unfolded, h_effhist_num, h_effhist_dem);


  MnvH2D *h_enu_data_effcor = NULL, *h_enu_mc_effcor = NULL;
  CorrectByEfficiency(h_enu_data_effcor, h_enu_data_unfolded,h_enu_effhist_num,h_enu_effhist_dem);
  CorrectByEfficiency(h_enu_mc_effcor, h_enu_mc_unfolded, h_enu_effhist_num, h_enu_effhist_dem);

  //cout << "Efficiency correcting the unfolding covariance matrix" << endl;
  //TMatrixD unfoldingCovMatrixEff;
  //TMatrixD tmpMat = divideCovByHists(unfoldingCovMatrixOrig,h_effhist_num[kQELike],h_effhist_dem[kQELike]);
  //int ncols = tmpMat.GetNcols();
  //unfoldingCovMatrixEff.ResizeTo(ncols,ncols);
  //unfoldingCovMatrixEff=tmpMat;
  //-----------------------------------------------
  // Special section to get the MnvH1D's we need for the signal component cross sections
  // Need QElike_qe, QELike_res, QELike_dis, QELike_oth
  //-----------------------------------------------
  //MnvH1D *h_truth_qelike_qe;
  //MnvH1D *h_truth_qelike_res;
  //MnvH1D *h_truth_qelike_dis;
  //MnvH1D *h_truth_qelike_2p2h;
  //h_truth_qelike_qe =   (MnvH1D*)h_effhist_dem[kQELike_QE]->Clone(Form("h_truth_q2qe_qelike_qe"));
  //h_truth_qelike_res =  (MnvH1D*)h_effhist_dem[kQELike_RES]->Clone(Form("h_truth_q2qe_qelike_res"));
  //h_truth_qelike_dis =  (MnvH1D*)h_effhist_dem[kQELike_DIS]->Clone(Form("h_truth_q2qe_qelike_dis"));
  //h_truth_qelike_2p2h = (MnvH1D*)h_effhist_dem[kQELike_2p2h]->Clone(Form("h_truth_q2qe_qelike_2p2h"));

  cout << "QELike_QE:\t" <<   h_effhist_dem[kQELike_QE]->GetTitle() << endl;
  cout << "QELike_RES:\t" <<  h_effhist_dem[kQELike_RES]->GetTitle() << endl;
  cout << "QELike_DIS:\t" <<  h_effhist_dem[kQELike_DIS]->GetTitle() << endl;
  cout << "QELike_OTH:\t" <<  h_effhist_dem[kQELike_OTH]->GetTitle() << endl;
  
  //-------------------------------------------------------
  // Normalize by Flux and number of nucleons on targets
  //-------------------------------------------------------
  MnvH1D *h_data_cross_section, *h_mc_cross_section;
  MnvH2D *h_enu_data_cross_section, *h_enu_mc_cross_section;
  //The components
    //MnvH1D *h_cross_section_qelike_qe;
    //MnvH1D *h_cross_section_qelike_res;
    //MnvH1D *h_cross_section_qelike_dis;
    //MnvH1D *h_cross_section_qelike_2p2h;
  cout<< "--------------------------------------------"  << endl;
  cout<< "Normalizing by Flux and Number of Targets"  << endl;
  cout<< "--------------------------------------------"  << endl;
  if(doenu){
    MnvH1D *h_flux_normalization = (MnvH1D*)h_data_effcor->Clone("h_flux_normalization");
    h_flux_normalization->ClearAllErrorBands();
    h_flux_normalization->Reset();
    
    //Calculate the flux integrated value for each universe and set those values in bins of pt-pz
    //Integrate from 0.0 to 100.0 GeV (JO:we're removing the neutrino energy cut)
    //double e_min = 0.0;//GeV
    //double e_max = 100.;//GeV
    //int b_min = h_flux->FindBin( e_min );
    //int b_max = h_flux->FindBin( e_max );
    
    //double flux_cv = h_flux->Integral( b_min, b_max, "width" );
    
    //const int lowBin = 0;
    const int highBinX = h_flux_normalization->GetNbinsX()+1;
    //const int highBinY = h_flux_normalization->GetNbinsY()+1;
    //const int highBin = h_flux_normalization->GetBin( highBinX);//, highBinY );
    cout << "Doing this by energy" << endl;
    //strategy is get 1D th1d and fill up h_flux_normalization
    TH1D* h_rebinned_flux = new TH1D(h_flux_normalization->GetCVHistoWithStatError());
    TH1D* tmphist = new TH1D(h_flux->GetCVHistoWithStatError());
    RebinFluxHist(tmphist, h_rebinned_flux);
    for(int i=0;i<highBinX;i++){
      //for(int j=0;j<highBinY;j++){
        //h_flux_normalization->SetBinContent(i,j,h_rebinned_flux->GetBinContent(i));
        h_flux_normalization->SetBinContent(i,h_rebinned_flux->GetBinContent(i));
      //}
    }
    cout << "Now for verts" << endl;
    //Do the same with the vertical error bands
    std::vector<std::string> vertNames = h_flux->GetVertErrorBandNames();
    for(unsigned int k=0; k<vertNames.size(); ++k ){
      MnvVertErrorBand *errBand = h_flux->GetVertErrorBand( vertNames[k] );
      int universes = errBand->GetNHists();
      if(vertNames[k]=="Flux") universes=100;
      std::vector<TH1D*> vert_hists;
      for(int u=0; u< universes; ++u)
      {
        TH1D *h_vert = new TH1D( (TH1D)*h_flux_normalization);
        h_vert->SetName( Form("h_vert_%s_%i", vertNames[k].c_str(), u) );
        
        //strategy is get 1D th1d and fill up h_flux_normalization
        TH1D* h_rebinned_flux = new TH1D(h_flux_normalization->GetCVHistoWithStatError());
        TH1D* tmphist = new TH1D(*errBand->GetHist( u ));
        RebinFluxHist(tmphist, h_rebinned_flux);
        for(int i=0;i<highBinX;i++){
          //for(int j=0;j<highBinY;j++){
            h_flux_normalization->SetBinContent(i,h_rebinned_flux->GetBinContent(i));
          //}
        }
        vert_hists.push_back( h_vert );
      }
      h_flux_normalization->AddVertErrorBand( vertNames[k], vert_hists );
      
      //cleaning
      for(std::vector<TH1D*>::iterator itHist = vert_hists.begin(); itHist != vert_hists.end(); ++itHist) 
        delete *itHist;
    }
    cout << "Now for lateral" << endl;
    //Do the same with the lateral error bands
    std::vector<std::string> latNames = h_flux->GetLatErrorBandNames();
    for(unsigned int k=0; k<latNames.size(); ++k ){
      cout << k << "\t" << latNames[k] << endl;
      MnvLatErrorBand *errBand = h_flux->GetLatErrorBand( latNames[k] );
      int universes = errBand->GetNHists();
      std::vector<TH1D*> lat_hists;
      cout << k << "\t" << latNames[k] << "\t" <<universes << endl;
      for(int u=0; u< universes; ++u)
      {
        TH1D *h_lat = new TH1D( (TH1D)*h_flux_normalization);
        h_lat->SetName( Form("h_lat_%s_%i", latNames[k].c_str(), u) );
        //strategy is get 1D th1d and fill up h_flux_normalization
        TH1D* h_rebinned_flux = new TH1D(h_flux_normalization->GetCVHistoWithStatError());
        TH1D* tmphist = new TH1D(*errBand->GetHist( u ));
        RebinFluxHist(tmphist, h_rebinned_flux);
        for(int i=0;i<highBinX;i++){
          //for(int j=0;j<highBinY;j++){
            h_flux_normalization->SetBinContent(i,h_rebinned_flux->GetBinContent(i));
          //}
        }
        lat_hists.push_back( h_lat );
      }
      h_flux_normalization->AddLatErrorBand( latNames[k], lat_hists );
      
      //cleaning
      for(std::vector<TH1D*>::iterator itHist = lat_hists.begin(); itHist != lat_hists.end(); ++itHist)
        delete *itHist;
    }
    
    
    //enu
    cout << "Enu Norm" << endl;
    NormalizeByFluxAndTargets( h_data_cross_section, h_data_effcor, h_flux, pot_data, applyFluxConstraint, utils, h_flux_normalization, false);
    NormalizeByFluxAndTargets( h_mc_cross_section, h_mc_effcor, h_flux, pot_data, applyFluxConstraint, utils, h_flux_normalization,true );



    MnvH2D* h_flux_normalization2D = (MnvH2D*) h_enu_data_effcor->Clone("h_flux_normalization2D");
    h_flux_normalization2D->ClearAllErrorBands();
    h_flux_normalization2D->Reset();
    ExpandHisto( h_flux_normalization, h_flux_normalization2D, 0 );

    NormalizeByFluxAndTargets( h_enu_data_cross_section, h_enu_data_effcor, h_flux, pot_data, applyFluxConstraint, utils, h_flux_normalization2D, false);
    NormalizeByFluxAndTargets( h_enu_mc_cross_section, h_enu_mc_effcor, h_flux, pot_data, applyFluxConstraint, utils, h_flux_normalization2D, true );
    //cout << "Enu Norm + comps" << endl;
    //NormalizeByFluxAndTargets( h_cross_section_qelike_qe, h_truth_qelike_qe, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff, h_flux_normalization );
    //cout << "res" << endl;
    //NormalizeByFluxAndTargets( h_cross_section_qelike_res, h_truth_qelike_res, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff, h_flux_normalization );
    //cout << "dis" << endl;
    //NormalizeByFluxAndTargets( h_cross_section_qelike_dis, h_truth_qelike_dis, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff, h_flux_normalization );
    //cout << "2p2h" << endl;
    //NormalizeByFluxAndTargets( h_cross_section_qelike_2p2h, h_truth_qelike_2p2h, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff, h_flux_normalization );
    //cout << "End enu comps" << endl;
  }
  else{

    cout<<"normalize q2"<<endl;
    NormalizeByFluxAndTargets( h_data_cross_section, h_data_effcor, h_flux, pot_data, applyFluxConstraint, utils , NULL, false);
    NormalizeByFluxAndTargets( h_mc_cross_section, h_mc_effcor, h_flux, pot_data, applyFluxConstraint, utils, NULL, true );

    cout<<"normalize enu-q2"<<endl;

    NormalizeByFluxAndTargets( h_enu_data_cross_section, h_enu_data_effcor, h_flux, pot_data, applyFluxConstraint, utils, NULL, false);
    NormalizeByFluxAndTargets( h_enu_mc_cross_section, h_enu_mc_effcor, h_flux, pot_data, applyFluxConstraint, utils, NULL, true);
    //cout << "doing specials" << endl;
    //NormalizeByFluxAndTargets( h_cross_section_qelike_qe, h_truth_qelike_qe, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff );
    //cout << "doing specials" << endl;
    //NormalizeByFluxAndTargets( h_cross_section_qelike_res, h_truth_qelike_res, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff );
    //cout << "doing specials" << endl;
    //NormalizeByFluxAndTargets( h_cross_section_qelike_dis, h_truth_qelike_dis, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff );
    //cout << "doing specials" << endl;
    //NormalizeByFluxAndTargets( h_cross_section_qelike_2p2h, h_truth_qelike_2p2h, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff );
    //cout << "done doing specials" << endl;

  }
  
  cout<< "----------------------------------------------------"  << endl;
  cout<< "Adding additional Systematics to DATA Cross Section"  << endl;
  cout<< "----------------------------------------------------"  << endl;
  cout<< "->Add Vertical Error non-reweightable GENIE: Xtalk+EFNUCR+FZONE+Hadronization_Alt1 (1 universe)"  << endl;
  //if(run_type==1){
  //  bool genieAdded = AddNonReweightableGENIESys(h_data_cross_section);
  //  if (!genieAdded)
  //    cout<<"CrossSectionsHists::AddNonReweightableGENIESys: Error adding GENIE non-reweightable systematic"<<endl;
  //}
  //==================================================================
  // Create ROOT file to store histograms
  // Write to file all the created histograms 
  //==================================================================
  TFile *f_output = new TFile( output_filename.c_str(), "RECREATE");
  f_output->cd();
  cout<< "--------------------------------------------"  << endl;
  cout<< "Writing histograms to file: "<<output_filename  << endl;
  cout<< "--------------------------------------------"  << endl;

  TVector2 *pot = (TVector2*) f_signal_regions->Get("pot");
  pot->Write("pot");
  h_data_nobck->Write();
  h_mc_nobck->Write();
  h_data_unfolded->Write();
  h_mc_unfolded->Write();
  h_data_effcor->Write();
  h_mc_effcor->Write();
  h_data_cross_section->Write();
  h_mc_cross_section->Write();
  h_effhist_num[kQELike_QE_H]->Write();
  h_effhist_dem[kQELike_QE_H]->Write();
  unfoldingCovMatrixData.Write("unfoldingMatrixData");

  h_enu_data_nobck->Write();
  h_enu_mc_nobck->Write();
  h_enu_data_unfolded->Write();
  h_enu_mc_unfolded->Write();
  h_enu_data_effcor->Write();
  h_enu_mc_effcor->Write();
  h_enu_data_cross_section->Write();
  h_enu_mc_cross_section->Write();
  h_enu_effhist_num[kQELike_QE_H]->Write();
  h_enu_effhist_dem[kQELike_QE_H]->Write();


  
  //h_cross_section_qelike_qe->Write(Form("%s_cross_section_qelike_qe",histnameprefix[run_type].c_str()));
  //h_cross_section_qelike_res->Write(Form("%s_cross_section_qelike_res",histnameprefix[run_type].c_str()));
  //h_cross_section_qelike_dis->Write(Form("%s_cross_section_qelike_dis",histnameprefix[run_type].c_str()));
  //h_cross_section_qelike_2p2h->Write(Form("%s_cross_section_qelike_2p2h",histnameprefix[run_type].c_str()));
  
  cout<< "All done..."  << endl;
  return 0;
};



//void SubtractBackground(MnvH2D** h, MnvH2D* h_bg_weights, MnvH2D* &h_data_nobck, MnvH2D* &h_mc_nobck, MnvH2D* &h_mc_tuned_bck, MnvH2D* &h_data_predicted_bck, std::string label)
//{
//  std::cout << "Subtracting background ..."<< std::endl;
//  //Tune MC background (For now, just bring the standard not tuned MC bakground)
//  h_mc_tuned_bck = (MnvH2D*) h[kQELikeNot]->Clone( Form( "_%s_tuned_bck_%s", h[kMC]->GetName(), label.c_str() ) );
//  ///h_mc_tuned_bck->Multiply( h_mc_tuned_bck, h_bg_weights );
//  //------------------------------------------------------
//  // Calculate data background = MC_bg x weight x DATA/MC
//  //-------------------------------------------------------
//  //Clone from MC so it keeps the error band universes...
//  h_data_predicted_bck = (MnvH2D*)h[kMC]->Clone( Form( "%s_predicted_bck_%s", h[kData]->GetName(), label.c_str() ) );
//  h_data_predicted_bck->Divide( h[kQELikeNot], h[kMC], 1.0, 1.0, "B");// MC background fraction
//  std::cout << "Multiply MC_bg/MC x bg_weights ..."<< std::endl;
//  h_data_predicted_bck->Multiply( h_data_predicted_bck, h_bg_weights); // weighted MC bkg frac = Data bkg frac`
//  //(Keep a different histogram up to this step for later)
//  MnvH2D *h_bg_fraction = (MnvH2D*)h_data_predicted_bck->Clone("h_bg_fraction");
//  std::cout << "Multiply DATA x MC_bg/MC x bg_weights ..."<< std::endl;
//  MnvH2D *h_data = (MnvH2D*)h[kData]->Clone("h_data");
//  h_data->AddMissingErrorBandsAndFillWithCV( *(h[kMC]) );
//  h_data_predicted_bck->Multiply( h_data_predicted_bck, h_data );//predicted background`
//  std::cout << "Done with multiplication ..."<< std::endl;
//
//  //-------------------------------------------------
//  // Subtract background
//  //--------------------------------------------------
// //-------------------------------------------------
//  //Note: Consider:
//  //Data_bg_predicted  = MC_bg x weight x DATA/MC
//  //DATA_bg_subtracted = DATA - DATA_bg_predicted
//  //                   = DATA x ( 1 - weight x MC_bg/MC)
//  // Instead of sutracting DATA - DATA_bg predicted, use
//  // the second derived option, because the first one
//  // double counts data stat uncertainties
//  // See docdb-8719, slide 3
//  //-------------------------------------------------
//
//  h_mc_nobck   = (MnvH2D*)h[kMC]->Clone( Form( "%s_nobck_%s", h[kMC]->GetName(), label.c_str() ) );
//  h_data_nobck = (MnvH2D*)h[kData]->Clone( Form( "%s_nobck_%s", h[kData]->GetName(), label.c_str() ) );
//  std::cout << "Add missing error bands ..."<< std::endl;
//  h_data_nobck->AddMissingErrorBandsAndFillWithCV( *h_mc_nobck );
//  std::cout << "Done adding missing error bands ..."<< std::endl;
//
//  h_mc_nobck->Add( h[kQELikeNot], -1.0);
//  const int highbin = h[kMC]->GetBin( h[kMC]->GetNbinsX()+1, h[kMC]->GetNbinsY()+1 ); // considering under/overflow
// 
//  //Get CV Data Background subtracted
//  //Create a new histogram tht will be = 1 - weight x MC_bg/MC
//  //Clone just the CV histogram, add error universes later 
//  MnvH2D *h_fraction_predicted = new MnvH2D(h_bg_fraction->GetCVHistoWithStatError() );
//  h_fraction_predicted->SetName("h_fraction_predicted");
//  //h_fraction_predicted->Reset();
//  for(int i=0; i<= highbin; ++i)
//  {
//    double bg_fraction_i       = h_bg_fraction->GetBinContent( i );
//    double bg_fraction_error_i = h_bg_fraction->GetBinError( i );
//    // if(bg_fraction_i < 0. || bg_fraction_i > 1.)
//    //   std::cout<<"WARNING: bkgd fraction is out of range for bin "<<i<<": "<<bg_fraction_i<<std::endl;
//    h_fraction_predicted->SetBinContent(i, 1. - bg_fraction_i );
//    h_fraction_predicted->SetBinError(i, bg_fraction_error_i );
//  }
//
//  //Do the same for vertical errors 
//  std::vector<std::string> vertNames = h_bg_fraction->GetVertErrorBandNames();
//  for(unsigned int i=0; i<vertNames.size(); ++i )
//  {
//    std::vector<TH2D*> vert_hists;
//    MnvVertErrorBand2D *errBand = h_bg_fraction->GetVertErrorBand(vertNames[i]);
//    int nUniverses = errBand->GetNHists();
//    for( int u = 0; u< nUniverses;  ++u)
//    {
//      TH2D* h_bg_fraction_u = errBand->GetHist( u );
//      TH2D* h_vert = (TH2D*)h_fraction_predicted->Clone("h_vert");
//      h_vert->Reset();
//      for(int j =0; j<=highbin;++j)
//      {
//        double bg_vert_fraction_j       = h_bg_fraction_u->GetBinContent( j );
//        double bg_vert_fraction_error_j = h_bg_fraction_u->GetBinError( j );
//        // if(bg_vert_fraction_j < 0. || bg_vert_fraction_j > 1.)
//        //   std::cout<<"WARNING: bkgd fraction is out of range for bin "<<j<<": "<<bg_vert_fraction_j<<std::endl;
//        h_vert->SetBinContent(j, 1. - bg_vert_fraction_j );
//        h_vert->SetBinError(j, bg_vert_fraction_error_j );
//      }
//      vert_hists.push_back( h_vert );
//    }
//    std::cout << "->Add Vertical Error ("<< vert_hists.size() << " Universes): " << vertNames[i].c_str() << std::endl;
//    h_fraction_predicted->AddVertErrorBand( vertNames[i], vert_hists );
//    
//    //cleaning
//    for(std::vector<TH2D*>::iterator itHist = vert_hists.begin(); itHist != vert_hists.end(); ++itHist)
//      delete *itHist;
//  }
//  
//  //Do the same for lateral errors
//  std::vector<std::string> latNames = h_bg_fraction->GetLatErrorBandNames();
//  for(unsigned int i=0; i<latNames.size(); ++i )
//  {
//    std::vector<TH2D*> lat_hists;
//    MnvLatErrorBand2D *errBand = h_bg_fraction->GetLatErrorBand(latNames[i]);
//    int nUniverses = errBand->GetNHists();
//    for( int u = 0; u< nUniverses;  ++u)
//    {
//      TH2D* h_bg_fraction_u = errBand->GetHist( u );
//      TH2D* h_lat = (TH2D*)h_fraction_predicted->Clone("h_lat");
//      h_lat->Reset();
//      for(int j =0; j<=highbin;++j)
//      {
//        double bg_lat_fraction_j       = h_bg_fraction_u->GetBinContent( j );
//        double bg_lat_fraction_error_j = h_bg_fraction_u->GetBinError( j );
//        // if(bg_lat_fraction_j < 0. || bg_lat_fraction_j > 1.)
//        //   std::cout<<"WARNING: bkgd fraction is out of range for bin "<<j<<": "<<bg_lat_fraction_j<<std::endl;
//        h_lat->SetBinContent(j, 1. - bg_lat_fraction_j );
//        h_lat->SetBinError(j, bg_lat_fraction_error_j );
//      }
//      lat_hists.push_back( h_lat );
//    }
//    std::cout << "->Add Lateral Error ("<< lat_hists.size() << " Universes): " << latNames[i].c_str() << std::endl;
//    h_fraction_predicted->AddLatErrorBand( latNames[i], lat_hists );
//    
//    //cleaning
//    for(std::vector<TH2D*>::iterator itHist = lat_hists.begin(); itHist != lat_hists.end(); ++itHist)
//      delete *itHist;
//  }
//  
//  //Background subtracted distribution is = DATA x h_fraction_predicted
//  h_data_nobck->Multiply( h_data_nobck, h_fraction_predicted );
//  
//  //cleaning
//  delete h_bg_fraction;
//  delete h_data;
//  delete h_fraction_predicted;
//
//}

//void SubtractBackground(MnvH2D** h, MnvH2D* h_bg_weights, MnvH2D* &h_data_nobck, MnvH2D* &h_mc_nobck){
//  MnvH2D* background_untuned = (MnvH2D*)h[kQELikeNot]->Clone( Form( "%s_untunedbck", h[kQELikeNot]->GetName()) );
//  MnvH2D* background_tuned = (MnvH2D*)h[kQELikeNot]->Clone( Form( "%s_tunedbck", h[kQELikeNot]->GetName()) );
//  h_mc_nobck   = (MnvH2D*)h[kMC]->Clone( Form( "%s_nobck", h[kMC]->GetName()) );
//  h_data_nobck = (MnvH2D*)h[kData]->Clone( Form( "%s_nobck", h[kData]->GetName()) );
//  h_data_nobck->AddMissingErrorBandsAndFillWithCV( *h_mc_nobck );//<---- make sure errors continue on!??
//  background_tuned->Multiply(background_untuned,h_bg_weights);//This is the data driven background!
//  h_data_nobck->Add(background_tuned,-1.0);
//  //  h_data_nobck->Add( h[kQELikeNot], -1.0);
//  h_mc_nobck->Add( h[kQELikeNot], -1.0);
//}

template<class MnvHnD> void SubtractBackground(MnvHnD** hists, MnvHnD* &h_data_nobck, MnvHnD* &h_data_predicted_bck, MnvHnD* &h_mc_nobck, MnvHnD*& h_mc_tuned_bck , std::string label)
{

  h_mc_nobck = (MnvHnD*) hists[kQELike_QE_H]->Clone( Form("%s_nobck_%s", hists[kQELike_QE_H]->GetName(), label.c_str()) );
  h_mc_tuned_bck = (MnvHnD*) hists[kMC]->Clone( Form("%s_tuned_bck_%s", hists[kMC]->GetName(), label.c_str()) );
  h_mc_tuned_bck->Add( h_mc_nobck , -1 );

  h_data_nobck = (MnvHnD*) hists[kData]->Clone( Form("%s_nobck_%s", hists[kData]->GetName(), label.c_str()) );
  h_data_nobck->AddMissingErrorBandsAndFillWithCV( *hists[kMC] );
  h_data_nobck->Add( h_mc_tuned_bck, -1 );

  h_data_predicted_bck = (MnvHnD*) hists[kData]->Clone( Form("%s_predicted_bck_%s", hists[kData]->GetName(), label.c_str()) );
  h_data_predicted_bck->AddMissingErrorBandsAndFillWithCV( *h_mc_tuned_bck );
  h_data_predicted_bck->Add( h_data_nobck, -1 );
}

template void SubtractBackground(MnvH1D** hists, MnvH1D* &h_data_nobck, MnvH1D* &h_data_predicted_bck, MnvH1D* &h_mc_nobck, MnvH1D*& h_mc_tuned_bck , std::string label);
template void SubtractBackground(MnvH2D** hists, MnvH2D* &h_data_nobck, MnvH2D* &h_data_predicted_bck, MnvH2D* &h_mc_nobck, MnvH2D*& h_mc_tuned_bck , std::string label);

template<class T>
void CorrectByEfficiency( T *&h_effcor, T* h_unfolded, T** h_effhist )
{
  //create efficiency corrected histo
  h_effcor = (T*)h_unfolded->Clone( Form( "%s_effcor", h_unfolded->GetName() ) );
 
  //efficiency doesn't have lateral error bands
  //so add to efficiency the lateral error bands found in unfolded data
  //and use the efficiency central value to fill the universes 
  h_effhist[kQELike_QE_H]->AddMissingErrorBandsAndFillWithCV( *h_unfolded );
  //  h_effhist[kQELike]->Scale(0.2);

  //divide by efficiency 
  h_effcor->Divide( h_unfolded, h_effhist[kQELike_QE_H] );

}


template<class T>
void CorrectByEfficiency( T *&h_effcor, T* h_unfolded, T** h_effhist_num, T** h_effhist_dem )
{
  //create efficiency corrected histo
  h_effcor = (T*)h_unfolded->Clone( Form( "%s_effcor", h_unfolded->GetName() ) );


  //T *h_eff_cor_temp[nHistos], *h_eff_cor_temp_dem[nHistos];
  //for( unsigned int j = 0; j < nHistos; ++j ){
  //  h_eff_cor_temp[j] = (T*)h_effhist_num[j]->Clone(Form("temp_eff_cor_%s-%d",h_eff_hist_num[j]->GetName(),j));
  //  h_eff_cor_temp_dem[j] = (T*)h_effhist_dem[j]->Clone(Form("temp_eff_cor_%s-%d-dem",h_eff_hist_dem[j]->GetName(),j));
  //  h_eff_cor_temp[j] = 
  //  h_eff_cor_temp_dem[j]

  //  h_eff_cor_temp[j]->Divide(h_effhist_dem[j],h_effhist_num[j],1.0,1.0,"B");
  //}

  T *h_eff_corr_temp = (T*) h_effhist_num[kQELike_QE_H]->Clone("h_eff_corr_temp");
  T *h_eff_corr_temp_dem = (T*) h_effhist_dem[kQELike_QE_H]->Clone("h_eff_corr_temp_dem");

  //h_eff_corr_temp->ClearAllErrorBands();
  //h_eff_corr_temp_dem->ClearAllErrorBands();

  //h_eff_corr_temp->AddMissingErrorBandsAndFillWithCV( *h_unfolded );
  //h_eff_corr_temp_dem->AddMissingErrorBandsAndFillWithCV( *h_unfolded );

  h_eff_corr_temp->Divide( h_eff_corr_temp, h_eff_corr_temp_dem, 1.0, 1.0, "B" );


  //efficiency doesn't have lateral error bands
  //so add to efficiency the lateral error bands found in unfolded data 
  //and use the efficiency central value to fill the universes 
  //This is for the efficiency denominator which is NOT affected by this... So CV is fine.
  //h_eff_cor_temp[kQELike]->AddMissingErrorBandsAndFillWithCV( *h_unfolded );
    
  //  h_effhist[kQELike]->Scale(0.2);

  //divide by efficiency 
  h_effcor->Divide( h_unfolded, h_eff_corr_temp );


}
template void CorrectByEfficiency( MnvH2D *&h_effcor, MnvH2D* h_unfolded, MnvH2D** h_effhist ); 
template void CorrectByEfficiency( MnvH2D *&h_effcor, MnvH2D* h_unfolded, MnvH2D** h_effhist_num, MnvH2D** h_effhist_dem ); 
template void CorrectByEfficiency( MnvH1D *&h_effcor, MnvH1D* h_unfolded, MnvH1D** h_effhist ); 
template void CorrectByEfficiency( MnvH1D *&h_effcor, MnvH1D* h_unfolded, MnvH1D** h_effhist_num, MnvH1D** h_effhist_dem ); 

void RebinFluxHist(TH1D* h_flux, TH1D *&h_rebinned_flux){

  //strategy is to recale orig by bin width (undo bin width normalization) then combine bins and then rescale by binwidth again
  TH1D* scaler = (TH1D*)h_flux->Clone("fluxscaler");
  TH1D* flux_cv = (TH1D*)h_flux->Clone("fluxcvtomod");
  for(int i=1;i<scaler->GetNbinsX();i++) scaler->SetBinContent(i,scaler->GetBinWidth(i));
  flux_cv->Multiply(scaler);//undid bin width normalization
  vector<double>rebinned_flux_bin_edges;
  for(int i=1;i<h_rebinned_flux->GetNbinsX()+2;i++){//need low edge of overflow (high edge of last bin)
    rebinned_flux_bin_edges.push_back(h_rebinned_flux->GetBinLowEdge(i));
  }
  h_rebinned_flux = (TH1D*)flux_cv->Rebin(rebinned_flux_bin_edges.size()-1,"Fluxrebinned",&rebinned_flux_bin_edges[0]);
  h_rebinned_flux->Scale(1.0,"width");//And redo bin width norm
  //  h_rebinned_flux->SaveAs("FluxRebinned_v3.root");

}


void NormalizeByFluxAndTargets(MnvH1D*&h_normalized,MnvH1D*h_effcor, MnvH1D* h_flux, double pot_scale, bool applyFluxConstraint, CCQENuUtils* utils, MnvH1D*flux_norm, bool isMC )
{
  //-------------------------------------------------
  //! Normalize by Flux
  //--------------------------------------------------
  //h_flux->AddMissingErrorBandsAndFillWithCV( *h_effcor );

  MnvH1D *h_flux_normalization = (MnvH1D*)h_effcor->Clone("h_flux_normalization");
  h_flux_normalization->ClearAllErrorBands();
  h_flux_normalization->Reset();

  //Calculate the flux integrated value for each universe and set those values in bins of pt-pz
  //Integrate from 0.0 to 100.0 GeV (JO:we're removing the neutrino energy cut)
  double e_min = 0.0;//GeV
  double e_max = 100.;//GeV
  int b_min = h_flux->FindBin( e_min );
  int b_max = h_flux->FindBin( e_max );

  double flux_cv = h_flux->Integral( b_min, b_max, "width" );

  const int lowBin = 0;
  const int highBinX = h_flux_normalization->GetNbinsX()+1;
  const int highBinY = h_flux_normalization->GetNbinsY()+1;
  const int highBin = h_flux_normalization->GetBin( highBinX, highBinY );
  if(flux_norm==NULL){
    cout <<"cv\t" << flux_cv << endl;
    for( int i=lowBin; i <= highBin; ++i ){
      h_flux_normalization->SetBinContent( i, flux_cv );
    }

    cout << "Now for verts" << endl;
    //Do the same with the vertical error bands
    std::vector<std::string> vertNames = h_flux->GetVertErrorBandNames();
    for(unsigned int k=0; k<vertNames.size(); ++k )
    {
      MnvVertErrorBand *errBand = h_flux->GetVertErrorBand( vertNames[k] );
      int universes = errBand->GetNHists();
      if(vertNames[k]=="Flux") universes=100;
      std::vector<TH1D*> vert_hists;
      for(int u=0; u< universes; ++u)
      {
        
        TH1D *h_vert = new TH1D( (TH1D)*h_flux_normalization);
        h_vert->SetName( Form("h_vert_%s_%i", vertNames[k].c_str(), u) );
        if(flux_norm==NULL){
          double flux_vert = errBand->GetHist( u )->Integral( b_min, b_max, "width");
          //	cout <<"vert\t" << u<<"\t" << flux_vert << endl;
          for( int i=0; i <= highBin; ++i ){
            h_vert->SetBinContent( i, flux_vert );
          }
        }
        vert_hists.push_back( h_vert );
      }
      h_flux_normalization->AddVertErrorBand( vertNames[k], vert_hists );
      
      //cleaning
      for(std::vector<TH1D*>::iterator itHist = vert_hists.begin(); itHist != vert_hists.end(); ++itHist)
        delete *itHist;
    }
    cout << "Now for lateral" << endl;
    //Do the same with the lateral error bands
    std::vector<std::string> latNames = h_flux->GetLatErrorBandNames();
    for(unsigned int k=0; k<latNames.size(); ++k )
    {
      cout << k << "\t" << latNames[k] << endl;
      MnvLatErrorBand *errBand = h_flux->GetLatErrorBand( latNames[k] );
      const int universes = errBand->GetNHists();
      std::vector<TH1D*> lat_hists;
      cout << k << "\t" << latNames[k] << "\t" <<universes << endl;
      for(int u=0; u< universes; ++u)
        {
          TH1D *h_lat = new TH1D( (TH1D)*h_flux_normalization);
          h_lat->SetName( Form("h_lat_%s_%i", latNames[k].c_str(), u) );
          if(flux_norm==NULL){
            double flux_lat = errBand->GetHist( u )->Integral( b_min, b_max, "width");
            for( int i=0; i <= highBin; ++i ){
              h_lat->SetBinContent( i, flux_lat );
            }
          }
          lat_hists.push_back( h_lat );
        }
      h_flux_normalization->AddLatErrorBand( latNames[k], lat_hists );
      
      //cleaning
      for(std::vector<TH1D*>::iterator itHist = lat_hists.begin(); itHist != lat_hists.end(); ++itHist)
        delete *itHist;
    }
  }
  else h_flux_normalization=(MnvH1D*)flux_norm->Clone("enuflux");


  /*Not done anymore as it is directly included in the FluxReweighter code
  //Apply the flux constraint if specified 
  if( applyFluxConstraint ) { 
    utils->loadFluxConstraint(); 
    utils->applyFluxConstraint<T, MnvVertErrorBand2D>( h_flux_normalization ); 
  }
  */
  cout << "Now the scales" << endl;
  //Convert flux units from m^2/POT to cm^2/POT
  h_flux_normalization->Scale( 1.0e-4 );
  h_normalized = (MnvH1D*)h_effcor->Clone( Form("%s_cross_section", h_effcor->GetName() ) );

  vector<string> normvert = h_normalized->GetVertErrorBandNames();
  vector<string> fluxnormvert = h_flux_normalization->GetVertErrorBandNames();

  for(UInt_t i=0;i<normvert.size();i++) cout << normvert[i] << "\t" << h_normalized->GetVertErrorBand(normvert[i])->GetNHists() << endl;
  for(UInt_t i=0;i<fluxnormvert.size();i++) cout << fluxnormvert[i] << "\t" << h_flux_normalization->GetVertErrorBand(fluxnormvert[i])->GetNHists() << endl;
  


  h_normalized->Divide( h_normalized, h_flux_normalization );

 
  //-------------------------------------------------
  //! Normalize by number of CarbonAtoms on targets
  //--------------------------------------------------
  const int nplanes = 2 * ( 80 - 27 + 1 );//fiducial volume -> modules 27-80
  TString name = h_normalized->GetName();

  double targets = 0.0;
  double apothem = 850;
  double mass = TargetUtils::Get().GetTrackerMass( nplanes, isMC, apothem );
  double massFrac = TargetUtils::Get().GetTrackerElementMassFraction(1, isMC);
  targets = mass*massFrac*TargetProp::AtomsPerGram::H  * TargetProp::ProtonsPerAtom::H;


  cout << TargetUtils::Get().GetTrackerNNeutrons(5980.,8422.,true,850.) << "\t" << TargetUtils::Get().GetNPlanes(5980,8422) <<"'\t This code." << endl;
  cout << TargetUtils::Get().GetTrackerNNeutrons(5990.,8340.,true,850.) << "\t" << TargetUtils::Get().GetNPlanes(5990,8340)<< "\t extractor. " << endl;
  cout << TargetUtils::Get().GetTrackerNCarbonAtoms( nplanes, /*isMC =*/ true, 850 ) << endl;
  cout << TargetUtils::Get().GetTrackerNProtons( nplanes, /*isMC =*/ true, 850 ) << endl;
  cout <<"MC\t" <<TargetUtils::Get().GetTrackerNNucleons( nplanes, /*isMC =*/ true, 850 ) << endl;//MC
  cout <<"Data\t" <<TargetUtils::Get().GetTrackerNNucleons( nplanes, /*isMC =*/ false, 850 ) << endl;//Data
  cout <<"MC\t" <<TargetUtils::Get().GetTrackerNNeutrons( nplanes, /*isMC =*/ true, 850 ) << endl;//MC
  cout <<"Data\t" <<TargetUtils::Get().GetTrackerNNeutrons( nplanes, /*isMC =*/ false, 850 ) << endl;//Data
  //Print out some useful information...
  double flux = h_flux_normalization->GetBinContent(1);
  double nC12 = TargetUtils::Get().GetTrackerNCarbonAtoms( nplanes, true, 850 );
  cout << " number of targets on H: " << targets << " protons" << endl;
  cout << " number of C12(MC): " << nC12 << endl;
  cout << " ratio of nH/nC12: " << targets/nC12 << endl;
  cout << " flux integrated: " << flux << endl;
  
  // Scale by targets and POT scale
  double scale = 1.0 / ( targets * pot_scale );
  h_normalized->Scale( scale );
  
  //Change units in Z axis
  h_normalized->GetYaxis()->SetTitle( "d#sigma_{QE}/dQ^{2} (cm^{2}/GeV^{2}/proton)" );
  cout << "Done" <<endl;
  //cout << "Scale the matrix" << endl;
  //scale,bin width norm, pushback

  //TMatrixD finalCovMatrix = divideCovByBinWidth(unfoldingMatrixToScale,h_normalized);
  // Jeremy tells me that the covariance matrix has the diagonal
  // errors on it, which are already included elsewhere, so we have to
  // subtract them off before adding the unfolding covariance matrix
  // back on
  //for(int i=0; i<finalCovMatrix.GetNrows(); ++i) finalCovMatrix(i,i)=0;

  //h_normalized->PushCovMatrix("unfoldingCov",finalCovMatrix);
  //finalCovMatrix.SaveAs("CovMat_tocheck.root");

  delete h_flux_normalization; 
}



void NormalizeByFluxAndTargets( MnvH2D *&h_normalized, MnvH2D *h_effcor, MnvH1D* h_flux, double pot_scale, bool applyFluxConstraint, CCQENuUtils* utils, MnvH2D *flux_norm, bool isMC )
{
  //-------------------------------------------------
  //! Normalize by Flux
  //--------------------------------------------------
  h_flux->AddMissingErrorBandsAndFillWithCV( *h_effcor->ProjectionX() );

  MnvH2D *h_flux_normalization = (MnvH2D*)h_effcor->Clone("h_flux_normalization");
  h_flux_normalization->ClearAllErrorBands();
  h_flux_normalization->Reset();

  TMatrixD unfoldCov = h_effcor->GetSysErrorMatrix("unfoldingCov");

  //Calculate the flux integrated value for each universe and set those values in bins of pt-pz
  //Integrate from 0.0 to 100.0 GeV (JO:we're removing the neutrino energy cut)
  double e_min = 0.0;//GeV
  double e_max = 100.;//GeV
  int b_min = h_flux->FindBin( e_min );
  int b_max = h_flux->FindBin( e_max );

  double flux_cv = h_flux->Integral( b_min, b_max, "width" );

  const int lowBin = 0;
  const int highBinX = h_flux_normalization->GetNbinsX()+1;
  const int highBinY = h_flux_normalization->GetNbinsY()+1;
  const int highBin = h_flux_normalization->GetBin( highBinX, highBinY );
  if(flux_norm==NULL){
    cout <<"cv\t" << flux_cv << endl;
    for( int i=lowBin; i <= highBin; ++i ){
      h_flux_normalization->SetBinContent( i, flux_cv );
    }

    cout << "Now for verts" << endl;
    //Do the same with the vertical error bands
    std::vector<std::string> vertNames = h_flux->GetVertErrorBandNames();
    for(unsigned int k=0; k<vertNames.size(); ++k )
      {
	MnvVertErrorBand *errBand = h_flux->GetVertErrorBand( vertNames[k] );
	int universes = errBand->GetNHists();
	if(vertNames[k]=="Flux") universes=100;
	std::vector<TH2D*> vert_hists;
	for(int u=0; u< universes; ++u)
	  {
	    
	    TH2D *h_vert = new TH2D( (TH2D)*h_flux_normalization);
	    h_vert->SetName( Form("h_vert_%s_%i", vertNames[k].c_str(), u) );
	    if(flux_norm==NULL){
	      double flux_vert = errBand->GetHist( u )->Integral( b_min, b_max, "width");
	      //	cout <<"vert\t" << u<<"\t" << flux_vert << endl;
	      for( int i=0; i <= highBin; ++i ){
		h_vert->SetBinContent( i, flux_vert );
	      }
	    }
	    vert_hists.push_back( h_vert );
	  }
	h_flux_normalization->AddVertErrorBand( vertNames[k], vert_hists );
	
	//cleaning
	for(std::vector<TH2D*>::iterator itHist = vert_hists.begin(); itHist != vert_hists.end(); ++itHist)
	  delete *itHist;
      }
    cout << "Now for lateral" << endl;
    //Do the same with the lateral error bands
    std::vector<std::string> latNames = h_flux->GetLatErrorBandNames();
    for(unsigned int k=0; k<latNames.size(); ++k )
      {
	cout << k << "\t" << latNames[k] << endl;
	MnvLatErrorBand *errBand = h_flux->GetLatErrorBand( latNames[k] );
	const int universes = errBand->GetNHists();
	std::vector<TH2D*> lat_hists;
	cout << k << "\t" << latNames[k] << "\t" <<universes << endl;
	for(int u=0; u< universes; ++u)
	  {
	    TH2D *h_lat = new TH2D( (TH2D)*h_flux_normalization);
	    h_lat->SetName( Form("h_lat_%s_%i", latNames[k].c_str(), u) );
	    if(flux_norm==NULL){
	      double flux_lat = errBand->GetHist( u )->Integral( b_min, b_max, "width");
	      for( int i=0; i <= highBin; ++i ){
		h_lat->SetBinContent( i, flux_lat );
	      }
	    }
	    lat_hists.push_back( h_lat );
	  }
	h_flux_normalization->AddLatErrorBand( latNames[k], lat_hists );
	
	//cleaning
	for(std::vector<TH2D*>::iterator itHist = lat_hists.begin(); itHist != lat_hists.end(); ++itHist)
	  delete *itHist;
      }
  }
  else
  { 
    h_flux_normalization=(MnvH2D*)flux_norm->Clone("enuflux");
  }


  /*Not done anymore as it is directly included in the FluxReweighter code
  //Apply the flux constraint if specified 
  if( applyFluxConstraint ) { 
    utils->loadFluxConstraint(); 
    utils->applyFluxConstraint<MnvH2D, MnvVertErrorBand2D>( h_flux_normalization ); 
  }
  */
  cout << "Now the scales" << endl;
  //Convert flux units from m^2/POT to cm^2/POT
  h_flux_normalization->Scale( 1.0e-4 );
  h_normalized = (MnvH2D*)h_effcor->Clone( Form("%s_cross_section", h_effcor->GetName() ) );

  vector<string> normvert = h_normalized->GetVertErrorBandNames();
  vector<string> fluxnormvert = h_flux_normalization->GetVertErrorBandNames();

  for(int i=0;i<normvert.size();i++) cout << normvert[i] << "\t" << h_normalized->GetVertErrorBand(normvert[i])->GetNHists() << endl;
  for(int i=0;i<fluxnormvert.size();i++) cout << fluxnormvert[i] << "\t" << h_flux_normalization->GetVertErrorBand(fluxnormvert[i])->GetNHists() << endl;
  


  h_normalized->Divide( h_normalized, h_flux_normalization );

  //flux correct cov matrix
  TMatrixD unfoldingMatrixToScale = divideCovByHist(unfoldCov,h_flux_normalization);
  
  //-------------------------------------------------
  //! Normalize by number of CarbonAtoms on targets
  //--------------------------------------------------
  const int nplanes = 2 * ( 80 - 27 + 1 );//fiducial volume -> modules 27-80
  TString name = h_normalized->GetName();


  double targets = 0.0;
  double apothem = 850;
  double mass = TargetUtils::Get().GetTrackerMass( nplanes, isMC, apothem );
  double massFrac = TargetUtils::Get().GetTrackerElementMassFraction(1, isMC);
  targets = mass*massFrac*TargetProp::AtomsPerGram::H  * TargetProp::ProtonsPerAtom::H;

  cout << TargetUtils::Get().GetTrackerNNeutrons(5980.,8422.,true,850.) << "\t" << TargetUtils::Get().GetNPlanes(5980,8422) <<"'\t This code." << endl;
  cout << TargetUtils::Get().GetTrackerNNeutrons(5990.,8340.,true,850.) << "\t" << TargetUtils::Get().GetNPlanes(5990,8340)<< "\t extractor. " << endl;
  cout << TargetUtils::Get().GetTrackerNCarbonAtoms( nplanes, /*isMC =*/ true, 850 ) << endl;
  cout << TargetUtils::Get().GetTrackerNProtons( nplanes, /*isMC =*/ true, 850 ) << endl;
  cout <<"MC\t" <<TargetUtils::Get().GetTrackerNNucleons( nplanes, /*isMC =*/ true, 850 ) << endl;//MC
  cout <<"Data\t" <<TargetUtils::Get().GetTrackerNNucleons( nplanes, /*isMC =*/ false, 850 ) << endl;//Data
  cout <<"MC\t" <<TargetUtils::Get().GetTrackerNNeutrons( nplanes, /*isMC =*/ true, 850 ) << endl;//MC
  cout <<"Data\t" <<TargetUtils::Get().GetTrackerNNeutrons( nplanes, /*isMC =*/ false, 850 ) << endl;//Data
  //Print out some useful information...
  double flux = h_flux_normalization->GetBinContent(1);
  double nC12 = TargetUtils::Get().GetTrackerNCarbonAtoms( nplanes, true, 850 );
  cout << " number of targets: " << targets << " neutrons" << endl;
  cout << " number of C12(MC): " << nC12 << endl;
  cout << " ratio of nC12/neutrons: " << nC12/targets << endl;
  cout << " flux integrated: " << flux << endl;
  
  // Scale by targets and POT scale
  double scale = 1.0 / ( targets * pot_scale );
  h_normalized->Scale( scale );
  
  //Change units in Z axis
  h_normalized->GetZaxis()->SetTitle( "d#sigma_{QELike}/d{#mu}_{P_Z}d{#mu}_{P_T} (cm^{2}/GeV^{2}/C^{12})" );
  cout << "Scale the matrix" << endl;
  //scale,bin width norm, pushback

  unfoldingMatrixToScale*= scale*scale;
  TMatrixD finalCovMatrix = divideCovByBinWidth(unfoldingMatrixToScale,h_normalized);
  // Jeremy tells me that the covariance matrix has the diagonal
  // errors on it, which are already included elsewhere, so we have to
  // subtract them off before adding the unfolding covariance matrix
  // back on
  for(int i=0; i<finalCovMatrix.GetNrows(); ++i) finalCovMatrix(i,i)=0;

  //OKAY ZERO OUT BINS (ONLY IN pzmu_ptmu
  string mytitle = h_normalized->GetName();
  cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
  cout << mytitle << endl;
  cout << mytitle.find("pzmu") << endl;
  if(mytitle.find("pzmu")!=string::npos){
    cout << "I'm going zero out bins" << endl;
    ZeroUnreported(h_normalized);
    ZeroUnreported(finalCovMatrix,h_normalized);
  }
  //h_normalized->PushCovMatrix("unfoldingCov",finalCovMatrix);
  

  finalCovMatrix.SaveAs("CovMat_tocheck_outside.root");
  TMatrixD finalCovMatrix2 = (TMatrixD) h_normalized->GetSysErrorMatrix("unfoldingCov");
  finalCovMatrix2.SaveAs("CovMat_tocheck_inside.root");


  delete h_flux_normalization; 
}


bool AddNonReweightableGENIESys(MnvH2D* &xs){
  TFile  *nwfile = new TFile("../rootfiles/NonReweightables.root");
  TH2D* h_xs_genie = (TH2D*)nwfile->Get("nrw_sys");
  TH2D* cv_xs = (TH2D*)xs->GetCVHistoWithStatError().Clone("tmp_cv_hist");
  h_xs_genie->Multiply(cv_xs);

  //Create only 1 universe
  vector<TH2D*> genie_nonre;
  genie_nonre.push_back(h_xs_genie);
  xs->AddVertErrorBand("GENIE_nonreweightable", genie_nonre);

  return true;

}

int main( int argc, char *argv[])
{
  ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if (argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t./CrossSectionHists  Output_file  Input_file_1TRACK  Input_file_2TRACK  Input_file_BKG  Input_file_MIGRATION  Input_file_EFFICIENCY  Input_file_FLUX  FluxConstraint Num_iterations  Run_Type DoENU \n\n"<<
      "\t-output_file\t =\t Cross Section output file\n"<<
      "\t-input_file\t =\t Muon Event Selection(Signal) for Neutron samples \n"<<
      "\t-input_file_angles\t =\t Muon Event Selection(Signal) for Neutron Region samples \n"<<
      "\t-input_file_bkg\t =\t Background weights filename\n"<<
      "\t-input_file_migration\t =\t Migration Histogram\n"<<
      "\t-input_file_efficiency\t =\t Efficiency Histograms\n"<<
      "\t-UseFluxConstraint\t =\t UseFluxConstraint\n"<< 
      "\t-Number of Iterations\t =\t Number of Iters\n" << 
      "\t-Run Type, unused\t =\t can be any integer\n" << 
      "\t-Do cross section with flux normalization\t = \t 0 or 1\n" << 
      "\t********************************************************************************************** \n"<<
      "\t Please see : MuonSelectionHists.cxx, BackgroundWeights.cxx, MigrationMatrixHists.cxx, EffPurityHists.cxx and Ana/Flux/python/compute_flux.py for getting the necessary input files"<< std::endl; 
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    return 0; 
  }

  //! Default parameters
  std::vector<std::string> par;
  par.push_back("CrossSectionPlots");//0
  par.push_back( Form("%s/ana/rootfiles/CrossSectionPlots.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/MuonSelectionHists_Signal_2track.root",getenv("CCQENUROOT") ) );//2
  par.push_back( Form("%s/ana/rootfiles/BackgroundWeights.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/MigrationHistogram.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/EffPurityHistogram.root",getenv("CCQENUROOT") ) );//5
  par.push_back("0"); 
  par.push_back("4"); 
  par.push_back("-1"); 
  par.push_back("1");//9

  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;


  string output_filename = par[1]+"_hydrogen.root";
  string filename_signal_regions = par[2];
  string filename_bg_weights = par[3];
  string filename_migration_histogram =  par[4];
  string filename_eff = par[5];
  bool fluxConstraint = (par[6] == "1")? true: false;
    bool makeFluxConstraintHisto = fluxConstraint;
    bool applyFluxConstraint = fluxConstraint;
  int num_iter = atoi( par[7].c_str() );
  int run_type = atoi( par[8].c_str());
  bool doenu = (par[9] == "1")? true: false;

 CrossSectionHists( output_filename, filename_signal_regions, filename_bg_weights, filename_migration_histogram, filename_eff, num_iter, makeFluxConstraintHisto, applyFluxConstraint, run_type, doenu);
  return 0;

}
