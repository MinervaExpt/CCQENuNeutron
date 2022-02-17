#include "include/CCQENuPlotUtils.h"
#include "include/CCQENuUtils.h"
#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/TargetUtils.h"


using namespace CCQENU_ANA;


void SubtractBackground(MnvH2D** h, MnvH2D* h_bg_weights, MnvH2D* &h_data_nobck, MnvH2D* &h_mc_nobck, MnvH2D* &h_mc_tuned_bck, MnvH2D* &h_data_predicted_bck, std::string label);
void SubtractBackground(MnvH2D** h, MnvH2D* h_bg_weights, MnvH2D* &h_data_nobck, MnvH2D* &h_mc_nobck);

void CorrectByEfficiency( MnvH2D *&h_effcor, MnvH2D* h_unfolded, MnvH2D** h_effhist ); 
void CorrectByEfficiency( MnvH2D *&h_effcor, MnvH2D* h_unfolded, MnvH2D** h_effhist_num, MnvH2D** h_effhist_dem ); 
void NormalizeByFluxAndTargets( MnvH2D *&h_normalized, MnvH2D *h_effcor, MnvH1D* h_flux, double pot_scale, bool applyFluxConstraint, CCQENuUtils* utils, TMatrixD covUnfold,MnvH2D *flux_norm= NULL );
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

TMatrixD divideCovByHist(TMatrixD& m, TH2D* h)
{
  TMatrixD ret(m);
  for(int i=0; i<h->fN; ++i){
    for(int j=0; j<h->fN; ++j){
      ret(i,j)=m(i,j)/(h->fArray[i]*h->fArray[j]);
    }
  }
  return ret;
}

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

//need to clean out result for unreportable bins
void SetBinZero(MnvH2D *&xs, int x, int y){

  xs->SetBinContent(x,y,0);
  vector<string> verterrnames = xs->GetVertErrorBandNames();
  vector<string> laterrnames = xs->GetLatErrorBandNames();
  //Vertical Errors
  for(int i=0;i<verterrnames.size();i++){
    MnvVertErrorBand2D *tmperr = xs->GetVertErrorBand(verterrnames[i]);
    int nhist = tmperr->GetNHists();
    for(int j=0;j<nhist;j++){
      tmperr->GetHist(j)->SetBinContent(x,y,0);
    }
  }
  //Lateral Errors
  for(int i=0;i<laterrnames.size();i++){
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


int CrossSectionHists( string output_filename,string filename_2track, string filename_bg_weights, string filename_migration_histogram, string filename_eff, int num_iter, bool makeFluxConstraintHisto, bool applyFluxConstraint,int run_type,bool doenu)
{

  cout <<"*********************************" << endl;
  cout << "Starting" << endl;
  cout << output_filename << endl;
  cout <<"*********************************" << endl;

  // read file to get histograms
  TFile *f_2track = new TFile( filename_2track.c_str(), "READ" );
  if (f_2track->IsZombie() || f_2track->GetListOfKeys()->IsEmpty()){
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
  utils->setPlaylist("minervame1D");
  utils->setFluxReweighterPlaylist();
  CCQENuPlotUtils *plotUtils = new CCQENuPlotUtils();
  MinervaUnfold::MnvUnfold unfold;

  //-------------------------------------------
  // Grab event counts before cuts for norm
  //-------------------------------------------
  double data_events,mc_events;
  TVector2 *evt = new TVector2();
  evt = (TVector2*)f_2track->Get("n_events");
  data_events = evt->X();
  mc_events   = evt->Y();
  
  cout<< "Number of Data Events (2-track) = " << data_events << endl;
  cout<< "Number of MC Events (2-track) = " <<   mc_events << endl;

  //--------------------------------------------
  // Load histos to plot 
  //--------------------------------------------
  //2-D
  MnvH2D *h_myanalysisvar_2track[nHistos];

  std::map<int,string> histnameprefix;
  histnameprefix[1]="h_pzmu_ptmu";
  histnameprefix[2]="h_q2_ptmu";
  histnameprefix[3]="h_enu_ptmu";
  histnameprefix[4]="h_pzptrecshort";
  histnameprefix[5]="h_pzptrec";
  histnameprefix[6]="h_dalphat_ptmu";
  histnameprefix[7]="h_dphit_ptmu";
  histnameprefix[8]="h_pn_ptmu";
  histnameprefix[9]="h_dpt_ptmu";
  histnameprefix[10]="h_dptx_ptmu";
  histnameprefix[11]="h_dpty_ptmu";
  histnameprefix[12]="h_signed_ptmu";
  histnameprefix[13]="h_signeddalphat_ptmu";
  histnameprefix[14]="h_signeddphit_ptmu";
  plotUtils->bookHistos( f_2track, h_myanalysisvar_2track, histnameprefix[run_type]);

  //combine event rate for correct efficiency
  MnvH2D *h_myanalysisvar_all[nHistos];
  for(int i=0;i<nHistos;i++){
      h_myanalysisvar_all[i] = (MnvH2D*)h_myanalysisvar_2track[i]->Clone(Form("eff_num_%s",h_myanalysisvar_2track[i]->GetName()));
  }
  
  //2D BG Weights
  MnvH2D *h_myanalysisvar_bg_weights_2track;

  std::map<int,string> bkgnames;
  bkgnames[1]="pzptbins";
  bkgnames[2]="q2ptbins";
  bkgnames[3]="enuptbins";
  bkgnames[4]="pzptrecshortbins";
  bkgnames[5]="pzptrecbins";
  bkgnames[6]="dalphatptbins";
  bkgnames[7]="dphitptbins";
  bkgnames[8]="pnptbins";
  bkgnames[9]="dptptbins";
  bkgnames[10]="dptxptbins";
  bkgnames[11]="dptyptbins";
  bkgnames[12]="signedptbins";
  bkgnames[13]="signeddalphatptbins";
  bkgnames[14]="signeddphitptbins";

  h_myanalysisvar_bg_weights_2track = (MnvH2D*) f_bg_weights->Get(Form("h_weights_2track_%s_qelikenot",bkgnames[run_type].c_str()));

  cout << "Starting migration" << endl;
  //2D Migration Histogram (and its corresponding reco/truth parts necessary for the unfolding)
  MnvH2D* h_myanalysisvar_migration, *h_myanalysisvar_reco, *h_myanalysisvar_generated;
  //pz
  h_myanalysisvar_reco = (MnvH2D*)f_migration->Get(Form("%s_qelike_reco",histnameprefix[run_type].c_str()));
  h_myanalysisvar_generated = (MnvH2D*)f_migration->Get(Form("%s_qelike_truth",histnameprefix[run_type].c_str()));
  utils->combineObject(h_myanalysisvar_all[kMC], h_myanalysisvar_migration, Form("%s_qelike_migration",histnameprefix[run_type].c_str()), f_migration);
  cout << h_myanalysisvar_migration->GetVertErrorBandNames().size() << endl;
  cout << "Done combining pzmu" << endl;

  cout << "Booking Efficiency " << endl;
  //Efficiency Histogram (These include signal components
  MnvH2D *h_myanalysisvar_effhist[nHistos], *h_myanalysisvar_effhist_num[nHistos], *h_myanalysisvar_effhist_dem[nHistos];
  plotUtils->bookHistos( f_efficiency, h_myanalysisvar_effhist_num,Form("%s",histnameprefix[run_type].c_str()));
  plotUtils->bookHistos( f_efficiency, h_myanalysisvar_effhist_dem,Form("%s_truth",histnameprefix[run_type].c_str()));
  cout << "DONE" << endl;
  //Getting the Flux
  
  //New method for flux
  MnvH1D* h_flux = utils->getReweightedFluxHist();
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
  
  double pot_data = plotUtils->getPOTData( f_2track );
  double pot_mc = plotUtils->getPOTMC( f_2track );
  double pot_scale = plotUtils->getPOTNormFactor( f_2track );
  
  if(pot_scale==0) pot_scale = 1.0;
  cout<< "---------------------------------------"  << endl;
  cout<< "POT INFORMATION:"  << endl;
  cout<< "POT Data = " << pot_data << endl;
  cout<< "POT MC   = " << pot_mc << endl;
  cout<< "POT Scale factor = " << pot_scale << endl;
  cout<< "---------------------------------------"  << endl;

  //2-D
  cout << "Scale Rate" << endl;
  plotUtils->scaleMCHistos( h_myanalysisvar_2track, pot_scale);
  cout << "Scale migration" << endl;
  h_myanalysisvar_reco->Scale(pot_scale);
  h_myanalysisvar_generated->Scale(pot_scale);
  h_myanalysisvar_migration->Scale(pot_scale);

  //--------------------------------------------
  // Background Subtraction
  //--------------------------------------------
  cout<< "---------------------------------------"  << endl;
  cout<< "Background Subtraction"  << endl;
  cout<< "---------------------------------------"  << endl;

  MnvH2D *h_pzmu_data_nobck_2track, *h_pzmu_data_predicted_bck_2track, *h_pzmu_mc_nobck_2track, *h_pzmu_mc_tuned_bck_2track;
  MnvH2D *h_myanalysisvar_data_nobck, *h_myanalysisvar_mc_nobck;

  cout<<"Background Subtraction of ==2 Track Sub-Sample"<<endl;
  SubtractBackground( h_myanalysisvar_2track, h_myanalysisvar_bg_weights_2track, h_pzmu_data_nobck_2track, h_pzmu_mc_nobck_2track);
  
  //-----------------------------------------------
  // Merge both background subtracted sub-samples
  //-----------------------------------------------
  cout<<"Merging 1 and >=2 Tracks Sub-Samples..."<<endl;
  h_myanalysisvar_data_nobck = (MnvH2D*)h_pzmu_data_nobck_2track->Clone(Form("%s_data_nobck",histnameprefix[run_type].c_str()));
    
  h_myanalysisvar_mc_nobck = (MnvH2D*)h_pzmu_mc_nobck_2track->Clone(Form("%s_mc_nobck",histnameprefix[run_type].c_str()));


  
  //--------------------------------------------
  // Unfold DATA and MC merged samples
  //--------------------------------------------
  MnvH2D *h_myanalysisvar_data_unfolded = NULL, *h_myanalysisvar_mc_unfolded = NULL;
  TMatrixD unfoldingCovMatrixOrig_myanalysisvar;
  cout<< "---------------------------------------"  << endl;
  cout<< "Now Unfolding with " << num_iter << " iterations" << endl;
  cout<< "---------------------------------------"  << endl;
  if (num_iter == 0) {
    //------------------------------------------------------------------
    // Special case code where there are no requested iterations
    // Copies the hist then changes the name to end in _unfold
    //------------------------------------------------------------------
    // DATA
    h_myanalysisvar_data_unfolded = new MnvH2D(*h_myanalysisvar_data_nobck);
    h_myanalysisvar_data_unfolded->SetName(Form( "%s_unfold", h_myanalysisvar_data_unfolded->GetName() ));
    // MC
    h_myanalysisvar_mc_unfolded = new MnvH2D(*h_myanalysisvar_mc_nobck);
    h_myanalysisvar_mc_unfolded->SetName(Form( "%s_unfold", h_myanalysisvar_mc_unfolded->GetName() ));

  } else {
    //---------------------------------------------------------------------------------
    // Last input variable (= true) ensures that we are unfolding each universe  
    // of the distribution with its corresponding universe in the migration matrix   
    //---------------------------------------------------------------------------------
    //pz
    cout << "Unfolding" << endl;
    bool data_unfolded = unfold.UnfoldHisto2D(h_myanalysisvar_data_unfolded, h_myanalysisvar_migration, h_myanalysisvar_reco, h_myanalysisvar_generated, h_myanalysisvar_data_nobck , num_iter, true, true);
    bool mc_unfolded = unfold.UnfoldHisto2D(h_myanalysisvar_mc_unfolded, h_myanalysisvar_migration, h_myanalysisvar_reco, h_myanalysisvar_generated, h_myanalysisvar_mc_nobck, num_iter, true, true);
    if(!data_unfolded || !mc_unfolded){
      std::cout << "Unfolding failed for either data or MC. Please check " << std::endl;
    }
    
    //NOW we have to get the damn covariance of the unfolding which means doing the unfolding AGAIN!@!
    //Only need this on the data as the MC is used as a CV only.
    //pz
    cout << "Getting the covariance of the unfolding for PZ" << endl;

    TH2D* hUnfoldedDummy=new TH2D(h_myanalysisvar_data_unfolded->GetCVHistoWithStatError());
    TH2D* hMigrationDummy=new TH2D(h_myanalysisvar_migration->GetCVHistoWithStatError());
    TH2D* hRecoDummy=new TH2D(h_myanalysisvar_reco->GetCVHistoWithStatError());
    TH2D* hTruthDummy=new TH2D(h_myanalysisvar_generated->GetCVHistoWithStatError());
    TH2D* hBGSubDataDummy=new TH2D(h_myanalysisvar_data_nobck->GetCVHistoWithStatError());

    unfold.UnfoldHisto2D(hUnfoldedDummy, unfoldingCovMatrixOrig_myanalysisvar, hMigrationDummy, hRecoDummy, hTruthDummy, hBGSubDataDummy, num_iter);
    // There's a bug in RooUnfold that's making it return covariance  // matrices with two extra bins. Kill them here, with a check.  
    // Conveniently, this bug was being hidden by an offsetting bug in  
    // MnvH2D, which is now fixed
    int correctNbins=hUnfoldedDummy->fN;
    int matrixRows=unfoldingCovMatrixOrig_myanalysisvar.GetNrows();
    if(correctNbins!=matrixRows){
      cout << "****************************************************************************" << endl;
      cout << "*  Fixing unfolding matrix size because of RooUnfold bug. From " << matrixRows << " to " << correctNbins << endl;
      cout << "****************************************************************************" << endl;
      // It looks like this DTRT, since the extra last two bins don't have any content
      unfoldingCovMatrixOrig_myanalysisvar.ResizeTo(correctNbins, correctNbins);
    }


    delete hUnfoldedDummy;
    delete hMigrationDummy;
    delete hRecoDummy;
    delete hTruthDummy;
    delete hBGSubDataDummy;
  }
    
    

  //--------------------------------------------
  // Efficiency Correction 
  //--------------------------------------------
  MnvH2D *h_myanalysisvar_data_effcor = NULL, *h_myanalysisvar_mc_effcor = NULL;
  cout<< "---------------------------------------"  << endl;
  cout<< "Correcting By Efficiency"  << endl;
  cout<< "---------------------------------------"  << endl;
  //Let's calculate the efficiency on the fly now
  //pz

  cout << "PZ EFF" << endl;
  CorrectByEfficiency(h_myanalysisvar_data_effcor, h_myanalysisvar_data_unfolded,h_myanalysisvar_effhist_num,h_myanalysisvar_effhist_dem);
  CorrectByEfficiency( h_myanalysisvar_mc_effcor, h_myanalysisvar_mc_unfolded, h_myanalysisvar_effhist_num, h_myanalysisvar_effhist_dem);
  

  cout << "Efficiency correcting the unfolding covariance matrix" << endl;
  TMatrixD unfoldingCovMatrixEff_myanalysisvar;
  TMatrixD tmpMat = divideCovByHists(unfoldingCovMatrixOrig_myanalysisvar,h_myanalysisvar_effhist_num[kQELike],h_myanalysisvar_effhist_dem[kQELike]);
  int ncols = tmpMat.GetNcols();
  unfoldingCovMatrixEff_myanalysisvar.ResizeTo(ncols,ncols);
  unfoldingCovMatrixEff_myanalysisvar=tmpMat;
  //-----------------------------------------------
  // Special section to get the MnvH2D's we need for the signal component cross sections
  // Need QElike_qe, QELike_res, QELike_dis, QELike_oth
  //-----------------------------------------------
  MnvH2D *h_myanalysisvar_truth_qelike_qe;
  MnvH2D *h_myanalysisvar_truth_qelike_res;
  MnvH2D *h_myanalysisvar_truth_qelike_dis;
  MnvH2D *h_myanalysisvar_truth_qelike_2p2h;
  h_myanalysisvar_truth_qelike_qe = (MnvH2D*)h_myanalysisvar_effhist_dem[kQELike_QE]->Clone(Form("%s_qecomp",histnameprefix[run_type].c_str()));
  h_myanalysisvar_truth_qelike_res = (MnvH2D*)h_myanalysisvar_effhist_dem[kQELike_RES]->Clone(Form("%s_rescomp",histnameprefix[run_type].c_str()));
  h_myanalysisvar_truth_qelike_dis = (MnvH2D*)h_myanalysisvar_effhist_dem[kQELike_DIS]->Clone(Form("%s_discomp",histnameprefix[run_type].c_str()));
  h_myanalysisvar_truth_qelike_2p2h = (MnvH2D*)h_myanalysisvar_effhist_dem[kQELike_2p2h]->Clone(Form("%s_2p2hcomp",histnameprefix[run_type].c_str()));

  cout << "QELike_QE:\t" <<  h_myanalysisvar_effhist_dem[kQELike_QE]->GetTitle() << endl;
  cout << "QELike_RES:\t" <<  h_myanalysisvar_effhist_dem[kQELike_RES]->GetTitle() << endl;
  cout << "QELike_DIS:\t" <<  h_myanalysisvar_effhist_dem[kQELike_DIS]->GetTitle() << endl;
  cout << "QELike_OTH:\t" <<  h_myanalysisvar_effhist_dem[kQELike_OTH]->GetTitle() << endl;
  
  //-------------------------------------------------------
  // Normalize by Flux and number of nucleons on targets
  //-------------------------------------------------------
  MnvH2D *h_myanalysisvar_data_cross_section, *h_myanalysisvar_mc_cross_section;
  //The components
  MnvH2D *h_myanalysisvar_cross_section_qelike_qe;
  MnvH2D *h_myanalysisvar_cross_section_qelike_res;
  MnvH2D *h_myanalysisvar_cross_section_qelike_dis;
  MnvH2D *h_myanalysisvar_cross_section_qelike_2p2h;
  cout<< "--------------------------------------------"  << endl;
  cout<< "Normalizing by Flux and Number of Targets"  << endl;
  cout<< "--------------------------------------------"  << endl;
  if(doenu){
    MnvH2D *h_flux_normalization = (MnvH2D*)h_myanalysisvar_data_effcor->Clone("h_flux_normalization");
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
    cout << "Doing this by energy" << endl;
    //strategy is get 1D th1d and fill up h_flux_normalization
    TH1D* h_rebinned_flux = new TH1D(h_flux_normalization->ProjectionX("proj_x",0,-1)->GetCVHistoWithStatError());
    TH1D* tmphist = new TH1D(h_flux->GetCVHistoWithStatError());
    RebinFluxHist(tmphist, h_rebinned_flux);
    for(int i=0;i<highBinX;i++){
      for(int j=0;j<highBinY;j++){
	h_flux_normalization->SetBinContent(i,j,h_rebinned_flux->GetBinContent(i));
      }
    }
    cout << "Now for verts" << endl;
    //Do the same with the vertical error bands
    std::vector<std::string> vertNames = h_flux->GetVertErrorBandNames();
    for(unsigned int k=0; k<vertNames.size(); ++k ){
      MnvVertErrorBand *errBand = h_flux->GetVertErrorBand( vertNames[k] );
      int universes = errBand->GetNHists();
      if(vertNames[k]=="Flux") universes=100;
      std::vector<TH2D*> vert_hists;
      for(int u=0; u< universes; ++u)
	{
	  TH2D *h_vert = new TH2D( (TH2D)*h_flux_normalization);
	  h_vert->SetName( Form("h_vert_%s_%i", vertNames[k].c_str(), u) );
	  
	  //strategy is get 1D th1d and fill up h_flux_normalization
	  TH1D* h_rebinned_flux = new TH1D(h_flux_normalization->ProjectionX("proj_x",0,-1)->GetCVHistoWithStatError());
	  TH1D* tmphist = new TH1D(*errBand->GetHist( u ));
	  RebinFluxHist(tmphist, h_rebinned_flux);
	  for(int i=0;i<highBinX;i++){
	    for(int j=0;j<highBinY;j++){
	      h_flux_normalization->SetBinContent(i,j,h_rebinned_flux->GetBinContent(i));
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
    for(unsigned int k=0; k<latNames.size(); ++k ){
      cout << k << "\t" << latNames[k] << endl;
      MnvLatErrorBand *errBand = h_flux->GetLatErrorBand( latNames[k] );
      int universes = errBand->GetNHists();
      std::vector<TH2D*> lat_hists;
      cout << k << "\t" << latNames[k] << "\t" <<universes << endl;
      for(int u=0; u< universes; ++u)
	{
	  TH2D *h_lat = new TH2D( (TH2D)*h_flux_normalization);
	  h_lat->SetName( Form("h_lat_%s_%i", latNames[k].c_str(), u) );
	  //strategy is get 1D th1d and fill up h_flux_normalization
	  TH1D* h_rebinned_flux = new TH1D(h_flux_normalization->ProjectionX("proj_x",0,-1)->GetCVHistoWithStatError());
	  TH1D* tmphist = new TH1D(*errBand->GetHist( u ));
	  RebinFluxHist(tmphist, h_rebinned_flux);
	  for(int i=0;i<highBinX;i++){
	    for(int j=0;j<highBinY;j++){
	      h_flux_normalization->SetBinContent(i,j,h_rebinned_flux->GetBinContent(i));
	    }
	  }
	  lat_hists.push_back( h_lat );
	}
      h_flux_normalization->AddLatErrorBand( latNames[k], lat_hists );
      
      //cleaning
      for(std::vector<TH2D*>::iterator itHist = lat_hists.begin(); itHist != lat_hists.end(); ++itHist)
	  delete *itHist;
    }
    
    
    //enu
    cout << "Enu Norm" << endl;
    NormalizeByFluxAndTargets( h_myanalysisvar_data_cross_section, h_myanalysisvar_data_effcor, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff_myanalysisvar, h_flux_normalization);
    NormalizeByFluxAndTargets( h_myanalysisvar_mc_cross_section, h_myanalysisvar_mc_effcor, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff_myanalysisvar, h_flux_normalization );
    cout << "Enu Norm + comps" << endl;
    NormalizeByFluxAndTargets( h_myanalysisvar_cross_section_qelike_qe, h_myanalysisvar_truth_qelike_qe, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff_myanalysisvar, h_flux_normalization );
    cout << "res" << endl;
    NormalizeByFluxAndTargets( h_myanalysisvar_cross_section_qelike_res, h_myanalysisvar_truth_qelike_res, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff_myanalysisvar, h_flux_normalization );
    cout << "dis" << endl;
    NormalizeByFluxAndTargets( h_myanalysisvar_cross_section_qelike_dis, h_myanalysisvar_truth_qelike_dis, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff_myanalysisvar, h_flux_normalization );
    cout << "2p2h" << endl;
    NormalizeByFluxAndTargets( h_myanalysisvar_cross_section_qelike_2p2h, h_myanalysisvar_truth_qelike_2p2h, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff_myanalysisvar, h_flux_normalization );
    cout << "End enu comps" << endl;
  }
  else{

    NormalizeByFluxAndTargets( h_myanalysisvar_data_cross_section, h_myanalysisvar_data_effcor, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff_myanalysisvar );
    NormalizeByFluxAndTargets( h_myanalysisvar_mc_cross_section, h_myanalysisvar_mc_effcor, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff_myanalysisvar );
    cout << "doing specials" << endl;
    NormalizeByFluxAndTargets( h_myanalysisvar_cross_section_qelike_qe, h_myanalysisvar_truth_qelike_qe, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff_myanalysisvar );
    cout << "doing specials" << endl;
    NormalizeByFluxAndTargets( h_myanalysisvar_cross_section_qelike_res, h_myanalysisvar_truth_qelike_res, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff_myanalysisvar );
    cout << "doing specials" << endl;
    NormalizeByFluxAndTargets( h_myanalysisvar_cross_section_qelike_dis, h_myanalysisvar_truth_qelike_dis, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff_myanalysisvar );
    cout << "doing specials" << endl;
    NormalizeByFluxAndTargets( h_myanalysisvar_cross_section_qelike_2p2h, h_myanalysisvar_truth_qelike_2p2h, h_flux, pot_data, applyFluxConstraint, utils, unfoldingCovMatrixEff_myanalysisvar );
    cout << "done doing specials" << endl;

  }
  
  cout<< "----------------------------------------------------"  << endl;
  cout<< "Adding additional Systematics to DATA Cross Section"  << endl;
  cout<< "----------------------------------------------------"  << endl;
  cout<< "->Add Vertical Error non-reweightable GENIE: Xtalk+EFNUCR+FZONE+Hadronization_Alt1 (1 universe)"  << endl;
  if(run_type==1){
    bool genieAdded = AddNonReweightableGENIESys(h_myanalysisvar_data_cross_section);
    if (!genieAdded)
      cout<<"CrossSectionsHists::AddNonReweightableGENIESys: Error adding GENIE non-reweightable systematic"<<endl;
  }
  //==================================================================
  // Create ROOT file to store histograms
  // Write to file all the created histograms 
  //==================================================================
  TFile *f_output = new TFile( output_filename.c_str(), "RECREATE");
  f_output->cd();
  cout<< "--------------------------------------------"  << endl;
  cout<< "Writing histograms to file: "<<output_filename  << endl;
  cout<< "--------------------------------------------"  << endl;

  h_myanalysisvar_data_nobck->Write();
  h_myanalysisvar_mc_nobck->Write();
  h_myanalysisvar_data_unfolded->Write();
  h_myanalysisvar_mc_unfolded->Write();
  h_myanalysisvar_data_effcor->Write();
  h_myanalysisvar_mc_effcor->Write();
  h_myanalysisvar_data_cross_section->Write();
  h_myanalysisvar_mc_cross_section->Write();
  
  h_myanalysisvar_cross_section_qelike_qe->Write(Form("%s_cross_section_qelike_qe",histnameprefix[run_type].c_str()));
  h_myanalysisvar_cross_section_qelike_res->Write(Form("%s_cross_section_qelike_res",histnameprefix[run_type].c_str()));
  h_myanalysisvar_cross_section_qelike_dis->Write(Form("%s_cross_section_qelike_dis",histnameprefix[run_type].c_str()));
  h_myanalysisvar_cross_section_qelike_2p2h->Write(Form("%s_cross_section_qelike_2p2h",histnameprefix[run_type].c_str()));
  
  cout<< "All done..."  << endl;
  return 0;
};



void SubtractBackground(MnvH2D** h, MnvH2D* h_bg_weights, MnvH2D* &h_data_nobck, MnvH2D* &h_mc_nobck, MnvH2D* &h_mc_tuned_bck, MnvH2D* &h_data_predicted_bck, std::string label)
{
  std::cout << "Subtracting background ..."<< std::endl;
  //Tune MC background (For now, just bring the standard not tuned MC bakground)
  h_mc_tuned_bck = (MnvH2D*) h[kQELikeNot]->Clone( Form( "_%s_tuned_bck_%s", h[kMC]->GetName(), label.c_str() ) );
  ///h_mc_tuned_bck->Multiply( h_mc_tuned_bck, h_bg_weights );
  //------------------------------------------------------
  // Calculate data background = MC_bg x weight x DATA/MC
  //-------------------------------------------------------
  //Clone from MC so it keeps the error band universes...
  h_data_predicted_bck = (MnvH2D*)h[kMC]->Clone( Form( "%s_predicted_bck_%s", h[kData]->GetName(), label.c_str() ) );
  h_data_predicted_bck->Divide( h[kQELikeNot], h[kMC], 1.0, 1.0, "B");// MC background fraction
  std::cout << "Multiply MC_bg/MC x bg_weights ..."<< std::endl;
  h_data_predicted_bck->Multiply( h_data_predicted_bck, h_bg_weights); // weighted MC bkg frac = Data bkg frac`
  //(Keep a different histogram up to this step for later)
  MnvH2D *h_bg_fraction = (MnvH2D*)h_data_predicted_bck->Clone("h_bg_fraction");
  std::cout << "Multiply DATA x MC_bg/MC x bg_weights ..."<< std::endl;
  MnvH2D *h_data = (MnvH2D*)h[kData]->Clone("h_data");
  h_data->AddMissingErrorBandsAndFillWithCV( *(h[kMC]) );
  h_data_predicted_bck->Multiply( h_data_predicted_bck, h_data );//predicted background`
  std::cout << "Done with multiplication ..."<< std::endl;

  //-------------------------------------------------
  // Subtract background
  //--------------------------------------------------
 //-------------------------------------------------
  //Note: Consider:
  //Data_bg_predicted  = MC_bg x weight x DATA/MC
  //DATA_bg_subtracted = DATA - DATA_bg_predicted
  //                   = DATA x ( 1 - weight x MC_bg/MC)
  // Instead of sutracting DATA - DATA_bg predicted, use
  // the second derived option, because the first one
  // double counts data stat uncertainties
  // See docdb-8719, slide 3
  //-------------------------------------------------

  h_mc_nobck   = (MnvH2D*)h[kMC]->Clone( Form( "%s_nobck_%s", h[kMC]->GetName(), label.c_str() ) );
  h_data_nobck = (MnvH2D*)h[kData]->Clone( Form( "%s_nobck_%s", h[kData]->GetName(), label.c_str() ) );
  std::cout << "Add missing error bands ..."<< std::endl;
  h_data_nobck->AddMissingErrorBandsAndFillWithCV( *h_mc_nobck );
  std::cout << "Done adding missing error bands ..."<< std::endl;

  h_mc_nobck->Add( h[kQELikeNot], -1.0);
  const int highbin = h[kMC]->GetBin( h[kMC]->GetNbinsX()+1, h[kMC]->GetNbinsY()+1 ); // considering under/overflow
 
  //Get CV Data Background subtracted
  //Create a new histogram tht will be = 1 - weight x MC_bg/MC
  //Clone just the CV histogram, add error universes later 
  MnvH2D *h_fraction_predicted = new MnvH2D(h_bg_fraction->GetCVHistoWithStatError() );
  h_fraction_predicted->SetName("h_fraction_predicted");
  //h_fraction_predicted->Reset();
  for(int i=0; i<= highbin; ++i)
  {
    double bg_fraction_i       = h_bg_fraction->GetBinContent( i );
    double bg_fraction_error_i = h_bg_fraction->GetBinError( i );
    // if(bg_fraction_i < 0. || bg_fraction_i > 1.)
    //   std::cout<<"WARNING: bkgd fraction is out of range for bin "<<i<<": "<<bg_fraction_i<<std::endl;
    h_fraction_predicted->SetBinContent(i, 1. - bg_fraction_i );
    h_fraction_predicted->SetBinError(i, bg_fraction_error_i );
  }

  //Do the same for vertical errors 
  std::vector<std::string> vertNames = h_bg_fraction->GetVertErrorBandNames();
  for(unsigned int i=0; i<vertNames.size(); ++i )
  {
    std::vector<TH2D*> vert_hists;
    MnvVertErrorBand2D *errBand = h_bg_fraction->GetVertErrorBand(vertNames[i]);
    int nUniverses = errBand->GetNHists();
    for( int u = 0; u< nUniverses;  ++u)
    {
      TH2D* h_bg_fraction_u = errBand->GetHist( u );
      TH2D* h_vert = (TH2D*)h_fraction_predicted->Clone("h_vert");
      h_vert->Reset();
      for(int j =0; j<=highbin;++j)
      {
        double bg_vert_fraction_j       = h_bg_fraction_u->GetBinContent( j );
        double bg_vert_fraction_error_j = h_bg_fraction_u->GetBinError( j );
        // if(bg_vert_fraction_j < 0. || bg_vert_fraction_j > 1.)
        //   std::cout<<"WARNING: bkgd fraction is out of range for bin "<<j<<": "<<bg_vert_fraction_j<<std::endl;
        h_vert->SetBinContent(j, 1. - bg_vert_fraction_j );
        h_vert->SetBinError(j, bg_vert_fraction_error_j );
      }
      vert_hists.push_back( h_vert );
    }
    std::cout << "->Add Vertical Error ("<< vert_hists.size() << " Universes): " << vertNames[i].c_str() << std::endl;
    h_fraction_predicted->AddVertErrorBand( vertNames[i], vert_hists );
    
    //cleaning
    for(std::vector<TH2D*>::iterator itHist = vert_hists.begin(); itHist != vert_hists.end(); ++itHist)
      delete *itHist;
  }
  
  //Do the same for lateral errors
  std::vector<std::string> latNames = h_bg_fraction->GetLatErrorBandNames();
  for(unsigned int i=0; i<latNames.size(); ++i )
  {
    std::vector<TH2D*> lat_hists;
    MnvLatErrorBand2D *errBand = h_bg_fraction->GetLatErrorBand(latNames[i]);
    int nUniverses = errBand->GetNHists();
    for( int u = 0; u< nUniverses;  ++u)
    {
      TH2D* h_bg_fraction_u = errBand->GetHist( u );
      TH2D* h_lat = (TH2D*)h_fraction_predicted->Clone("h_lat");
      h_lat->Reset();
      for(int j =0; j<=highbin;++j)
      {
        double bg_lat_fraction_j       = h_bg_fraction_u->GetBinContent( j );
        double bg_lat_fraction_error_j = h_bg_fraction_u->GetBinError( j );
        // if(bg_lat_fraction_j < 0. || bg_lat_fraction_j > 1.)
        //   std::cout<<"WARNING: bkgd fraction is out of range for bin "<<j<<": "<<bg_lat_fraction_j<<std::endl;
        h_lat->SetBinContent(j, 1. - bg_lat_fraction_j );
        h_lat->SetBinError(j, bg_lat_fraction_error_j );
      }
      lat_hists.push_back( h_lat );
    }
    std::cout << "->Add Lateral Error ("<< lat_hists.size() << " Universes): " << latNames[i].c_str() << std::endl;
    h_fraction_predicted->AddLatErrorBand( latNames[i], lat_hists );
    
    //cleaning
    for(std::vector<TH2D*>::iterator itHist = lat_hists.begin(); itHist != lat_hists.end(); ++itHist)
      delete *itHist;
  }
  
  //Background subtracted distribution is = DATA x h_fraction_predicted
  h_data_nobck->Multiply( h_data_nobck, h_fraction_predicted );
  
  //cleaning
  delete h_bg_fraction;
  delete h_data;
  delete h_fraction_predicted;

}

void SubtractBackground(MnvH2D** h, MnvH2D* h_bg_weights, MnvH2D* &h_data_nobck, MnvH2D* &h_mc_nobck){
  MnvH2D* background_untuned = (MnvH2D*)h[kQELikeNot]->Clone( Form( "%s_untunedbck", h[kQELikeNot]->GetName()) );
  MnvH2D* background_tuned = (MnvH2D*)h[kQELikeNot]->Clone( Form( "%s_tunedbck", h[kQELikeNot]->GetName()) );
  h_mc_nobck   = (MnvH2D*)h[kMC]->Clone( Form( "%s_nobck", h[kMC]->GetName()) );
  h_data_nobck = (MnvH2D*)h[kData]->Clone( Form( "%s_nobck", h[kData]->GetName()) );
  h_data_nobck->AddMissingErrorBandsAndFillWithCV( *h_mc_nobck );//<---- make sure errors continue on!??
  background_tuned->Multiply(background_untuned,h_bg_weights);//This is the data driven background!
  h_data_nobck->Add(background_tuned,-1.0);
  //  h_data_nobck->Add( h[kQELikeNot], -1.0);
  h_mc_nobck->Add( h[kQELikeNot], -1.0);
}
  

void CorrectByEfficiency( MnvH2D *&h_effcor, MnvH2D* h_unfolded, MnvH2D** h_effhist )
{
  //create efficiency corrected histo
  h_effcor = (MnvH2D*)h_unfolded->Clone( Form( "%s_effcor", h_unfolded->GetName() ) );
 
  //efficiency doesn't have lateral error bands
  //so add to efficiency the lateral error bands found in unfolded data
  //and use the efficiency central value to fill the universes 
  h_effhist[kQELike]->AddMissingErrorBandsAndFillWithCV( *h_unfolded );
  //  h_effhist[kQELike]->Scale(0.2);

  //divide by efficiency 
  h_effcor->Divide( h_unfolded, h_effhist[kQELike] );

}


void CorrectByEfficiency( MnvH2D *&h_effcor, MnvH2D* h_unfolded, MnvH2D** h_effhist_num, MnvH2D** h_effhist_dem )
{
  //create efficiency corrected histo
  h_effcor = (MnvH2D*)h_unfolded->Clone( Form( "%s_effcor", h_unfolded->GetName() ) );


  MnvH2D *h_eff_cor_temp[nHistos];
  for( unsigned int j = 0; j < nHistos; ++j ){
    h_effhist_dem[j]->AddMissingErrorBandsAndFillWithCV( *h_unfolded );
    h_eff_cor_temp[j] = (MnvH2D*)h_effhist_num[j]->Clone(Form("temp_eff_cor_%s-%d",h_unfolded->GetName(),j));
    h_eff_cor_temp[j]->Divide(h_effhist_num[j],h_effhist_dem[j],1.0,1.0,"B");
  }
  //efficiency doesn't have lateral error bands
  //so add to efficiency the lateral error bands found in unfolded data 
  //and use the efficiency central value to fill the universes 
  //This is for the efficiency denominator which is NOT affected by this... So CV is fine.
  //h_eff_cor_temp[kQELike]->AddMissingErrorBandsAndFillWithCV( *h_unfolded );
    
  //  h_effhist[kQELike]->Scale(0.2);

  //divide by efficiency 
  h_effcor->Divide( h_unfolded, h_eff_cor_temp[kQELike] );


}

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


void NormalizeByFluxAndTargets( MnvH2D *&h_normalized, MnvH2D *h_effcor, MnvH1D* h_flux, double pot_scale, bool applyFluxConstraint, CCQENuUtils* utils, TMatrixD unfoldCov, MnvH2D *flux_norm )
{
  //-------------------------------------------------
  //! Normalize by Flux
  //--------------------------------------------------
  h_flux->AddMissingErrorBandsAndFillWithCV( *h_effcor->ProjectionX() );

  MnvH2D *h_flux_normalization = (MnvH2D*)h_effcor->Clone("h_flux_normalization");
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
  else h_flux_normalization=(MnvH2D*)flux_norm->Clone("enuflux");


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
  //const int nplanes = 2 * ( 75 - 30 + 1 );//fiducial volume -> modules 30-75 ...
  TString name = h_normalized->GetName();

  double targets = 0.0;
  if( name.Contains( "mc" ) ){
    //targets = TargetUtils::Get().GetTrackerNCarbonAtoms( nplanes, /*isMC =*/ true, 850 );
    targets = TargetUtils::Get().GetTrackerNNucleons( nplanes, /*isMC =*/ true, 850 );
  }
  else {
    //targets = TargetUtils::Get().GetTrackerNCarbonAtoms( nplanes, /*isMC =*/ false, 850 );
    targets = TargetUtils::Get().GetTrackerNNucleons( nplanes, /*isMC =*/ false, 850 );
  }

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
  h_normalized->PushCovMatrix("unfoldingCov",finalCovMatrix);
  finalCovMatrix.SaveAs("CovMat_tocheck.root");

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
      "\t./CrossSectionHists  Output_file  Input_file_1TRACK  Input_file_2TRACK  Input_file_BKG  Input_file_MIGRATION  Input_file_EFFICIENCY  Input_file_FLUX  Num_iterations  Make_Flux_Constraint_Histo  Apply_Flux_Constraint \n\n"<<
      "\t-output_file\t =\t Cross Section output file\n"<<
      "\t-input_file_1track\t =\t Muon Event Selection(Signal) for 1 track sub-sample\n"<<
      "\t-input_file_2track\t =\t Muon Event Selection(Signal) for >=2 track sub-sample\n"<<
      "\t-input_file_bkg\t =\t Background weights filename\n"<<
      "\t-input_file_migration\t =\t Migration Histogram\n"<<
      "\t-input_file_efficiency\t =\t Efficiency Histograms\n"<<
      "\t-num_iterations\t =\t No. of iterations\n"<< 
      "\t-Make_Flux_Constraint_Histo\t =\t If TRUE Enter 1; If FALSE Enter 0\n" << 
      "\t-Apply_Flux_Constraint\t =\t If TRUE Enter 1; If FALSE Enter 0\n" << 
      "\t********************************************************************************************** \n"<<
      "\t Please see : MuonSelectionHists.cxx, BackgroundWeights.cxx, MigrationMatrixHists.cxx, EffPurityHists.cxx and Ana/Flux/python/compute_flux.py for getting the necessary input files"<< std::endl; 
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    return 0; 
  }

  //! Default parameters
  std::vector<std::string> par;
  par.push_back("CrossSectionPlots");//0
  par.push_back( Form("%s/ana/rootfiles/CrossSectionPlots.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/MuonSelectionHists_Signal_2track.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/BackgroundWeights.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/MigrationHistogram.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/EffPurityHistogram.root",getenv("CCQENUROOT") ) );//5
  par.push_back("1"); 
  par.push_back("0"); 
  par.push_back("0"); 
  par.push_back("-1");//9

  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;

  bool fluxConstraintHistoNeeded = ( *(par.end()-2) == "1" ) ? true : false; 

  bool applyFluxConstraint = ( par.back() == "1" ) ? true : false; 


  string pzmuout = par[1]+"_pzmu.root"; 
  string q2out = par[1]+"_q2.root"; 
  string enuout = par[1]+"_enu.root";
  string short3out = par[1]+"_short3d.root";
  string big3out = par[1]+"_big3d.root";
  string dalphatout = par[1]+"_dalphat.root";
  string dphitout = par[1]+"_dphit.root";
  string pnout = par[1]+"_pn.root";
  string dptout = par[1]+"_dpt.root";
  string dptxout = par[1]+"_dptx.root";
  string dptyout = par[1]+"_dpty.root";
  string signedout = par[1]+"_signed.root";
  string signeddalphatout = par[1]+"_signeddalphat.root";
  string signeddphitout = par[1]+"_signeddphit.root";
  int run = atoi(par[9].c_str());

  if(run==1) CrossSectionHists( pzmuout, par[2], par[3], par[4], par[5], atoi(par[6].c_str()), fluxConstraintHistoNeeded, applyFluxConstraint,run,false);
  else if(run==2) CrossSectionHists( q2out, par[2], par[3], par[4], par[5], atoi(par[6].c_str()), fluxConstraintHistoNeeded, applyFluxConstraint,run,false);
  else if (run==3) CrossSectionHists( enuout, par[2], par[3], par[4], par[5], atoi(par[6].c_str()), fluxConstraintHistoNeeded, applyFluxConstraint,run,true);
  else if (run==4) CrossSectionHists( short3out, par[2], par[3], par[4], par[5], atoi(par[6].c_str()), fluxConstraintHistoNeeded, applyFluxConstraint,run,false);//short3d
  else if (run==5) CrossSectionHists( big3out, par[2], par[3], par[4], par[5], atoi(par[6].c_str()), fluxConstraintHistoNeeded, applyFluxConstraint,run,false);
  else if (run==6) CrossSectionHists( dalphatout, par[2], par[3], par[4], par[5], atoi(par[6].c_str()), fluxConstraintHistoNeeded, applyFluxConstraint,run,false);
  else if (run==7) CrossSectionHists( dphitout, par[2], par[3], par[4], par[5], atoi(par[6].c_str()), fluxConstraintHistoNeeded, applyFluxConstraint,run,false);
  else if (run==8) CrossSectionHists( pnout, par[2], par[3], par[4], par[5], atoi(par[6].c_str()), fluxConstraintHistoNeeded, applyFluxConstraint,run,false);
  else if (run==9) CrossSectionHists( dptout, par[2], par[3], par[4], par[5], atoi(par[6].c_str()), fluxConstraintHistoNeeded, applyFluxConstraint,run,false);
  else if (run==10) CrossSectionHists( dptxout, par[2], par[3], par[4], par[5], atoi(par[6].c_str()), fluxConstraintHistoNeeded, applyFluxConstraint,run,false);
  else if (run==11) CrossSectionHists( dptyout, par[2], par[3], par[4], par[5], atoi(par[6].c_str()), fluxConstraintHistoNeeded, applyFluxConstraint,run,false);
  else if (run==12) CrossSectionHists( signedout, par[2], par[3], par[4], par[5], atoi(par[6].c_str()), fluxConstraintHistoNeeded, applyFluxConstraint,run,false);
  else if (run==13) CrossSectionHists( signeddalphatout, par[2], par[3], par[4], par[5], atoi(par[6].c_str()), fluxConstraintHistoNeeded, applyFluxConstraint,run,false);
  else if (run==14) CrossSectionHists( signeddphitout, par[2], par[3], par[4], par[5], atoi(par[6].c_str()), fluxConstraintHistoNeeded, applyFluxConstraint,run,false);
  else cout << "Please provide a run type to run with current mapping run==1 = pzmu; run==2 = q2; run==3 = enu" << endl;

  return 0;

}
