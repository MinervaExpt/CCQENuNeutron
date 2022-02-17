#include "include/CCQENuUtils.h"
#include "PlotUtils/HyperDimLinearizer.h"
#include "PlotUtils/modifier_Eb.h"
#include "TParameter.h"


#include "include/Parameters.h"
#include "include/GeneralFunc.h" //in CCQENuNeutron/include/
#include "include/EventCuts.h"
#include "include/EventSelector.h"
#include "include/CommonBins.h"

using namespace CCQENU_ANA;


int testEvent = -1;
bool testMode = false;

void FillHisto( CCQENuUtils *utils, string namebase, double reco_xvar, double true_xvar, CCQENuEvent* &mc, double wgt, bool pass_cv=true )
{
  bool isData = false;
  if( pass_cv )
  {
    utils->fillHistosV3( utils->histos2D[namebase+"_migration"], reco_xvar, true_xvar, isData, mc, wgt );
    utils->fillHistosV3( utils->histos1D[namebase+"_reco"], reco_xvar, isData, mc, wgt );
    utils->fillHistosV3( utils->histos1D[namebase+"_truth"], true_xvar, isData, mc, wgt );
    if(RunCodeWithSystematics)
    {
      utils->fillVertErrorBands( utils->histos2D[namebase+"_migration"], reco_xvar, true_xvar, mc );
      utils->fillVertErrorBands( utils->histos1D[namebase+"_reco"], reco_xvar,  mc );
      utils->fillVertErrorBands( utils->histos1D[namebase+"_truth"], true_xvar,  mc );
    }
  }

  if(RunCodeWithSystematics)
  {
    utils->fillLatErrorBandsMigration1D( utils->histos2D[namebase+"_migration"], "q2", reco_xvar, true_xvar, mc ,pass_cv);
    utils->fillLatErrorBands( utils->histos1D[namebase+"_reco"], "q2", reco_xvar,  mc ,pass_cv);
    utils->fillLatErrorBands( utils->histos1D[namebase+"_truth"], "q2", true_xvar,  mc,pass_cv );
  }
}

void WriteHistograms(TFile* f, CCQENuUtils *utils, MnvH2D* &h_xvar_migration, MnvH1D* &h_xvar_reco, MnvH1D* &h_xvar_truth)
{
  f->cd();
  h_xvar_migration->Write();
  h_xvar_reco->Write();
  h_xvar_truth->Write();

  MnvH2D *h_X_xvar_migration[200]={NULL};
  utils->splitObject( h_xvar_migration, h_X_xvar_migration);
  for(int i=0;i<200;i++){
    cout << "Migration " << i << endl;
    if(h_X_xvar_migration[i]==NULL) break;
    MnvH1D* hreco = (MnvH1D*) h_X_xvar_migration[i]->ProjectionX( Form("%s_reco",h_X_xvar_migration[i]->GetName() ) );
    MnvH1D* htrue = (MnvH1D*) h_X_xvar_migration[i]->ProjectionY( Form("%s_truth",h_X_xvar_migration[i]->GetName() ) );

    h_X_xvar_migration[i]->Write();
    hreco->Write();
    htrue->Write();

    delete hreco;
    delete htrue;
    delete h_X_xvar_migration[i];
  }
}

void WriteHistograms(TFile* f, CCQENuUtils *utils, MnvH2D* &h_xvar_migration, MnvH2D* &h_xvar_reco, MnvH2D* &h_xvar_truth)
{

  f->cd();cout<<"f->cd()"<<endl;
  h_xvar_migration->Write();cout<<"h_xvar_migration->Write()"<<endl;
  h_xvar_reco->Write();cout<<"h_xvar_reco->Write();"<<endl;
  h_xvar_truth->Write();cout<<"h_xvar_truth->Write()"<<endl;
  cout<<"Integral: "<<h_xvar_migration->Integral()<<endl;
  //cout<<h_xvar_migration->GetName()<<endl;

  MnvH2D *h_X_xvar_migration[200]={NULL};
  utils->splitObject( h_xvar_migration, h_X_xvar_migration);
  for(int i=0;i<200;i++){
    if(h_X_xvar_migration[i]==NULL) break;
    cout << "Migration " << i << endl;

    h_X_xvar_migration[i]->Write();
  }
  for( int i = 0; i<200;i++) delete h_X_xvar_migration[i];
}


int MigrationMatrixHists( string filename, string playlist, bool makeFluxConstraintHisto, bool applyFluxConstraint, int multiplicity = 0, string sample = "hydrogen", int n_mcfiles = -1 , int npMode = 1, string xvarname = "Q2") 
{
  //---------------------------------------------
  // create file to store histograms
  //---------------------------------------------
  CCQENuUtils  *utils  = new CCQENuUtils( false, makeFluxConstraintHisto );
  utils->setPlaylist(playlist);
  utils->setFluxReweighterPlaylist();

  CCQENuCuts   *cutter = new CCQENuCuts();
  NeutronBlobCuts *ncutter = utils->GetNeutronCutter();

  PhysicsUtils* physicsUtils = new PhysicsUtils();
  GeoUtils *geoUtils = new GeoUtils();

  KinematicsModifier* kineMod = new KinematicsModifier();
  //--------------------------------------------
  // Xvars
  //--------------------------------------------
  axis_binning xbins = muonPtbins;
  string xvartitle = "ptmu";
  if (xvarname == "PT") xbins = muonPtbins;
  if (xvarname == "Q2") { xbins = Q2bins;xvartitle = "q2qe";}



  //---------------------------------------------
  // Book Histograms
  //---------------------------------------------
  //
  // 2-D Smearing histograms for 1-D unfolding
  // This is only used to check the smearing in Reco vs Truth from 1-D variables.
  MinervaUnfold::MnvResponse *response_h_angle_10 = NULL;
  MnvH2D *h_xvar_migration_10[nHistos];
  MnvH1D *h_xvar_reco_10[nHistos];
  MnvH1D *h_xvar_truth_10[nHistos];

  MinervaUnfold::MnvResponse *response_h_angle_01 = NULL;
  MnvH2D *h_xvar_migration_01[nHistos];
  MnvH1D *h_xvar_reco_01[nHistos];
  MnvH1D *h_xvar_truth_01[nHistos];

  MinervaUnfold::MnvResponse *response_h_angle_02 = NULL;
  MnvH2D *h_xvar_migration_02[nHistos];
  MnvH1D *h_xvar_reco_02[nHistos];
  MnvH1D *h_xvar_truth_02[nHistos];

  MinervaUnfold::MnvResponse *response_h_angle_03 = NULL;
  MnvH2D *h_xvar_migration_03[nHistos];
  MnvH1D *h_xvar_reco_03[nHistos];
  MnvH1D *h_xvar_truth_03[nHistos];

  MinervaUnfold::MnvResponse *response_h_angle_04 = NULL;
  MnvH2D *h_xvar_migration_04[nHistos];
  MnvH1D *h_xvar_reco_04[nHistos];
  MnvH1D *h_xvar_truth_04[nHistos];

  MinervaUnfold::MnvResponse *response_h_angle_all = NULL;
  MnvH2D *h_xvar_migration_all[nHistos];
  MnvH1D *h_xvar_reco_all[nHistos];
  MnvH1D *h_xvar_truth_all[nHistos];

  if(sample=="hydrogen") 
  {
    response_h_angle_10  = new MinervaUnfold::MnvResponse(("h_"+xvartitle+"_angle_10").c_str(), string("Q2").c_str(), xbins, xbins );
    response_h_angle_01  = new MinervaUnfold::MnvResponse(("h_"+xvartitle+"_angle_01").c_str(), string("Q2").c_str(), xbins, xbins );
    response_h_angle_02  = new MinervaUnfold::MnvResponse(("h_"+xvartitle+"_angle_02").c_str(), string("Q2").c_str(), xbins, xbins );
    response_h_angle_03  = new MinervaUnfold::MnvResponse(("h_"+xvartitle+"_angle_03").c_str(), string("Q2").c_str(), xbins, xbins );
    response_h_angle_04  = new MinervaUnfold::MnvResponse(("h_"+xvartitle+"_angle_04").c_str(), string("Q2").c_str(), xbins, xbins );
    response_h_angle_all  = new MinervaUnfold::MnvResponse(("h_"+xvartitle+"_angle_all").c_str(), string("Q2").c_str(), xbins, xbins );

    utils->bookHistos( h_xvar_migration_10, "h_q2qe_angle_10_migration", "reco vs true", xbins, xbins );
    utils->bookHistos( h_xvar_reco_10, "h_q2qe_angle_10_reco", "reco", xbins);
    utils->bookHistos( h_xvar_truth_10, "h_q2qe_angle_10_truth", "truth", xbins);

    utils->bookHistos( h_xvar_migration_01, "h_q2qe_angle_01_migration", "reco vs true", xbins, xbins );
    utils->bookHistos( h_xvar_reco_01, "h_q2qe_angle_01_reco", "reco", xbins);
    utils->bookHistos( h_xvar_truth_01, "h_q2qe_angle_01_truth", "truth", xbins);

    utils->bookHistos( h_xvar_migration_02, "h_q2qe_angle_02_migration", "reco vs true", xbins, xbins );
    utils->bookHistos( h_xvar_reco_02, "h_q2qe_angle_02_reco", "reco", xbins);
    utils->bookHistos( h_xvar_truth_02, "h_q2qe_angle_02_truth", "truth", xbins);

    utils->bookHistos( h_xvar_migration_03, "h_q2qe_angle_03_migration", "reco vs true", xbins, xbins );
    utils->bookHistos( h_xvar_reco_03, "h_q2qe_angle_03_reco", "reco", xbins);
    utils->bookHistos( h_xvar_truth_03, "h_q2qe_angle_03_truth", "truth", xbins);

    utils->bookHistos( h_xvar_migration_04, "h_q2qe_angle_04_migration", "reco vs true", xbins, xbins );
    utils->bookHistos( h_xvar_reco_04, "h_q2qe_angle_04_reco", "reco", xbins);
    utils->bookHistos( h_xvar_truth_04, "h_q2qe_angle_04_truth", "truth", xbins);

    utils->bookHistos( h_xvar_migration_all, "h_q2qe_angle_all_migration", "reco vs true", xbins, xbins );
    utils->bookHistos( h_xvar_reco_all, "h_q2qe_angle_all_reco", "reco", xbins);
    utils->bookHistos( h_xvar_truth_all, "h_q2qe_angle_all_truth", "truth", xbins);
  }

  MinervaUnfold::MnvResponse *response_dthetaP = NULL;
  MnvH2D *h_dthetaP_xvar_migration = NULL;
  MnvH2D *h_dthetaP_xvar_reco = NULL;
  MnvH2D *h_dthetaP_xvar_true = NULL;
  utils->setupResponse(response_dthetaP, "h_dthetaP_"+xvartitle, "Delta ThetaP vs "+xvarname,dthetaPerpbins,xbins,dthetaPerpbins,xbins );

  MinervaUnfold::MnvResponse *response_dthetaR= NULL;
  MnvH2D *h_dthetaR_xvar_migration = NULL;
  MnvH2D *h_dthetaR_xvar_reco = NULL;
  MnvH2D *h_dthetaR_xvar_true = NULL;
  utils->setupResponse(response_dthetaR, "h_dthetaR_"+xvartitle, "Delta ThetaP vs "+xvarname,dthetaReactbins,xbins,dthetaReactbins,xbins );

  MinervaUnfold::MnvResponse *response_muontheta= NULL;
  MnvH2D *h_muontheta_xvar_migration = NULL;
  MnvH2D *h_muontheta_xvar_reco = NULL;
  MnvH2D *h_muontheta_xvar_true = NULL;
  utils->setupResponse(response_muontheta, "h_muontheta_"+xvartitle, "Muon Theta vs "+xvarname,thetabins,xbins,thetabins,xbins );

  MinervaUnfold::MnvResponse *response_muonmomentum= NULL;
  MnvH2D *h_muonmomentum_xvar_migration = NULL;
  MnvH2D *h_muonmomentum_xvar_reco = NULL;
  MnvH2D *h_muonmomentum_xvar_true = NULL;
  utils->setupResponse(response_muonmomentum, "h_muonmomentum_"+xvartitle, "Muon Momentum vs "+xvarname,leptonmomentumbins,xbins,leptonmomentumbins,xbins );

  MinervaUnfold::MnvResponse *response_enu= NULL;
  MnvH2D *h_enu_xvar_migration = NULL;
  MnvH2D *h_enu_xvar_reco = NULL;
  MnvH2D *h_enu_xvar_true = NULL;
  utils->setupResponse(response_enu, "h_enu_"+xvartitle, "Enu vs "+xvarname,leptonmomentumbins,xbins,leptonmomentumbins,xbins );





  //  cout << response_q2 << endl;
  //---------------------------------------------
  // Add Error Bands
  //---------------------------------------------
  if(RunCodeWithSystematics){
    utils->addVertErrorBands( response_h_angle_10 );
    utils->addLatErrorBands(  response_h_angle_10 );

    utils->addVertErrorBands( response_h_angle_01 );
    utils->addLatErrorBands(  response_h_angle_01 );

    utils->addVertErrorBands( response_h_angle_02 );
    utils->addLatErrorBands(  response_h_angle_02 );
    
    utils->addVertErrorBands( response_h_angle_03 );
    utils->addLatErrorBands(  response_h_angle_03 );

    utils->addVertErrorBands( response_h_angle_04 );
    utils->addLatErrorBands(  response_h_angle_04 );

    utils->addVertErrorBands( response_h_angle_all );
    utils->addLatErrorBands(  response_h_angle_all );

    utils->addVertErrorBands( h_xvar_migration_10 );
    utils->addVertErrorBands( h_xvar_migration_01 );
    utils->addVertErrorBands( h_xvar_migration_02 );
    utils->addVertErrorBands( h_xvar_migration_03 );
    utils->addVertErrorBands( h_xvar_migration_04 );
    utils->addVertErrorBands( h_xvar_reco_10 );
    utils->addVertErrorBands( h_xvar_reco_01 );
    utils->addVertErrorBands( h_xvar_reco_02 );
    utils->addVertErrorBands( h_xvar_reco_03 );
    utils->addVertErrorBands( h_xvar_reco_04 );
    utils->addVertErrorBands( h_xvar_reco_all );
    utils->addVertErrorBands( h_xvar_truth_10 );
    utils->addVertErrorBands( h_xvar_truth_01 );
    utils->addVertErrorBands( h_xvar_truth_02 );
    utils->addVertErrorBands( h_xvar_truth_03 );
    utils->addVertErrorBands( h_xvar_truth_04 );
    utils->addVertErrorBands( h_xvar_truth_all );

    utils->addLatErrorBands( h_xvar_migration_10 );
    utils->addLatErrorBands( h_xvar_migration_01 );
    utils->addLatErrorBands( h_xvar_migration_02 );
    utils->addLatErrorBands( h_xvar_migration_03 );
    utils->addLatErrorBands( h_xvar_migration_04 );
    utils->addLatErrorBands( h_xvar_migration_all );
    utils->addLatErrorBands( h_xvar_reco_10 );
    utils->addLatErrorBands( h_xvar_reco_01 );
    utils->addLatErrorBands( h_xvar_reco_02 );
    utils->addLatErrorBands( h_xvar_reco_03 );
    utils->addLatErrorBands( h_xvar_reco_04 );
    utils->addLatErrorBands( h_xvar_reco_all );
    utils->addLatErrorBands( h_xvar_truth_10 );
    utils->addLatErrorBands( h_xvar_truth_01 );
    utils->addLatErrorBands( h_xvar_truth_02 );
    utils->addLatErrorBands( h_xvar_truth_03 );
    utils->addLatErrorBands( h_xvar_truth_04 );
    utils->addLatErrorBands( h_xvar_truth_all );

    utils->addVertErrorBands( response_dthetaP );
    utils->addLatErrorBands(  response_dthetaP );

    utils->addVertErrorBands( response_dthetaR );
    utils->addLatErrorBands(  response_dthetaR );


    utils->addVertErrorBands( response_muontheta );
    utils->addLatErrorBands(  response_muontheta );


    utils->addVertErrorBands( response_muonmomentum );
    utils->addLatErrorBands(  response_muonmomentum );

    utils->addVertErrorBands( response_enu );
    utils->addLatErrorBands(  response_enu );
  }//end run with sysetmatics




  //---------------------------------------------
  // Get MC Tree
  //---------------------------------------------
  //TChain *truth_mc = utils->getMCTree("Truth", n_mcfiles );
  //utils->setmnvHadronReweightTruthTree(truth_mc);

  TChain* tree_mc  = utils->getMCTree("CCQENu", n_mcfiles );
  int entries_mc   = tree_mc->GetEntries();

  cout << "MC entries: " << entries_mc << endl;
  CCQENuEvent* mc = new CCQENuEvent( tree_mc );
  utils->setmnvHadronReweightDataTree(tree_mc);

  //---------------------------------------------
  // Fill MC Histograms
  //---------------------------------------------

  bool isData = false;
  bool isMC = true;

  double mc_events = 0., data_events = 0.;
  map<int,int> mc_counter, data_counter;
  for (int i=0; i<entries_mc ; ++i) {
    mc->GetEntry(i);
    if( (i+1)%1000 == 0 ) 
      cout << "Reading MC entry : " << i+1 << " - " << 100*(i+1)/entries_mc << " % " << endl; 

    if( IgnoreEvent( mc, isData )) continue;
    if( !cutter->passTrueCCQELike( mc ) ) continue;
    if(testMode)
    {
      if( !cutter->passTrueCCQE( mc ) ) continue;
      if( mc->mc_targetZ != 1 ) continue;
    }
   
    //==============Reco Cuts======================

    if( !dropClusters ) ncutter->UpdateCandidates( mc, isMC );
    else ncutter->UpdateCandidates( mc, isMC , discardFraction, discardThreshold );

    bool PassEventCut = (EventCut( mc, cutter, ncutter, mc_counter, multiplicity, "Signal" ,  isData, npMode, useRecoilCut, orderPar ));

    if (!PassEventCut) continue;
    std::vector< XYZTVector > particles_XYZT;
    bool PassRecoSel = PassReco( mc, ncutter, kineMod, physicsUtils,geoUtils, isMC, npMode, particles_XYZT );
    if(!PassRecoSel) continue;
    std::vector< TLorentzVector > particles_lorentz = Convert<TLorentzVector, XYZTVector>( particles_XYZT );

    bool PassRecoilBlobCutLoose = PassRecoilBlobQsqCutsLoose(mc, cutter, ncutter );
    if(!PassRecoilBlobCutLoose) continue;

    bool pass_cv = PassRecoilBlobQsqCutsCV(mc, cutter, ncutter );

    std::vector<bool> commonNeutronCuts = CommonNeutronCuts(mc,ncutter,orderPar);
    bool hasNC = commonNeutronCuts[0] && commonNeutronCuts[1];
    bool isForward = commonNeutronCuts[2] && hasNC;
    bool isBackward = commonNeutronCuts[3] && hasNC;


    if(npMode == 10 && !hasNC ) continue;

    if(testMode) cout<<i<<endl;

    XYZTVector reco_neutrino_XYZT = particles_XYZT[0];
    XYZTVector reco_muon_XYZT = particles_XYZT[1];   XYZVector reco_muon_XYZ = (XYZVector) reco_muon_XYZT.Vect();
    XYZTVector reco_proton_XYZT = particles_XYZT[2]; XYZVector reco_proton_XYZ = (XYZVector) reco_proton_XYZT.Vect();
    XYZTVector reco_neutron_XYZT = particles_XYZT[3]; XYZVector reco_neutron_XYZ = (XYZVector) reco_neutron_XYZT.Vect();


    TVector3 reco_proton_V3 = Convert( reco_proton_XYZ );
    TVector3 reco_neutron_V3 = Convert( reco_neutron_XYZ );

    TVector3 primary_nucleon_V3 = (npMode ==0)? reco_neutron_V3 : reco_proton_V3;


    //Muon Variables
    //________________________________________________________________________________________________
    //________________________________________________________________________________________________
    //Reco________________________________________________
    double reco_muon_theta     = mc->muon_theta_allNodes[20]; // get theta from 20th node, note NX kludge fixes this...
    double reco_cos_muon_theta = cos(reco_muon_theta);


    double reco_muon_p         = reco_muon_XYZT.P();
    double reco_muon_T         = reco_muon_XYZT.E() - reco_muon_XYZT.M();
    double reco_muon_pt_beam = sin(reco_muon_theta)*reco_muon_p;
    double reco_muon_pz_beam = cos(reco_muon_theta)*reco_muon_p;
    reco_muon_theta*=180/TMath::Pi();

    //True________________________________________________
    vector<XYZTVector> true_particles_XYZT;
    bool passTrue = GetTrueParticles(  mc, utils, kineMod, physicsUtils, geoUtils, isMC, npMode, true_particles_XYZT );
    XYZTVector true_neutrino_XYZT = true_particles_XYZT[0];
    XYZTVector true_muon_XYZT = true_particles_XYZT[1];   XYZVector true_muon_XYZ = (XYZVector) true_muon_XYZT.Vect();
    XYZTVector true_proton_XYZT = true_particles_XYZT[2]; XYZVector true_proton_XYZ = (XYZVector) true_proton_XYZT.Vect();
    XYZTVector true_neutron_XYZT = true_particles_XYZT[3]; XYZVector true_neutron_XYZ = (XYZVector) true_neutron_XYZT.Vect();


    TVector3 true_proton_V3 = Convert( true_proton_XYZ );
    TVector3 true_neutron_V3 = Convert( true_neutron_XYZ );


    double true_muon_p         = true_muon_XYZT.P();
    double true_muon_theta     = utils->GetTheta( true_muon_XYZ.X(),true_muon_XYZ.Y(),true_muon_XYZ.Z() );
    double true_muon_pt_beam   = utils->GetTransverseMomentumWRTBeam( true_muon_XYZ.X(),true_muon_XYZ.Y(),true_muon_XYZ.Z() );
    double true_muon_pz_beam   = utils->GetLongitudinalMomentumWRTBeam( true_muon_XYZ.X(),true_muon_XYZ.Y(),true_muon_XYZ.Z() );
    true_muon_theta*=180/TMath::Pi();


    XYZTVector reco_vtx = XYZTVector( mc->CCQENu_vtx[0],mc->CCQENu_vtx[1],mc->CCQENu_vtx[2],mc->CCQENu_vtx[3]);
    XYZTVector true_vtx = XYZTVector( mc->mc_vtx[0],mc->mc_vtx[1],mc->mc_vtx[2],mc->mc_vtx[3]);

    int charge = (neutrinoMode)? -1:1;
    double binding_e = (npMode==0)? 0 : bindingE;
    double is_part_mass = (neutrinoMode)? 0.9395654133: 0.93827231; //neutron/proton
    double fs_part_mass = (neutrinoMode)? 0.93827231: 0.9395654133;


    double beam_bias = 0;
    XYZVector beam = geoUtils->BeamAxis( beam_bias ); 






    //Neutron variables
    //________________________________________________________________________________________________
    //________________________________________________________________________________________________

    XYZTVector expected_nucleon_XYZT = geoUtils->ComputeExpectedNucleon( beam, reco_muon_XYZT, is_part_mass, fs_part_mass, binding_e);
    XYZVector expected_nucleon_XYZ = (XYZVector) expected_nucleon_XYZT.Vect();

    XYZTVector expected_nucleon_XYZT_true = geoUtils->ComputeExpectedNucleon( beam, true_muon_XYZT, is_part_mass, fs_part_mass, binding_e);
    XYZVector expected_nucleon_XYZ_true = (XYZVector) expected_nucleon_XYZT.Vect();

    XYZVector blobVtxDir = ncutter->MainNeutronCandidate()->BlobVtxDir;
    XYZVector true_blobVtxDir = (XYZVector) (ncutter->MainNeutronCandidate()->MCTrackPos - true_vtx ).Vect();

    std::vector<double> reco_neutron_vars = geoUtils->ComputeNeutronAngularVars( beam, expected_nucleon_XYZ, blobVtxDir );
    std::vector<double> true_neutron_vars = geoUtils->ComputeNeutronAngularVars( beam, expected_nucleon_XYZ_true, true_blobVtxDir );


    double reco_q2qe = physicsUtils->qsqCCQE( beam, reco_muon_XYZT, charge, binding_e ); //GeV^2
    double true_q2qe = physicsUtils->qsqCCQE( beam, true_muon_XYZT, charge, binding_e ); //GeV^2
    double reco_enu = physicsUtils->nuEnergyCCQE( beam, reco_muon_XYZT, charge, binding_e ); //GeV^2
    double true_enu = physicsUtils->nuEnergyCCQE( beam, true_muon_XYZT, charge, binding_e ); //GeV^2

    double reco_dthetaP = reco_neutron_vars[vdthetaP];
    double reco_dthetaR = reco_neutron_vars[vdthetaR];
    double reco_dtheta = reco_neutron_vars[vdtheta];
    reco_dthetaP *= 180/TMath::Pi();
    reco_dthetaR *= 180/TMath::Pi();
    reco_dtheta *= 180/TMath::Pi();

    ncutter->SetAngularVar( reco_dthetaP, reco_dthetaR );


    double true_dthetaP = true_neutron_vars[vdthetaP];
    double true_dthetaR = true_neutron_vars[vdthetaR];
    double true_dtheta = true_neutron_vars[vdtheta];
    true_dthetaP *= 180/TMath::Pi();
    true_dthetaR *= 180/TMath::Pi();
    true_dtheta *= 180/TMath::Pi();



    //cout<<reco_dthetaR<<", "<<true_dthetaR<<endl;

    bool reco_01 = (abs(reco_dthetaP)<=1 && abs(reco_dthetaR)<=1);
    bool reco_02 = (abs(reco_dthetaP)<=2 && abs(reco_dthetaR)<=2);
    bool reco_03 = (abs(reco_dthetaP)<=3 && abs(reco_dthetaR)<=3);
    bool reco_04 = (abs(reco_dthetaP)<=4 && abs(reco_dthetaR)<=4);
    bool reco_10 = (abs(reco_dthetaP)<=10 && abs(reco_dthetaR)<=10);
    bool reco_15 = (abs(reco_dthetaP)<=15 && abs(reco_dthetaR)<=15);



    //Kinematics Variables
    //________________________________________________________________________________________________
    //________________________________________________________________________________________________
    double reco_xvar = -99999;
    double true_xvar = -99999;

    if( xvarname == "PT" )
    {
      reco_xvar = reco_muon_pt_beam;
      true_xvar = true_muon_pt_beam;
    }
    if( xvarname == "Q2" )
    {
      reco_xvar = reco_q2qe;
      true_xvar = true_q2qe;
    }


    //Get cvweight

    TVector3 reco_primary_nucleon_V3 = (npMode ==0)? reco_neutron_V3.Unit()*expected_nucleon_XYZ.R() : reco_proton_V3;
    utils->setTransverseParticle( reco_primary_nucleon_V3, 0, false );
    double wgt = utils->GetCVWeight( mc );
    mc_events+=wgt;
    
    //Filling Histos

    //cout<<wgt<<endl;


    if (pass_cv)
    {
      if(reco_01) FillHisto( utils, "h_q2qe_angle_01",reco_xvar, true_xvar, mc, wgt );
      if(reco_02) FillHisto( utils, "h_q2qe_angle_02",reco_xvar, true_xvar, mc, wgt );
      if(reco_03) FillHisto( utils, "h_q2qe_angle_03",reco_xvar, true_xvar, mc, wgt );
      if(reco_04) FillHisto( utils, "h_q2qe_angle_04",reco_xvar, true_xvar, mc, wgt );
      FillHisto( utils, "h_q2qe_angle_all",reco_xvar, true_xvar, mc, wgt );

      bool fillMuon = (mc->mc_targetZ == 1) && reco_10 && cutter->passTrueCCQE(mc);
      if( cutter->passTrueCCQELike( mc ) && passQ2( reco_xvar) )
      {
        if(fillMuon) utils->fillResponse( response_muontheta, reco_muon_theta, reco_xvar, true_muon_theta, true_xvar, mc,wgt );
        if(fillMuon) utils->fillResponse( response_muonmomentum, reco_muon_p, reco_xvar, true_muon_p, true_xvar, mc,wgt );
        if(fillMuon) utils->fillResponse( response_enu, reco_enu, reco_xvar, true_enu, true_xvar, mc,wgt );
        utils->fillResponse( response_dthetaR, reco_dthetaR, reco_xvar, true_dthetaR, true_xvar, mc,wgt );
        utils->fillResponse( response_dthetaP, reco_dthetaP, reco_xvar, true_dthetaP, true_xvar, mc,wgt );

        if(RunCodeWithSystematics)
        {
          if(fillMuon) utils->fillVertErrorBands( response_muontheta, reco_muon_theta, reco_xvar, true_muon_theta, true_xvar, mc );
          if(fillMuon) utils->fillVertErrorBands( response_muonmomentum, reco_muon_p, reco_xvar, true_muon_p, true_xvar, mc );
          if(fillMuon) utils->fillVertErrorBands( response_enu, reco_enu, reco_xvar, true_enu, true_xvar, mc );
          utils->fillVertErrorBands( response_dthetaR, reco_dthetaR, reco_xvar, true_dthetaR, true_xvar, mc );
          utils->fillVertErrorBands( response_dthetaP, reco_dthetaP, reco_xvar, true_dthetaP, true_xvar, mc );
        }
      }
    }

    ncutter->SetLateralRegion(0);
    ncutter->SetAngularVar( reco_dthetaP, reco_dthetaR );
    double wgt_lat = utils->GetCVWeight( mc );
    FillHisto( utils, "h_q2qe_angle_10",reco_xvar, true_xvar, mc, wgt_lat, reco_10&&pass_cv );
    ncutter->UnsetFillLateralMode();


    //



    //Filling for dthetaR_q2qe etc
  } //End of for loop over MC entries 
 
  cout << "Finished looping over all MC entries " << endl; 
  delete mc; 
  delete tree_mc; 

  //-----------------------------------------------
  // Getting h_migration, h_reco, h_truth
  //-----------------------------------------------
  //Get Data POT for this playlist
  //Get data POT for particular playlist
  TChain* tree_data  = utils->getDataTree( "CCQENu", -1 );
  delete tree_data;


  //==================================================================
  // Create ROOT file to store histograms
  // Write to file all the created histograms 
  //==================================================================
  TFile *f = new TFile( filename.c_str(), "RECREATE" );  
  cout << "write POT" << endl;
  //Write MC and DATA POT to file 
  utils->writePOT( f );
  cout << "write Multiplicity" << endl;
  //Write multiplicity and sample definitions
  TMap *sampleMap = new TMap(2,0);
  string multiplicity_label = (multiplicity==0)? "alltracks": Form("%i_track",multiplicity);
  sampleMap->Add(new TObjString("multiplicity"), new TObjString(multiplicity_label.c_str()) );
  f->WriteTObject( sampleMap, "sample");
  cout << "write n_events" << endl;
  //Write DATA and MC number of events
  TVector2 *evt = new TVector2( data_events, mc_events );
  f->WriteTObject( evt, "n_events" );
    cout << "write fluxParam" << endl;
  //Write out if output ROOT file has the constrained flux or not 
  TParameter<bool> *fluxParam = new TParameter<bool>( "appliedFluxConstraint", applyFluxConstraint ); 
  f->WriteTObject( fluxParam ); 

  f->cd();
  cout << "write hists q2qe hydrogen" << endl;
  WriteHistos( utils,f);

  cout<<"Write others"<<endl;

  utils->getResponseObjects( response_dthetaR, h_dthetaR_xvar_migration, h_dthetaR_xvar_reco, h_dthetaR_xvar_true);
  utils->getResponseObjects( response_dthetaP, h_dthetaP_xvar_migration, h_dthetaP_xvar_reco, h_dthetaP_xvar_true);
  utils->getResponseObjects( response_muontheta, h_muontheta_xvar_migration, h_muontheta_xvar_reco, h_muontheta_xvar_true);
  utils->getResponseObjects( response_muonmomentum, h_muonmomentum_xvar_migration, h_muonmomentum_xvar_reco, h_muonmomentum_xvar_true);
  utils->getResponseObjects( response_enu, h_enu_xvar_migration, h_enu_xvar_reco, h_enu_xvar_true);
  cout<<h_muontheta_xvar_migration->GetName()<<endl;

  WriteHistograms( f, utils, h_muontheta_xvar_migration, h_muontheta_xvar_reco, h_muontheta_xvar_true );
  WriteHistograms( f, utils, h_muonmomentum_xvar_migration, h_muonmomentum_xvar_reco, h_muonmomentum_xvar_true );
  WriteHistograms( f, utils, h_enu_xvar_migration, h_enu_xvar_reco, h_enu_xvar_true );
  WriteHistograms( f, utils, h_dthetaR_xvar_migration, h_dthetaR_xvar_reco, h_dthetaR_xvar_true );
  WriteHistograms( f, utils, h_dthetaP_xvar_migration, h_dthetaP_xvar_reco, h_dthetaP_xvar_true );

  cout<<"Done"<<endl;


  
  f->Close();
  delete f;
  delete utils;
  delete cutter; 
  delete binner;
  delete minmodbinner; 

  return 0;
}

int main( int argc, char *argv[])
{
  ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);
  if (argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t-./MigrationMatrixHists  Name_and_path_to_Output_file  Make_Flux_Constraint_Histo  Apply_Flux_Constraint  Multiplicity  Number_MC_files \n\n"<<
      "\t-Name_and_path_to_Output_file\t =\t Name of and path to the Output ROOT file that will be created \n"<<
      "\t-Playlist\t =\t Name of the playlist you want to run over e.g. minerva1 \n"<<
      "\t-Make_Flux_Constraint_Histo\t =\t If TRUE Enter 1; If FALSE Enter 0 \n"<< 
      "\t-Apply_Flux_Constraint\t =\t If TRUE Enter 1; If FALSE Enter 0 \n"<< 
      "\t-Multiplicity\t =\t Enter 0, 1 or 2. Here 0: 1 or 2 tracks; 1: 1-track events only; 2: 2-track events only \n"<<
      "\t-Sample to run\t =\t Enter pzmu,enu,q2,transverseVar,neutronVars \n"<<
      "\t-Number_MC_files\t =\t Number of MonteCarlo files. To use all files, set this to -1\n" <<
      "\t-xType \n"<<
      "\t-yType \n"<<
      "\t-npMode \n"<<
      "\t-Dependent Var \n" << std::endl; 
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    return 0; 
  }
  
  //! Default parameters
  std::vector<std::string> par;
  par.push_back("MigrationMatrixHists");
  par.push_back( Form("%s/ana/rootfiles/MigrationMatrixHists.root",getenv("CCQENUROOT") ) );
  par.push_back("minerva1");
  par.push_back("0");//flux contrain histo needed
  par.push_back("0");//apply flux constrain
  par.push_back("1");//multiplicity
  par.push_back("hydrogen");//sample
  par.push_back("-1");//n_mcfiles
  par.push_back("1");//proton, neutron,proton/neutron mode
  par.push_back("PT");//dependent var
  par.push_back("0");//drop clusters?

  
  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }
  

  if(par.at(9)!="PT" &&
     par.at(9)!="Q2"  )
  {
    cout << "Invalid xvar choice. Please choose PT or Q2. You used: "<<par.at(11) << endl;
    return 1;
  }
  if(par.at(6)!="hydrogen" ){
    cout << "Invalid sample choice. Please choose pzmu, enu, or q2. You used: "<<par.at(6) << endl;
    return 1;
  }



  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;
  
  string filename = par.at(1);
  string playlist = par.at(2);
  bool fluxConstraintHistoNeeded = ( par.at(3) == "1" ) ? true : false; 
  bool applyFluxConstraint = ( par.at(4) == "1" ) ? true : false; 
  int multiplicity = atoi(par[5].c_str());
  string var = par[6];
  int n_mcfiles = atoi(par[7].c_str());
  int npMode = atoi(par[8].c_str());
  string xvarname = par[9];
  dropClusters = bool(atoi(par[10].c_str() ) );


  return MigrationMatrixHists( filename, playlist, fluxConstraintHistoNeeded, applyFluxConstraint, multiplicity, var, n_mcfiles, npMode, xvarname );
}
