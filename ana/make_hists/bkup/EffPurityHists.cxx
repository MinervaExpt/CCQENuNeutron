#include "include/CCQENuUtils.h"
#include "PlotUtils/HyperDimLinearizer.h"
#include "TParameter.h" 

using namespace CCQENU_ANA;
using namespace Neutron_ANA;

void FillEfficiencyHisto( const MnvH1D* num, MnvH1D* den, MnvH1D* eff );
void FillEfficiencyHisto( const MnvH2D* num, MnvH2D* den, MnvH2D* eff );

int EffPurityHists(string filename, string playlist, bool makeFluxConstraintHisto, bool applyFluxConstraint, int n_mcfiles = -1, int xtype = -1, int ytype = -1)
{
  //---------------------------------------------
  // create file to store histograms
  //---------------------------------------------
  //TFile *f = new TFile( filename.c_str(), "RECREATE" );
  bool typeCheck = xtype==-1 && ytype==-1;//Do normal
  CCQENuUtils  *utils          = new CCQENuUtils( false, makeFluxConstraintHisto );
  utils->setPlaylist(playlist);
  utils->setFluxReweighterPlaylist();
  gStyle->SetNumberContours(500);
  CCQENuCuts   *cutter         = new CCQENuCuts();
  AnaBinning *binner                  = new AnaBinning();
  CCQENuBinning *minmodbinner  = new CCQENuBinning();

  NeutronBlobBinning *neutbinner = new NeutronBlobBinning();
  GeoUtils *geoUtils = new GeoUtils();

  axis_binning q2bins           = binner->GetQ2BinsGeV();

  //! Picking this option for 'theta' binning for J. Wolcott's plots (June 2015)  
  axis_binning thetabins        = minmodbinner->GetMuonThetaBins();
  //axis_binning thetabins        = binner->GetMuonAngleBinsDegrees();

  //! Picking this option for 'tmu' binning for J. Wolcott's plots (June 2015) 
  axis_binning tmubins          = binner->GetMuonEnergyBinsGeV();
  //axis_binning tmubins          = binner->GetMuonEnergyUniformBinsGeV();
  
  axis_binning muonpTbins         = minmodbinner->GetMuonPtBinsGeV();
  axis_binning muonPzbins         = minmodbinner->GetMuonPzBinsGeV();

  axis_binning Q2bins         = minmodbinner->GetQ2BinsGeV();
  axis_binning enubins         = minmodbinner->GetMuonPzBinsGeV();
  
  //! This option is used for the efficiency study in each pZ-pT bin  
  //axis_binning muonpTbins         = minmodbinner->GetMuonPt_FinerBinsGeV();
  //axis_binning muonPzbins         = minmodbinner->GetMuonPz_FinerBinsGeV();
  
  //! Reminder: "nu" here signifies the recoil (non-leptonic) energy 
  axis_binning nubins           = minmodbinner->GetNuBinsGeV(); 
  axis_binning muonPbins        = binner->GetMuonEnergyBinsGeV();

  axis_binning dalphatbins    = minmodbinner->GetDalphaT();
  axis_binning dphitbins    = minmodbinner->GetDphiT();
  axis_binning pnbins    = minmodbinner->GetPn();
  axis_binning dptbins    = minmodbinner->GetDeltaPt();
  axis_binning dptxbins    = minmodbinner->GetDeltaPtx();
  axis_binning dptybins    = minmodbinner->GetDeltaPty();
  axis_binning protonThetabins = minmodbinner->GetProtonThetaBins();
  axis_binning protonKineticbins = minmodbinner->GetProtonEnergyBinsGeV();  
  axis_binning signbins          = minmodbinner->GetSignedBins();
  axis_binning signeddalphatbins = minmodbinner->GetSignedDeltaAlphaT();
  axis_binning signeddphitbins = minmodbinner->GetSignedDeltaPhiT();

  axis_binning dalphatbins_TKI    = minmodbinner->GetDalphaT_TKI();
  axis_binning dphitbins_TKI    = minmodbinner->GetDphiT_TKI();
  axis_binning pnbins_TKI    = minmodbinner->GetPn_TKI();
  axis_binning dptbins_TKI    = minmodbinner->GetDeltaPt_TKI();
  axis_binning dptxbins_TKI    = minmodbinner->GetDeltaPtx_TKI();
  axis_binning dptybins_TKI    = minmodbinner->GetDeltaPty_TKI();
  axis_binning protonThetabins_TKI = minmodbinner->GetProtonThetaBins_TKI();
  axis_binning protonKineticbins_TKI = minmodbinner->GetProtonEnergyBinsGeV_TKI();
  axis_binning signeddalphatbins_TKI = minmodbinner->GetSignedDeltaAlphaT_TKI(); 
  axis_binning signeddphitbins_TKI = minmodbinner->GetSignedDeltaPhiT_TKI();


  axis_binning dthetaPerpbins = neutbinner->GetThetaPerpBinsDegree();
  axis_binning dthetaReactbins = neutbinner->GetThetaReactionPlaneBinsDegree();


  //------------------------------------------------------------------------------------------
  // Book 2D Efficiency Histogram for PZ vs PT (for the efficiency correction) 
  // Note: These are "NOT for cut-after-cut", it is after all the cuts for "Signal" events  
  //------------------------------------------------------------------------------------------
  MnvH2D *h_dalphat_ptmu_truth[nHistos],*h_dphit_ptmu_truth[nHistos],*h_pn_ptmu_truth[nHistos],*h_dpt_ptmu_truth[nHistos],*h_dptx_ptmu_truth[nHistos],*h_dpty_ptmu_truth[nHistos],*h_tp_ptmu_truth[nHistos],*h_ptheta_ptmu_truth[nHistos],*h_signed_ptmu_truth[nHistos],*h_signeddalphat_ptmu_truth[nHistos],*h_signeddphit_ptmu_truth[nHistos];
  MnvH2D *h_dalphat_ptmu[nHistos],*h_dphit_ptmu[nHistos],*h_pn_ptmu[nHistos],*h_dpt_ptmu[nHistos],*h_dptx_ptmu[nHistos],*h_dpty_ptmu[nHistos],*h_tp_ptmu[nHistos],*h_ptheta_ptmu[nHistos],*h_signed_ptmu[nHistos],*h_signeddalphat_ptmu[nHistos],*h_signeddphit_ptmu[nHistos];

  MnvH2D *h_dthetaR_ptmu[nHistos], *h_dthetaP_ptmu[nHistos], *h_dPP_ptmu[nHistos], *h_dPR_ptmu[nHistos], *h_dPPi_ptmu[nHistos], *h_dPRi_ptmu[nHistos];
  MnvH2D *h_dthetaR_ptmu_truth[nHistos], *h_dthetaP_ptmu_truth[nHistos], *h_dPP_ptmu_truth[nHistos], *h_dPR_ptmu_truth[nHistos], *h_dPPi_ptmu_truth[nHistos], *h_dPRi_ptmu_truth[nHistos];
          
  utils->bookHistos( h_dalphat_ptmu_truth, "h_dalphat_ptmu_truth","dat:ptmu",dalphatbins,muonpTbins);
  utils->bookHistos( h_dphit_ptmu_truth, "h_dphit_ptmu_truth","dphit:ptmu",dphitbins,muonpTbins);
  utils->bookHistos( h_pn_ptmu_truth, "h_pn_ptmu_truth","pn:ptmu",pnbins,muonpTbins);
  utils->bookHistos( h_dpt_ptmu_truth, "h_dpt_ptmu_truth","dpt:ptmu",dptbins,muonpTbins);
  utils->bookHistos( h_dptx_ptmu_truth, "h_dptx_ptmu_truth","dptx:ptmu",dptxbins,muonpTbins);
  utils->bookHistos( h_dpty_ptmu_truth, "h_dpty_ptmu_truth","dpty:ptmu",dptybins,muonpTbins);
  utils->bookHistos( h_tp_ptmu_truth, "h_tp_ptmu_truth","tp:ptmu",protonKineticbins,muonpTbins);
  utils->bookHistos( h_ptheta_ptmu_truth, "h_ptheta_ptmu_truth","ptheta:ptmu",protonThetabins,muonpTbins);
  utils->bookHistos( h_signed_ptmu_truth, "h_signed_ptmu_truth","signed:ptmu",signbins,muonpTbins);
  utils->bookHistos( h_signeddalphat_ptmu_truth, "h_signeddalphat_ptmu_truth","signeddalphat:ptmu",signeddalphatbins,muonpTbins);
  utils->bookHistos( h_signeddphit_ptmu_truth, "h_signeddphit_ptmu_truth","signeddphit:ptmu",signeddphitbins,muonpTbins);

  utils->bookHistos( h_dthetaR_ptmu_truth, "h_dthetaR_ptmu_truth", "dthetaR:ptmu", dthetaReactbins, muonpTbins);
  utils->bookHistos( h_dthetaP_ptmu_truth, "h_dthetaP_ptmu_truth", "dthetaP:ptmu", dthetaPerpbins, muonpTbins);
  utils->bookHistos( h_dPR_ptmu_truth, "h_dPR_ptmu_truth", "dPR:ptmu", dptybins, muonpTbins);
  utils->bookHistos( h_dPP_ptmu_truth, "h_dPP_ptmu_truth", "dPP:ptmu", dptxbins, muonpTbins);
  utils->bookHistos( h_dPRi_ptmu_truth, "h_dPRi_ptmu_truth", "dPRi:ptmu", dptybins, muonpTbins);
  utils->bookHistos( h_dPPi_ptmu_truth, "h_dPPi_ptmu_truth", "dPPi:ptmu", dptxbins, muonpTbins);



  utils->bookHistos( h_dalphat_ptmu, "h_dalphat_ptmu","dat:ptmu",dalphatbins,muonpTbins);
  utils->bookHistos( h_dphit_ptmu, "h_dphit_ptmu","dphit:ptmu",dphitbins,muonpTbins);
  utils->bookHistos( h_pn_ptmu, "h_pn_ptmu","pn:ptmu",pnbins,muonpTbins);
  utils->bookHistos( h_dpt_ptmu, "h_dpt_ptmu","dpt:ptmu",dptbins,muonpTbins);
  utils->bookHistos( h_dptx_ptmu, "h_dptx_ptmu","dptx:ptmu",dptxbins,muonpTbins);
  utils->bookHistos( h_dpty_ptmu, "h_dpty_ptmu","dpty:ptmu",dptybins,muonpTbins);
  utils->bookHistos( h_tp_ptmu, "h_tp_ptmu","tp:ptmu",protonKineticbins,muonpTbins);
  utils->bookHistos( h_ptheta_ptmu, "h_ptheta_ptmu","ptheta:ptmu",protonThetabins,muonpTbins);
  utils->bookHistos( h_signed_ptmu, "h_signed_ptmu","signed:ptmu",signbins,muonpTbins);
  utils->bookHistos( h_signeddalphat_ptmu, "h_signeddalphat_ptmu","signeddalphat:ptmu",signeddalphatbins,muonpTbins);
  utils->bookHistos( h_signeddphit_ptmu, "h_signeddphit_ptmu","signeddphit:ptmu",signeddphitbins,muonpTbins);

  utils->bookHistos( h_dthetaR_ptmu, "h_dthetaR_ptmu", "dthetaR:ptmu", dthetaReactbins, muonpTbins);
  utils->bookHistos( h_dthetaP_ptmu, "h_dthetaP_ptmu", "dthetaP:ptmu", dthetaPerpbins, muonpTbins);
  utils->bookHistos( h_dPR_ptmu, "h_dPR_ptmu", "dPR:ptmu", dptybins, muonpTbins);
  utils->bookHistos( h_dPP_ptmu, "h_dPP_ptmu", "dPP:ptmu", dptxbins, muonpTbins);
  utils->bookHistos( h_dPRi_ptmu, "h_dPRi_ptmu", "dPRi:ptmu", dptybins, muonpTbins);
  utils->bookHistos( h_dPPi_ptmu, "h_dPPi_ptmu", "dPPi:ptmu", dptxbins, muonpTbins);




  //special case
 //Special 3D for TKI vs TKI
  string spec_title = "";
  string spec_title_truth = "";
  vector<vector<double> > full3D;
  //Now to pick the variables
  map<int,axis_binning > axis_bin_map;
  map<int,string > axis_name;
  axis_bin_map[0] = muonPzbins;
  axis_bin_map[1] = dalphatbins_TKI;
  axis_bin_map[2] = dphitbins_TKI;
  axis_bin_map[3] = pnbins_TKI;
  axis_bin_map[4] = dptbins_TKI;
  axis_bin_map[5] = dptxbins_TKI;
  axis_bin_map[6] = dptybins_TKI;
  axis_bin_map[7] = protonKineticbins_TKI;
  axis_bin_map[8] = protonThetabins_TKI;
  axis_bin_map[9] = signbins;
  axis_bin_map[10] = signeddalphatbins_TKI;
  axis_bin_map[11] = signeddphitbins_TKI;

  axis_name[0] = "muonPz";
  axis_name[1] = "dalphat";
  axis_name[2] = "dphit";
  axis_name[3] = "pn";
  axis_name[4] = "dpt";
  axis_name[5] = "dptx";
  axis_name[6] = "dpty";
  axis_name[7] = "tp";
  axis_name[8] = "prottheta";
  axis_name[9] = "sign";
  axis_name[10] = "signeddalphat";
  axis_name[11] = "signeddphit";

  MnvH2D *h_tki_tki[nHistos],*h_tki_tki_truth[nHistos];
  if(!typeCheck){
    spec_title = "h_"+axis_name[xtype]+"_"+axis_name[ytype];
    spec_title_truth = "h_"+axis_name[xtype]+"_"+axis_name[ytype]+"_truth";
    cout << "I will be running " << spec_title << endl;
    utils->bookHistos(h_tki_tki,spec_title.c_str(),"tki:tki",axis_bin_map[xtype],axis_bin_map[ytype]);
    utils->bookHistos(h_tki_tki_truth,spec_title_truth.c_str(),"tki:tki",axis_bin_map[xtype],axis_bin_map[ytype]);
  }







  //---------------------------------------------------------------------
  // Add Vertical Error Bands for the 2D PZ vs PT distributions 
  // They will be used for the efficiency correction 
  //---------------------------------------------------------------------
  if(RunCodeWithSystematics){
    if(typeCheck){
      utils->addVertErrorBands( h_dalphat_ptmu );
      utils->addVertErrorBands( h_dalphat_ptmu_truth );
      
      utils->addVertErrorBands( h_dphit_ptmu );
      utils->addVertErrorBands( h_dphit_ptmu_truth );
      
      utils->addVertErrorBands( h_pn_ptmu );
      utils->addVertErrorBands( h_pn_ptmu_truth );
      
      
      utils->addVertErrorBands( h_dpt_ptmu );
      utils->addVertErrorBands( h_dpt_ptmu_truth );
      
      utils->addVertErrorBands( h_dptx_ptmu );
      utils->addVertErrorBands( h_dptx_ptmu_truth );
      
      utils->addVertErrorBands( h_dpty_ptmu );
      utils->addVertErrorBands( h_dpty_ptmu_truth );
      
      utils->addVertErrorBands( h_tp_ptmu );
      utils->addVertErrorBands( h_tp_ptmu_truth );
      
      utils->addVertErrorBands( h_ptheta_ptmu );
      utils->addVertErrorBands( h_ptheta_ptmu_truth );
      
      utils->addVertErrorBands( h_signed_ptmu );
      utils->addVertErrorBands( h_signed_ptmu_truth );
      
      utils->addVertErrorBands( h_signeddalphat_ptmu );
      utils->addVertErrorBands( h_signeddalphat_ptmu_truth );
      
      utils->addVertErrorBands( h_signeddphit_ptmu );
      utils->addVertErrorBands( h_signeddphit_ptmu_truth );

      utils->addVertErrorBands( h_dthetaR_ptmu );
      utils->addVertErrorBands( h_dthetaP_ptmu );
      utils->addVertErrorBands( h_dPP_ptmu );
      utils->addVertErrorBands( h_dPR_ptmu );
      utils->addVertErrorBands( h_dPPi_ptmu );
      utils->addVertErrorBands( h_dPRi_ptmu );
      utils->addVertErrorBands( h_dthetaR_ptmu_truth );
      utils->addVertErrorBands( h_dthetaP_ptmu_truth );
      utils->addVertErrorBands( h_dPP_ptmu_truth );
      utils->addVertErrorBands( h_dPR_ptmu_truth );
      utils->addVertErrorBands( h_dPPi_ptmu_truth );
      utils->addVertErrorBands( h_dPRi_ptmu_truth );
 
 


      utils->addLatErrorBands( h_dalphat_ptmu );
      utils->addLatErrorBands( h_dphit_ptmu );
      utils->addLatErrorBands( h_pn_ptmu );
      utils->addLatErrorBands( h_dpt_ptmu );
      utils->addLatErrorBands( h_dptx_ptmu );
      utils->addLatErrorBands( h_dpty_ptmu );
      utils->addLatErrorBands( h_tp_ptmu );
      utils->addLatErrorBands( h_ptheta_ptmu );
      utils->addLatErrorBands( h_signed_ptmu );
      utils->addLatErrorBands( h_signeddalphat_ptmu );
      utils->addLatErrorBands( h_signeddphit_ptmu );
      
      utils->addLatErrorBands( h_dthetaR_ptmu );
      utils->addLatErrorBands( h_dthetaP_ptmu );
      utils->addLatErrorBands( h_dPP_ptmu );
      utils->addLatErrorBands( h_dPR_ptmu );
      utils->addLatErrorBands( h_dPPi_ptmu );
      utils->addLatErrorBands( h_dPRi_ptmu );



    }
    else{
      utils->addVertErrorBands( h_tki_tki );
      utils->addVertErrorBands( h_tki_tki_truth );

      utils->addLatErrorBands( h_tki_tki );
    }

  }//end if run with systematics
  //---------------------------------------------
  // Get Truth Tree
  //---------------------------------------------
  TChain* tree_truth  = utils->getMCTree( "Truth", n_mcfiles );
  int entries_truth = tree_truth->GetEntries();
  cout << "Truth entries: " << entries_truth << endl;
  utils->setmnvHadronReweightTruthTree(tree_truth);

  double truth_events = 0.;
  CCQENuTruth* truth = new CCQENuTruth( tree_truth );
  //-------------------------------------------
  // Fill Truth histograms 
  //-------------------------------------------
  cout << "********************************************" << endl;
  cout << "********************************************" << endl;
  cout << "WILL BE CUTTING ON TRUE ANGLE AND REMOVING FROM THE DENOMINATOR ALL TRACKS WITH THETA GREATER THAN 20 DEGREES. COMMENT OUT THE ANGLE CUT IN THE CODE IF YOU DON'T DO AN ANGLE CUT IN THE NUMERATOR. SEE passCCQESelection TO MAKE SURE" << endl;
  cout << "********************************************" << endl;
  cout << "********************************************" << endl;
  for (int i=0; i<entries_truth ; ++i) {
    if( (i+1)%2000000 == 0 )
      cout << "Reading Truth entry : " << i+1 << " - " << int(100.0*(i+1)/entries_truth) << " %" << endl;
    
    truth->GetEntry(i);


    //cout << "True event true_nu " << true_nu << endl;  
    if( !truth->truth_is_fiducial ) continue;

    
    //Get CVWeight
    double wgt = utils->GetCVWeight( truth );

    //Efficiency denominator modification due to POT counting but. POT_Used/POT_Total is an additional weight needed to bring the "exposure" to the same value.

    double potwgt = utils->GetPOTWeight();//global_used_pot_mc/global_total_pot_mc
    wgt*=potwgt;

    const double bindingE = 34.0;


    TVector3 true_proton_p_comp = utils->GetHighestTrueProtonMomentumComp(truth);
    double true_proton_px = true_proton_p_comp.X();
    double true_proton_py = true_proton_p_comp.Y();
    double true_proton_pz = true_proton_p_comp.Z();

    double true_proton_P = TMath::Sqrt(true_proton_px*true_proton_px+
				       true_proton_py*true_proton_py+
				       true_proton_pz*true_proton_pz);

    double true_protonTn = TMath::Sqrt(true_proton_P*true_proton_P+(utils->M_p)*(utils->M_p))-utils->M_p;
    double true_protonAngle = TMath::ACos( TMath::Sqrt(true_proton_px*true_proton_px+true_proton_py*true_proton_py)/(true_proton_pz))*180.0/3.14159;

    double true_muon_px   = truth->mc_primFSLepton[0];
    double true_muon_py   = truth->mc_primFSLepton[1];
    double true_muon_pz   = truth->mc_primFSLepton[2];
    double true_muon_E    = truth->mc_primFSLepton[3];
    double true_muon_p    = utils->GetTotalMomentum( true_muon_px, true_muon_py, true_muon_pz ); //GeV
    double true_theta     = utils->GetTheta( true_muon_px, true_muon_py, true_muon_pz );
    double true_muon_pt_beam   = utils->GetTransverseMomentumWRTBeam( true_muon_px, true_muon_py, true_muon_pz ); //GeV 
    double true_muon_pz_beam   = utils->GetLongitudinalMomentumWRTBeam( true_muon_px, true_muon_py, true_muon_pz ); //GeV
  

    //Scales
    true_proton_px /= pow(10,3);
    true_proton_py /= pow(10,3);
    true_proton_pz /= pow(10,3);
    true_proton_P  /= pow(10,3);
    true_protonTn        /= pow(10,3);
    true_muon_px   /= pow(10,3);
    true_muon_py   /= pow(10,3);
    true_muon_pz   /= pow(10,3);
    true_muon_E    /= pow(10,3);
    true_muon_p    /= pow(10,3);
    true_muon_pt_beam /= pow(10,3);
    true_muon_pz_beam /= pow(10,3);


    double true_dalphat   = utils->GetDeltaAlphaT(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_dphit     = utils->GetDeltaPhiT(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_pn        = utils->GetPn(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_dpt       = utils->GetDeltaPt(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_dptx      = utils->GetDeltaPtx(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_dpty      = utils->GetDeltaPty(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_sign      = utils->GetMuonProtonSign(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_signed_dalphat = utils->GetSignedDeltaAlphaT(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_signed_dphit   = utils->GetSignedDeltaPhiT(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);


  //Directional "transverse" variables, Tejin
  // beam_bias default to 0, using this method allows calculation of angle uncertainties. 
  // 1. Define the beam direction, and muon 4P
  double beam_bias = 0;
  XYZVector beam = geoUtils->BeamAxis( beam_bias ); 
  //XYZVector beam_true = XYZVector( evt->mc_incomingPartVec[0],evt->mc_incomingPartVec[1],evt->mc_incomingPartVec[2]).Unit();
  XYZTVector true_muon_4P( true_muon_px, true_muon_py, true_muon_pz ,true_muon_E);
  
  // 2. Compute expected stationary nucleon scattering momentum
  // XYZTVector ComputeExpectedNucleon ( XYZVector &nu, XYZTVector &muon, double ISMass = 0.93827231,double FSMass =   0.93956563, double BindingE = 0.00 ); // vbar p --> mu n, Initial Mass, Final state mass
  double binding_e = 0; // theoretically we should include binding E for carbon.. but I'm comparing to hydrogen in anti-nu, so probably set it to 0.
  double is_neutron_mass = 0.9395654133;
  double fs_proton_mass = 0.93827231;

  XYZTVector expectedProton4P = geoUtils->ComputeExpectedNucleon( beam, true_muon_4P, is_neutron_mass, fs_proton_mass, binding_e);
  XYZVector expectedProton3P = (XYZVector) expectedProton4P.Vect();

    // 3. 
    // ComputeNeutronAngularVars( XYZVector &nu, XYZVector &expVec, XYZVector &targetVec ), variables are : dThetaP, dThetaR, dThetaTotal,  dPp, dPr, inferred dPp and dPr
    XYZVector protonVect_XYZ( true_proton_px, true_proton_py, true_proton_pz );
    vector<double> true_vars = geoUtils->ComputeNeutronAngularVars( beam, expectedProton3P, protonVect_XYZ );
    double true_dThetaP   = true_vars[0]*180/TMath::Pi();
    double true_dThetaR   = true_vars[1]*180/TMath::Pi();
    double true_dTheta    = true_vars[2]*180/TMath::Pi();
    double true_dPp       = true_vars[3];
    double true_dPr       = true_vars[4];
    double true_dPp_infer = true_vars[5];
    double true_dPr_infer = true_vars[6];

 

    map<int,double> true_varmap;
    true_varmap[0] = true_muon_pz_beam;//muonPzbins;
    true_varmap[1] = true_dalphat;//dalphatbins;
    true_varmap[2] = true_dphit;//dphitbins;
    true_varmap[3] = true_pn;//pnbins;
    true_varmap[4] = true_dpt;//dptbins;
    true_varmap[5] = true_dptx;//dptxbins;
    true_varmap[6] = true_dpty;//dptybins;
    true_varmap[7] = true_protonTn;//protonKineticbins;
    true_varmap[8] = true_protonAngle;//protonThetabins;
    true_varmap[9] = true_sign;//signbins;
    true_varmap[10] = true_signed_dalphat;//signeddalphatbins;
    true_varmap[11] = true_signed_dphit;//signeddphitbins;
    true_varmap[12] = true_dThetaP;  
    true_varmap[13] = true_dThetaR;  
    true_varmap[14] = true_dTheta;
    true_varmap[15] = true_dPp;      
    true_varmap[16] = true_dPr;     
    true_varmap[17] = true_dPp_infer;
    true_varmap[18] = true_dPr_infer;

    //Fiducial signal definition
    
    if(true_theta > 20.0) continue; //This is a cut you should use ONLY if you are also applying a 20 degree reco cut... outputting warning
    if(true_proton_p_comp.Mag()<450 || true_proton_p_comp.Mag()>1000) continue; //tracking threshold
    if(true_muon_pz<1.5) continue; //matching threshold

    //No. of True events passing the fiducial volume cut 
    truth_events+=wgt;
    if(typeCheck){
      utils->fillHistos(h_dalphat_ptmu_truth, true_dalphat, true_muon_pt_beam,truth,wgt);
      utils->fillHistos(h_dphit_ptmu_truth, true_dphit, true_muon_pt_beam,truth,wgt);
      utils->fillHistos(h_pn_ptmu_truth, true_pn, true_muon_pt_beam,truth,wgt);
      utils->fillHistos(h_dpt_ptmu_truth, true_dpt, true_muon_pt_beam,truth,wgt);
      utils->fillHistos(h_dptx_ptmu_truth, true_dptx, true_muon_pt_beam,truth,wgt);
      utils->fillHistos(h_dpty_ptmu_truth, true_dpty, true_muon_pt_beam,truth,wgt);
      utils->fillHistos(h_tp_ptmu_truth, true_protonTn, true_muon_pt_beam,truth,wgt);
      utils->fillHistos(h_ptheta_ptmu_truth, true_protonAngle, true_muon_pt_beam,truth,wgt);
      utils->fillHistos(h_signed_ptmu_truth, true_sign, true_muon_pt_beam,truth,wgt);
      utils->fillHistos(h_signeddalphat_ptmu_truth, true_signed_dalphat, true_muon_pt_beam,truth,wgt);
      utils->fillHistos(h_signeddphit_ptmu_truth, true_signed_dphit, true_muon_pt_beam,truth,wgt);
      //Tejin
      utils->fillHistos( h_dthetaP_ptmu_truth, true_dThetaP, true_muon_pt_beam, truth, wgt );
      utils->fillHistos( h_dPP_ptmu_truth, true_dPp, true_muon_pt_beam, truth, wgt );
      utils->fillHistos( h_dPPi_ptmu_truth, true_dPp_infer, true_muon_pt_beam, truth, wgt );
      utils->fillHistos( h_dthetaR_ptmu_truth, true_dThetaR, true_muon_pt_beam, truth, wgt );
      utils->fillHistos( h_dPR_ptmu_truth, true_dPr, true_muon_pt_beam, truth, wgt );
      utils->fillHistos( h_dPRi_ptmu_truth, true_dPr_infer, true_muon_pt_beam, truth, wgt );

    }
    else{
      utils->fillHistos(h_tki_tki_truth, true_varmap[xtype],true_varmap[ytype],truth,wgt);
    }
    //-------------------------------------------------------------------
    // Fill Vertical Error Bands for the 2D PZ vs PT distributions
    // They will be used for the efficiency correction 
    // Also fill in the bands for the 1D pmu and theta distributions 
    //------------------------------------------------------------------- 
    if(RunCodeWithSystematics){
      if(typeCheck){
	utils->fillVertErrorBands( h_dalphat_ptmu_truth, true_dalphat, true_muon_pt_beam, truth);
	utils->fillVertErrorBands( h_dphit_ptmu_truth, true_dphit, true_muon_pt_beam, truth);
	utils->fillVertErrorBands( h_pn_ptmu_truth, true_pn, true_muon_pt_beam, truth);
	utils->fillVertErrorBands( h_dpt_ptmu_truth, true_dpt, true_muon_pt_beam, truth);
	utils->fillVertErrorBands( h_dptx_ptmu_truth, true_dptx, true_muon_pt_beam, truth);
	utils->fillVertErrorBands( h_dpty_ptmu_truth, true_dpty, true_muon_pt_beam, truth);
	utils->fillVertErrorBands( h_tp_ptmu_truth, true_protonTn, true_muon_pt_beam, truth);
	utils->fillVertErrorBands( h_ptheta_ptmu_truth, true_protonAngle, true_muon_pt_beam, truth);
	utils->fillVertErrorBands( h_signed_ptmu_truth, true_sign, true_muon_pt_beam, truth);
	utils->fillVertErrorBands( h_signeddalphat_ptmu_truth, true_signed_dalphat, true_muon_pt_beam, truth);
	utils->fillVertErrorBands( h_signeddphit_ptmu_truth, true_signed_dphit, true_muon_pt_beam, truth);

  //Tejin
  utils->fillVertErrorBands( h_dthetaP_ptmu_truth, true_dThetaP, true_muon_pt_beam,truth);
  utils->fillVertErrorBands( h_dPP_ptmu_truth, true_dPp, true_muon_pt_beam,truth);
  utils->fillVertErrorBands( h_dPPi_ptmu_truth, true_dPp_infer, true_muon_pt_beam,truth);
  utils->fillVertErrorBands( h_dthetaR_ptmu_truth, true_dThetaR, true_muon_pt_beam,truth);
  utils->fillVertErrorBands( h_dPR_ptmu_truth, true_dPr, true_muon_pt_beam,truth);
  utils->fillVertErrorBands( h_dPRi_ptmu_truth, true_dPr_infer, true_muon_pt_beam,truth);


      }
      else{
	utils->fillVertErrorBands(h_tki_tki_truth, true_varmap[xtype],true_varmap[ytype],truth);
      }
    }//end run with systematics
  } //End of for loop over Truth entries 
  
  cout << "Finished looping over all TRUTH entries. Number of weighted TRUE events satisfying all selection criteria = " << truth_events << endl; 
  delete truth;
  delete tree_truth; 

  //---------------------------------------------
  // Get MC Tree
  //---------------------------------------------
  TChain* tree_mc  = utils->getMCTree("CCQENu", n_mcfiles );
  int entries_mc   = tree_mc->GetEntries();
  cout << "MC entries: " << entries_mc << endl;

  
  double mc_events = 0.;
  CCQENuEvent* mc = new CCQENuEvent( tree_mc );
  utils->setmnvHadronReweightDataTree(tree_mc);  
  //---------------------------------------------
  // Fill MC Histograms
  //---------------------------------------------
  for (int i=0; i<entries_mc ; ++i) {
    if( (i+1)%100000 == 0 )
      cout << "Reading MC entry : " << i+1 << " - " << int(100.0*(i+1)/entries_mc) << " %" << endl;
    mc->GetEntry(i);

    if(!cutter->passInteractionVertex(mc))          continue;
    if(!cutter->passExtraTracksProtonsCut(mc))      continue;
    if(mc->multiplicity!=2)                         continue;
    if( !cutter->passInteractionVertex( mc ) )      continue;	
    if( !cutter->passDeadTimeCuts( mc ) )           continue;
    if( !cutter->passNuHelicityCut( mc ) )   	    continue;  
    if( !cutter->passImprovedMichelCut( mc ) )      continue;
    if( !cutter->passExtraTracksProtonsCut ( mc ) ) continue;
    if( !cutter->passNBlobs( mc ) )                 continue;
    if( !cutter->passTrackAngleCut( mc ) )          continue;

    //Get CVWeight
    double wgt = utils->GetCVWeight( mc, "Signal", true );

    const double bindingE = 34.0;

    TVector3 true_proton_p_comp = utils->GetHighestTrueProtonMomentumComp(mc);
    double true_proton_px = true_proton_p_comp.X();
    double true_proton_py = true_proton_p_comp.Y();
    double true_proton_pz = true_proton_p_comp.Z();

    double true_proton_P = TMath::Sqrt(true_proton_px*true_proton_px+
				       true_proton_py*true_proton_py+
				       true_proton_pz*true_proton_pz);

    double true_protonTn = TMath::Sqrt(true_proton_P*true_proton_P+(utils->M_p)*(utils->M_p))-utils->M_p;
    double true_protonAngle = TMath::ACos( TMath::Sqrt(true_proton_px*true_proton_px+true_proton_py*true_proton_py)/(true_proton_pz))*180.0/3.14159;

    double true_muon_px   = mc->mc_primFSLepton[0];
    double true_muon_py   = mc->mc_primFSLepton[1];
    double true_muon_pz   = mc->mc_primFSLepton[2];
    double true_muon_p    = utils->GetTotalMomentum( true_muon_px, true_muon_py, true_muon_pz ); //GeV
    double true_theta     = utils->GetTheta( true_muon_px, true_muon_py, true_muon_pz );
    double true_muon_pt_beam   = utils->GetTransverseMomentumWRTBeam( true_muon_px, true_muon_py, true_muon_pz ); //GeV 
    double true_muon_pz_beam   = utils->GetLongitudinalMomentumWRTBeam( true_muon_px, true_muon_py, true_muon_pz ); //GeV

    //Scales
    true_proton_px /= pow(10,3);
    true_proton_py /= pow(10,3);
    true_proton_pz /= pow(10,3);
    true_proton_P  /= pow(10,3);
    true_protonTn        /= pow(10,3);
    true_muon_px   /= pow(10,3);
    true_muon_py   /= pow(10,3);
    true_muon_pz   /= pow(10,3);
    true_muon_p    /= pow(10,3);
    true_muon_pt_beam /= pow(10,3);
    true_muon_pz_beam /= pow(10,3);



    double true_dalphat   = utils->GetDeltaAlphaT(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_dphit     = utils->GetDeltaPhiT(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_pn        = utils->GetPn(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_dpt       = utils->GetDeltaPt(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_dptx      = utils->GetDeltaPtx(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_dpty      = utils->GetDeltaPty(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_sign      = utils->GetMuonProtonSign(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_signed_dalphat = utils->GetSignedDeltaAlphaT(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_signed_dphit   = utils->GetSignedDeltaPhiT(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);

  //Directional "transverse" variables, Tejin
  // beam_bias default to 0, using this method allows calculation of angle uncertainties. 
  // 1. Define the beam direction, and muon 4P
  double beam_bias = 0;
  XYZVector beam = geoUtils->BeamAxis( beam_bias ); 
  XYZTVector true_muon_4P( true_muon_px, true_muon_py, true_muon_pz ,mc->CCQENu_muon_E/1e3);
  
  // 2. Compute expected stationary nucleon scattering momentum
  // XYZTVector ComputeExpectedNucleon ( XYZVector &nu, XYZTVector &muon, double ISMass = 0.93827231,double FSMass =   0.93956563, double BindingE = 0.00 ); // vbar p --> mu n, Initial Mass, Final state mass
  double binding_e = 0; // theoretically we should include binding E for carbon.. but I'm comparing to hydrogen in anti-nu, so probably set it to 0.
  double is_neutron_mass = 0.9395654133;
  double fs_proton_mass = 0.93827231;

  XYZTVector expectedProton4P = geoUtils->ComputeExpectedNucleon( beam, true_muon_4P, is_neutron_mass, fs_proton_mass, binding_e);
  XYZVector expectedProton3P = (XYZVector) expectedProton4P.Vect();

    // 3. 
    // ComputeNeutronAngularVars( XYZVector &nu, XYZVector &expVec, XYZVector &targetVec ), variables are : dThetaP, dThetaR, dThetaTotal,  dPp, dPr, inferred dPp and dPr
    XYZVector protonVect_XYZ( true_proton_px, true_proton_py, true_proton_pz );
    vector<double> true_vars = geoUtils->ComputeNeutronAngularVars( beam, expectedProton3P, protonVect_XYZ );
    double true_dThetaP   = true_vars[0]*180/TMath::Pi();
    double true_dThetaR   = true_vars[1]*180/TMath::Pi();
    double true_dTheta    = true_vars[2]*180/TMath::Pi();
    double true_dPp       = true_vars[3];
    double true_dPr       = true_vars[4];
    double true_dPp_infer = true_vars[5];
    double true_dPr_infer = true_vars[6];

 

    map<int,double> true_varmap;
    true_varmap[0] = true_muon_pz_beam;//muonPzbins;
    true_varmap[1] = true_dalphat;//dalphatbins;
    true_varmap[2] = true_dphit;//dphitbins;
    true_varmap[3] = true_pn;//pnbins;
    true_varmap[4] = true_dpt;//dptbins;
    true_varmap[5] = true_dptx;//dptxbins;
    true_varmap[6] = true_dpty;//dptybins;
    true_varmap[7] = true_protonTn;//protonKineticbins;
    true_varmap[8] = true_protonAngle;//protonThetabins;
    true_varmap[9] = true_sign;//signbins;
    true_varmap[10] = true_signed_dalphat;//signeddalphatbins;
    true_varmap[11] = true_signed_dphit;//signeddphitbins;
    true_varmap[12] = true_dThetaP;  
    true_varmap[13] = true_dThetaR;  
    true_varmap[14] = true_dTheta;
    true_varmap[15] = true_dPp;      
    true_varmap[16] = true_dPr;     
    true_varmap[17] = true_dPp_infer;
    true_varmap[18] = true_dPr_infer;

    //---------------------------------------------------------------------------------------------------------------
    // These histograms are not filled cut-after-cut. These are filled up for the efficiency correction procedure. 
    // Hence the events are required to pass all cuts before a histo of a type (e.g. QELike, RES, QE) gets filled.  
    //---------------------------------------------------------------------------------------------------------------

    //-------------------------------------------------------------------
    // Check that the CV of events pass the proton and recoil cuts 
    // This is because we have systematics due to these cuts 
    // and a change in their values can affect the composition 
    // of the selected sample. 
    //-------------------------------------------------------------------
    bool doesCVpass_SingleProtonCut = (mc->multiplicity==1)? true : cutter->passSingleProtonCut( mc, 0, 0 );
    bool doesCVpass_ExtraProtonsCut = (mc->multiplicity==1)? true : cutter->passExtraProtonsCut( mc, NULL, 0 );
    bool doesCVpass_RecoilCut = cutter->passSignalFunction( mc, 0, 0 ); 

    if( doesCVpass_SingleProtonCut && doesCVpass_ExtraProtonsCut && doesCVpass_RecoilCut ) {
      
      mc_events+=wgt; 
      if(typeCheck){
	utils->fillHistos( h_dalphat_ptmu, true_dalphat, true_muon_pt_beam,mc, wgt);
	utils->fillHistos( h_dphit_ptmu, true_dphit, true_muon_pt_beam,mc, wgt);
	utils->fillHistos( h_pn_ptmu, true_pn, true_muon_pt_beam,mc, wgt);
	utils->fillHistos( h_dpt_ptmu, true_dpt, true_muon_pt_beam,mc, wgt);
	utils->fillHistos( h_dptx_ptmu, true_dptx, true_muon_pt_beam,mc, wgt);
	utils->fillHistos( h_dpty_ptmu, true_dpty, true_muon_pt_beam,mc, wgt);
	utils->fillHistos( h_tp_ptmu, true_protonTn, true_muon_pt_beam,mc, wgt);
	utils->fillHistos( h_ptheta_ptmu, true_protonAngle, true_muon_pt_beam,mc, wgt);
	utils->fillHistos( h_signed_ptmu, true_sign, true_muon_pt_beam,mc,wgt);
	utils->fillHistos( h_signeddalphat_ptmu, true_signed_dalphat, true_muon_pt_beam,mc,wgt);
	utils->fillHistos( h_signeddphit_ptmu, true_signed_dphit, true_muon_pt_beam,mc,wgt);

  //Tejin
  utils->fillHistos( h_dthetaP_ptmu, true_dThetaP, true_muon_pt_beam, mc, wgt );
  utils->fillHistos( h_dPP_ptmu, true_dPp, true_muon_pt_beam, mc, wgt );
  utils->fillHistos( h_dPPi_ptmu, true_dPp_infer, true_muon_pt_beam, mc, wgt );
  utils->fillHistos( h_dthetaR_ptmu, true_dThetaR, true_muon_pt_beam, mc, wgt );
  utils->fillHistos( h_dPR_ptmu, true_dPr, true_muon_pt_beam, mc, wgt );
  utils->fillHistos( h_dPRi_ptmu, true_dPr_infer, true_muon_pt_beam, mc, wgt );

      }
      else{
	utils->fillHistos( h_tki_tki, true_varmap[xtype],true_varmap[ytype],mc,wgt);
      }

      //Fill Vertical Errors for PZ vs PT
      if(RunCodeWithSystematics){
	if(typeCheck){
	  utils->fillVertErrorBands( h_dalphat_ptmu, true_dalphat, true_muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_dphit_ptmu, true_dphit, true_muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_pn_ptmu, true_pn, true_muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_dpt_ptmu, true_dpt, true_muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_dptx_ptmu, true_dptx, true_muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_dpty_ptmu, true_dpty, true_muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_tp_ptmu, true_protonTn, true_muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_ptheta_ptmu, true_protonAngle, true_muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_signed_ptmu, true_sign, true_muon_pt_beam, mc);
	  utils->fillVertErrorBands( h_signeddalphat_ptmu, true_signed_dalphat, true_muon_pt_beam, mc);
	  utils->fillVertErrorBands( h_signeddphit_ptmu, true_signed_dphit, true_muon_pt_beam, mc);

    //Tejin
    utils->fillVertErrorBands( h_dthetaP_ptmu, true_dThetaP, true_muon_pt_beam, mc);
    utils->fillVertErrorBands( h_dPP_ptmu, true_dPp, true_muon_pt_beam, mc);
    utils->fillVertErrorBands( h_dPPi_ptmu, true_dPp_infer, true_muon_pt_beam, mc);
    utils->fillVertErrorBands( h_dthetaR_ptmu, true_dThetaR, true_muon_pt_beam, mc);
    utils->fillVertErrorBands( h_dPR_ptmu, true_dPr, true_muon_pt_beam, mc);
    utils->fillVertErrorBands( h_dPRi_ptmu, true_dPr_infer, true_muon_pt_beam, mc);


	}
	else{
	  utils->fillVertErrorBands( h_tki_tki, true_varmap[xtype],true_varmap[ytype],mc);
	}
      }//end run with systematics
    }
    if(RunCodeWithSystematics){

      //ADD LAT ERRORS FOR tranvserse    

    }
    
  } // End of for loop over MC entries
  
  cout << "Finished looping over all MC entries. Number of weighted MC events satisfying all selection criteria = " << mc_events << endl; 
  delete mc; 
  delete tree_mc; 
  //Get data POT for particular playlist
  TChain* tree_data  = utils->getDataTree( "CCQENu", -1 );

  //==================================================================
  // Create ROOT file to store histograms
  // Write to file all the created histograms 
  //==================================================================
  if(!typeCheck) filename = filename+"_"+axis_name[xtype]+"_"+axis_name[ytype]+".root";
  TFile *f = new TFile( filename.c_str(), "RECREATE" );

  //Write MC and DATA POT to file 
  utils->writePOT( f );
  
  //Write truth and reconstructed MC number of events 
  TVector2 *evt = new TVector2( truth_events, mc_events );
  f->WriteTObject( evt, "n_events" );

  //Write out if output ROOT file has the constrained flux or not 
  TParameter<bool> *fluxParam = new TParameter<bool>( "appliedFluxConstraint", applyFluxConstraint ); 
  f->WriteTObject( fluxParam ); 
  
  f->cd();
  for( unsigned int j=kData; j<nHistos; ++j ) {
    if(typeCheck){
      h_dalphat_ptmu_truth[j]->Write();
      h_dphit_ptmu_truth[j]->Write();
      h_pn_ptmu_truth[j]->Write();
      h_dpt_ptmu_truth[j]->Write();
      h_dptx_ptmu_truth[j]->Write();
      h_dpty_ptmu_truth[j]->Write();
      h_tp_ptmu_truth[j]->Write();
      h_ptheta_ptmu_truth[j]->Write();
      h_signed_ptmu_truth[j]->Write();
      h_signeddalphat_ptmu_truth[j]->Write();
      h_signeddphit_ptmu_truth[j]->Write();
      
      h_dalphat_ptmu[j]->Write();
      h_dphit_ptmu[j]->Write();
      h_pn_ptmu[j]->Write();
      h_dpt_ptmu[j]->Write();
      h_dptx_ptmu[j]->Write();
      h_dpty_ptmu[j]->Write();
      h_tp_ptmu[j]->Write();
      h_ptheta_ptmu[j]->Write();
      h_signed_ptmu[j]->Write();
      h_signeddalphat_ptmu[j]->Write();
      h_signeddphit_ptmu[j]->Write();

      h_dthetaR_ptmu_truth[j]->Write();
      h_dthetaP_ptmu_truth[j]->Write();
      h_dPP_ptmu_truth[j]->Write();
      h_dPR_ptmu_truth[j]->Write();
      h_dPPi_ptmu_truth[j]->Write();
      h_dPRi_ptmu_truth[j]->Write();

      h_dthetaR_ptmu[j]->Write();
      h_dthetaP_ptmu[j]->Write();
      h_dPP_ptmu[j]->Write();
      h_dPR_ptmu[j]->Write();
      h_dPPi_ptmu[j]->Write();
      h_dPRi_ptmu[j]->Write();

    }
    else{
      h_tki_tki_truth[j]->Write();
      h_tki_tki[j]->Write();
    }
  }

  f->Close();
  delete f;
  delete utils;
  delete cutter; 
  delete binner;
  delete minmodbinner; 

  return 0;
}

void FillEfficiencyHisto( const MnvH2D* num, MnvH2D* den, MnvH2D* eff ){

  //  den->AddMissingErrorBandsAndFillWithCV( *num );
  eff->Divide( num, den, 1.0, 1.0, "B" );

}

void FillEfficiencyHisto( const MnvH1D* num, MnvH1D* den, MnvH1D* eff ){

  //  den->AddMissingErrorBandsAndFillWithCV( *num );
  eff->Divide( num, den, 1.0, 1.0, "B" );

}

int main( int argc, char *argv[])
{

  ROOT::Cintex::Cintex::Enable();

  if (argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t-./EffPurityHists  Name_and_path_to_Output_file  Make_Flux_Constraint_Histo  Apply_Flux_Constraint  Number_MC_files \n\n"<<
      "\t-Name_and_path_to_Output_file\t =\t Name of and path to the Output ROOT file that will be created \n"<< 
      "\t-Playlist\t =\t Name of the playlist you want to run over e.g. minerva1 \n"<<
      "\t-Make_Flux_Constraint_Histo\t =\t If TRUE Enter 1; If FALSE Enter 0 \n"<< 
      "\t-Apply_Flux_Constraint\t =\t If TRUE Enter 1; If FALSE Enter 0 \n"<< 
      "\t-Number_MC_files\t =\t Number of MonteCarlo files. To use all files, set this to: -1"<< std::endl; 
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    return 0; 
  }

  //! Default parameters
  std::vector<std::string> par;
  par.push_back("EffPurityHists");
  par.push_back( Form("%s/ana/rootfiles/EffPurityHists.root",getenv("CCQENUROOT") ) );
  par.push_back("minerva1");
  par.push_back("0");
  par.push_back("0");
  par.push_back("-1");
  par.push_back("-1");//xtype
  par.push_back("-1");//ytype

  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }

  bool fluxConstraintHistoNeeded = ( par.at(3) == "1" ) ? true : false; 

  bool applyFluxConstraint = ( par.at(4) == "1" ) ? true : false; 

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;

  return EffPurityHists( par[1], par[2], fluxConstraintHistoNeeded, applyFluxConstraint, atoi(par[5].c_str()),atoi(par[6].c_str()),atoi(par[7].c_str()) );
}
