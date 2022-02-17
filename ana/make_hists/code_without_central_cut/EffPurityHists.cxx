#include "include/CCQENuUtils.h"
#include "PlotUtils/HyperDimLinearizer.h"
#include "PlotUtils/modifier_Eb.h"
#include "TParameter.h" 
#include "PlotUtils/modifier_Eb.h"

#include "include/Parameters.h"
#include "include/GeneralFunc.h" //in CCQENuNeutron/include/
#include "include/EventCuts.h"
#include "include/CommonBins.h"

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

using namespace CCQENU_ANA;
using namespace ROOT::Math;

void FillEfficiencyHisto( const MnvH1D* num, MnvH1D* den, MnvH1D* eff );
void FillEfficiencyHisto( const MnvH2D* num, MnvH2D* den, MnvH2D* eff );

bool RunSystematics = true;
//RunCodeWithSystematics = true;
bool DoFullKinematicRegion = false;

double EfficiencyLooperMC(CCQENuEvent* evt, CCQENuUtils* utils, CCQENuCuts* cutter, NeutronBlobCuts* ncutter, string sample, int multiplicity, int n_entries, map<int,int> &cuts_counter, int npMode = 0)
{
  GeoUtils* geoUtils = new GeoUtils();
  PhysicsUtils* physicsUtils = new PhysicsUtils();
  KinematicsModifier* kineMod = new KinematicsModifier();
  bool isData = false;
  bool isMC = true;

  double n_evts=0;

  for( int i = 0; i< n_entries; ++i)
  {
    evt->GetEntry(i);
    if (i%1000 == 0) cout<<"At Event "<<i<<endl;
    if( IgnoreEvent( evt, isData )) continue;
    if( !cutter->passTrueCCQELike( evt ) ) continue;

    //==============Reco Cuts======================
    //
    if( !dropClusters ) ncutter->UpdateCandidates( evt, isMC );
    else ncutter->UpdateCandidates( evt, isMC, discardFraction, discardThreshold );
    bool PassEventCut = (EventCut( evt, cutter, ncutter, cuts_counter, multiplicity, sample ,  isData, npMode, useRecoilCut, orderPar));
    //now add in neutron-->
    std::vector<bool> commonNeutronCuts = CommonNeutronCuts(evt,ncutter,orderPar);
    bool hasNC = commonNeutronCuts[0] && commonNeutronCuts[1];
    bool isForward = commonNeutronCuts[2] && hasNC;
    bool isBackward = commonNeutronCuts[3] && hasNC;

    if(npMode == 10 && !hasNC ) continue;
    if (!PassEventCut) continue;


    //cout<<"passed cut"<<endl;
    //==============MC-Neutron-Reco Selection=======================
    // all reco cuts should be applied before
    //-----------
    //
    //

    XYZTVector reco_muon_4P( evt->CCQENu_leptonE[0]/1000,
                             evt->CCQENu_leptonE[1]/1000,
                             evt->CCQENu_leptonE[2]/1000,
                             evt->CCQENu_leptonE[3]/1000);
    if( modify_Eb  ) reco_muon_4P = kineMod->GetCorrectedMuon( evt, reco_muon_4P, mode_eb );
    XYZVector reco_muon_3P = (XYZVector) reco_muon_4P.Vect();

    double reco_muon_theta     = evt->muon_theta_allNodes[20]; // get theta from 20th node, note NX kludge fixes this...
    double reco_cos_muon_theta = cos(reco_muon_theta);
    double reco_muon_T         =  reco_muon_4P.E() - reco_muon_4P.M();
    double reco_muon_px        =  reco_muon_4P.X(); 
    double reco_muon_py        =  reco_muon_4P.Y(); 
    double reco_muon_pz        =  reco_muon_4P.Z(); 
    double reco_muon_e         =  reco_muon_4P.E(); 
    double reco_muon_p         =  reco_muon_4P.R();
    reco_muon_theta*=180/TMath::Pi();


    if(reco_muon_p<MuonPMin && reco_muon_p > MuonPMax && useMuonPCut) continue; // a muon momentum cut


    //=========================================================================
    // Other than selection, all variables should be calculated using truth
    //cout<<"Other than selection, all variables should be calculated using truth"<<endl;


    //Convert to degrees, GeV
    //
    int charge = (neutrinoMode)? -1:1;
    double binding_e = 0; 
    double is_part_mass = (neutrinoMode)? 0.9395654133: 0.93827231; //neutron/proton
    double fs_part_mass = (neutrinoMode)? 0.93827231: 0.9395654133;

    //Beam
    double beam_bias = 0;
    XYZVector beam = geoUtils->BeamAxis( beam_bias ); 



    //cout<<q2qe<<", "<<wgt<<endl;


    XYZVector blobVtxDir = ncutter->MainNeutronCandidate()->BlobVtxDir;

    XYZTVector reco_expectedNucleon4P = geoUtils->ComputeExpectedNucleon( beam, reco_muon_4P, is_part_mass, fs_part_mass, binding_e);
    XYZVector reco_expectedNucleon3P = (XYZVector) reco_expectedNucleon4P.Vect();

    std::vector<double> reco_neutron_vars = geoUtils->ComputeNeutronAngularVars( beam, reco_expectedNucleon3P, blobVtxDir );

    double reco_dthetaP = reco_neutron_vars[vdthetaP];
    double reco_dthetaR = reco_neutron_vars[vdthetaR];
    double reco_dtheta = reco_neutron_vars[vdtheta];
    reco_dthetaP *= 180/3.14159;
    reco_dthetaR *= 180/3.14159;
    reco_dtheta *= 180/3.14159;


    bool reco_10 = (abs(reco_dthetaP)<=10 && abs(reco_dthetaR)<=10);
    bool reco_01 = (abs(reco_dthetaP)<=1 && abs(reco_dthetaR)<=1);
    bool reco_02 = (abs(reco_dthetaP)<=2 && abs(reco_dthetaR)<=2);
    bool reco_03 = (abs(reco_dthetaP)<=3 && abs(reco_dthetaR)<=3);
    bool reco_04 = (abs(reco_dthetaP)<=4 && abs(reco_dthetaR)<=4);

    //==============MC Event Selection========================

    //cout<<"MC Event Selection"<<endl;
    NeutronCandidates* ptrNeutronCandidates; 
    if (npMode != 1 || hasNC) ptrNeutronCandidates = ncutter->GetCandidates();
    bool fill_common_histos = true;

    //neutron candidate energy cut
    double Mn=0.9395654133;
    double blobEnergy = ncutter->MainNeutronCandidate()->TotalE;
    //if(npMode==0 && q2qe/2/Mn < blobEnergy/1000. ) continue;
    //==============MC/TRUTH weight ========================
    double wgt = utils->GetCVWeight( evt );
    n_evts+=wgt;

    //==============MC-True Variables=======================
    // all reco cuts should be applied before
    //-----------

    //cout<<"MC-True Variables"<<endl;
    //Lepton

    XYZTVector muon( evt->mc_primFSLepton[0]/1000,
                     evt->mc_primFSLepton[1]/1000,
                     evt->mc_primFSLepton[2]/1000,
                     evt->mc_primFSLepton[3]/1000);
    if( modify_Eb ) muon = kineMod->GetCorrectedMuon( evt, muon, mode_eb );


    double true_muon_px   = muon.X(); 
    double true_muon_py   = muon.Y(); 
    double true_muon_pz   = muon.Z(); 
    double true_muon_E    = muon.E(); 
    double true_muon_p    = muon.R(); 
    double true_muon_theta     = utils->GetTheta( true_muon_px, true_muon_py, true_muon_pz );
    double true_muon_pt_beam   = utils->GetTransverseMomentumWRTBeam( true_muon_px, true_muon_py, true_muon_pz ); //GeV 
    double true_muon_pz_beam   = utils->GetLongitudinalMomentumWRTBeam( true_muon_px, true_muon_py, true_muon_pz ); //GeV
    true_muon_theta*=180/TMath::Pi();

    //if( true_muon_theta > 20 ) continue;
    //---> true muon cut should only appear in efficiency denominator


    //Calculating Q2qe:


    XYZTVector nu( evt->mc_incomingPartVec[0]/1000,evt->mc_incomingPartVec[1]/1000,evt->mc_incomingPartVec[2]/1000,evt->mc_incomingPartVec[3]/1000);


    double q2qe = physicsUtils->qsqCCQE( beam, muon, charge, binding_e ); //GeV^2
    double enu = physicsUtils->nuEnergyCCQE( beam, muon, charge, binding_e ); //GeV^2

    //Nucleons
    XYZTVector proton = Convert(  utils->GetHighestTruePart4P( evt, 2212 ) )/1000;
    XYZTVector neutron = Convert( utils->GetHighestTruePart4P( evt, 2112 ) )/1000;
    bool hasProton = proton.E() > 0;
    bool hasNeutron = neutron.E() > 0;
    if( modify_Eb && hasProton ) proton = kineMod->GetCorrectedNucleon( evt, proton, mode_eb );
    if( modify_Eb && hasNeutron ) neutron = kineMod->GetCorrectedNucleon( evt, neutron, mode_eb );


    //double true_muon_theta = ACos( nu.Vect().Unit().Dot( muon.Vect().Unit() ) )*180/TMath::Pi();
    //if( !evt->truth_is_fiducial ) continue;

    XYZTVector primaryNucleon = (npMode == 0)? neutron : proton;
    XYZTVector secondaryNucleon = (npMode == 0)? proton : neutron;
    bool hasPrimary = primaryNucleon.E()>0;
    bool hasSecondary = secondaryNucleon.E()>0;



    //cout<<"hasPrimary:"<<hasPrimary<<endl;
    //if(npMode!=0 && !hasPrimary) continue;

    
    //cout<<"npMode == 0 and hasProton:"<<(npMode == 0 && hasProton)<<endl;
    //if(npMode == 0 && hasProton) continue;
    //if ( npMode !=0 && !NucleonKinematicsCut( beam, (XYZVector) proton.Vect() ) ) continue;

    //Calculations

    double cos_muon_theta =  beam.Dot( muon.Vect().Unit() );
    double muon_theta= ACos( cos_muon_theta );
    double muon_T    = (muon.E()-muon.M());
    double muon_px   = muon.Px();
    double muon_py   = muon.Py();
    double muon_pz   = muon.Pz();
    double muon_p    = muon.P();

    //New method (maybe temporary)
    double muon_pt_beam = sin(muon_theta)*muon_p;
    double muon_pz_beam = cos(muon_theta)*muon_p;

    //Convert to degrees
    muon_theta = true_muon_theta;

    // 3. expected nucleon
    XYZTVector expectedNucleon4P = geoUtils->ComputeExpectedNucleon( beam, muon, is_part_mass, fs_part_mass, binding_e);
    XYZVector expectedNucleon3P = (XYZVector) expectedNucleon4P.Vect();


    //preparing for STKI and angular vars:
    double Mcarbon = 6*fs_part_mass+6*is_part_mass-92.163/1E3;
    double Enuc=primaryNucleon.E();


    //calculating transverse variables
    std::vector<double> transverseVars = geoUtils->GetTransverseVariables( beam, muon, primaryNucleon, Mcarbon, 27.13/1E3);
    //calculating angular variables
    XYZVector primaryNucleon3D = primaryNucleon.Vect();
    std::vector<double> neutron_vars = geoUtils->ComputeNeutronAngularVars( beam, expectedNucleon3P, primaryNucleon3D );
    // test if the expected and muon cancel, they should

    double En      = transverseVars[vEnu];
    double pn      = transverseVars[vPn];
    double dpt     = transverseVars[vdpT];
    double dptx    = transverseVars[vdpTx];
    double dpty    = transverseVars[vdpTy];
    double dalphat = transverseVars[vdalphaT];
    double dphit   = transverseVars[vdphiT];
    double sign    = transverseVars[vsign];
    double signed_dalphat = sign*dalphat ;
    double signed_dphit = dphit *sign;

    double dthetaP = neutron_vars[vdthetaP];
    double dthetaR = neutron_vars[vdthetaR];
    double dtheta = neutron_vars[vdtheta];
    //if (isMC) cout<< evt->mc_intType<<endl;
    dthetaP *= 180/TMath::Pi();
    dthetaR *= 180/TMath::Pi();
    dtheta *= 180/TMath::Pi();

    //cout<<dthetaP<<endl;
    double dPp = neutron_vars[vdPp];
    double dPr = neutron_vars[vdPr];
    double dPp_infer = neutron_vars[vdPpi];
    double dPr_infer = neutron_vars[vdPri];

    double protonTn = ( proton.E()- proton.M());
    double protonAngle = ACos( proton.Vect().Unit().Dot( beam ) ) * 180/TMath::Pi();

    string quadrant="";
    if( hasNeutron && hasProton ) 
    {
      XYZVector mu = muon.Vect();
      XYZVector prot = proton.Vect();
      XYZVector neut = neutron.Vect();
      quadrant = GetQuadrantCode( beam, mu, prot, neut );
      //cout<<quadrant<<endl;
    }

    if (fill_common_histos)
    {
      //---------------------------------------------
      //Fill 1-D Plots
      //---------------------------------------------
     // cout<<"Fill start"<<endl;
      utils->fillHistosV3( utils->histos1D["h_ptmu"], muon_pt_beam,  isData, evt, wgt );
      utils->fillHistosV3( utils->histos1D["h_q2qe"], q2qe,  isData, evt, wgt );
      utils->fillHistosV3( utils->histos1D["h_enu"], enu,  isData, evt, wgt );
      //cout<<abs(reco_dthetaR)<<", "<<abs(reco_dthetaP)<<endl;
      if(reco_10) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_10"], q2qe,  isData, evt, wgt );
      if(reco_01) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_01"], q2qe,  isData, evt, wgt );
      if(reco_02) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_02"], q2qe,  isData, evt, wgt );
      if(reco_03) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_03"], q2qe,  isData, evt, wgt );
      if(reco_04) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_04"], q2qe,  isData, evt, wgt );
      utils->fillHistosV3( utils->histos2D["h_q2qe_ptmu"], q2qe, muon_pt_beam, isData, evt, wgt );
      

      //=============Possible conflict in tki analysis ==========================
      if(hasNeutron && hasProton && hasNC)
      {
        utils->fillHistosV3( utils->histos2D[ Form("h_pn_q2qe_%s", quadrant.c_str())], pn, q2qe, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_dphit_q2qe_%s", quadrant.c_str())], dphit, q2qe, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_dalphat_q2qe_%s", quadrant.c_str())], dalphat, q2qe, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_dpt_q2qe_%s", quadrant.c_str())], dpt, q2qe, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_dptx_q2qe_%s", quadrant.c_str())], dptx, q2qe, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_dpty_q2qe_%s", quadrant.c_str())], dpty, q2qe, isData,evt,wgt);
      }

      FillHDL2D( "h_dthetaPdthetaR_q2qe", utils, "dthetaP", "q2", "dthetaR", dthetaP, q2qe, dthetaR, isData, evt, wgt, true);
      if(hasNeutron) 
      {
        FillHDL2D( "h_dthetaPdthetaR_q2qe_nc", utils, "dthetaP", "q2", "dthetaR", dthetaP, q2qe, dthetaR, isData, evt, wgt, true);
      }

      utils->fillHistosV3( utils->histos2D[ "h_dthetaR_q2qe"], dthetaR, q2qe, isData, evt, wgt );
      utils->fillHistosV3( utils->histos2D[ "h_dthetaP_q2qe"], dthetaP, q2qe, isData, evt, wgt );
      if(reco_10) utils->fillHistosV3( utils->histos2D[ "h_muontheta_q2qe"], muon_theta, q2qe, isData, evt, wgt );
      if(reco_10) utils->fillHistosV3( utils->histos2D[ "h_muonmomentum_q2qe"], muon_p, q2qe, isData, evt, wgt );
      if(reco_10) utils->fillHistosV3( utils->histos2D[ "h_enu_q2qe"], enu, q2qe, isData, evt, wgt );

      //------------------------------------------------------
      //Fill 2D Vertical Error Bands
      //Error Bands are not filled up for tmu_costheta 
      //and pmu_ptmu histos to avoid code slow down 
      //------------------------------------------------------
      if(RunSystematics && !isData){

        utils->fillVertErrorBands( utils->histos1D["h_ptmu"],  muon_pt_beam, evt );
        utils->fillVertErrorBands( utils->histos1D["h_q2qe"],  q2qe, evt );
        utils->fillVertErrorBands( utils->histos1D["h_enu"],  enu, evt );
        utils->fillVertErrorBands( utils->histos2D["h_q2qe_ptmu"],  q2qe, muon_pt_beam, evt );
        if(reco_10) utils->fillVertErrorBands( utils->histos2D[ "h_muontheta_q2qe"], muon_theta, q2qe,evt);
        if(reco_10) utils->fillVertErrorBands( utils->histos2D[ "h_muonmomentum_q2qe"], muon_p, q2qe, evt);
        if(reco_10) utils->fillVertErrorBands( utils->histos2D[ "h_enu_q2qe"], enu, q2qe, evt );

       if(reco_10) utils->fillVertErrorBands( utils->histos1D["h_q2qe_angle_10"], q2qe, evt );
       if(reco_01) utils->fillVertErrorBands( utils->histos1D["h_q2qe_angle_01"], q2qe, evt );
       if(reco_02) utils->fillVertErrorBands( utils->histos1D["h_q2qe_angle_02"], q2qe, evt );
       if(reco_03) utils->fillVertErrorBands( utils->histos1D["h_q2qe_angle_03"], q2qe, evt );
       if(reco_04) utils->fillVertErrorBands( utils->histos1D["h_q2qe_angle_04"], q2qe, evt );
 

        //---------------------------------------------------------------------------------------
        //Now fillLatErrorBands for events whose CV did not pass the proton and recoil cuts  
        //What if the shifted proton score and recoil pass the selection cuts ? 
        //Those events should also be included in the LatErrorBands
        //fillLatErrorBands have their own checks for shifted proton scores and recoils  
        //---------------------------------------------------------------------------------------

        //----------------------------------------------------
        //Fill 1D Lateral Error Bands
        //JO has turned them OFF since they're already 
        //being filled in 2D below, avoids code slow down. 
        //----------------------------------------------------
        utils->fillLatErrorBands( utils->histos1D["h_ptmu"],  "pT",    muon_pt_beam, evt, true );
        utils->fillLatErrorBands( utils->histos1D["h_q2qe"],  "q2",    q2qe, evt, true );
        utils->fillLatErrorBands( utils->histos2D["h_q2qe_ptmu"],  "q2", "pT",  q2qe, muon_pt_beam, evt );


        if(reco_10) utils->fillLatErrorBands( utils->histos2D[ "h_enu_q2qe"], "enu", "q2",enu, q2qe, evt );
        if(reco_10) utils->fillLatErrorBands( utils->histos2D[ "h_muontheta_q2qe"],"theta","q2", muon_theta, q2qe,evt);
        if(reco_10) utils->fillLatErrorBands( utils->histos2D[ "h_muonmomentum_q2qe"], "emu", "q2",muon_p, q2qe, evt);

        if(reco_10) utils->fillLatErrorBands( utils->histos1D["h_q2qe_angle_10"], "q2", q2qe, evt, true );
        if(reco_01) utils->fillLatErrorBands( utils->histos1D["h_q2qe_angle_01"], "q2", q2qe, evt, true );
        if(reco_02) utils->fillLatErrorBands( utils->histos1D["h_q2qe_angle_02"], "q2", q2qe, evt, true );
        if(reco_03) utils->fillLatErrorBands( utils->histos1D["h_q2qe_angle_03"], "q2", q2qe, evt, true );
        if(reco_04) utils->fillLatErrorBands( utils->histos1D["h_q2qe_angle_04"], "q2", q2qe, evt, true );
 

      }//end run with systematics        
    } // End of fillHistosV3 and fillVertErrorBands where CV of events pass the proton and recoil cuts 
  }
  return n_evts;
}


double EfficiencyLooperTruth(CCQENuTruth* evt, CCQENuUtils* utils, CCQENuCuts* cutter, NeutronBlobCuts* ncutter, string sample, int multiplicity, int n_entries, map<int,int> &cuts_counter, int npMode = 0)
{
  GeoUtils* geoUtils = new GeoUtils();
  PhysicsUtils* physicsUtils = new PhysicsUtils();
  KinematicsModifier* kineMod = new KinematicsModifier();
  bool isData = false;
  bool isMC = true;

  double n_evts=0;

  for( int i = 0; i< n_entries; ++i)
  {
    evt->GetEntry(i);
    if( IgnoreEvent( evt, isData )) continue;
    if( !cutter->passTrueCCQELike( evt ) ) continue;
    if (i%1000 == 0) cout<<"At Event "<<i<<endl;
    //==============MC Event Selection========================
    bool fill_common_histos = true;
    // fiducial?


    //==============MC/TRUTH weight ========================
    double wgt = utils->GetCVWeight( evt );
    n_evts+=wgt;

    //==============MC-True Variables=======================
    // all reco cuts should be applied before
    //-----------
    //Beam
    double beam_bias = 0;
    XYZVector beam = geoUtils->BeamAxis( beam_bias ); 

    //vtx
    XYZVector vertex( evt->mc_vtx[0],evt->mc_vtx[1],evt->mc_vtx[2] );

    //Lepton
    XYZTVector muon( evt->mc_primFSLepton[0]/1000.,evt->mc_primFSLepton[1]/1000.,evt->mc_primFSLepton[2]/1000.,evt->mc_primFSLepton[3]/1000.);
    if( modify_Eb ) muon = kineMod->GetCorrectedMuon( evt, muon, mode_eb );
    XYZTVector nu( evt->mc_incomingPartVec[0]/1000.,evt->mc_incomingPartVec[1]/1000.,evt->mc_incomingPartVec[2]/1000.,evt->mc_incomingPartVec[3]/1000.);
    
    //Fiducial and Muon Angle Cut
    double true_muon_theta = ACos( nu.Vect().Unit().Dot( muon.Vect().Unit() ) )*180/3.1415927;
    if( !evt->truth_is_fiducial ) continue;
    //cout<<"DoFullKinematicRegion?"<<DoFullKinematicRegion<<endl;
    if( !DoFullKinematicRegion && true_muon_theta > 20 ) continue;
    //cout<<true_muon_theta<<endl;
    //Nucleons
    XYZTVector proton = Convert(  utils->GetHighestTruePart4P( evt, 2212 ) )/1000.;
    XYZTVector neutron = Convert( utils->GetHighestTruePart4P( evt, 2112 ) )/1000.;
    bool hasProton = proton.E() > 0;
    bool hasNeutron = neutron.E() > 0;

    //correct for binding e
    if( modify_Eb && hasProton ) proton = kineMod->GetCorrectedNucleon( evt, proton, mode_eb );
    if( modify_Eb && hasNeutron ) neutron = kineMod->GetCorrectedNucleon( evt, neutron, mode_eb );

    if(npMode == 0 && hasProton) continue;
    if ( npMode !=0 && !NucleonKinematicsCut( beam, (XYZVector) proton.Vect() ) ) continue;

    XYZTVector primaryNucleon = (npMode == 0)? neutron : proton;
    XYZTVector secondaryNucleon = (npMode == 0)? proton : neutron;
    bool hasPrimary = primaryNucleon.E()>0;
    bool hasSecondary = secondaryNucleon.E()>0;

    if(!hasPrimary) continue;

    //Calculations

    double cos_muon_theta =  beam.Dot( muon.Vect().Unit() );
    double muon_theta= ACos( cos_muon_theta );
    double muon_T    = (muon.E()-muon.M());
    double muon_px   = muon.Px();
    double muon_py   = muon.Py();
    double muon_pz   = muon.Pz();
    double muon_p    = muon.P();
    if(muon_p<MuonPMin && muon_p > MuonPMax && useMuonPCut) continue; // a muon momentum cut

    //New method (maybe temporary)
    double muon_pt_beam = sin(muon_theta)*muon_p;
    double muon_pz_beam = cos(muon_theta)*muon_p;

    //Convert to degrees
    muon_theta*= 180. / 3.14159;

    // Binding and is/fs mass
    double binding_e = 0; 
    double is_part_mass = (neutrinoMode)? 0.9395654133: 0.93827231; //neutron/proton
    double fs_part_mass = (neutrinoMode)? 0.93827231: 0.9395654133;

    // 3. expected nucleon
    XYZTVector expectedNucleon4P = geoUtils->ComputeExpectedNucleon( beam, muon, is_part_mass, fs_part_mass, binding_e);
    XYZVector expectedNucleon3P = (XYZVector) expectedNucleon4P.Vect();

    //Calculating Q2qe:
    int charge = (neutrinoMode)? -1:1;
    double q2qe = physicsUtils->qsqCCQE( beam, muon, charge, binding_e ); //GeV^2
    double enu = physicsUtils->nuEnergyCCQE( beam, muon, charge, binding_e ); //GeV^2
    //cout<<q2qe<<", "<<wgt<<endl;

    //preparing for STKI and angular vars:
    double Mcarbon = 6*fs_part_mass+6*is_part_mass-92.163/1E3;
    double Enuc=primaryNucleon.E();


    //calculating transverse variables
    std::vector<double> transverseVars = geoUtils->GetTransverseVariables( beam, muon, primaryNucleon, Mcarbon, 27.13/1E3);
    //calculating angular variables
    XYZVector primaryNucleon3D = primaryNucleon.Vect();
    std::vector<double> neutron_vars = geoUtils->ComputeNeutronAngularVars( beam, expectedNucleon3P, primaryNucleon3D );
    // test if the expected and muon cancel, they should

    double En      = transverseVars[vEnu];
    double pn      = transverseVars[vPn];
    double dpt     = transverseVars[vdpT];
    double dptx    = transverseVars[vdpTx];
    double dpty    = transverseVars[vdpTy];
    double dalphat = transverseVars[vdalphaT];
    double dphit   = transverseVars[vdphiT];
    double sign    = transverseVars[vsign];
    double signed_dalphat = sign*dalphat ;
    double signed_dphit = dphit *sign;

    double dthetaP = neutron_vars[vdthetaP];
    double dthetaR = neutron_vars[vdthetaR];
    double dtheta = neutron_vars[vdtheta];
    //if (isMC) cout<< evt->mc_intType<<endl;
    dthetaP *= 180/TMath::Pi();
    dthetaR *= 180/TMath::Pi();
    dtheta *= 180/TMath::Pi();
    //cout<<dthetaP<<endl;
    double dPp = neutron_vars[vdPp];
    double dPr = neutron_vars[vdPr];
    double dPp_infer = neutron_vars[vdPpi];
    double dPr_infer = neutron_vars[vdPri];

    double protonTn = ( proton.E()- proton.M());
    double protonAngle = ACos( proton.Vect().Unit().Dot( beam ) ) * 180/TMath::Pi();

    bool q2Cut = passQ2( q2qe );

    string quadrant="";
    if( hasNeutron && hasProton ) 
    {
      XYZVector mu = muon.Vect();
      XYZVector prot = proton.Vect();
      XYZVector neut = neutron.Vect();
      quadrant = GetQuadrantCode( beam, mu, prot, neut );
      //cout<<quadrant<<endl;
    }

    if (fill_common_histos)
    {
      //---------------------------------------------
      //Fill 1-D Plots
      //---------------------------------------------
      //cout<<"Fill start"<<endl;
      utils->fillHistosV3( utils->histos1D["h_truth_ptmu"], muon_pt_beam,  isData, evt, wgt );
      utils->fillHistosV3( utils->histos1D["h_truth_q2qe"], q2qe,  isData, evt, wgt );
      utils->fillHistosV3( utils->histos1D["h_truth_enu"], enu,  isData, evt, wgt );
      utils->fillHistosV3( utils->histos2D["h_truth_q2qe_ptmu"], q2qe, muon_pt_beam, isData, evt, wgt );
      
        //cout<<utils->histos1D["h_truth_q2qe"][kMC]->Integral()<<endl;

    
      if(hasNeutron && hasProton)
      {
        utils->fillHistosV3( utils->histos2D[ Form("h_truth_pn_q2qe_%s", quadrant.c_str())], pn, q2qe, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_truth_dphit_q2qe_%s", quadrant.c_str())], dphit, q2qe, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_truth_dalphat_q2qe_%s", quadrant.c_str())], dalphat, q2qe, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_truth_dpt_q2qe_%s", quadrant.c_str())], dpt, q2qe, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_truth_dptx_q2qe_%s", quadrant.c_str())], dptx, q2qe, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_truth_dpty_q2qe_%s", quadrant.c_str())], dpty, q2qe, isData,evt,wgt);
      }


      utils->fillHistosV3( utils->histos2D["h_truth_dthetaR_q2qe"], dthetaR, q2qe, isData, evt, wgt );
      utils->fillHistosV3( utils->histos2D["h_truth_dthetaP_q2qe"], dthetaP, q2qe, isData, evt, wgt );
      utils->fillHistosV3( utils->histos2D["h_truth_muontheta_q2qe"], muon_theta, q2qe, isData, evt, wgt );
      utils->fillHistosV3( utils->histos2D["h_truth_muonmomentum_q2qe"], muon_p, q2qe, isData, evt, wgt );
      utils->fillHistosV3( utils->histos2D["h_truth_enu_q2qe"], enu, q2qe, isData, evt, wgt );

      FillHDL2D( "h_truth_dthetaPdthetaR_q2qe", utils, dthetaP, q2qe, dthetaR, isData, evt, wgt, true);

      //------------------------------------------------------
      //Fill 2D Vertical Error Bands
      //Error Bands are not filled up for tmu_costheta 
      //and pmu_ptmu histos to avoid code slow down 
      //------------------------------------------------------
      if(RunSystematics && !isData){

        utils->fillVertErrorBands( utils->histos1D["h_truth_ptmu"],  muon_pt_beam, evt );
        utils->fillVertErrorBands( utils->histos1D["h_truth_q2qe"],  q2qe, evt );
        utils->fillVertErrorBands( utils->histos1D["h_truth_enu"],  enu, evt );
        utils->fillVertErrorBands( utils->histos2D["h_truth_q2qe_ptmu"],  q2qe, muon_pt_beam, evt );


        utils->fillVertErrorBands( utils->histos2D["h_truth_dthetaR_q2qe"], dthetaR, q2qe, evt );
        utils->fillVertErrorBands( utils->histos2D["h_truth_dthetaP_q2qe"], dthetaP, q2qe, evt );
        utils->fillVertErrorBands( utils->histos2D["h_truth_muontheta_q2qe"], muon_theta, q2qe, evt );
        utils->fillVertErrorBands( utils->histos2D["h_truth_muonmomentum_q2qe"], muon_p, q2qe, evt );
      }//end run with systematics        
    } // End of fillHistosV3 and fillVertErrorBands where CV of events pass the proton and recoil cuts 
  }
  return n_evts;
}






int EffPurityHists(string filename, string playlist, bool makeFluxConstraintHisto, bool applyFluxConstraint, int multiplicity, int n_mcfiles = -1,  int npMode = 0)
{
  //---------------------------------------------
  // create file to store histograms
  //---------------------------------------------
  //TFile *f = new TFile( filename.c_str(), "RECREATE" );
  CCQENuUtils  *utils          = new CCQENuUtils( false, makeFluxConstraintHisto );
  utils->setPlaylist(playlist);
  utils->setFluxReweighterPlaylist();

  gStyle->SetNumberContours(500);
  AnaBinning *binner                  = new AnaBinning();

  CCQENuCuts   *cutter         = new CCQENuCuts();
  NeutronBlobCuts  *ncutter = new NeutronBlobCuts();
  PhysicsUtils* physicsUtils = new PhysicsUtils();
  GeoUtils *geoUtils = new GeoUtils();

  cout<<"NPMode = "<<npMode<<endl;


  //! Picking this option for 'theta' binning for J. Wolcott's plots (June 2015)  
  
  //! This option is used for the efficiency study in each pZ-pT bin  
  //axis_binning muonPtbins         = minmodbinner->GetMuonPt_FinerBinsGeV();
  //axis_binning muonPzbins         = minmodbinner->GetMuonPz_FinerBinsGeV();
  
  //! Reminder: "nu" here signifies the recoil (non-leptonic) energy 
  //axis_binning nubins           = minmodbinner->GetNuBinsGeV(); 


  MnvH1D *h_ptmu[nHistos],*h_q2qe[nHistos], *h_enu[nHistos];
  MnvH1D *h_q2qe_angle_10[nHistos];
  MnvH1D *h_q2qe_angle_01[nHistos],*h_q2qe_angle_02[nHistos],*h_q2qe_angle_03[nHistos],*h_q2qe_angle_04[nHistos];

  utils->bookHistos( h_ptmu, "h_ptmu", Form( "Muon Pt;Reconstructed Muon p_{T} (GeV);Events / %.3f GeV", muonPtbins.default_width ), muonPtbins );
  utils->bookHistos( h_q2qe, "h_q2qe", Form( "Q^{2}_{QE};Reconstructed Q^{2}_{QE} (GeV^{2});Events / %.3f GeV^{2}", Q2bins.default_width ), Q2bins );
  utils->bookHistos( h_enu, "h_enu", Form( "E_{#nu};Reconstructed E_{#nu} (GeV);Events / %.3f GeV", Enubins.default_width ), Enubins );
  utils->bookHistos( h_q2qe_angle_10, "h_q2qe_angle_10", Form( "Q^{2}_{QE};Reconstructed Q^{2}_{QE} (GeV^{2});Events / %.3f GeV^{2}", Q2bins.default_width ), Q2bins );
  utils->bookHistos( h_q2qe_angle_01, "h_q2qe_angle_01", Form( "Q^{2}_{QE};Reconstructed Q^{2}_{QE} (GeV^{2});Events / %.3f GeV^{2}", Q2bins.default_width ), Q2bins );
  utils->bookHistos( h_q2qe_angle_02, "h_q2qe_angle_02", Form( "Q^{2}_{QE};Reconstructed Q^{2}_{QE} (GeV^{2});Events / %.3f GeV^{2}", Q2bins.default_width ), Q2bins );
  utils->bookHistos( h_q2qe_angle_03, "h_q2qe_angle_03", Form( "Q^{2}_{QE};Reconstructed Q^{2}_{QE} (GeV^{2});Events / %.3f GeV^{2}", Q2bins.default_width ), Q2bins );
  utils->bookHistos( h_q2qe_angle_04, "h_q2qe_angle_04", Form( "Q^{2}_{QE};Reconstructed Q^{2}_{QE} (GeV^{2});Events / %.3f GeV^{2}", Q2bins.default_width ), Q2bins );

  MnvH2D *h_q2qe_ptmu[nHistos];
  utils->bookHistos( h_q2qe_ptmu, "h_q2qe_ptmu","q2qe:ptmu",Q2bins,muonPtbins);

  axis_binning   dthetaPdthetaRBins;
  MnvH2D   *h_dthetaPdthetaR_q2qe[nHistos];
  MnvH2D   *h_dthetaPdthetaR_q2qe_nc[nHistos];

  HyperDimLinearizer *hdl_dthetaPdthetaR_q2qe = DeclareHDLHistos2D( utils, h_dthetaPdthetaR_q2qe, "h_dthetaPdthetaR_q2qe", "q2qe vs dthetaPdthetaR: dthetaP : dthetaR", dthetaPerpbins, Q2bins, dthetaReactbins, dthetaPdthetaRBins );
  HyperDimLinearizer *hdl_dthetaPdthetaR_q2qe_nc = DeclareHDLHistos2D( utils, h_dthetaPdthetaR_q2qe_nc, "h_dthetaPdthetaR_q2qe_nc", "q2qe vs dthetaPdthetaR with nc: dthetaP : dthetaR", dthetaPerpbins, Q2bins, dthetaReactbins, dthetaPdthetaRBins );

  MnvH2D *h_pn_q2qe_ll[nHistos],*h_pn_q2qe_lr[nHistos],*h_pn_q2qe_rl[nHistos],*h_pn_q2qe_rr[nHistos];
  MnvH2D *h_dpt_q2qe_ll[nHistos],*h_dpt_q2qe_lr[nHistos],*h_dpt_q2qe_rl[nHistos],*h_dpt_q2qe_rr[nHistos];
  MnvH2D *h_dptx_q2qe_ll[nHistos],*h_dptx_q2qe_lr[nHistos],*h_dptx_q2qe_rl[nHistos],*h_dptx_q2qe_rr[nHistos];
  MnvH2D *h_dpty_q2qe_ll[nHistos],*h_dpty_q2qe_lr[nHistos],*h_dpty_q2qe_rl[nHistos],*h_dpty_q2qe_rr[nHistos];
  MnvH2D *h_dalphat_q2qe_ll[nHistos],*h_dalphat_q2qe_lr[nHistos],*h_dalphat_q2qe_rl[nHistos],*h_dalphat_q2qe_rr[nHistos];
  MnvH2D *h_dphit_q2qe_ll[nHistos],*h_dphit_q2qe_lr[nHistos],*h_dphit_q2qe_rl[nHistos],*h_dphit_q2qe_rr[nHistos];


  utils->bookHistos( h_pn_q2qe_ll, "h_pn_q2qe_ll", "pn_ll:q2qe",pnbins,Q2bins); utils->bookHistos( h_pn_q2qe_lr, "h_pn_q2qe_lr", "pn_lr:q2qe",pnbins,Q2bins); utils->bookHistos( h_pn_q2qe_rl, "h_pn_q2qe_rl", "pn_rl:q2qe",pnbins,Q2bins);utils->bookHistos( h_pn_q2qe_rr, "h_pn_q2qe_rr", "pn_rr:q2qe",pnbins,Q2bins);
  utils->bookHistos( h_dpt_q2qe_ll, "h_dpt_q2qe_ll", "dpt_ll:q2qe",dptbins,Q2bins); utils->bookHistos( h_dpt_q2qe_lr, "h_dpt_q2qe_lr", "dpt_lr:q2qe",dptbins,Q2bins); utils->bookHistos( h_dpt_q2qe_rl, "h_dpt_q2qe_rl", "dpt_rl:q2qe",dptbins,Q2bins);utils->bookHistos( h_dpt_q2qe_rr, "h_dpt_q2qe_rr", "dpt_rr:q2qe",dptbins,Q2bins);
  utils->bookHistos( h_dptx_q2qe_ll, "h_dptx_q2qe_ll", "dptx_ll:q2qe",dptxbins,Q2bins); utils->bookHistos( h_dptx_q2qe_lr, "h_dptx_q2qe_lr", "dptx_lr:q2qe",dptxbins,Q2bins); utils->bookHistos( h_dptx_q2qe_rl, "h_dptx_q2qe_rl", "dptx_rl:q2qe",dptxbins,Q2bins);utils->bookHistos( h_dptx_q2qe_rr, "h_dptx_q2qe_rr", "dptx_rr:q2qe",dptxbins,Q2bins);
  utils->bookHistos( h_dpty_q2qe_ll, "h_dpty_q2qe_ll", "dpty_ll:q2qe",dptybins,Q2bins); utils->bookHistos( h_dpty_q2qe_lr, "h_dpty_q2qe_lr", "dpty_lr:q2qe",dptybins,Q2bins); utils->bookHistos( h_dpty_q2qe_rl, "h_dpty_q2qe_rl", "dpty_rl:q2qe",dptybins,Q2bins);utils->bookHistos( h_dpty_q2qe_rr, "h_dpty_q2qe_rr", "dpty_rr:q2qe",dptybins,Q2bins);
  utils->bookHistos( h_dphit_q2qe_ll, "h_dphit_q2qe_ll", "dphit_ll:q2qe",dphitbins,Q2bins); utils->bookHistos( h_dphit_q2qe_lr, "h_dphit_q2qe_lr", "dphit_lr:q2qe",dphitbins,Q2bins); utils->bookHistos( h_dphit_q2qe_rl, "h_dphit_q2qe_rl", "dphit_rl:q2qe",dphitbins,Q2bins);utils->bookHistos( h_dphit_q2qe_rr, "h_dphit_q2qe_rr", "dphit_rr:q2qe",dphitbins,Q2bins);
  utils->bookHistos( h_dalphat_q2qe_ll, "h_dalphat_q2qe_ll", "dalphat_ll:q2qe",dalphatbins,Q2bins); utils->bookHistos( h_dalphat_q2qe_lr, "h_dalphat_q2qe_lr", "dalphat_lr:q2qe",dalphatbins,Q2bins); utils->bookHistos( h_dalphat_q2qe_rl, "h_dalphat_q2qe_rl", "dalphat_rl:q2qe",dalphatbins,Q2bins);utils->bookHistos( h_dalphat_q2qe_rr, "h_dalphat_q2qe_rr", "dalphat_rr:q2qe",dalphatbins,Q2bins);

  


  cout<<"booked quandrants"<<endl;

  MnvH2D* h_dthetaR_q2qe[nHistos], *h_dthetaP_q2qe[nHistos];
  MnvH2D* h_muontheta_q2qe[nHistos], *h_muonmomentum_q2qe[nHistos];
  MnvH2D* h_enu_q2qe[nHistos];
  utils->bookHistos( h_dthetaR_q2qe, "h_dthetaR_q2qe", "dthetaR:q2qe", dthetaReactbins, Q2bins );
  utils->bookHistos( h_dthetaP_q2qe, "h_dthetaP_q2qe", "dthetaP:q2qe", dthetaPerpbins, Q2bins );
  utils->bookHistos( h_muontheta_q2qe, "h_muontheta_q2qe", "muontheta:q2qe", thetabins, Q2bins );
  utils->bookHistos( h_muonmomentum_q2qe, "h_muonmomentum_q2qe", "muonmomentum:q2qe", leptonmomentumbins, Q2bins );
  utils->bookHistos( h_enu_q2qe, "h_enu_q2qe", "enu:q2qe", leptonmomentumbins, Q2bins );

  //--------------------------------------------------------------
  if(RunSystematics){
      utils->addVertErrorBands( h_ptmu );
      utils->addLatErrorBands ( h_ptmu );
      cout<<1<<endl;
      
      utils->addVertErrorBands( h_q2qe );
      utils->addLatErrorBands ( h_q2qe );

      utils->addVertErrorBands( h_enu );
      utils->addLatErrorBands ( h_enu );

      cout<<2<<endl;
      utils->addVertErrorBands( h_q2qe_angle_10 );
      utils->addLatErrorBands ( h_q2qe_angle_10 );

      cout<<3<<endl;
      utils->addVertErrorBands( h_q2qe_angle_01 );
      utils->addLatErrorBands ( h_q2qe_angle_01 );
      utils->addVertErrorBands( h_q2qe_angle_02 );
      utils->addLatErrorBands ( h_q2qe_angle_02 );
      utils->addVertErrorBands( h_q2qe_angle_03 );
      utils->addLatErrorBands ( h_q2qe_angle_03 );
      utils->addVertErrorBands( h_q2qe_angle_04 );
      utils->addLatErrorBands ( h_q2qe_angle_04 );

      cout<<4<<endl;
      utils->addVertErrorBands( h_q2qe_ptmu );
      utils->addLatErrorBands ( h_q2qe_ptmu );

      utils->addVertErrorBands( h_dthetaR_q2qe );
      utils->addVertErrorBands( h_dthetaP_q2qe );
      utils->addVertErrorBands( h_muontheta_q2qe );
      utils->addVertErrorBands( h_muonmomentum_q2qe );
      utils->addVertErrorBands( h_enu_q2qe );
      utils->addLatErrorBands( h_dthetaR_q2qe );
      utils->addLatErrorBands( h_dthetaP_q2qe );
      utils->addLatErrorBands( h_muontheta_q2qe );
      utils->addLatErrorBands( h_muonmomentum_q2qe );
      utils->addLatErrorBands( h_enu_q2qe );

  }//end run with systematics

  cout<<"booked reco systematics"<<endl;

  //==================== Define Truth Histos =======================

  MnvH1D *h_truth_ptmu[nHistos],*h_truth_q2qe[nHistos], *h_truth_enu[nHistos];

  utils->bookHistos( h_truth_ptmu, "h_truth_ptmu", Form( "Muon Pt;Reconstructed Muon p_{T} (GeV);Events / %.3f GeV", muonPtbins.default_width ), muonPtbins );
  utils->bookHistos( h_truth_q2qe, "h_truth_q2qe", Form( "Q^{2}_{QE};Reconstructed Q^{2}_{QE} (GeV^{2});Events / %.3f GeV^{2}", Q2bins.default_width ), Q2bins );
  utils->bookHistos( h_truth_enu, "h_truth_enu", Form( "E_{#nu};Reconstructed E_{#nu} (GeV);Events / %.3f GeV", Enubins.default_width ), Enubins );

  MnvH2D *h_truth_q2qe_ptmu[nHistos];
  utils->bookHistos( h_truth_q2qe_ptmu, "h_truth_q2qe_ptmu","q2qe:ptmu",Q2bins,muonPtbins);

  MnvH2D *h_truth_pn_q2qe_ll[nHistos],*h_truth_pn_q2qe_lr[nHistos],*h_truth_pn_q2qe_rl[nHistos],*h_truth_pn_q2qe_rr[nHistos];
  MnvH2D *h_truth_dpt_q2qe_ll[nHistos],*h_truth_dpt_q2qe_lr[nHistos],*h_truth_dpt_q2qe_rl[nHistos],*h_truth_dpt_q2qe_rr[nHistos];
  MnvH2D *h_truth_dptx_q2qe_ll[nHistos],*h_truth_dptx_q2qe_lr[nHistos],*h_truth_dptx_q2qe_rl[nHistos],*h_truth_dptx_q2qe_rr[nHistos];
  MnvH2D *h_truth_dpty_q2qe_ll[nHistos],*h_truth_dpty_q2qe_lr[nHistos],*h_truth_dpty_q2qe_rl[nHistos],*h_truth_dpty_q2qe_rr[nHistos];
  MnvH2D *h_truth_dalphat_q2qe_ll[nHistos],*h_truth_dalphat_q2qe_lr[nHistos],*h_truth_dalphat_q2qe_rl[nHistos],*h_truth_dalphat_q2qe_rr[nHistos];
  MnvH2D *h_truth_dphit_q2qe_ll[nHistos],*h_truth_dphit_q2qe_lr[nHistos],*h_truth_dphit_q2qe_rl[nHistos],*h_truth_dphit_q2qe_rr[nHistos];

  utils->bookHistos( h_truth_pn_q2qe_ll, "h_truth_pn_q2qe_ll", "pn_ll:q2qe",pnbins,Q2bins); utils->bookHistos( h_truth_pn_q2qe_lr, "h_truth_pn_q2qe_lr", "pn_lr:q2qe",pnbins,Q2bins); utils->bookHistos( h_truth_pn_q2qe_rl, "h_truth_pn_q2qe_rl", "pn_rl:q2qe",pnbins,Q2bins);utils->bookHistos( h_truth_pn_q2qe_rr, "h_truth_pn_q2qe_rr", "pn_rr:q2qe",pnbins,Q2bins);
  utils->bookHistos( h_truth_dpt_q2qe_ll, "h_truth_dpt_q2qe_ll", "dpt_ll:q2qe",dptbins,Q2bins); utils->bookHistos( h_truth_dpt_q2qe_lr, "h_truth_dpt_q2qe_lr", "dpt_lr:q2qe",dptbins,Q2bins); utils->bookHistos( h_truth_dpt_q2qe_rl, "h_truth_dpt_q2qe_rl", "dpt_rl:q2qe",dptbins,Q2bins);utils->bookHistos( h_truth_dpt_q2qe_rr, "h_truth_dpt_q2qe_rr", "dpt_rr:q2qe",dptbins,Q2bins);
  utils->bookHistos( h_truth_dptx_q2qe_ll, "h_truth_dptx_q2qe_ll", "dptx_ll:q2qe",dptxbins,Q2bins); utils->bookHistos( h_truth_dptx_q2qe_lr, "h_truth_dptx_q2qe_lr", "dptx_lr:q2qe",dptxbins,Q2bins); utils->bookHistos( h_truth_dptx_q2qe_rl, "h_truth_dptx_q2qe_rl", "dptx_rl:q2qe",dptxbins,Q2bins);utils->bookHistos( h_truth_dptx_q2qe_rr, "h_truth_dptx_q2qe_rr", "dptx_rr:q2qe",dptxbins,Q2bins);
  utils->bookHistos( h_truth_dpty_q2qe_ll, "h_truth_dpty_q2qe_ll", "dpty_ll:q2qe",dptybins,Q2bins); utils->bookHistos( h_truth_dpty_q2qe_lr, "h_truth_dpty_q2qe_lr", "dpty_lr:q2qe",dptybins,Q2bins); utils->bookHistos( h_truth_dpty_q2qe_rl, "h_truth_dpty_q2qe_rl", "dpty_rl:q2qe",dptybins,Q2bins);utils->bookHistos( h_truth_dpty_q2qe_rr, "h_truth_dpty_q2qe_rr", "dpty_rr:q2qe",dptybins,Q2bins);
  utils->bookHistos( h_truth_dphit_q2qe_ll, "h_truth_dphit_q2qe_ll", "dphit_ll:q2qe",dphitbins,Q2bins); utils->bookHistos( h_truth_dphit_q2qe_lr, "h_truth_dphit_q2qe_lr", "dphit_lr:q2qe",dphitbins,Q2bins); utils->bookHistos( h_truth_dphit_q2qe_rl, "h_truth_dphit_q2qe_rl", "dphit_rl:q2qe",dphitbins,Q2bins);utils->bookHistos( h_truth_dphit_q2qe_rr, "h_truth_dphit_q2qe_rr", "dphit_rr:q2qe",dphitbins,Q2bins);
  utils->bookHistos( h_truth_dalphat_q2qe_ll, "h_truth_dalphat_q2qe_ll", "dalphat_ll:q2qe",dalphatbins,Q2bins); utils->bookHistos( h_truth_dalphat_q2qe_lr, "h_truth_dalphat_q2qe_lr", "dalphat_lr:q2qe",dalphatbins,Q2bins); utils->bookHistos( h_truth_dalphat_q2qe_rl, "h_truth_dalphat_q2qe_rl", "dalphat_rl:q2qe",dalphatbins,Q2bins);utils->bookHistos( h_truth_dalphat_q2qe_rr, "h_truth_dalphat_q2qe_rr", "dalphat_rr:q2qe",dalphatbins,Q2bins);


  cout<<"booked quandrants"<<endl;

  MnvH2D* h_truth_dthetaPdthetaR_q2qe[nHistos];
  MnvH2D* h_truth_dthetaPdthetaR_q2qe_nc[nHistos];

  HyperDimLinearizer *hdl_truth_dthetaPdthetaR_q2qe = DeclareHDLHistos2D( utils, h_dthetaPdthetaR_q2qe, "h_truth_dthetaPdthetaR_q2qe", "q2qe vs dthetaPdthetaR: dthetaP : dthetaR", dthetaPerpbins, Q2bins, dthetaReactbins, dthetaPdthetaRBins );
  HyperDimLinearizer *hdl_truth_dthetaPdthetaR_q2qe_nc = DeclareHDLHistos2D( utils, h_dthetaPdthetaR_q2qe_nc, "h_truth_dthetaPdthetaR_q2qe_nc", "q2qe vs dthetaPdthetaR with nc: dthetaP : dthetaR", dthetaPerpbins, Q2bins, dthetaReactbins, dthetaPdthetaRBins );


  MnvH2D* h_truth_dthetaR_q2qe[nHistos], *h_truth_dthetaP_q2qe[nHistos];
  MnvH2D* h_truth_muontheta_q2qe[nHistos], *h_truth_muonmomentum_q2qe[nHistos];
  MnvH2D* h_truth_enu_q2qe[nHistos];
  utils->bookHistos( h_truth_dthetaR_q2qe, "h_truth_dthetaR_q2qe", "dthetaR:q2qe", dthetaReactbins, Q2bins );
  utils->bookHistos( h_truth_dthetaP_q2qe, "h_truth_dthetaP_q2qe", "dthetaP:q2qe", dthetaPerpbins, Q2bins );
  utils->bookHistos( h_truth_muontheta_q2qe, "h_truth_muontheta_q2qe", "muontheta:q2qe", thetabins, Q2bins );
  utils->bookHistos( h_truth_muonmomentum_q2qe, "h_truth_muonmomentum_q2qe", "muonmomentum:q2qe", leptonmomentumbins, Q2bins );
  utils->bookHistos( h_truth_enu_q2qe, "h_truth_enu_q2qe", "enu:q2qe", leptonmomentumbins, Q2bins );

  //--------------------------------------------------------------
  if(RunSystematics){
      utils->addVertErrorBands( h_truth_ptmu );
      utils->addVertErrorBands( h_truth_q2qe );
      utils->addVertErrorBands( h_truth_q2qe_ptmu );  

      utils->addVertErrorBands( h_truth_dthetaR_q2qe );
      utils->addVertErrorBands( h_truth_dthetaP_q2qe );
      utils->addVertErrorBands( h_truth_muontheta_q2qe );
      utils->addVertErrorBands( h_truth_muonmomentum_q2qe );
  }//end run with systematics

  cout<<"booked truth systematics"<<endl;







  map<int,int> cuts_counter;


  //---------------------------------------------
  // Get Truth Tree
  //---------------------------------------------
  TChain* tree_truth  = utils->getMCTree( "Truth", n_mcfiles );
  int entries_truth = tree_truth->GetEntries();
  cout << "Truth entries: " << entries_truth << endl;
  utils->setmnvHadronReweightTruthTree(tree_truth);

  CCQENuTruth* truth = new CCQENuTruth( tree_truth );

  //-------------------------------------------
  // Fill Truth histograms 
  //-------------------------------------------
  //entries_truth=0;
  double truth_events = EfficiencyLooperTruth(truth, utils,cutter, ncutter, "Signal", multiplicity, entries_truth, cuts_counter, npMode);

  cout << "Finished looping over all TRUTH entries. Number of weighted TRUE events satisfying all selection criteria = " << truth_events << endl; 
  delete truth;
  delete tree_truth; 

  //---------------------------------------------
  // Get MC Tree
  //---------------------------------------------
  TChain* tree_mc  = utils->getMCTree("CCQENu", n_mcfiles );
  int entries_mc   = tree_mc->GetEntries();
  cout << "MC entries: " << entries_mc << endl;

  
  CCQENuEvent* mc = new CCQENuEvent( tree_mc );
  utils->setmnvHadronReweightDataTree(tree_mc);  
  //---------------------------------------------
  // Fill MC Histograms
  //---------------------------------------------
  double mc_events = EfficiencyLooperMC(mc, utils,cutter, ncutter, "Signal", multiplicity, entries_mc, cuts_counter, npMode);
  cout << "Finished looping over all MC entries. Number of weighted MC events satisfying all selection criteria = " << mc_events << endl; 
  delete mc; 
  delete tree_mc; 







  //Get data POT for particular playlist
  TChain* tree_data  = utils->getDataTree( "CCQENu", -1 );

  //==================================================================
  // Create ROOT file to store histograms
  // Write to file all the created histograms 
  //==================================================================
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

  WriteHistos( utils, f );

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
  par.push_back("1");
  par.push_back("-1");
  par.push_back("0");//proton, neutron,proton/neutron mode, 1,0,10
  par.push_back("0");//drop clusters?
  par.push_back("0");//do full kinematic region?

  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }

  string filename = par[1];
  string playlist = par[2];
  bool fluxConstraintHistoNeeded = ( par.at(3) == "1" ) ? true : false; 
  bool applyFluxConstraint = ( par.at(4) == "1" ) ? true : false; 
  int multiplicity = atoi(par[5].c_str());
  int nMCFiles = atoi(par[6].c_str());
  int NPMode = atoi(par[7].c_str());
  dropClusters = bool(atoi(par[8].c_str() ) );
  DoFullKinematicRegion = bool( atoi(par[9].c_str() ) );

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;

  return EffPurityHists( filename, playlist, fluxConstraintHistoNeeded, applyFluxConstraint, multiplicity,nMCFiles,NPMode );
}
