#include "include/GlobalParameters.h"
#include "include/CCQENuUtils.h"
#include "TParameter.h" 
#include "PlotUtils/HyperDimLinearizer.h"

#include "include/GeneralFunc.h" //in CCQENuNeutron/include/
#include "include/EventCuts.h"

using namespace CCQENU_ANA;


//forward declaration
//HyperDimLinearizer* DeclareHDLHistos2D( CCQENuUtils* utils, MnvH2D **h, std::string name, std::string title, axis_binning &xbins, axis_binning &ybins, axis_binning &zbins, axis_binning &global_x );
//void FillHDL2D( std::string name, CCQENuUtils *utils, double x, double y, double z, bool isData, CCQENuEvent* evt, double wgt );
//================================================

//int npMode = 1;//0:neutron, 1:proton, 10:neutron+proton

template<class T>
double AnalysisLooperBase(T* evt, CCQENuUtils* utils, CCQENuCuts* cutter, NeutronBlobCuts* ncutter, string sample, int multiplicity, int n_entries, map<int,int> &cuts_counter, bool isData=false, bool isTruth=false, int npMode = 1)
{
  GeoUtils* geoUtils = new GeoUtils();
  PhysicsUtils* physicsUtils = new PhysicsUtils();

  double n_evts=0;

  for( int i = 0; i< n_entries; ++i)
  {
    evt->GetEntry(i);
    if (i%1000 == 0) cout<<"At Event "<<i<<endl;
    //==============Fill Neutron Blobs======================
    bool isMC = !isData;
    ncutter->UpdateCandidates( evt, isMC );
    if(isDebug) cout<<"Updated ncutter"<<endl;
    //==============Event Selection========================
    bool PassEventCut = (EventCut( evt, cutter, ncutter, cuts_counter, multiplicity, sample ,  isData, npMode));
    if(isDebug) cout<<PassEventCut<<endl;
    //bool PassCCQESelection = (evt->multiplicity == multiplicity && cutter->passCCQESelection(evt, sample) );
    //if (PassEventCut != PassCCQESelection)
    //{
    //  cout<<multiplicity<<", "<<evt->multiplicity<<endl;
    //  cout<<"event cutters don't agree at "<<i<<endl;
    //  return n_evts;
    //}
    if (!PassEventCut) continue;
    bool fill_common_histos = true;
    //bool pass_signalFunc = cutter->passSignalFunction( evt, 0, 0 ) ;
    //bool pass_singleProton = cutter->passSingleProtonCut( evt, 0, 0 );
    //bool pass_extraProtons = cutter->passExtraProtonsCut( evt, 0, 0 );
    //fill_common_histos =( pass_signalFunc && pass_singleProton && pass_extraProtons );
    //===============Reconstruction ========================
    double wgt = (isData)? 1: utils->GetCVWeight( evt, sample );
    n_evts+=wgt;
    if(isDebug) cout<<wgt<<endl;

    //-----------
    //Muon Variables
    double muon_theta     = evt->muon_theta_allNodes[20]; // get theta from 20th node, note NX kludge fixes this...
    double cos_muon_theta =  cos(muon_theta);
    double muon_T         =  evt->CCQENu_muon_T / pow(10,3); //GeV
    double reco_muon_px   =  evt->CCQENu_leptonE[0];
    double reco_muon_py   =  evt->CCQENu_leptonE[1];
    double reco_muon_pz   =  evt->CCQENu_leptonE[2];
    double muon_p         = utils->GetTotalMomentum( reco_muon_px, reco_muon_py, reco_muon_pz ) / pow(10,3); //GeV

    //New method (maybe temporary)
    double muon_pt_beam = sin(muon_theta)*muon_p;
    double muon_pz_beam = cos(muon_theta)*muon_p;

    //Convert to degrees, GeV
    muon_theta*= 180. / 3.14159;
    reco_muon_px /= pow(10,3);
    reco_muon_py /= pow(10,3);
    reco_muon_pz /= pow(10,3);

    // 1. beam and muon
    double beam_bias = 0;
    XYZVector beam = geoUtils->BeamAxis( beam_bias ); 
    TVector3 beamV3(beam.X(), beam.Y(), beam.Z() );

    XYZTVector reco_muon_4P( reco_muon_px, reco_muon_py, reco_muon_pz ,evt->CCQENu_muon_E/1e3);
    XYZVector reco_muon_3P( reco_muon_px, reco_muon_py, reco_muon_pz);



    //Proton Variables
    TVector3 protonVect;
    double protonAngle, protonMom, protonTn;

    if (npMode != 0 )
    {
      protonVect.SetXYZ(evt->CCQENu_proton_Px_fromdEdx,evt->CCQENu_proton_Py_fromdEdx,evt->CCQENu_proton_Pz_fromdEdx);//MeV
      protonAngle = evt->CCQENu_proton_theta*180./3.14159;
      protonMom = protonVect.Mag();//MeV
      protonTn = TMath::Sqrt(protonMom*protonMom+utils->M_p*utils->M_p)-utils->M_p;//Done in MeV    }

      //Convert to degrees, GeV
      protonVect *= pow(10,-3);
      protonTn /= pow(10,3);
    }

    //Truth Proton Variables
    //TVector3 true_proton_p_comp = utils->GetHighestTrueProtonMomentumComp(evt);
    TVector3 true_proton_p_comp(evt->proton_prong_4p[0],evt->proton_prong_4p[1],evt->proton_prong_4p[2]);
    double true_protonMom = true_proton_p_comp.Mag();
    double proton_res = 1-protonMom/true_protonMom;
    
    //Neutron Variables
    TVector3 neutronDir(-1,-1,-1);
    XYZVector blobVtxDir(-1,-1,-1);
    if (npMode !=1)
    {
      blobVtxDir = ncutter->MainNeutronCandidate()->BlobVtxDir;
      neutronDir.SetXYZ(blobVtxDir.X(),blobVtxDir.Y(), blobVtxDir.Z() ) ; 
    }
    TVector3 true_neutron_p_comp;
    double true_neutronMom=-999; 
    double neutron_res_angle=-999;
    double neutron_res_angle_x=-999;
    double neutron_res_angle_y=-999;
    //1. look for the most energetic fs neutron
    if( isMC )
    {
      std::vector<TLorentzVector> TrueFSNeutrons = GetAllParticlesMomenta( evt, 2112, 0 );
      if(TrueFSNeutrons.size() > 0 )
      {
        true_neutron_p_comp = TrueFSNeutrons[0].Vect();
        TVector3 hx = true_neutron_p_comp.Cross( beamV3 ).Unit();
        TVector3 hy = true_neutron_p_comp.Cross( hx ).Unit();

        true_neutronMom = true_neutron_p_comp.Mag();
        neutron_res_angle = acos( true_neutron_p_comp.Unit().Dot( neutronDir.Unit() ) );
        neutron_res_angle_x = acos( true_neutron_p_comp.Unit().Dot( (neutronDir.Dot( hx )*hx).Unit() ) );
        neutron_res_angle_y = acos( true_neutron_p_comp.Unit().Dot( (neutronDir.Dot( hy )*hy).Unit() ) );
      }
    }


    //Expected nucleon kinematics

    // 2. binding and is/fs mass
    double binding_e = 0; 
    double is_part_mass = (GlobalParameters::Get().neutrinoMode)? 0.9395654133: 0.93827231; //neutron/proton
    double fs_part_mass = (GlobalParameters::Get().neutrinoMode)? 0.93827231: 0.9395654133;

    // 3. expected nucleon

    XYZTVector expectedNucleon4P = geoUtils->ComputeExpectedNucleon( beam, reco_muon_4P, is_part_mass, fs_part_mass, binding_e);
    XYZVector expectedNucleon3P = (XYZVector) expectedNucleon4P.Vect();


    // Set the primary nucleon momentum vec
    // proton/pn mode : proton. neutron mode: neutron direction
    TVector3 primaryNucleonVect = ( npMode == 0 )? neutronDir.Unit()*expectedNucleon3P.R() : protonVect ;
    //cout<<"Is proton primary? "<< (primaryNucleonVect==protonVect)<<endl;

    XYZVector primaryNucleonVect_XYZ( primaryNucleonVect.X(), primaryNucleonVect.Y(), primaryNucleonVect.Z() );

    //Calculating Q2qe:
    int charge = (GlobalParameters::Get().neutrinoMode)? -1:1;
    double q2qe = physicsUtils->qsqCCQE( beam, reco_muon_4P, charge, binding_e ); //GeV^2



    //preparing for STKI and angular vars:
    double Mcarbon = 6*fs_part_mass+6*is_part_mass-92.163/1E3;
    double Enuc=sqrt(primaryNucleonVect_XYZ.Mag2()+fs_part_mass*fs_part_mass);
    //XYZTVector primaryNucleonVect_XYZT( protonVect.X(), protonVect.Y(), protonVect.Z(), Enuc );
    XYZTVector primaryNucleonVect_XYZT( primaryNucleonVect_XYZ.X(), primaryNucleonVect_XYZ.Y(), primaryNucleonVect_XYZ.Z(), Enuc );

    //calculating transverse variables
    std::vector<double> transverseVars = geoUtils->GetTransverseVariables( beam, reco_muon_4P, primaryNucleonVect_XYZT, Mcarbon, 27.13/1E3);
    //calculating angular variables
    std::vector<double> reco_neutron_vars = geoUtils->ComputeNeutronAngularVars( beam, expectedNucleon3P, primaryNucleonVect_XYZ );
    // test if the expected and muon cancel, they should

    double reco_En      = transverseVars[0];
    double reco_pn      = transverseVars[1];
    double reco_dpt     = transverseVars[2];
    double reco_dptx    = transverseVars[3];
    double reco_dpty    = transverseVars[4];
    double reco_dalphat = transverseVars[5];
    double reco_dphit   = transverseVars[6];
    double reco_sign    = transverseVars[7];
    double reco_signed_dalphat = reco_sign*reco_dalphat ;
    double reco_signed_dphit = reco_dphit *reco_sign;

    double reco_dThetaP = reco_neutron_vars[0];
    double reco_dThetaR = reco_neutron_vars[1];
    double reco_dTheta = reco_neutron_vars[2];
    //if (isMC) cout<< evt->mc_intType<<endl;
    reco_dThetaP *= 180/3.14159;
    reco_dThetaR *= 180/3.14159;
    reco_dTheta *= 180/3.14159;
    //cout<<reco_dThetaP<<endl;
    double reco_dPp = reco_neutron_vars[3];
    double reco_dPr = reco_neutron_vars[4];
    double reco_dPp_infer = reco_neutron_vars[5];
    double reco_dPr_infer = reco_neutron_vars[6];


    //cout<<"dpTx - dPp = "<<abs(reco_dptx) - abs(reco_dPp)<<endl;

    //I feel there is also a need  to include the left-right uncertainty of the beam, since LR-direction is pretty important....

    if (fill_common_histos)
    {
      //---------------------------------------------
      //Fill 1-D Plots
      //---------------------------------------------
      //cout<<"Fill start"<<endl;
      utils->fillHistosV3( utils->histos1D["h_multiplicity"], evt->multiplicity,  isData, evt, wgt ); 
      utils->fillHistosV3( utils->histos1D["h_ptmu"], muon_pt_beam,  isData, evt, wgt );
      utils->fillHistosV3( utils->histos1D["h_q2qe"], q2qe,  isData, evt, wgt );
      utils->fillHistosV3( utils->histos2D["h_q2qe_ptmu"], q2qe, muon_pt_beam, isData, evt, wgt );
      
      int n_nodes = evt->CCQENu_proton_nodes_nodesNormE_sz;
      int pattrec = evt->CCQENu_proton_patternRec;
      //cout<<"Fill 1"<<endl;
      utils->fillHistosV3(utils->histos1D["h_nodes"],n_nodes, isData, evt,wgt);
      if(pattrec==1 || pattrec==2)utils->fillHistosV3(utils->histos1D["h_nodes_long"],n_nodes, isData, evt,wgt);
      else if(pattrec==3)         utils->fillHistosV3(utils->histos1D["h_nodes_anchored"],n_nodes, isData, evt,wgt);
      else if (pattrec==6)        utils->fillHistosV3(utils->histos1D["h_nodes_vest"],n_nodes, isData, evt,wgt);
      
      double summed_01 = 0.0;
      for(int node=0;node<n_nodes;node++)
      {
        utils->fillHistosV3(utils->histos2D["h_node_nodeenergy"],node, evt->CCQENu_proton_nodes_nodesNormE[node], isData, evt,wgt);
        if(pattrec==1||pattrec==2)
          utils->fillHistosV3(utils->histos2D["h_node_nodeenergy_long"],node, evt->CCQENu_proton_nodes_nodesNormE[node], isData, evt,wgt);
        else if(pattrec==3)
          utils->fillHistosV3(utils->histos2D["h_node_nodeenergy_anchored"],node, evt->CCQENu_proton_nodes_nodesNormE[node], isData, evt,wgt);
        else if(pattrec==6) 
          utils->fillHistosV3(utils->histos2D["h_node_nodeenergy_vest"],node, evt->CCQENu_proton_nodes_nodesNormE[node], isData, evt,wgt);
        
        if(node==0) summed_01+= evt->CCQENu_proton_nodes_nodesNormE[node];
        if(node==1)
        {
          summed_01+= evt->CCQENu_proton_nodes_nodesNormE[node];
          utils->fillHistosV3(utils->histos2D["h_protonRes_nodeenergy_node01"],summed_01,proton_res, isData, evt,wgt);
          summed_01=0;
        }
        if(node==2) utils->fillHistosV3(utils->histos2D["h_protonRes_nodeenergy_node2"],  evt->CCQENu_proton_nodes_nodesNormE[node],proton_res, isData, evt,wgt);
        if(node==3) utils->fillHistosV3(utils->histos2D["h_protonRes_nodeenergy_node3"],  evt->CCQENu_proton_nodes_nodesNormE[node],proton_res, isData, evt,wgt);
        if(node==4) utils->fillHistosV3(utils->histos2D["h_protonRes_nodeenergy_node4"],  evt->CCQENu_proton_nodes_nodesNormE[node],proton_res, isData, evt,wgt);
        if(node==5) utils->fillHistosV3(utils->histos2D["h_protonRes_nodeenergy_node5"],  evt->CCQENu_proton_nodes_nodesNormE[node],proton_res, isData, evt,wgt);
        
	  
      }
      
      //-------------------------------------------------------
      //Fill 2-D Plots
      //-------------------------------------------------------
      utils->fillHistosV3( utils->histos2D["h_dalphat_ptmu"],reco_dalphat,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dphit_ptmu"],reco_dphit,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_pn_ptmu"],reco_pn,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dpt_ptmu"],reco_dpt,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dptx_ptmu"],reco_dptx,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dpty_ptmu"],reco_dpty,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_tp_ptmu"],protonTn, muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_ptheta_ptmu"],protonAngle, muon_pt_beam, isData, evt,wgt);
      //utils->fillHistosV3( utils->histos2D["h_sign_ptmu"],reco_sign,muon_pt_beam, isData, evt,wgt);
      //utils->fillHistosV3( utils->histos2D["h_signdalphat_ptmu"],reco_signed_dalphat,muon_pt_beam, isData, evt,wgt);
      //utils->fillHistosV3( utils->histos2D["h_signdphit_ptmu"],reco_signed_dphit,muon_pt_beam, isData, evt,wgt);

      utils->fillHistosV3( utils->histos2D["h_dthetaP_ptmu"],reco_dThetaP,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dthetaR_ptmu"],reco_dThetaR,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPP_ptmu"],reco_dPp,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPR_ptmu"],reco_dPr,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPPi_ptmu"],reco_dPp_infer,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPRi_ptmu"],reco_dPr_infer,muon_pt_beam, isData, evt,wgt);

      utils->fillHistosV3( utils->histos2D["h_dthetaP_q2qe"],reco_dThetaP,  q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dthetaR_q2qe"],reco_dThetaR,  q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPP_q2qe"],    reco_dPp,      q2qe, isData,evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPR_q2qe"],    reco_dPr,      q2qe, isData,evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPPi_q2qe"],   reco_dPp_infer,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPRi_q2qe"],   reco_dPr_infer,q2qe, isData, evt,wgt);

      FillHDL2D( "h_dthetaRq2qe_ptmu", utils, reco_dThetaR, muon_pt_beam, q2qe, isData, evt, wgt );
      FillHDL2D( "h_dthetaPq2qe_ptmu", utils, reco_dThetaP, muon_pt_beam, q2qe, isData, evt, wgt );
      //------------------------------------------------------
      //Fill 2D Vertical Error Bands
      //Error Bands are not filled up for tmu_costheta 
      //and pmu_ptmu histos to avoid code slow down 
      //------------------------------------------------------
      if(RunCodeWithSystematics && !isData){
        utils->fillVertErrorBands( utils->histos1D["h_ptmu"],  muon_pt_beam, evt );
        utils->fillVertErrorBands( utils->histos1D["h_q2qe"],  q2qe, evt );
        utils->fillVertErrorBands( utils->histos2D["h_q2qe_ptmu"],  q2qe, muon_pt_beam, evt );
        
        utils->fillVertErrorBands( utils->histos2D["h_dalphat_ptmu"],reco_dalphat,muon_pt_beam,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dphit_ptmu"],reco_dphit,muon_pt_beam,evt);
        utils->fillVertErrorBands( utils->histos2D["h_pn_ptmu"],reco_pn,muon_pt_beam,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dpt_ptmu"],reco_dpt,muon_pt_beam,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dptx_ptmu"],reco_dptx,muon_pt_beam,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dpty_ptmu"],reco_dpty,muon_pt_beam,evt);
        utils->fillVertErrorBands( utils->histos2D["h_tp_ptmu"],protonTn, muon_pt_beam,evt);
        utils->fillVertErrorBands( utils->histos2D["h_ptheta_ptmu"],protonAngle, muon_pt_beam,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_sign_ptmu"],reco_sign,muon_pt_beam,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_signdalphat_ptmu"],reco_signed_dalphat,muon_pt_beam,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_signdphit_ptmu"],reco_signed_dphit,muon_pt_beam,evt);

        utils->fillVertErrorBands( utils->histos2D["h_dthetaP_ptmu"],reco_dThetaP,muon_pt_beam,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dthetaR_ptmu"],reco_dThetaR,muon_pt_beam,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPP_ptmu"],reco_dPp,muon_pt_beam,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPR_ptmu"],reco_dPr,muon_pt_beam,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPPi_ptmu"],reco_dPp_infer,muon_pt_beam,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPRi_ptmu"],reco_dPr_infer,muon_pt_beam,evt);

        utils->fillVertErrorBands( utils->histos2D["h_dthetaP_q2qe"],reco_dThetaP,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dthetaR_q2qe"],reco_dThetaR,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPP_q2qe"],reco_dPp,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPR_q2qe"],reco_dPr,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPPi_q2qe"],reco_dPp_infer,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPRi_q2qe"],reco_dPr_infer,q2qe,evt);


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
      }//end run with systematics        
    } // End of fillHistosV3 and fillVertErrorBands where CV of events pass the proton and recoil cuts 
  }
  return n_evts;
}


double AnalysisLooper(CCQENuEvent* evt, CCQENuUtils* utils, CCQENuCuts* cutter, NeutronBlobCuts* ncutter, string sample, int multiplicity, int n_entries, map<int,int>&cutsC, bool isData=false, bool isTruth=false, int npMode = 1 )
{ 
  return AnalysisLooperBase(evt, utils, cutter, ncutter, sample, multiplicity,n_entries, cutsC, isData, isTruth, npMode );
}

int MuonSelectionHists(string filename, string sample, bool makeFluxConstraintHisto, int multiplicity, string playlist, int n_mcfiles = -1, int n_datafiles = -1, int npMode = 1)
{
  double mc_events, data_events;
  map<int,int> mc_counter, data_counter;
  //---------------------------------------------
  // create file to store histograms
  //---------------------------------------------
  //TFile *f = new TFile( filename.c_str(), "RECREATE" );
  CCQENuUtils  *utils  = new CCQENuUtils( false, makeFluxConstraintHisto );
  utils->setPlaylist(playlist);
  utils->setFluxReweighterPlaylist();
  CCQENuCuts   *cutter = new CCQENuCuts();
  NeutronBlobCuts  *ncutter = new NeutronBlobCuts();


  //--------------------------------------------
  // Define histograms
  //--------------------------------------------
  AnaBinning *binner = new AnaBinning();
  CCQENuBinning *minmodbinner = new CCQENuBinning();

  NeutronBlobBinning *neutbinner = new NeutronBlobBinning();
  GeoUtils *geoUtils = new GeoUtils();

  axis_binning thetabins      = minmodbinner->GetMuonThetaBins();
  axis_binning costhetabins   = minmodbinner->GetMuonCosThetaBinsMiniBoone();
  axis_binning muonTbins      = binner->GetMuonEnergyBinsGeV();
  axis_binning muonPbins      = binner->GetMuonEnergyBinsGeV();
  //axis_binning muonPbins      = binner->GetMuonEnergyUniformBinsGeV();
  //axis_binning muonPtbins      = minmodbinner->GetMuonPtUniformBinsGeV();
  axis_binning muonPtbins     = minmodbinner->GetMuonPtBinsGeV();
  axis_binning muonPzbins     = minmodbinner->GetMuonPzBinsGeV();
  axis_binning Q2bins         = minmodbinner->GetQ2BinsGeV();
  axis_binning Enubins        = minmodbinner->GetMuonPzBinsGeV();
  axis_binning blobEnergyBinsMeV= minmodbinner->GetBlobEnergyBinsMeV();// up to 500 MeV: Events / 10 MeV
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


  axis_binning dthetaPerpbins = neutbinner->GetThetaPerpBinsDegreeFine();
  axis_binning dthetaReactbins = neutbinner->GetThetaPerpBinsDegreeFine();
  
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




  //---------------------------------------------
  // Book Histograms
  //---------------------------------------------
  //

  //1D

  MnvH2D *h_dalphat_ptmu[nHistos],*h_dphit_ptmu[nHistos],*h_pn_ptmu[nHistos],*h_dpt_ptmu[nHistos],*h_dptx_ptmu[nHistos],*h_dpty_ptmu[nHistos],*h_tp_ptmu[nHistos],*h_ptheta_ptmu[nHistos],*h_sign_ptmu[nHistos],*h_signdalphat_ptmu[nHistos],*h_signdphit_ptmu[nHistos];
  MnvH2D *h_node_nodeenergy[nHistos],*h_node_nodeenergy_long[nHistos],*h_node_nodeenergy_anchored[nHistos],*h_node_nodeenergy_vest[nHistos];
  MnvH2D *h_protonRes_nodeenergy_node01[nHistos],*h_protonRes_nodeenergy_node2[nHistos],*h_protonRes_nodeenergy_node3[nHistos],*h_protonRes_nodeenergy_node4[nHistos],*h_protonRes_nodeenergy_node5[nHistos];


  utils->bookHistos( h_dalphat_ptmu, "h_dalphat_ptmu","dat:ptmu",dalphatbins,muonPtbins);
  utils->bookHistos( h_dphit_ptmu, "h_dphit_ptmu","dpt:ptmu",dphitbins,muonPtbins);
  utils->bookHistos( h_pn_ptmu, "h_pn_ptmu","pn:ptmu",pnbins,muonPtbins);
  utils->bookHistos( h_dpt_ptmu, "h_dpt_ptmu","dtp:ptmu",dptbins,muonPtbins);
  utils->bookHistos( h_dptx_ptmu, "h_dptx_ptmu","dtpx:ptmu",dptxbins,muonPtbins);
  utils->bookHistos( h_dpty_ptmu, "h_dpty_ptmu","dtpy:ptmu",dptybins,muonPtbins);
  utils->bookHistos( h_tp_ptmu, "h_tp_ptmu","proton kin:ptmu",protonKineticbins,muonPtbins);
  utils->bookHistos( h_ptheta_ptmu, "h_ptheta_ptmu", "protonang:ptmu",protonThetabins,muonPtbins);
  //utils->bookHistos( h_sign_ptmu, "h_signed_ptmu","sign:ptmu",signbins,muonPtbins);
  //utils->bookHistos( h_signdalphat_ptmu,"h_signeddalphat_ptmu","sdalphat:ptm",signeddalphatbins,muonPtbins);
  //utils->bookHistos( h_signdphit_ptmu,"h_signeddphit_ptmu","sdphit:ptm",signeddphitbins,muonPtbins);

  MnvH2D *h_dthetaR_ptmu[nHistos], *h_dthetaP_ptmu[nHistos], *h_dPP_ptmu[nHistos], *h_dPR_ptmu[nHistos], *h_dPPi_ptmu[nHistos], *h_dPRi_ptmu[nHistos];
  MnvH2D *h_dthetaR_q2qe[nHistos], *h_dthetaP_q2qe[nHistos], *h_dPP_q2qe[nHistos], *h_dPR_q2qe[nHistos], *h_dPPi_q2qe[nHistos], *h_dPRi_q2qe[nHistos];

  utils->bookHistos( h_dthetaR_ptmu, "h_dthetaR_ptmu", "dthetaR:ptmu", dthetaReactbins, muonPtbins);
  utils->bookHistos( h_dthetaP_ptmu, "h_dthetaP_ptmu", "dthetaP:ptmu", dthetaPerpbins, muonPtbins);
  utils->bookHistos( h_dPR_ptmu, "h_dPR_ptmu", "dPR:ptmu", dptybins, muonPtbins);
  utils->bookHistos( h_dPP_ptmu, "h_dPP_ptmu", "dPP:ptmu", dptxbins, muonPtbins);
  utils->bookHistos( h_dPRi_ptmu, "h_dPRi_ptmu", "dPRi:ptmu", dptybins, muonPtbins);
  utils->bookHistos( h_dPPi_ptmu, "h_dPPi_ptmu", "dPPi:ptmu", dptxbins, muonPtbins);

  utils->bookHistos( h_dthetaR_q2qe, "h_dthetaR_q2qe", "dthetaR:q2qe", dthetaReactbins, Q2bins);
  utils->bookHistos( h_dthetaP_q2qe, "h_dthetaP_q2qe", "dthetaP:q2qe", dthetaPerpbins, Q2bins);
  utils->bookHistos( h_dPR_q2qe, "h_dPR_q2qe", "dPR:q2qe", dptybins, Q2bins);
  utils->bookHistos( h_dPP_q2qe, "h_dPP_q2qe", "dPP:q2qe", dptxbins, Q2bins);
  utils->bookHistos( h_dPRi_q2qe, "h_dPRi_q2qe", "dPRi:q2qe", dptybins, Q2bins);
  utils->bookHistos( h_dPPi_q2qe, "h_dPPi_q2qe", "dPPi:q2qe", dptxbins, Q2bins);

  utils->bookHistos( h_node_nodeenergy, "h_node_nodeenergy","node:energy",50,0,50,500,0,50);
  utils->bookHistos( h_node_nodeenergy_long, "h_node_nodeenergy_long","node:energy",50,0,50,500,0,50);
  utils->bookHistos( h_node_nodeenergy_anchored, "h_node_nodeenergy_anchored","node:energy",50,0,50,500,0,50);
  utils->bookHistos( h_node_nodeenergy_vest, "h_node_nodeenergy_vest","node:energy",50,0,50,500,0,50);

  utils->bookHistos( h_protonRes_nodeenergy_node01, "h_protonRes_nodeenergy_node01","res:node",500,0,500,100,-1,1);
  utils->bookHistos( h_protonRes_nodeenergy_node2, "h_protonRes_nodeenergy_node2","res:node",500,0,500,100,-1,1);
  utils->bookHistos( h_protonRes_nodeenergy_node3, "h_protonRes_nodeenergy_node3","res:node",500,0,500,100,-1,1);
  utils->bookHistos( h_protonRes_nodeenergy_node4, "h_protonRes_nodeenergy_node4","res:node",500,0,500,100,-1,1);
  utils->bookHistos( h_protonRes_nodeenergy_node5, "h_protonRes_nodeenergy_node5","res:node",500,0,500,100,-1,1);


  axis_binning dthetaRq2qeBins, dthetaPq2qeBins;
  MnvH2D *h_dthetaRq2qe_ptmu[nHistos], *h_dthetaPq2qe_ptmu[nHistos];

  HyperDimLinearizer *hdl_dthetaRq2qe_ptmu = DeclareHDLHistos2D( utils, h_dthetaRq2qe_ptmu, "h_dthetaRq2qe_ptmu", "pt vs dthetaRq2qe: dthetaR : q2qe", dthetaReactbins, muonPtbins, Q2bins, dthetaRq2qeBins );
  HyperDimLinearizer *hdl_dthetaPq2qe_ptmu = DeclareHDLHistos2D( utils, h_dthetaPq2qe_ptmu, "h_dthetaPq2qe_ptmu", "pt vs dthetaPq2qe: dthetaP : q2qe", dthetaPerpbins, muonPtbins, Q2bins, dthetaPq2qeBins );




  //utils->SetHyperDim(my3d);  
 
  // Add Vertical and Lateral Error Bands
  // JO is removing error bands from all 1D because they are all 
  // present in 2D. Plz make projections directly from the 2D. 
  //--------------------------------------------------------------



  TChain *truth_mc = utils->getMCTree("Truth", n_mcfiles );
  
  //---------------------------------------------
  // Get MC Tree
  //---------------------------------------------
  TChain* tree_mc   =  utils->getMCTree("CCQENu", n_mcfiles );
  int entries_mc   = tree_mc->GetEntries();
  utils->setmnvHadronReweightTruthTree(truth_mc);
  utils->setmnvHadronReweightDataTree(tree_mc);
  cout << "MC entries: " << entries_mc << endl;
  CCQENuEvent* mc = new CCQENuEvent( tree_mc );


  //-----------------------------------------------------------
  // Running over event selection criteria for each MC event 
  //-----------------------------------------------------------
  bool isData = false;
  bool isTruth= false;
  cout<< "AnalysisLooper will run in NPMode = "<<npMode<<endl;
  mc_events=AnalysisLooper( mc, utils, cutter, ncutter, sample,  multiplicity, entries_mc, mc_counter, isData, isTruth, npMode );
  cout << "Finished looping over all MC entries. Number of weighted MC events satisfying all selection criteria = " << mc_events << endl; 
  delete mc; 
  delete tree_mc; 
 
  //==================================================================
  // Create ROOT file to store histograms
  // Write to file all the created histograms 
  //==================================================================
  TFile *f = new TFile( filename.c_str(), "RECREATE" );
   
  //Write MC and DATA POT to file
  utils->writePOT( f );
  
  //Write multiplicity and sample definitions
  TMap *sampleMap = new TMap(2,0);
  string multiplicity_label = (multiplicity==0)? "alltracks": Form("%i_track",multiplicity);
  sampleMap->Add(new TObjString("multiplicity"), new TObjString(multiplicity_label.c_str()) );
  sampleMap->Add(new TObjString("sample"), new TObjString(sample.c_str()) );
  f->WriteTObject( sampleMap, "sample");
  
  //Write DATA and MC number of events
  TVector2 *evt = new TVector2( data_events, mc_events );
  f->WriteTObject( evt, "n_events" );

  //Write out whether this ROOT file contains histograms that will be used for the flux constraint  
  TParameter<bool> *fluxConstrain = new TParameter<bool>( "fluxConstraintHisto", makeFluxConstraintHisto ); 
  f->WriteTObject( fluxConstrain ); 
  
  WriteHistos( utils, f );
  //f->cd();
  //for( unsigned int i=kData; i<nHistos; ++i ){
  //    h_multiplicity[i]->Write();
  //    h_node_nodeenergy[i]->Write();
  //    h_nodes[i]->Write();
  //    h_ptmu[i]->Write();
  //    
  //    h_node_nodeenergy_long[i]->Write();
  //    h_nodes_long[i]->Write();
  //    h_node_nodeenergy_anchored[i]->Write();
  //    h_nodes_anchored[i]->Write();
  //    h_node_nodeenergy_vest[i]->Write();
  //    h_nodes_vest[i]->Write();
  //    
  //    h_protonRes_nodeenergy_node01[i]->Write();
  //    h_protonRes_nodeenergy_node2[i]->Write();
  //    h_protonRes_nodeenergy_node3[i]->Write();
  //    h_protonRes_nodeenergy_node4[i]->Write();
  //    h_protonRes_nodeenergy_node5[i]->Write();
  //    
  //    
  //    h_dalphat_ptmu[i]->Write();
  //    h_dphit_ptmu[i]->Write();
  //    h_pn_ptmu[i]->Write();
  //    h_dpt_ptmu[i]->Write();
  //    h_dptx_ptmu[i]->Write();
  //    h_dpty_ptmu[i]->Write();
  //    h_tp_ptmu[i]->Write();
  //    h_ptheta_ptmu[i]->Write();
  //    h_sign_ptmu[i]->Write();
  //    h_signdalphat_ptmu[i]->Write();
  //    h_signdphit_ptmu[i]->Write();

  //    h_dthetaR_ptmu[i]->Write();
  //    h_dthetaP_ptmu[i]->Write();
  //    h_dPP_ptmu[i]->Write();
  //    h_dPR_ptmu[i]->Write();
  //    h_dPPi_ptmu[i]->Write();
  //    h_dPRi_ptmu[i]->Write();
  //  }
  //}
   
  f->Close();
  delete f;
  delete utils;
  delete cutter; 

  return 0;
}

//HyperDimLinearizer* DeclareHDLHistos2D( CCQENuUtils* utils, MnvH2D **h, std::string name, std::string title, axis_binning &xbins, axis_binning &ybins, axis_binning &zbins, axis_binning &global_x )
//{
//  std::vector< std::vector<double> > bins3D({xbins.bin_edges, ybins.bin_edges, zbins.bin_edges});
//
//  HyperDimLinearizer *hdl = new HyperDimLinearizer(bins3D,0);
//
//  int n_bins = (xbins.bin_edges.size()+1)*(zbins.bin_edges.size()+1 );
//  vector<double> bins;
//  for( int i = 0; i<n_bins; ++i ) bins.push_back( i );
//
//  global_x.bin_edges = bins;
//  global_x.nbins = bins.size()-1;
//  global_x.min = bins.front();
//  global_x.max = bins.back();
//  utils->bookHistos( h, name.c_str(), title.c_str(), global_x, ybins );
//  utils->HDLMap[ name.c_str() ] = hdl;
//  return hdl;
//}
//
//void FillHDL2D( std::string name, CCQENuUtils *utils, double x, double y, double z, bool isData, CCQENuEvent* evt, double wgt )
//{
//  //std::vector<double> values({x,y,z});
//  HyperDimLinearizer* hdl = utils->HDLMap[ name.c_str() ];
//  double globalx = hdl->GetBin( std::vector<double>({x,y,z}) ).first+0.0001;
//  utils->fillHistosV3( utils->histos2D[ name.c_str() ], globalx, y, isData, evt, wgt );
//
//  if( RunCodeWithSystematics && !isData)
//  {
//	  utils->fillVertErrorBands( utils->histos2D[ name.c_str() ],globalx, y,evt);
//  }
//}


int main( int argc, char *argv[])
{
  ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if (argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t-./MuonSelectionHists  Name_and_path_to_Output_file  Selection_Signal  Make_Flux_Constraint_Histo  Multiplicity  Number_MC_files  Number_DATA_files \n\n"<<
      "\t-Name_and_path_to_Output_file\t =\t Name of and path to the Output ROOT file that will be created \n"<<
      "\t-Playlist\t =\t Name of the playlist you want to run over e.g. minerva1 \n"<<
      "\t-Selection_Signal\t =\t Can be: Signal, MichelSideBand, SideBand (this is:non-vtx vs Q2 sideband) \n"<<
      "\t-Make_Flux_Constraint_Histo\t =\t If TRUE Enter 1; If FALSE Enter 0 \n"<<
      "\t-Multiplicity\t =\t Enter 0, 1 or 2. Here 0: 1 or 2 tracks; 1: 1-track events only; 2: 2-track events only \n"<<
      "\t-Number_MC_files\t =\t Number of MonteCarlo files. To use all files, set this to: -1 \n"<<
      "\t-Number_DATA_files\t =\t Number of Data files. To use all files, set this to: -1\n" <<
      "\t-Variable combo\t = \tChoose TKI vs TKI combinations, default is off (-1)" << 
      "\t-Proton Neutron Mode\t = \t1: proton, 0: neutron, 10: proton+neutron, default is (1)" << std::endl;
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    return 0; 
  }
  
  //! Default parameters
  std::vector<std::string> par;
  par.push_back("MuonSelectionHists");
  par.push_back( Form("%s/ana/rootfiles/MuonSelectionHists.root",getenv("CCQENUROOT") ) );
  par.push_back("minerva1");
  par.push_back("Signal");
  par.push_back("0");
  par.push_back("0");
  par.push_back("-1");
  par.push_back("-1");
  par.push_back("1");//proton, neutron,proton/neutron mode

  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }

  bool fluxConstraintHistoNeeded = ( par.at(4) == "1" ) ? true : false; 

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;

  return MuonSelectionHists(par[1], par[3], fluxConstraintHistoNeeded, atoi(par[5].c_str()), par[2],atoi(par[6].c_str()), atoi(par[7].c_str()), atoi(par[8].c_str()) );
}
