#include "include/GlobalParameters.h"
#include "include/CCQENuUtils.h"
#include "TParameter.h" 
#include "PlotUtils/HyperDimLinearizer.h"

#include "include/GeneralFunc.h" //in CCQENuNeutron/include/
#include "include/EventCuts.h"

using namespace CCQENU_ANA;

bool myRunCodeWithSystematics = false;

//forward declaration
//HyperDimLinearizer* DeclareHDLHistos2D( CCQENuUtils* utils, MnvH2D **h, std::string name, std::string title, axis_binning &xbins, axis_binning &ybins, axis_binning &zbins, axis_binning &global_x );
//void FillHDL2D( std::string name, CCQENuUtils *utils, double x, double y, double z, bool isData, CCQENuEvent* evt, double wgt );
//================================================

//int npMode = 1;//0:neutron, 1:proton, 10:neutron+proton


double SumTotalE( std::vector<NeutronBlob*> &blobs )
{
  double e = 0;
  for( std::vector<NeutronBlob*>::iterator it = blobs.begin(); it!= blobs.end(); ++it ) e+= (*it)->TotalE;
  return e;
}


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
    XYZVector true_neutron_p_comp_XYZ;
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
        true_neutron_p_comp_XYZ.SetXYZ( true_neutron_p_comp.X(),true_neutron_p_comp.Y(),true_neutron_p_comp.Z());
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
    TVector3 primaryNucleonVect = ( npMode == 0 )? protonVect : neutronDir.Unit()*expectedNucleon3P.R() ;

    XYZVector primaryNucleonVect_XYZ( primaryNucleonVect.X(), primaryNucleonVect.Y(), primaryNucleonVect.Z() );

    //Calculating Q2qe:
    int charge = (GlobalParameters::Get().neutrinoMode)? -1:1;
    double q2qe = physicsUtils->qsqCCQE( beam, reco_muon_4P, charge, binding_e ); //GeV^2



    //preparing for STKI and angular vars:
    double Mcarbon = 6*fs_part_mass+6*is_part_mass-92.163/1E3;
    double Enuc=sqrt(primaryNucleonVect_XYZ.Mag2()+fs_part_mass*fs_part_mass);
    XYZTVector primaryNucleonVect_XYZT( protonVect.X(), protonVect.Y(), protonVect.Z(), Enuc );

    //calculating transverse variables
    std::vector<double> transverseVars = geoUtils->GetTransverseVariables( beam, reco_muon_4P, primaryNucleonVect_XYZT, Mcarbon, 27.13/1E3);
    //calculating angular variables
    std::vector<double> reco_neutron_vars = geoUtils->ComputeNeutronAngularVars( beam, expectedNucleon3P, primaryNucleonVect_XYZ );
    std::vector<double> true_ish_neutron_vars = geoUtils->ComputeNeutronAngularVars( beam, expectedNucleon3P, true_neutron_p_comp_XYZ );

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
    reco_dThetaP *= 180/3.14159;
    reco_dThetaR *= 180/3.14159;
    reco_dTheta *= 180/3.14159;
    double reco_dPp = reco_neutron_vars[3];
    double reco_dPr = reco_neutron_vars[4];
    double reco_dPp_infer = reco_neutron_vars[5];
    double reco_dPr_infer = reco_neutron_vars[6];


    double true_dThetaP = true_ish_neutron_vars[0];
    double true_dThetaR = true_ish_neutron_vars[1];
    double true_dTheta = true_ish_neutron_vars[2];
    true_dThetaP *= 180/3.14159;
    true_dThetaR *= 180/3.14159;
    true_dTheta *= 180/3.14159;
    double true_dPp = true_ish_neutron_vars[3];
    double true_dPr = true_ish_neutron_vars[4];
    double true_dPp_infer = true_ish_neutron_vars[5];
    double true_dPr_infer = true_ish_neutron_vars[6];

    //=======================================
    // Different recoil energy
    //=======================================

    //1 Search main candidate --- done automatically
    // form cone about main candidate
    //Cone Neutron energy

    XYZVector muonPWRTDetector = XYZVector( reco_muon_px, reco_muon_py, reco_muon_pz );
    XYZVector vtx(evt->vtx[0], evt->vtx[1], evt->vtx[2]);

    NeutronBlob mainNeutronCandidate = *ncutter->MainNeutronCandidate();
    XYZVector mcPos =  (XYZVector) mainNeutronCandidate.BlobPos.Vect();
    XYZVector blobVtxVect = mcPos - vtx;
    double blobVtxLength = blobVtxVect.R();


    Radiator muonRadiator;
    muonRadiator.AddCone( vtx, muonPWRTDetector, 15.);
    muonRadiator.AddCylinder( vtx, muonPWRTDetector, 100/*rho mm*/ );
    muonRadiator.AddBall( vtx, 100 );

    Radiator neutronRadiator1;
    neutronRadiator1.AddCylinder( vtx, blobVtxVect, 0, blobVtxLength );
    neutronRadiator1.AddCone( mcPos, blobVtxVect, 15 );




    vector<NeutronBlob*>::iterator itBlob;
    vector<NeutronBlob*> all2D_blobs = ncutter->GetCandidates()->TwoDBlobs;
    vector<NeutronBlob*> all3D_blobs = ncutter->GetCandidates()->ThreeDBlobs;
    vector<NeutronBlob*> all_blobs = ncutter->GetCandidates()->AllBlobs;
    std::vector<NeutronBlob*> blobs_in_r1;
    std::vector<NeutronBlob*> blobs_in_muon;
    std::vector<NeutronBlob*> blobs_unused, blobs_unused3D;
    for ( itBlob = all_blobs.begin(); itBlob != all_blobs.end(); ++itBlob )
    {
      if ( muonRadiator.PassRadiator( (**itBlob) ) )
      {
        blobs_in_muon.push_back( *itBlob );
      } else if ( neutronRadiator1.PassRadiator( (**itBlob) ) )
      {
        blobs_in_r1.push_back( *itBlob );
      } else
      {
        blobs_unused.push_back( *itBlob );
        if( (*itBlob)->Is3D) blobs_unused3D.push_back( *itBlob );
      }
    }
    double e_muon = SumTotalE( blobs_in_muon );
    double e_neut = SumTotalE( blobs_in_r1 );
    double e_recoil = SumTotalE( blobs_unused );
    double nblobs_ununsed = blobs_unused.size();
    double nblobs_ununsed_3D = blobs_unused3D.size();



    //I feel there is also a need  to include the left-right uncertainty of the beam, since LR-direction is pretty important....

    if (fill_common_histos)
    {
      //---------------------------------------------
      //Fill 1-D Plots
      //---------------------------------------------
      //cout<<"Fill start"<<endl;
      utils->fillHistosV3( utils->histos1D["h_neutronRes"],neutron_res_angle ,  isData, evt, wgt ); 
      utils->fillHistosV3( utils->histos1D["h_neutronResX"],neutron_res_angle_x ,  isData, evt, wgt ); 
      utils->fillHistosV3( utils->histos1D["h_neutronResY"],neutron_res_angle_y ,  isData, evt, wgt ); 
          
      //-------------------------------------------------------
      //Fill 2-D Plots
      //-------------------------------------------------------
      //h_true_vs_reco_dthetaR
      utils->fillHistosV3( utils->histos2D["h_true_vs_reco_dthetaP"], true_dThetaP, reco_dThetaP, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_true_vs_reco_dthetaR"], true_dThetaP, reco_dThetaR, isData, evt,wgt);

      utils->fillHistosV3( utils->histos2D["h_study_recoilNew_qsq"], e_recoil, q2qe, isData, evt, wgt );
      utils->fillHistosV3( utils->histos2D["h_study_recoilNew_neutronE"], e_recoil, e_neut, isData, evt, wgt );
      utils->fillHistosV3( utils->histos2D["h_study_recoilNew_dThetaP"], e_recoil, reco_dThetaP, isData, evt, wgt );
      utils->fillHistosV3( utils->histos2D["h_study_recoilNew_nblobs"], e_recoil, nblobs_ununsed, isData, evt, wgt );
      utils->fillHistosV3( utils->histos2D["h_study_recoilNew_nblobs3D"], e_recoil, nblobs_ununsed_3D, isData, evt, wgt );
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


  axis_binning dthetaPerpbins = neutbinner->GetThetaPerpBinsDegree();
  axis_binning dthetaReactbins = neutbinner->GetThetaReactionPlaneBinsDegree();
  
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


  std::vector<double> bins;
  bins.clear();
  axis_binning uniformbins_0_3rad;
  for(float i=0;i<=31;i+=1) bins.push_back(i*0.1);
  uniformbins_0_3rad.bin_edges = bins;
  uniformbins_0_3rad.nbins = bins.size()-1;
  uniformbins_0_3rad.min = bins.front();
  uniformbins_0_3rad.max = bins.back();

  //1D
  MnvH1D *h_neutronRes[nHistos], *h_neutronResX[nHistos], *h_neutronResY[nHistos];
  utils->bookHistos( h_neutronRes, "h_neutronRes", "NeutRes",uniformbins_0_3rad );
  utils->bookHistos( h_neutronResX, "h_neutronResX", "NeutResX",uniformbins_0_3rad );
  utils->bookHistos( h_neutronResY, "h_neutronResY", "NeutResY",uniformbins_0_3rad );
  //2-D
  MnvH2D *h_true_vs_reco_dthetaR[nHistos], *h_true_vs_reco_dthetaP[nHistos];
  MnvH2D *h_dperp_dreactplane[nHistos], *h_dtheta_qsq[nHistos], *h_dperp_qsq[nHistos], *h_dreactplane_qsq[nHistos];

  utils->bookHistos( h_true_vs_reco_dthetaR, "h_true_vs_reco_dthetaR","h_true_vs_reco_dthetaR", dthetaReactbins, dthetaReactbins );
  utils->bookHistos( h_true_vs_reco_dthetaP, "h_true_vs_reco_dthetaP","h_true_vs_reco_dthetaP", dthetaPerpbins, dthetaPerpbins );



  axis_binning uniformbins_0_500MeV;
  bins.clear();
  for(float i=0;i<=500;i+=5) bins.push_back(i);
  uniformbins_0_500MeV.bin_edges = bins;
  uniformbins_0_500MeV.nbins = bins.size()-1;
  uniformbins_0_500MeV.min = bins.front();
  uniformbins_0_500MeV.max = bins.back();

  axis_binning uniformbins_0_10;
  bins.clear();
  for(float i=0;i<=11;i+=1) bins.push_back(i-.5);
  uniformbins_0_10.bin_edges = bins;
  uniformbins_0_10.nbins = bins.size()-1;
  uniformbins_0_10.min = bins.front();
  uniformbins_0_10.max = bins.back();


  //histograms looking at all recoil energy:
  MnvH2D *h_study_recoilNew_qsq[nHistos];
  MnvH2D *h_study_recoilNew_dThetaP[nHistos];
  MnvH2D *h_study_recoilNew_neutronE[nHistos];
  MnvH2D *h_study_recoilNew_nblobs[nHistos];
  MnvH2D *h_study_recoilNew_nblobs3D[nHistos];
  utils->bookHistos( h_study_recoilNew_qsq, "h_study_recoilNew_qsq","recoil,qsq: recoil :qsq", uniformbins_0_500MeV, Q2bins );
  utils->bookHistos( h_study_recoilNew_neutronE, "h_study_recoilNew_neutronE","recoil,neutronE: recoil :neutronE", uniformbins_0_500MeV, uniformbins_0_500MeV );
  utils->bookHistos( h_study_recoilNew_dThetaP, "h_study_recoilNew_dThetaP","recoil,dThetaP: recoil :dThetaP", uniformbins_0_500MeV, dthetaPerpbins );
  utils->bookHistos( h_study_recoilNew_nblobs, "h_study_recoilNew_nblobs","recoil,nblobs: recoil :nblobs", uniformbins_0_500MeV, uniformbins_0_10 );
  utils->bookHistos( h_study_recoilNew_nblobs3D, "h_study_recoilNew_nblobs3D","recoil,nblobs: recoil :nblobs", uniformbins_0_500MeV, uniformbins_0_10 );
  

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
 
  //Print out some information
  cout<<" --------------------------------------------- " << endl;
  cout<<" MC " << endl;
  cout<<" --------------------------------------------- " << endl;
  cout<<" ****************** " << endl;
  cout<<" MC 1 Track Events " << endl;
  cout<<" ****************** " << endl;
  cout<< Form("nocuts = %i; cut1 = %i, cut2 = %i, cut3 = %i, cut4 = %i, cut5 = %i, cut6 = %i, cut6_sideband = %i",mc_counter[10],mc_counter[11], mc_counter[12], mc_counter[13],mc_counter[14],mc_counter[15],mc_counter[16],mc_counter[17]) << endl;
  cout<<" ****************** " << endl;
  cout<<" MC 2 Track Events " << endl;
  cout<<" ****************** " << endl;
  cout<< Form("nocuts = %i; cut1 = %i, cut2 = %i, cut3 = %i, cut4 = %i, cut5 = %i, cut6 = %i, cut7 = %i, cut7_sideband = %i",mc_counter[20],mc_counter[21], mc_counter[22], mc_counter[23],mc_counter[24],mc_counter[25],mc_counter[26],mc_counter[27]) << endl;
  
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
//  if( myRunCodeWithSystematics && !isData)
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
