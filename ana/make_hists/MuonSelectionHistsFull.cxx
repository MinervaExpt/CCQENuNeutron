#include "include/GlobalParameters.h"
#include "include/CCQENuUtils.h"
#include "TParameter.h" 
#include "PlotUtils/HyperDimLinearizer.h"
#include "PlotUtils/modifier_Eb.h"

#include "include/Parameters.h"
#include "include/GeneralFunc.h" //in CCQENuNeutron/include/
#include "include/EventCuts.h"
#include "include/CommonBins.h"
#include "include/EventSelector.h"

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
using namespace CCQENU_ANA;



//forward declaration
//HyperDimLinearizer* DeclareHDLHistos2D( CCQENuUtils* utils, MnvH2D **h, std::string name, std::string title, axis_binning &xbins, axis_binning &ybins, axis_binning &zbins, axis_binning &global_x );
//void FillHDL2D( std::string name, CCQENuUtils *utils, double x, double y, double z, bool isData, CCQENuEvent* evt, double wgt );
//================================================

//int npMode = 1;//0:neutron, 1:proton, 10:neutron+proton


double SolveProtonP( double a )
{
  if (a>0.5) return (a-0.5);
  double A = 0.935673;
  double B = 4.77731/pow(  -247040. + 438615. *a + 662.28*pow( 139264. - 494080. *a + 438615. *a*a        ,0.5 )             ,1/3.);
  double C = 0.012599*pow(-247040. + 438615. *a + 662.28 *pow(139264. - 494080. *a + 438615.*a*a,0.5), 1/3.);
  return (A+B)/C;
}


double SolveNeutronP( XYZVector &reco_neutron_V3, XYZTVector &protonP )
{
  //XYZVector baseAxis(0,0,1);
  //AxisAngle R1( reco_neutron_V3.Cross( baseAxis ), ACos( baseAxis.Dot( reco_neutron_V3.Unit() ) ) );
  //XYZTVector protonR1 = (R1*protonP );
  double Mn = 0.9395654133, Mp=0.93827231;
  double Mn2 = Mn*Mn, Mp2 = Mp*Mp;
  double Ep = protonP.T(), Pz = protonP.Vect().Dot( reco_neutron_V3.Unit() );
  double Ep2 = Ep*Ep, Pz2 = Pz*Pz;
  double P = protonP.P();
  double P2 = P*P;
  double Pt = pow( P*P-Pz*Pz,0.5);
  double Pt2 = Pt*Pt;

  double A = Ep*Mp*Pz-Mp2*Pz;
  double B0 = 2*Ep*Mn2+Ep2*Mp2-2*Mn2*Mp2-2*Ep*Mp2*Mp+Mp2*Mp2-Mn2*Pt2;
  double B1 = 2*Ep*Mp-Mp2-Ep2;
  double B = pow(-B0*B1, 0.5);
  double C = 2*Ep*Mp-2*Mp2-Pt2;
  double ans0 = (A+B)/C;
  double ans1 = (A-B)/C;
  if (ans0>0) return ans0;
  return -9999;
}


int GetMultiplicity( NeutronCandidates* candidates, double degree, 
    vector< vector<NeutronBlob*>>& groupedBlobs, 
    vector< NeutronBlob* > & ungroupedBlobs,
    vector< Cone > &cones,
    double &ungroupedEnergy, int &sizeOfUngroupedBlobs )
{
  bool doPrint = false;
  if(doPrint) cout<<"___ GetMultiplicity ___"<<endl;
  groupedBlobs.clear();
  cones.clear();

  vector<NeutronBlob*> AllBlobs= candidates->ThreeDBlobs;
  if( AllBlobs.size() == 0 ) 
  {
    ungroupedEnergy = candidates->RecoilEnergy(1);
    ungroupedBlobs = candidates->TwoDBlobs;
    sizeOfUngroupedBlobs = ungroupedBlobs.size();
    return -1;
  }


  XYZVector vtx = (XYZVector) candidates->EventVtx.Vect();

  //3D loop
  vector<NeutronBlob*>::iterator itblob = AllBlobs.begin();
  while( itblob != AllBlobs.end() )
  {
    //loop through all cones to check if blob is contained.
    // at first blob, cones.begin() == cones.end(), so loop automatially ends.

    XYZVector axis = (XYZVector) ( (*itblob)->BlobBegPos - candidates->EventVtx ).Vect();
    vector<Cone>::iterator itcone = cones.begin();
    for( ; itcone != cones.end(); itcone++ )
    {
      if( itcone->InsideCone( *(*itblob) )  )
      {
        if( doPrint ) cout<<ACos( itcone->GetAxis().Unit().Dot( axis.Unit() ) )*180/TMath::Pi() <<endl;
        groupedBlobs[itcone-cones.begin()].push_back( *itblob );
        AllBlobs.erase( itblob );
        break;
      }
    }
    //create a new cone if the old cones don't contain this blob
    if( itcone == cones.end() )
    {
      cones.push_back( Cone( vtx, axis, degree ) );
      if(doPrint)cout<<"push back cone: "<<cones.back().GetDefaultAngleRadian()<<", with degree: "<<degree<<endl;
      groupedBlobs.push_back( vector<NeutronBlob*>({*itblob}) );
      AllBlobs.erase( itblob );
    }
    //now if Blobs not empty, restart the process
  }

  //Now start 2D loop
  for( vector<NeutronBlob*>::iterator itblob = candidates->TwoDBlobs.begin(); itblob!=candidates->TwoDBlobs.end(); itblob++ )
  {
    bool inCone = false;
    for( uint i = 0; i< cones.size(); i++ )
    {
      if( cones[i].InsideCone( *(*itblob) ) )
      {
        groupedBlobs[i].push_back( *itblob );
        inCone = true;
        break;
      }
    }
    if( !inCone )
    {
    }
  }

  sizeOfUngroupedBlobs = ungroupedBlobs.size();
  if(doPrint) cout<<"end______Multiplicity"<<endl;
  return groupedBlobs.size();
}

template<class T>
double AnalysisLooperBase(T* evt, CCQENuUtils* &utils, CCQENuCuts* cutter, NeutronBlobCuts* ncutter, string sample, int multiplicity, int n_entries, map<int,int> &cuts_counter, bool isData=false, bool isTruth=false, int npMode = 1)
{
  GeoUtils* geoUtils = new GeoUtils();
  PhysicsUtils* physicsUtils = new PhysicsUtils();
  KinematicsModifier* kineMod = new KinematicsModifier();

  double n_evts=0;
  //return 0;

  for( int i = 0; i< n_entries; ++i)
  {
    evt->GetEntry(i);
    if( IgnoreEvent( evt, isData )) continue;
    if (i%1000 == 0) cout<<"At Event "<<i<<" ("<<(i+1)*100./float(n_entries)<<"%)"<<endl;
    //==============Fill Neutron Blobs======================
    bool isMC = !isData;
    
    if( !dropClusters ) ncutter->UpdateCandidates( evt, isMC );
    else ncutter->UpdateCandidates( evt, isMC, discardFraction, discardThreshold );
    if(isDebug) cout<<"Updated ncutter"<<endl;
    //==============Event Selection========================

    //--------------------- fill low q2qe plots ---------------------

    //if(npMode == 0)
    //{
    //  if (EventCutWithoutNeutron( evt, cutter, ncutter, cuts_counter, multiplicity, sample ,  isData, npMode))
    //  {
    //    double vtxE = evt->recoil_energy_nonmuon_vtx100mm;
    //    double recoilE = evt->recoil_energy_nonmuon_nonvtx100mm;
    //    double wgt = (isData)? 1: utils->GetCVWeight( evt, sample );
    //    if (recoilE<1)
    //    {
    //      //cout<<recoilE<<endl;

    //      XYZTVector reco_muon_XYZT( evt->CCQENu_leptonE[0]/1000., evt->CCQENu_leptonE[1]/1000., evt->CCQENu_leptonE[2]/1000. ,evt->CCQENu_muon_E/1e3);
    //      XYZVector decayPoint(0.231135*1000, 45.368069*1000, -766.384058*1000);
    //      XYZVector vtx( evt->vtx[0], evt->vtx[1], evt->vtx[2] );
    //      ///XYZVector beam = geoUtils->BeamAxis( beam_bias );
    //      XYZVector beam = (vtx-decayPoint).Unit();
    //      int charge = (neutrinoMode)? -1:1;
    //      double binding_e = 0;
    //      double q2qe = physicsUtils->qsqCCQE( beam, reco_muon_XYZT, charge, binding_e ); //GeV^2
    //      //cout<<vtxE<<", "<<q2qe<<endl;
    //      utils->fillHistosV3( utils->histos2D["h_vtxEnergy_q2qe_lowRecoil"], vtxE, q2qe, isData,evt,wgt);
    //      //cout<<"filled"<<endl;
    //    }
    //  }
    //}

    bool PassEventCut = (EventCut( evt, cutter, ncutter, cuts_counter, multiplicity, sample ,  isData, npMode, useRecoilCut, orderPar));
    if(isDebug) cout<<PassEventCut<<endl;




    if (!PassEventCut) continue;

    std::vector< XYZTVector > particles_XYZT;
    bool PassRecoSel = PassReco(evt, ncutter, kineMod, physicsUtils,geoUtils, isMC, npMode, particles_XYZT );
    if(!PassRecoSel) continue;
    std::vector< TLorentzVector > particles_lorentz = Convert<TLorentzVector, XYZTVector>( particles_XYZT );

    //now add in neutron-->
    std::vector<bool> commonNeutronCuts = CommonNeutronCuts(evt,ncutter, orderPar);
    bool hasNC = commonNeutronCuts[0] && commonNeutronCuts[1];
    bool isForward = commonNeutronCuts[2] && hasNC;
    bool isBackward = commonNeutronCuts[3] && hasNC;
    bool hasProton = (npMode!=0);

    NeutronCandidates* ptrNeutronCandidates; 
    NeutronBlob* mainCandidate; 
    if (npMode != 1 || hasNC) 
    {
      ptrNeutronCandidates = ncutter->GetCandidates();
      mainCandidate = ncutter->MainNeutronCandidate();
    }



    bool fill_common_histos = true;
    //===============Reconstruction ========================
    //-----------
    //Vertex
    XYZVector vertex( evt->vtx[0],evt->vtx[1],evt->vtx[2]);
    
    XYZTVector neutrino_XYZT = particles_XYZT[0];
    XYZTVector reco_muon_XYZT = particles_XYZT[1]; XYZVector reco_muon_XYZ = (XYZVector) reco_muon_XYZT.Vect();
    XYZTVector reco_proton_XYZT = particles_XYZT[2]; XYZVector reco_proton_XYZ = (XYZVector) reco_proton_XYZT.Vect();
    XYZTVector reco_neutron_XYZT = particles_XYZT[3]; XYZVector reco_neutron_XYZ = (XYZVector) reco_neutron_XYZT.Vect();


    

    TVector3 reco_proton_V3 = Convert( reco_proton_XYZ );


    TVector3 reco_neutron_V3 = Convert( reco_neutron_XYZ );
    
    
    //double muon_T         =  evt->CCQENu_muon_T/1000.; //GeV
    double muon_T         =  reco_muon_XYZT.E() - reco_muon_XYZT.M();
    double reco_muon_px   =  reco_muon_XYZT.X(); 
    double reco_muon_py   =  reco_muon_XYZT.Y(); 
    double reco_muon_pz   =  reco_muon_XYZT.Z(); 
    double reco_muon_e    =  reco_muon_XYZT.E(); 
    double muon_p         =  reco_muon_XYZT.R();

    double muon_theta     = evt->muon_theta_allNodes[20]; // get theta from 20th node, note NX kludge fixes this...
    double ptmu = sin(muon_theta)*muon_p;
    double pzmu = cos(muon_theta)*muon_p;

    //Convert to degrees
    muon_theta*= 180. / TMath::Pi();

    // 1. beam and muon
    double beam_bias = 0;
    XYZVector beam = geoUtils->BeamAxis( beam_bias ); 
    TVector3 beamV3(beam.X(), beam.Y(), beam.Z() );


    XYZVector xcoord(1,0,0);
    XYZVector ycoord = beam.Unit().Cross(xcoord);
    double muon_phi_detector = ATan2( reco_muon_py, reco_muon_px )*180./TMath::Pi();




    double protonAngle, protonMom, protonTn;
    if (npMode != 0 )
    {
      protonAngle = evt->CCQENu_proton_theta*180./TMath::Pi();
      protonMom = reco_proton_XYZT.P();//GeV
      protonTn = reco_proton_XYZT.E() - reco_proton_XYZT.M();
    }


    //Truth Proton Variables
    //TVector3 true_proton_p_comp = utils->GetHighestTrueProtonMomentumComp(evt);
    TVector3 true_proton_p_comp(evt->proton_prong_4p[0],evt->proton_prong_4p[1],evt->proton_prong_4p[2]);
    double true_protonMom = true_proton_p_comp.Mag();
    double proton_res = 1-protonMom/true_protonMom;
    

    TVector3 true_fs_nucleon_p_comp;
    double true_FSPartMom=-999; 
    double fspart_res_angle=-999;
    double fspart_res_angle_x=-999;
    double fspart_res_angle_y=-999;
    double fspart_res_dthetaR = -999;
    double fspart_res_dthetaP = -999;
    //1. look for the most energetic fs neutron
    if( isMC )
    {
      int pdg = (npMode == 0)? 2112:2212;
      std::vector<TLorentzVector> TrueFSParts = GetAllParticlesMomenta( evt, pdg, 0 );
      //cout<<"Number of FS Part = "<<TrueFSParts.size()<<endl;
      if(TrueFSParts.size() > 0 )
      {
        true_fs_nucleon_p_comp = TrueFSParts[0].Vect();
        TVector3 hx = true_fs_nucleon_p_comp.Cross( beamV3 ).Unit();
        TVector3 hy = true_fs_nucleon_p_comp.Cross( hx ).Unit();

        true_FSPartMom = true_fs_nucleon_p_comp.Mag();
        fspart_res_angle = acos( true_fs_nucleon_p_comp.Unit().Dot( reco_neutron_V3.Unit() ) )*180/TMath::Pi();
        fspart_res_angle_x = acos( true_fs_nucleon_p_comp.Unit().Dot( (reco_neutron_V3.Dot( hx )*hx).Unit() ) )*180/TMath::Pi();
        fspart_res_angle_y = acos( true_fs_nucleon_p_comp.Unit().Dot( (reco_neutron_V3.Dot( hy )*hy).Unit() ) )*180/TMath::Pi();

      }
    }


    //Expected nucleon kinematics

    // 2. binding and is/fs mass
    int charge = (neutrinoMode)? -1:1;
    double binding_e = (npMode==0)? 0 : bindingE;
    double is_part_mass = (neutrinoMode)? 0.9395654133: 0.93827231; //neutron/proton
    double fs_part_mass = (neutrinoMode)? 0.93827231: 0.9395654133;

    // 3. expected nucleon

    XYZTVector expected_nucleon_XYZT = geoUtils->ComputeExpectedNucleon( beam, reco_muon_XYZT, is_part_mass, fs_part_mass, binding_e);
    XYZVector expected_nucleon_XYZ = (XYZVector) expected_nucleon_XYZT.Vect();


    // Set the primary nucleon momentum vec
    // proton/pn mode : proton. neutron mode: neutron direction
    TVector3 primary_nucleon_V3 = ( npMode == 0 )? reco_neutron_V3.Unit()*expected_nucleon_XYZ.R() : reco_proton_V3 ;
    XYZVector primary_nucleon_XYZ= Convert( primary_nucleon_V3 );

    //cout<<"Is proton primary? "<< (primary_nucleon_V3==reco_proton_V3)<<endl;


    //Calculating Q2qe:
    double q2qe = physicsUtils->qsqCCQE( beam, reco_muon_XYZT, charge, binding_e ); //GeV^2
    double enu = physicsUtils->nuEnergyCCQE( beam, reco_muon_XYZT, charge, binding_e );
    double q0 = enu - reco_muon_XYZT.E();
    double q3 = pow(q0*q0+q2qe,.5);


    //preparing for STKI and angular vars:
    double Mcarbon = 6*fs_part_mass+6*is_part_mass-92.163/1E3;
    double Enuc=sqrt(primary_nucleon_XYZ.Mag2()+fs_part_mass*fs_part_mass);

    XYZTVector primary_nucleon_XYZT( primary_nucleon_XYZ.X(), primary_nucleon_XYZ.Y(), primary_nucleon_XYZ.Z(), Enuc );

    //calculating transverse variables
    std::vector<double> reco_transverse_vars = geoUtils->GetTransverseVariables( beam, reco_muon_XYZT, primary_nucleon_XYZT, Mcarbon, 27.13/1E3);
    //calculating angular variables
    std::vector<double> reco_neutron_vars = geoUtils->ComputeNeutronAngularVars( beam, expected_nucleon_XYZ, primary_nucleon_XYZ );
    // test if the expected and muon cancel, they should

    double reco_En      = reco_transverse_vars[vEnu];
    double reco_pn      = reco_transverse_vars[vPn];
    double reco_dpt     = reco_transverse_vars[vdpT];
    double reco_dptx    = reco_transverse_vars[vdpTx];
    double reco_dpty    = reco_transverse_vars[vdpTy];
    double reco_dalphat = reco_transverse_vars[vdalphaT];
    double reco_dphit   = reco_transverse_vars[vdphiT];
    double reco_sign    = reco_transverse_vars[vsign];
    double reco_signed_dalphat = reco_sign*reco_dalphat ;
    double reco_signed_dphit = reco_dphit *reco_sign;

    double reco_dthetaP = reco_neutron_vars[vdthetaP];
    double reco_dthetaR = reco_neutron_vars[vdthetaR];
    double reco_dtheta = reco_neutron_vars[vdtheta];
    //if (isMC) cout<< evt->mc_intType<<endl;
    reco_dthetaP *= 180/TMath::Pi();
    reco_dthetaR *= 180/TMath::Pi();
    reco_dtheta *= 180/TMath::Pi();
    //cout<<reco_dthetaP<<endl;
    double reco_dPp = reco_neutron_vars[vdPp];
    double reco_dPr = reco_neutron_vars[vdPr];
    double reco_dPp_infer = reco_neutron_vars[vdPpi];
    double reco_dPr_infer = reco_neutron_vars[vdPri];

    int region = getRegions( reco_dthetaP, reco_dthetaR );
    int regionCode = (region == -1)? 99: region;
    //cout<<regionCode<<endl;

    //Get 2D angulars
    double dtheta2D = -99999999999;
    if(npMode == 0)
    {
      int view = mainCandidate->View;
      int view2;

      XYZVector posTClus= ptrNeutronCandidates->GetClusterTPosClosestToVtx( mainCandidate, view2 );
      XYZVector vtx2D = (ncutter->DetectorUtilsPtr()->ViewProjection2D( vertex, view2 ) );

      XYZVector blobDir2D = (mainCandidate->BlobPosT - vtx2D).Unit();
      
      XYZVector posT= (ncutter->DetectorUtilsPtr()->ViewProjection2D( mainCandidate->BlobBegPos.Vect(), view2 ));

      XYZVector nu2D =( ncutter->DetectorUtilsPtr()->ViewProjection2D( beam, view2 ) ).Unit();
      XYZVector expected2D = ( ncutter->DetectorUtilsPtr()->ViewProjection2D( expected_nucleon_XYZ, view2 ) ).Unit();


      double dtheta2D = geoUtils->ComputeNeutronAngle2D( nu2D, expected2D, blobDir2D)*180/TMath::Pi();
    }

    //=================SetTransverseParticles=========================================
    utils->setTransverseParticle( primary_nucleon_V3, 0, false );

    //=================Finally Get CV Weight==========================================
    if(isDebug) cout<<"Getting CV Weight:"<<endl;
    double wgt = (isData)? 1: utils->GetCVWeight( evt, sample );
    n_evts+=wgt;
    if(isDebug) cout<<wgt<<endl;




    //======================Fit Neutron ==================================
    string quadrant="";
    double dphitt = -9999;
    double dthetaRt = -9999, dthetaPt = -9999;
    double pneutron = -99999, pneutronErr = -99999;
    double pneutronRes = -99999;
    double EnuFit = -99999;
    double EnuResolution = -9999;
    double EnuQ2QEResolution = -9999;
    XYZVector fitNeutron = reco_neutron_XYZ;
    if( hasNC && hasProton  && npMode == 10) 
    {
      quadrant = GetQuadrantCode( beam, reco_muon_XYZ, reco_proton_XYZ, reco_neutron_XYZ );
      //cout<<quadrant<<endl;

      vector<double> fitMomentum = geoUtils->FitMomentum( beam, reco_muon_XYZ, reco_proton_XYZ, (XYZVector) reco_neutron_XYZ.Unit(), geoUtils->MnGeV() );
      pneutron = fitMomentum[0];
      if (pneutron>0)
      {
        pneutronErr = fitMomentum[1];
        XYZVector fitNeutron = reco_neutron_XYZ.Unit()*abs(fitMomentum[0]);

        double M = geoUtils->MnGeV();
        double E = pow(pneutron*pneutron+M*M,0.5);
        XYZTVector fitNeutron_XYZT( fitNeutron.X(),fitNeutron.Y(),fitNeutron.Z(),E);

        EnuFit = beam.Dot( reco_muon_XYZ+reco_proton_XYZ+fitNeutron );

        std::vector<double> ddreco_transverse_vars = geoUtils->GetTransverseVariables( beam, reco_muon_XYZT+primary_nucleon_XYZT, fitNeutron_XYZT, Mcarbon, 27.13/1E3);
        dphitt = ddreco_transverse_vars[vdphiT];

        if(isMC)
        {
          double enuTrue = evt->mc_incomingPartVec[3]/1000;
          EnuResolution = (EnuFit - enuTrue)/enuTrue;
          EnuQ2QEResolution = (enu - enuTrue)/enuTrue;
          TLorentzVector NeutronPart = utils->GetHighestTruePart4P( evt, 2112 );
          double neutronPTrue = NeutronPart.P()/1000;
          pneutronRes = (pneutron - neutronPTrue)/neutronPTrue;
        }
      }
    }


    //Recoil Characteristic
    double recoil =ncutter->GetDefaultRecoilEnergy(evt);
    double recoil_inc = evt->recoil_summed_energy[0]; //MeV

    double vtxE = evt->recoil_energy_nonmuon_vtx100mm;
    double recoilE = evt->recoil_energy_nonmuon_nonvtx100mm;

    //Blob Characteristics
    double blobDist, blobEnergy, nBlobs, n2DBlobs, n3DBlobs;
    double nonBlobEnergy;
    double proton_score=-1, pion_score=-1;
    bool inVtx = false;
    if( hasNC )
    {
      blobDist = mainCandidate->BlobVtxVec.R();
      blobEnergy = mainCandidate->TotalE;
      nonBlobEnergy = recoil - blobEnergy;
      n3DBlobs = ptrNeutronCandidates->ThreeDBlobs.size();
      n2DBlobs = ptrNeutronCandidates->TwoDBlobs.size();
      nBlobs = ptrNeutronCandidates->AllBlobs.size();
      inVtx = blobDist<200;
      if( mainCandidate->hasTrack )
      {
        proton_score = mainCandidate->dEdX_proton;
        pion_score = mainCandidate->dEdX_pion;

      }

      double Mn=0.9395654133;
      //cout<<"in code: "<< q2qe<<"\t"<<q2qe/2/Mn<<"\t"<<blobEnergy/1000.<<endl;
      //if(npMode==0 && q2qe/2/Mn < blobEnergy/1000. ) continue;
    }

  double nHits = 0; 
  double totalE = 0; 
  int size = evt->nonvtx_iso_blobs_distance_in_prong_sz;

  for(int i=0;i<size;i++){
    if(evt->nonvtx_iso_blobs_start_position_z_in_prong[i] > 4750){
      nHits+=evt->nonvtx_iso_blobs_n_hits_in_prong[i];
      totalE+=evt->nonvtx_iso_blobs_energy_in_prong[i];
    }
  }





    vector< vector<NeutronBlob*>> blobMult_015;
    vector< NeutronBlob* > blobUngrouped_015;
    vector< Cone > cones;
    double ungroupedE_015;
    int nUngroupedBlobs_015;

    double angle = 30;
    int blob_multiplicity_015=GetMultiplicity( ptrNeutronCandidates, angle, 
        blobMult_015, blobUngrouped_015, cones,
        ungroupedE_015, nUngroupedBlobs_015 );

    if (isMC)
    {
      //get dthetaR/P resolution
      double true_fsp_p = true_fs_nucleon_p_comp.Mag(); 
      XYZVector true_fsp_3p( true_fs_nucleon_p_comp.X(), true_fs_nucleon_p_comp.Y(), true_fs_nucleon_p_comp.Z()  );
      //calculating angular variables
      std::vector<double> mc_neutron_vars = geoUtils->ComputeNeutronAngularVars( beam, expected_nucleon_XYZ, true_fsp_3p );
      // test if the expected and muon cancel, they should

      double mc_dThetaP = mc_neutron_vars[vdthetaP];
      double mc_dThetaR = mc_neutron_vars[vdthetaR];
      //if (isMC) cout<< evt->mc_intType<<endl;
      mc_dThetaP *= 180/TMath::Pi();
      mc_dThetaR *= 180/TMath::Pi();
      fspart_res_dthetaR = reco_dthetaR - mc_dThetaR;
      fspart_res_dthetaP = reco_dthetaP - mc_dThetaP;


      if(fill_common_histos )
      {
        utils->fillHistosV3( utils->histos2D["h_true_vs_reco_dthetaP"], mc_dThetaP, reco_dthetaP, isData, evt, wgt);
        utils->fillHistosV3( utils->histos2D["h_true_vs_reco_dthetaR"], mc_dThetaR, reco_dthetaR, isData, evt, wgt);
      }

      int main_type = parsePartType( mainCandidate->MCPID );
      utils->fillHistosV3( utils->histos2D["h_mainBlobE_PartType"], mainCandidate->TotalE, main_type, isData, evt, wgt );
      utils->fillHistosV3( utils->histos2D["h_mainBlobMaxE_PartType"], mainCandidate->ClusMaxE, main_type, isData, evt, wgt );
      for( uint i = 0; i<ptrNeutronCandidates->AllBlobs.size(); i++ )
      {
        NeutronBlob* blob = ptrNeutronCandidates->AllBlobs[i];
        double blobE = blob->TotalE;
        double blobMaxE = blob->ClusMaxE;
        int part_type = parsePartType( blob->MCPID );
        utils->fillHistosV3( utils->histos2D["h_blobE_PartType"], blobE, part_type, isData, evt, wgt );
        utils->fillHistosV3( utils->histos2D["h_blobMaxE_PartType"], blobMaxE, part_type, isData, evt, wgt );
        if( blob->Is3D )
        {
          utils->fillHistosV3( utils->histos2D["h_blobE3D_PartType"], blobE, part_type, isData, evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_blobMaxE3D_PartType"], blobMaxE, part_type, isData, evt, wgt );
        }
      }
    }


    bool angle_03 = (abs(reco_dthetaP) < 3 && abs(reco_dthetaR) < 3 );
    bool angle_04 = (abs(reco_dthetaP) < 4 && abs(reco_dthetaR) < 4 );
    bool angle_05 = (abs(reco_dthetaP) < 5 && abs(reco_dthetaR) < 5 );
    bool angle_06 = (abs(reco_dthetaP) < 6 && abs(reco_dthetaR) < 6 );
    bool angle_07 = (abs(reco_dthetaP) < 7 && abs(reco_dthetaR) < 7 );
    bool angle_08 = (abs(reco_dthetaP) < 8 && abs(reco_dthetaR) < 8 );
    bool angle_09 = (abs(reco_dthetaP) < 9 && abs(reco_dthetaR) < 9 );
    bool angle_10 = (abs(reco_dthetaP) < 10 && abs(reco_dthetaR) < 10 );



    //========================== Histo Filling Proper ============================

    if (fill_common_histos)
    {
      //---------------------------------------------
      //Fill 1-D Plots
      //---------------------------------------------
      //cout<<"Fill start"<<endl;
      utils->fillHistosV3( utils->histos1D["h_multiplicity"], evt->multiplicity,  isData, evt, wgt ); 
      utils->fillHistosV3( utils->histos1D["h_ptmu"], ptmu,  isData, evt, wgt );
      utils->fillHistosV3( utils->histos1D["h_q2qe"], q2qe,  isData, evt, wgt );
      utils->fillHistosV3( utils->histos2D["h_q2qe_ptmu"], q2qe, ptmu, isData, evt, wgt );
      utils->fillHistosV3( utils->histos1D["h_enu"], reco_En,  isData, evt, wgt );
      utils->fillHistosV3( utils->histos1D["h_pmu"], muon_p, isData, evt, wgt );

      utils->fillHistosV3( utils->histos2D["h_nHits_totalE"], nHits, totalE, isData, evt, wgt );

      utils->fillHistosV3( utils->histos1D["h_dthetaP"], reco_dthetaP, isData, evt, wgt );
      utils->fillHistosV3( utils->histos1D["h_dthetaR"], reco_dthetaR, isData, evt, wgt );
      utils->fillHistosV3( utils->histos1D["h_pn"], reco_pn, isData, evt, wgt );

      utils->fillHistosV3( utils->histos2D["h_multiplicity_q2qe"], blob_multiplicity_015, q2qe, isData,evt,wgt);

      //cout<<"fill region :"<<regionCode<<endl;
      //utils->fillHistosV3( utils->histos2D[Form("h_vtxEnergy_q2qe_region_%02d",regionCode)],vtxE, q2qe,isData,evt,wgt);
      //cout<<"filled"<<endl;

      utils->fillHistosV3( utils->histos2D[Form( "h_q2qe_ptmu_region_%02d",regionCode) ], q2qe, ptmu, isData, evt, wgt );
      if( region != -1 )utils->fillHistosV3( utils->histos2D["h_q2qe_ptmu_region_non99"], q2qe, ptmu, isData, evt, wgt );


      utils->fillHistosV3( utils->histos2D[Form( "h_enu_q2qe_region_%02d",regionCode) ], enu, q2qe, isData, evt, wgt );
      if( region != -1 )utils->fillHistosV3( utils->histos2D["h_enu_q2qe_region_non99"], enu, q2qe, isData, evt, wgt );


      utils->fillHistosV3( utils->histos2D[Form( "h_muontheta_q2qe_region_%02d",regionCode) ], muon_theta, q2qe, isData, evt, wgt );
      if( region != -1 )utils->fillHistosV3( utils->histos2D["h_muontheta_q2qe_region_non99"], muon_theta, q2qe, isData, evt, wgt );

      if(isDebug) cout<<"Filled q2qe_ptmu_region"<<endl;


      //if( inVtx )
      //{
      //  if( region == 0) utils->fillHistosV3( utils->histos1D["h_q2qe_vtx_region_00"], q2qe, isData, evt, wgt );
      //  if( region == 1) utils->fillHistosV3( utils->histos1D["h_q2qe_vtx_region_01"], q2qe, isData, evt, wgt );
      //  if( region == 2) utils->fillHistosV3( utils->histos1D["h_q2qe_vtx_region_02"], q2qe, isData, evt, wgt );
      //  if( region == 3) utils->fillHistosV3( utils->histos1D["h_q2qe_vtx_region_03"], q2qe, isData, evt, wgt );
      //  if( region == 4) utils->fillHistosV3( utils->histos1D["h_q2qe_vtx_region_04"], q2qe, isData, evt, wgt );
      //  if( region == 5) utils->fillHistosV3( utils->histos1D["h_q2qe_vtx_region_05"], q2qe, isData, evt, wgt );
      //  if( region ==-1) utils->fillHistosV3( utils->histos1D["h_q2qe_vtx_region_99"], q2qe, isData, evt, wgt );
      //}
      //else
      //{
      //  if( region == 0) utils->fillHistosV3( utils->histos1D["h_q2qe_nonvtx_region_00"], q2qe, isData, evt, wgt );
      //  if( region == 1) utils->fillHistosV3( utils->histos1D["h_q2qe_nonvtx_region_01"], q2qe, isData, evt, wgt );
      //  if( region == 2) utils->fillHistosV3( utils->histos1D["h_q2qe_nonvtx_region_02"], q2qe, isData, evt, wgt );
      //  if( region == 3) utils->fillHistosV3( utils->histos1D["h_q2qe_nonvtx_region_03"], q2qe, isData, evt, wgt );
      //  if( region == 4) utils->fillHistosV3( utils->histos1D["h_q2qe_nonvtx_region_04"], q2qe, isData, evt, wgt );
      //  if( region == 5) utils->fillHistosV3( utils->histos1D["h_q2qe_nonvtx_region_05"], q2qe, isData, evt, wgt );
      //  if( region ==-1) utils->fillHistosV3( utils->histos1D["h_q2qe_nonvtx_region_99"], q2qe, isData, evt, wgt );

      //}
      
      int n_nodes = evt->CCQENu_proton_nodes_nodesNormE_sz;
      int pattrec = evt->CCQENu_proton_patternRec;
      //cout<<"Fill 1"<<endl;
      
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
      utils->fillHistosV3( utils->histos2D["h_dalphat_ptmu"],reco_dalphat,ptmu, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dphit_ptmu"],reco_dphit,ptmu, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_pn_ptmu"],reco_pn,ptmu, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dpt_ptmu"],reco_dpt,ptmu, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dptx_ptmu"],reco_dptx,ptmu, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dpty_ptmu"],reco_dpty,ptmu, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_tp_ptmu"],protonTn, ptmu, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_ptheta_ptmu"],protonAngle, ptmu, isData, evt,wgt);

      utils->fillHistosV3( utils->histos2D["h_dalphat_q2qe"],reco_dalphat,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dphit_q2qe"],reco_dphit,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_pn_q2qe"],reco_pn,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dpt_q2qe"],reco_dpt,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dptx_q2qe"],reco_dptx,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dpty_q2qe"],reco_dpty,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_tp_q2qe"],protonTn, q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_ptheta_q2qe"],protonAngle, q2qe, isData, evt,wgt);

      utils->fillHistosV3( utils->histos2D["h_muontheta_q2qe"], muon_theta, q2qe, isData, evt, wgt);
      utils->fillHistosV3( utils->histos2D["h_muontheta_ptmu"], muon_theta, ptmu, isData, evt, wgt);
      utils->fillHistosV3( utils->histos2D["h_muonmomentum_q2qe"], reco_muon_e, q2qe, isData, evt, wgt);
      utils->fillHistosV3( utils->histos2D["h_muonmomentum_ptmu"], reco_muon_e, ptmu, isData, evt, wgt);
      if(angle_10) utils->fillHistosV3( utils->histos2D["h_enu_q2qe"], enu, q2qe, isData, evt, wgt);
      if(angle_10) utils->fillHistosV3( utils->histos2D["h_enu_ptmu"], enu, ptmu, isData, evt, wgt);



      //utils->fillHistosV3( utils->histos2D["h_sign_ptmu"],reco_sign,ptmu, isData, evt,wgt);
      //utils->fillHistosV3( utils->histos2D["h_signdalphat_ptmu"],reco_signed_dalphat,ptmu, isData, evt,wgt);
      //utils->fillHistosV3( utils->histos2D["h_signdphit_ptmu"],reco_signed_dphit,ptmu, isData, evt,wgt);

      utils->fillHistosV3( utils->histos2D["h_dthetaP_ptmu"],reco_dthetaP,ptmu, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dthetaR_ptmu"],reco_dthetaR,ptmu, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPP_ptmu"],reco_dPp,ptmu, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPR_ptmu"],reco_dPr,ptmu, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPPi_ptmu"],reco_dPp_infer,ptmu, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPRi_ptmu"],reco_dPr_infer,ptmu, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dtheta2D_ptmu"],   dtheta2D,ptmu, isData, evt,wgt);

      utils->fillHistosV3( utils->histos2D["h_dthetaP_q2qe"],reco_dthetaP,  q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dthetaR_q2qe"],reco_dthetaR,  q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPP_q2qe"],    reco_dPp,      q2qe, isData,evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPR_q2qe"],    reco_dPr,      q2qe, isData,evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPPi_q2qe"],   reco_dPp_infer,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPRi_q2qe"],   reco_dPr_infer,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dtheta2D_q2qe"],   dtheta2D,q2qe, isData, evt,wgt);


      utils->fillHistosV3( utils->histos2D["h_dthetaP_dthetaR"],reco_dthetaP,  reco_dthetaR, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dthetaP_dthetaR_Fine"],reco_dthetaP,  reco_dthetaR, isData, evt,wgt);

      utils->fillHistosV3( utils->histos2D["h_recoil_inc_q2qe"], recoil_inc, q2qe, isData, evt, wgt );
      utils->fillHistosV3( utils->histos2D["h_recoil_inc_q3"], recoil_inc, q3, isData, evt, wgt );

      utils->fillHistosV3( utils->histos2D["h_res_dthetaP_ptmu"],   fspart_res_dthetaP ,ptmu, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dthetaR_ptmu"],   fspart_res_dthetaR ,ptmu, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dtheta_ptmu"],   fspart_res_angle ,ptmu, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dthetaX_ptmu"],   fspart_res_angle_x ,ptmu, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dthetaY_ptmu"],   fspart_res_angle_y ,ptmu, isData, evt,wgt);

      utils->fillHistosV3( utils->histos2D["h_res_dthetaP_q2qe"],   fspart_res_dthetaP ,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dthetaR_q2qe"],   fspart_res_dthetaR ,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dtheta_q2qe"],   fspart_res_angle ,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dthetaX_q2qe"],   fspart_res_angle_x ,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dthetaY_q2qe"],   fspart_res_angle_y ,q2qe, isData, evt,wgt);

      utils->fillHistosV3( utils->histos2D["h_res_dthetaP"],   fspart_res_dthetaP ,reco_dthetaP, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dthetaR"],   fspart_res_dthetaR ,reco_dthetaR, isData, evt,wgt);

      if(hasNC && hasProton)
      {
        utils->fillHistosV3( utils->histos2D[ Form("h_pn_q2qe_%s", quadrant.c_str())], reco_pn, q2qe, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_dphit_q2qe_%s", quadrant.c_str())], reco_dphit, q2qe, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_dalphat_q2qe_%s", quadrant.c_str())], reco_dalphat, q2qe, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_dpt_q2qe_%s", quadrant.c_str())], reco_dpt, q2qe, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_dptx_q2qe_%s", quadrant.c_str())], reco_dptx, q2qe, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_dpty_q2qe_%s", quadrant.c_str())], reco_dpty, q2qe, isData,evt,wgt);


        utils->fillHistosV3( utils->histos2D[ Form("h_pn_ptmu_%s", quadrant.c_str())], reco_pn, ptmu, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_dphit_ptmu_%s", quadrant.c_str())], reco_dphit, ptmu, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_dalphat_ptmu_%s", quadrant.c_str())], reco_dalphat, ptmu, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_dpt_ptmu_%s", quadrant.c_str())], reco_dpt, ptmu, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_dptx_ptmu_%s", quadrant.c_str())], reco_dptx, ptmu, isData,evt,wgt);
        utils->fillHistosV3( utils->histos2D[ Form("h_dpty_ptmu_%s", quadrant.c_str())], reco_dpty, ptmu, isData,evt,wgt);

        utils->fillHistosV3( utils->histos2D[ "h_fitP_q2qe" ], pneutron, q2qe, isData, evt, wgt );
        utils->fillHistosV3( utils->histos2D[ "h_fitPRes_q2qe" ], pneutronErr, q2qe, isData, evt, wgt );//reco
        utils->fillHistosV3( utils->histos2D[ "h_fitEnu_q2qe" ],EnuFit, q2qe, isData, evt, wgt );

        utils->fillHistosV3( utils->histos2D[ "h_dFitP_q2qe" ], pneutronRes, q2qe, isData, evt, wgt );//true
        utils->fillHistosV3( utils->histos2D[ "h_dFitEnu_q2qe" ],EnuResolution, q2qe, isData, evt, wgt );//true quantity
        utils->fillHistosV3( utils->histos2D[ "h_dEnu_q2qe" ],EnuQ2QEResolution, q2qe, isData, evt, wgt );//true quantity

        utils->fillHistosV3( utils->histos2D[ "h_fitP_ptmu" ], pneutron, ptmu, isData, evt, wgt );
        utils->fillHistosV3( utils->histos2D[ "h_fitPRes_ptmu" ], pneutronErr, ptmu, isData, evt, wgt );//reco
        utils->fillHistosV3( utils->histos2D[ "h_fitEnu_ptmu" ],EnuFit, ptmu, isData, evt, wgt );

        utils->fillHistosV3( utils->histos2D[ "h_dFitP_ptmu" ], pneutronRes, ptmu, isData, evt, wgt );//true
        utils->fillHistosV3( utils->histos2D[ "h_dFitEnu_ptmu" ],EnuResolution, ptmu, isData, evt, wgt );//true quantity
        utils->fillHistosV3( utils->histos2D[ "h_dEnu_ptmu" ],EnuQ2QEResolution, ptmu, isData, evt, wgt );//true quantity

        utils->fillHistosV3( utils->histos2D[ "h_dphitt_q2qe" ], dphitt, q2qe, isData, evt, wgt );
        utils->fillHistosV3( utils->histos2D[ "h_dphitt_ptmu" ], dphitt, ptmu, isData, evt, wgt );


      }


      FillHDL2D( "h_dthetaPdthetaR_q2qe", utils, "dthetaP", "q2", "dthetaR", reco_dthetaP, q2qe, reco_dthetaR, isData, evt, wgt, true);//no systematics
      if(hasNC) FillHDL2D( "h_dthetaPdthetaR_q2qe_nc", utils, reco_dthetaP, q2qe, reco_dthetaR, isData, evt, wgt, true);//no systematics
      if( abs(reco_dthetaR)<10 && abs(reco_dthetaP) < 10 )
      {
        FillHDL2D( "h_dthetaPdthetaR_q2qe_fine", utils,"dthetaP", "q2", "dthetaR",  reco_dthetaP, q2qe, reco_dthetaR, isData, evt, wgt, true);//no systematics
        if(hasNC) FillHDL2D( "h_dthetaPdthetaR_q2qe_nc_fine", utils,"dthetaP", "q2", "dthetaR",  reco_dthetaP, q2qe, reco_dthetaR, isData, evt, wgt, true);//no systematics
      }

      //if( abs(reco_dthetaR) < 10 ) FillHDL2D( "h_dthetaPq2qe_ptmu_center", utils, reco_dthetaP, ptmu, q2qe, isData, evt, wgt, RunCodeWithSystematics );
      //else FillHDL2D( "h_dthetaPq2qe_ptmu_side", utils, reco_dthetaP, ptmu, q2qe, isData, evt, wgt, RunCodeWithSystematics );
      //
      //
      //
      //
      //
      //
      //
      //
      //
          //utils->fillHistosV3( utils->histos2D["h_blobEnergy_q2qe"], blobEnergy, q2qe, isData,evt, wgt );
          //utils->fillHistosV3( utils->histos2D["h_nonBlobEnergy_q2qe"], nonBlobEnergy, q2qe, isData,evt, wgt );
          //utils->fillHistosV3( utils->histos2D["h_blobDist_q2qe"], blobDist, q2qe, isData,evt, wgt );
          //utils->fillHistosV3( utils->histos2D["h_nBlobs_q2qe"], nBlobs, q2qe, isData,evt, wgt );
          //utils->fillHistosV3( utils->histos2D["h_n2DBlobs_q2qe"], n2DBlobs, q2qe, isData,evt, wgt );
          //utils->fillHistosV3( utils->histos2D["h_n3DBlobs_q2qe"], n3DBlobs, q2qe, isData,evt, wgt );
          //utils->fillHistosV3( utils->histos2D["h_muonTheta_q2qe"], muon_theta, q2qe, isData,evt, wgt );
          //utils->fillHistosV3( utils->histos2D["h_muonPhi_q2qe"], muon_phi_detector, q2qe, isData,evt, wgt );
          //utils->fillHistosV3( utils->histos2D["h_muonE_q2qe"], muon_p, q2qe, isData,evt, wgt );

          //utils->fillHistosV3( utils->histos2D["h_protonDEDX_q2qe"], proton_score, q2qe, isData,evt, wgt );
          //utils->fillHistosV3( utils->histos2D["h_pionDEDX_q2qe"], pion_score, q2qe, isData,evt, wgt );
          //utils->fillHistosV3( utils->histos2D["h_hasTrack_q2qe"], mainCandidate->hasTrack, q2qe, isData,evt, wgt );

      if( abs(reco_dthetaR) < 10 ) 
      {
        utils->fillHistosV3( utils->histos2D["h_dthetaP_center_q2qe"], reco_dthetaP, q2qe, isData, evt, wgt);
        utils->fillHistosV3( utils->histos2D["h_dthetaPpos_center_q2qe"], abs(reco_dthetaP), q2qe, isData, evt, wgt);
        utils->fillHistosV3( utils->histos2D["h_dthetaP_center_ptmu"], reco_dthetaP, ptmu, isData, evt, wgt);

        int reco_dthetaP_pos = abs(reco_dthetaP);
        //bool angle_01 = (abs(reco_dthetaP) < 1 && abs(reco_dthetaR) < 1 );
        //bool angle_02 = (abs(reco_dthetaP) < 2 && abs(reco_dthetaR) < 2 );
        //if( angle_01 ) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_01"], q2qe, isData, evt, wgt );
        //else if( angle_02 ) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_02"], q2qe, isData, evt, wgt );
        if( angle_03 ) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_03"], q2qe, isData, evt, wgt );
        //else if( angle_04 ) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_04"], q2qe, isData, evt, wgt );
        else if( angle_05 ) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_05"], q2qe, isData, evt, wgt );
        //else if( angle_06 ) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_06"], q2qe, isData, evt, wgt );
        else if( angle_07 ) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_07"], q2qe, isData, evt, wgt );
        //else if( angle_08 ) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_08"], q2qe, isData, evt, wgt );
        //else if( angle_09 ) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_09"], q2qe, isData, evt, wgt );
        else if( angle_10 ) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_10"], q2qe, isData, evt, wgt );



        if(RunCodeWithSystematics && !isData)
        {
          utils->fillVertErrorBands( utils->histos2D["h_dthetaP_center_ptmu"],reco_dthetaP,ptmu,evt);
          utils->fillVertErrorBands( utils->histos2D["h_dthetaP_center_q2qe"],reco_dthetaP,q2qe,evt);
          utils->fillVertErrorBands( utils->histos2D["h_dthetaPpos_center_q2qe"],abs(reco_dthetaP),q2qe,evt);





          //if(angle_01) 
          //{
          //  utils->fillVertErrorBands( utils->histos1D["h_q2qe_angle_01"], q2qe,  evt );
          //  utils->fillLatErrorBands(  utils->histos1D["h_q2qe_angle_01"],  "q2", q2qe, evt, true );
          //}
          //else if(angle_02) 
          //{
          //  utils->fillVertErrorBands( utils->histos1D["h_q2qe_angle_02"], q2qe,  evt );
          //  utils->fillLatErrorBands(  utils->histos1D["h_q2qe_angle_02"],  "q2", q2qe, evt, true );
          //}
          if(angle_03) 
          {
            utils->fillVertErrorBands( utils->histos1D["h_q2qe_angle_03"], q2qe,  evt );
            utils->fillLatErrorBands(  utils->histos1D["h_q2qe_angle_03"],  "q2", q2qe, evt, true );
          }
          //else if(angle_04) 
          //{
          //  utils->fillVertErrorBands( utils->histos1D["h_q2qe_angle_04"], q2qe,  evt );
          //  utils->fillLatErrorBands(  utils->histos1D["h_q2qe_angle_04"],  "q2", q2qe, evt, true );
          //}
          else if(angle_05) 
          {
            utils->fillVertErrorBands( utils->histos1D["h_q2qe_angle_05"], q2qe,  evt );
            utils->fillLatErrorBands(  utils->histos1D["h_q2qe_angle_05"],  "q2", q2qe, evt, true );
          }
          //else if(angle_06) 
          //{
          //  utils->fillVertErrorBands( utils->histos1D["h_q2qe_angle_06"], q2qe,  evt );
          //  utils->fillLatErrorBands(  utils->histos1D["h_q2qe_angle_06"],  "q2", q2qe, evt, true );
          //}
          else if(angle_07) 
          {
            utils->fillVertErrorBands( utils->histos1D["h_q2qe_angle_07"], q2qe,  evt );
            utils->fillLatErrorBands(  utils->histos1D["h_q2qe_angle_07"],  "q2", q2qe, evt, true );
          }
          else if(angle_10) 
          {
            utils->fillVertErrorBands( utils->histos1D["h_q2qe_angle_10"], q2qe,  evt );
            utils->fillLatErrorBands(  utils->histos1D["h_q2qe_angle_10"],  "q2", q2qe, evt, true );
          }
        }

        if(angle_10)
        {
          //fill neutron diagnositic plots
          utils->fillHistosV3( utils->histos2D["h_blobEnergy_q2qe"], blobEnergy, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_nonBlobEnergy_q2qe"], nonBlobEnergy, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_blobDist_q2qe"], blobDist, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_nBlobs_q2qe"], nBlobs, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_n2DBlobs_q2qe"], n2DBlobs, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_n3DBlobs_q2qe"], n3DBlobs, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_muonTheta_q2qe"], muon_theta, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_muonPhi_q2qe"], muon_phi_detector, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_muonE_q2qe"], muon_p, q2qe, isData,evt, wgt );

          utils->fillHistosV3( utils->histos2D["h_protonDEDX_q2qe"], proton_score, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_pionDEDX_q2qe"], pion_score, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_hasTrack_q2qe"], mainCandidate->hasTrack, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_blobMaxE_q2qe"], mainCandidate->ClusMaxE, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_nClus_q2qe"], mainCandidate->NClusters, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_blobEnergyLow_q2qe"], blobEnergy, q2qe, isData,evt, wgt );

          utils->fillHistosV3( utils->histos2D["h_blobDist_blobEnergy"], blobDist, blobEnergy, isData,evt,wgt);
          utils->fillHistosV3( utils->histos2D["h_blobDist_blobEnergyLow"], blobDist, blobEnergy, isData,evt,wgt);
        }

        if( region == 1 )
        {
          utils->fillHistosV3( utils->histos2D["h_blobDist_q2qe_r1"], blobDist, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_blobEnergy_q2qe_r1"], blobEnergy, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_nonBlobEnergy_q2qe_r1"], nonBlobEnergy, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_nBlobs_q2qe_r1"], nBlobs, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_n2DBlobs_q2qe_r1"], n2DBlobs, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_n3DBlobs_q2qe_r1"], n3DBlobs, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_blobMaxE_q2qe_r1"], mainCandidate->ClusMaxE, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_nClus_q2qe_r1"], mainCandidate->NClusters, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_blobEnergyLow_q2qe_r1"], blobEnergy, q2qe, isData,evt, wgt );
        }
        if( region == 5 )
        {
          utils->fillHistosV3( utils->histos2D["h_blobDist_q2qe_r5"], blobDist, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_blobEnergy_q2qe_r5"], blobEnergy, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_nonBlobEnergy_q2qe_r5"], nonBlobEnergy, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_nBlobs_q2qe_r5"], nBlobs, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_n2DBlobs_q2qe_r5"], n2DBlobs, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_n3DBlobs_q2qe_r5"], n3DBlobs, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_blobMaxE_q2qe_r5"], mainCandidate->ClusMaxE, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_nClus_q2qe_r5"], mainCandidate->NClusters, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_blobEnergyLow_q2qe_r5"], blobEnergy, q2qe, isData,evt, wgt );
        }
        if( region == -1 )
        {
          utils->fillHistosV3( utils->histos2D["h_blobDist_q2qe_r99"], blobDist, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_blobEnergy_q2qe_r99"], blobEnergy, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_nonBlobEnergy_q2qe_r99"], nonBlobEnergy, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_nBlobs_q2qe_r99"], nBlobs, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_n2DBlobs_q2qe_r99"], n2DBlobs, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_n3DBlobs_q2qe_r99"], n3DBlobs, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_blobMaxE_q2qe_r99"], mainCandidate->ClusMaxE, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_nClus_q2qe_r99"], mainCandidate->NClusters, q2qe, isData,evt, wgt );
          utils->fillHistosV3( utils->histos2D["h_blobEnergyLow_q2qe_r99"], blobEnergy, q2qe, isData,evt, wgt );
        }




      }



      //------------------------------------------------------
      //Fill 2D Vertical Error Bands
      //Error Bands are not filled up for tmu_costheta 
      //and pmu_ptmu histos to avoid code slow down 
      //------------------------------------------------------
      if(RunCodeWithSystematics && !isData){

        utils->fillVertErrorBands( utils->histos1D["h_ptmu"],  ptmu, evt );
        utils->fillVertErrorBands( utils->histos1D["h_q2qe"],  q2qe, evt );

        //if( inVtx )
        //{
        //  if( region == 0) utils->fillVertErrorBands( utils->histos1D["h_q2qe_vtx_region_00"], q2qe, evt );
        //  if( region == 1) utils->fillVertErrorBands( utils->histos1D["h_q2qe_vtx_region_01"], q2qe, evt );
        //  if( region == 2) utils->fillVertErrorBands( utils->histos1D["h_q2qe_vtx_region_02"], q2qe, evt );
        //  if( region == 3) utils->fillVertErrorBands( utils->histos1D["h_q2qe_vtx_region_03"], q2qe, evt );
        //  if( region == 4) utils->fillVertErrorBands( utils->histos1D["h_q2qe_vtx_region_04"], q2qe, evt );
        //  if( region == 5) utils->fillVertErrorBands( utils->histos1D["h_q2qe_vtx_region_05"], q2qe, evt );
        //  if( region ==-1) utils->fillVertErrorBands( utils->histos1D["h_q2qe_vtx_region_99"], q2qe, evt );
        //}
        //else
        //{
        //  if( region == 0) utils->fillVertErrorBands( utils->histos1D["h_q2qe_nonvtx_region_00"], q2qe, evt );
        //  if( region == 1) utils->fillVertErrorBands( utils->histos1D["h_q2qe_nonvtx_region_01"], q2qe, evt );
        //  if( region == 2) utils->fillVertErrorBands( utils->histos1D["h_q2qe_nonvtx_region_02"], q2qe, evt );
        //  if( region == 3) utils->fillVertErrorBands( utils->histos1D["h_q2qe_nonvtx_region_03"], q2qe, evt );
        //  if( region == 4) utils->fillVertErrorBands( utils->histos1D["h_q2qe_nonvtx_region_04"], q2qe, evt );
        //  if( region == 5) utils->fillVertErrorBands( utils->histos1D["h_q2qe_nonvtx_region_05"], q2qe, evt );
        //  if( region ==-1) utils->fillVertErrorBands( utils->histos1D["h_q2qe_nonvtx_region_99"], q2qe, evt );
        //}
   


        utils->fillVertErrorBands( utils->histos2D[Form( "h_q2qe_ptmu_region_%02d",regionCode) ], q2qe, ptmu, evt );
        if( region != -1 )utils->fillVertErrorBands( utils->histos2D["h_q2qe_ptmu_region_non99"], q2qe, ptmu, evt );

        utils->fillVertErrorBands( utils->histos2D[Form( "h_enu_q2qe_region_%02d",regionCode) ], enu, q2qe, evt );
        if( region != -1 )utils->fillVertErrorBands( utils->histos2D["h_enu_q2qe_region_non99"], enu, q2qe, evt );


        utils->fillVertErrorBands( utils->histos2D[Form( "h_muontheta_q2qe_region_%02d",regionCode) ], muon_theta, q2qe, evt );
        if( region != -1 )utils->fillVertErrorBands( utils->histos2D["h_muontheta_q2qe_region_non99"], muon_theta, q2qe, evt );


        if(isDebug) cout<<"Filled q2qe_ptmu_region Vert"<<endl;



        utils->fillVertErrorBands( utils->histos1D["h_enu"],  enu, evt );
        utils->fillVertErrorBands( utils->histos2D["h_q2qe_ptmu"],  q2qe, ptmu, evt );
        
        utils->fillVertErrorBands( utils->histos2D["h_dalphat_ptmu"],reco_dalphat,ptmu,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dphit_ptmu"],reco_dphit,ptmu,evt);
        utils->fillVertErrorBands( utils->histos2D["h_pn_ptmu"],reco_pn,ptmu,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dpt_ptmu"],reco_dpt,ptmu,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dptx_ptmu"],reco_dptx,ptmu,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dpty_ptmu"],reco_dpty,ptmu,evt);
        utils->fillVertErrorBands( utils->histos2D["h_tp_ptmu"],protonTn, ptmu,evt);
        utils->fillVertErrorBands( utils->histos2D["h_ptheta_ptmu"],protonAngle, ptmu,evt);

        utils->fillVertErrorBands( utils->histos2D["h_dalphat_q2qe"],reco_dalphat,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dphit_q2qe"],reco_dphit,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_pn_q2qe"],reco_pn,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dpt_q2qe"],reco_dpt,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dptx_q2qe"],reco_dptx,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dpty_q2qe"],reco_dpty,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_tp_q2qe"],protonTn, q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_ptheta_q2qe"],protonAngle, q2qe,evt);


        //utils->fillVertErrorBands( utils->histos2D["h_sign_ptmu"],reco_sign,ptmu,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_signdalphat_ptmu"],reco_signed_dalphat,ptmu,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_signdphit_ptmu"],reco_signed_dphit,ptmu,evt);

        utils->fillVertErrorBands( utils->histos2D["h_dthetaP_ptmu"],reco_dthetaP,ptmu,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dthetaR_ptmu"],reco_dthetaR,ptmu,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPP_ptmu"],reco_dPp,ptmu,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPR_ptmu"],reco_dPr,ptmu,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPPi_ptmu"],reco_dPp_infer,ptmu,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPRi_ptmu"],reco_dPr_infer,ptmu,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dtheta2D_ptmu"],dtheta2D,ptmu,evt);

        utils->fillVertErrorBands( utils->histos2D["h_dthetaP_q2qe"],reco_dthetaP,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dthetaR_q2qe"],reco_dthetaR,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPP_q2qe"],reco_dPp,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPR_q2qe"],reco_dPr,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPPi_q2qe"],reco_dPp_infer,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPRi_q2qe"],reco_dPr_infer,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dtheta2D_q2qe"],dtheta2D,q2qe,evt);

        utils->fillVertErrorBands( utils->histos2D["h_muontheta_q2qe"], muon_theta, q2qe, evt );
        utils->fillVertErrorBands( utils->histos2D["h_muontheta_ptmu"], muon_theta, ptmu, evt );
        utils->fillVertErrorBands( utils->histos2D["h_muonmomentum_q2qe"], reco_muon_e, q2qe, evt );
        utils->fillVertErrorBands( utils->histos2D["h_muonmomentum_ptmu"], reco_muon_e, ptmu, evt );
        if(angle_10) utils->fillVertErrorBands( utils->histos2D["h_enu_q2qe"], enu, q2qe, evt );
        if(angle_10) utils->fillVertErrorBands( utils->histos2D["h_enu_ptmu"], enu, ptmu, evt );

        utils->fillVertErrorBands( utils->histos2D["h_dthetaP_dthetaR"],reco_dthetaP,reco_dthetaR,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dthetaP_dthetaR_Fine"],reco_dthetaP,reco_dthetaR,evt);

        utils->fillVertErrorBands( utils->histos2D["h_recoil_inc_q2qe"], recoil_inc, q2qe, evt );
        utils->fillVertErrorBands( utils->histos2D["h_recoil_inc_q3"], recoil_inc, q3, evt );

        if(hasNC && hasProton && false)
        {
          utils->fillVertErrorBands( utils->histos2D[ Form("h_pn_q2qe_%s", quadrant.c_str())], reco_pn, q2qe, evt);
          utils->fillVertErrorBands( utils->histos2D[ Form("h_dphit_q2qe_%s", quadrant.c_str())], reco_dphit, q2qe,evt);
          utils->fillVertErrorBands( utils->histos2D[ Form("h_dalphat_q2qe_%s", quadrant.c_str())], reco_dalphat, q2qe,evt);
          utils->fillVertErrorBands( utils->histos2D[ Form("h_dpt_q2qe_%s", quadrant.c_str())], reco_dpt, q2qe,evt);
          utils->fillVertErrorBands( utils->histos2D[ Form("h_dptx_q2qe_%s", quadrant.c_str())], reco_dptx, q2qe,evt);
          utils->fillVertErrorBands( utils->histos2D[ Form("h_dpty_q2qe_%s", quadrant.c_str())], reco_dpty, q2qe,evt);

          utils->fillLatErrorBands( utils->histos2D[ Form("h_pn_q2qe_%s", quadrant.c_str())], "pn","q2",reco_pn, q2qe, evt);
          utils->fillLatErrorBands( utils->histos2D[ Form("h_dphit_q2qe_%s", quadrant.c_str())],"dphiT","q2", reco_dphit, q2qe,evt);
          utils->fillLatErrorBands( utils->histos2D[ Form("h_dalphat_q2qe_%s", quadrant.c_str())],"dalphaT","q2", reco_dalphat, q2qe,evt);
          utils->fillLatErrorBands( utils->histos2D[ Form("h_dpt_q2qe_%s", quadrant.c_str())],"dpT ","q2", reco_dpt, q2qe,evt);
          utils->fillLatErrorBands( utils->histos2D[ Form("h_dptx_q2qe_%s", quadrant.c_str())],"dpTx","q2", reco_dptx, q2qe,evt);
          utils->fillLatErrorBands( utils->histos2D[ Form("h_dpty_q2qe_%s", quadrant.c_str())],"dpTy","q2", reco_dpty, q2qe,evt);


          utils->fillVertErrorBands( utils->histos2D[ Form("h_pn_ptmu_%s", quadrant.c_str())], reco_pn, ptmu, evt);
          utils->fillVertErrorBands( utils->histos2D[ Form("h_dphit_ptmu_%s", quadrant.c_str())], reco_dphit, ptmu,evt);
          utils->fillVertErrorBands( utils->histos2D[ Form("h_dalphat_ptmu_%s", quadrant.c_str())], reco_dalphat, ptmu,evt);
          utils->fillVertErrorBands( utils->histos2D[ Form("h_dpt_ptmu_%s", quadrant.c_str())], reco_dpt, ptmu,evt);
          utils->fillVertErrorBands( utils->histos2D[ Form("h_dptx_ptmu_%s", quadrant.c_str())], reco_dptx, ptmu,evt);
          utils->fillVertErrorBands( utils->histos2D[ Form("h_dpty_ptmu_%s", quadrant.c_str())], reco_dpty, ptmu,evt);

          utils->fillLatErrorBands( utils->histos2D[ Form("h_pn_ptmu_%s", quadrant.c_str())], "pn","pT",reco_pn, ptmu, evt);
          utils->fillLatErrorBands( utils->histos2D[ Form("h_dphit_ptmu_%s", quadrant.c_str())],"dphiT","pT", reco_dphit, ptmu,evt);
          utils->fillLatErrorBands( utils->histos2D[ Form("h_dalphat_ptmu_%s", quadrant.c_str())],"dalphaT","pT", reco_dalphat, ptmu,evt);
          utils->fillLatErrorBands( utils->histos2D[ Form("h_dpt_ptmu_%s", quadrant.c_str())],"dpT ","pT", reco_dpt, ptmu,evt);
          utils->fillLatErrorBands( utils->histos2D[ Form("h_dptx_ptmu_%s", quadrant.c_str())],"dpTx","pT", reco_dptx, ptmu,evt);
          utils->fillLatErrorBands( utils->histos2D[ Form("h_dpty_ptmu_%s", quadrant.c_str())],"dpTy","pT", reco_dpty, ptmu,evt);

        }



       

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
        utils->fillLatErrorBands( utils->histos1D["h_ptmu"],  "pT",    ptmu, evt, true );
        utils->fillLatErrorBands( utils->histos1D["h_q2qe"],  "q2",    q2qe, evt, true );
        utils->fillLatErrorBands( utils->histos2D["h_q2qe_ptmu"],  "q2", "pT",  q2qe, ptmu, evt );
        utils->fillLatErrorBands( utils->histos1D["h_dthetaP"],  "dthetaP", reco_dthetaP, evt, true );
        utils->fillLatErrorBands( utils->histos1D["h_dthetaR"],  "dthetaR", reco_dthetaR, evt, true );
        utils->fillLatErrorBands( utils->histos1D["h_pn"],  "pn", reco_pn, evt, true );
        //This is a test
        utils->fillLatErrorBands( utils->histos2D["h_dthetaP_q2qe"], "dthetaP", "q2", reco_dthetaP, q2qe, evt );
        utils->fillLatErrorBands( utils->histos2D["h_dthetaR_q2qe"], "dthetaR", "q2", reco_dthetaR, q2qe, evt );
        utils->fillLatErrorBands( utils->histos2D["h_pn_q2qe"], "pn", "q2", reco_pn, q2qe, evt );
        utils->fillLatErrorBands( utils->histos2D["h_dptx_q2qe"], "dpTx", "q2", reco_dptx, q2qe, evt );
        utils->fillLatErrorBands( utils->histos2D["h_dpty_q2qe"], "dpTy", "q2", reco_dpty, q2qe, evt );
        utils->fillLatErrorBands( utils->histos2D["h_dpt_q2qe"], "dpT", "q2", reco_dpt, q2qe, evt );
        utils->fillLatErrorBands( utils->histos2D["h_dalphat_q2qe"], "dalphaT", "q2", reco_dalphat, q2qe, evt );
        utils->fillLatErrorBands( utils->histos2D["h_dphit_q2qe"], "dphiT", "q2", reco_dphit, q2qe, evt );
        utils->fillLatErrorBands( utils->histos2D["h_muontheta_q2qe"], "theta", "q2", muon_theta, q2qe, evt );
        utils->fillLatErrorBands( utils->histos2D["h_muonmomentum_q2qe"], "emu", "q2", reco_muon_e, q2qe, evt );
        if(angle_10) utils->fillLatErrorBands( utils->histos2D["h_enu_q2qe"], "enu", "q2", enu, q2qe, evt );

        utils->fillLatErrorBands( utils->histos2D["h_dthetaP_ptmu"], "dthetaP", "pT", reco_dthetaP, ptmu, evt );
        utils->fillLatErrorBands( utils->histos2D["h_dthetaR_ptmu"], "dthetaR", "pT", reco_dthetaR, ptmu, evt );
        utils->fillLatErrorBands( utils->histos2D["h_pn_ptmu"], "pn", "pT", reco_pn, ptmu, evt );
        utils->fillLatErrorBands( utils->histos2D["h_dptx_ptmu"], "dpTx", "pT", reco_dptx, ptmu, evt );
        utils->fillLatErrorBands( utils->histos2D["h_dpty_ptmu"], "dpTy", "pT", reco_dpty, ptmu, evt );
        utils->fillLatErrorBands( utils->histos2D["h_dpt_ptmu"], "dpT", "pT", reco_dpt, ptmu, evt );
        utils->fillLatErrorBands( utils->histos2D["h_dalphat_ptmu"], "dalphaT", "pT", reco_dalphat, ptmu, evt );
        utils->fillLatErrorBands( utils->histos2D["h_dphit_ptmu"], "dphiT", "pT", reco_dphit, ptmu, evt );
        utils->fillLatErrorBands( utils->histos2D["h_muontheta_ptmu"], "theta", "pT", muon_theta, ptmu, evt );
        utils->fillLatErrorBands( utils->histos2D["h_muonmomentum_ptmu"], "emu", "pT", reco_muon_e, ptmu, evt );
        if(angle_10) utils->fillLatErrorBands( utils->histos2D["h_enu_ptmu"], "enu", "pT", enu, ptmu, evt );

        utils->fillLatErrorBands( utils->histos2D[Form( "h_q2qe_ptmu_region_%02d",regionCode) ], "q2qe","pT", q2qe, ptmu, evt );
        if( region != -1 )utils->fillLatErrorBands( utils->histos2D["h_q2qe_ptmu_region_non99"], "q2qe","pT", q2qe, ptmu, evt );

        utils->fillLatErrorBands( utils->histos2D[Form( "h_enu_q2qe_region_%02d",regionCode) ], "enu","q2qe", enu, q2qe, evt );
        if( region != -1 )utils->fillLatErrorBands( utils->histos2D["h_enu_q2qe_region_non99"], "enu","q2qe", enu, q2qe, evt );

        utils->fillLatErrorBands( utils->histos2D[Form( "h_muontheta_q2qe_region_%02d",regionCode) ], "theta","q2qe", muon_theta, q2qe, evt );
        if( region != -1 )utils->fillLatErrorBands( utils->histos2D["h_muontheta_q2qe_region_non99"], "theta","q2qe", muon_theta, q2qe, evt );


        if(isDebug) cout<<"Filled q2qe_ptmu_region Lat"<<endl;

        //if( inVtx )
        //{
        //  if( region == 0) utils->fillLatErrorBands( utils->histos1D["h_q2qe_vtx_region_00"], "q2", q2qe, evt, true );
        //  if( region == 1) utils->fillLatErrorBands( utils->histos1D["h_q2qe_vtx_region_01"], "q2", q2qe, evt, true );
        //  if( region == 2) utils->fillLatErrorBands( utils->histos1D["h_q2qe_vtx_region_02"], "q2", q2qe, evt, true );
        //  if( region == 3) utils->fillLatErrorBands( utils->histos1D["h_q2qe_vtx_region_03"], "q2", q2qe, evt, true );
        //  if( region == 4) utils->fillLatErrorBands( utils->histos1D["h_q2qe_vtx_region_04"], "q2", q2qe, evt, true );
        //  if( region == 5) utils->fillLatErrorBands( utils->histos1D["h_q2qe_vtx_region_05"], "q2", q2qe, evt, true );
        //  if( region ==-1) utils->fillLatErrorBands( utils->histos1D["h_q2qe_vtx_region_99"], "q2", q2qe, evt, true );
        //}
        //else
        //{
        //  if( region == 0) utils->fillLatErrorBands( utils->histos1D["h_q2qe_nonvtx_region_00"], "q2", q2qe, evt, true );
        //  if( region == 1) utils->fillLatErrorBands( utils->histos1D["h_q2qe_nonvtx_region_01"], "q2", q2qe, evt, true );
        //  if( region == 2) utils->fillLatErrorBands( utils->histos1D["h_q2qe_nonvtx_region_02"], "q2", q2qe, evt, true );
        //  if( region == 3) utils->fillLatErrorBands( utils->histos1D["h_q2qe_nonvtx_region_03"], "q2", q2qe, evt, true );
        //  if( region == 4) utils->fillLatErrorBands( utils->histos1D["h_q2qe_nonvtx_region_04"], "q2", q2qe, evt, true );
        //  if( region == 5) utils->fillLatErrorBands( utils->histos1D["h_q2qe_nonvtx_region_05"], "q2", q2qe, evt, true );
        //  if( region ==-1) utils->fillLatErrorBands( utils->histos1D["h_q2qe_nonvtx_region_99"], "q2", q2qe, evt, true );
        //}





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

  GeoUtils *geoUtils = new GeoUtils();

  //---------------------------------------------
  // Book Histograms
  //---------------------------------------------
  //

  //1D
  MnvH1D *h_multiplicity[nHistos],*h_ptmu[nHistos];
  MnvH1D *h_q2qe[nHistos], *h_enu[nHistos];
  MnvH1D *h_pmu[nHistos];
  MnvH1D *h_dthetaP[nHistos];
  MnvH1D *h_dthetaR[nHistos];
  MnvH1D *h_pn[nHistos];

  utils->bookHistos( h_multiplicity, "h_multiplicity", "Multiplicity of Outgoing Tracks; Number of Outgoing Tracks; Number of Events", 9, 0., 9.);
  utils->bookHistos( h_ptmu, "h_ptmu", Form( "Muon Pt;Reconstructed Muon p_{T} (GeV);Events / %.3f GeV", muonPtbins.default_width ), muonPtbins );
  utils->bookHistos( h_q2qe, "h_q2qe", Form( "Q^{2}_{QE};Reconstructed Q^{2}_{QE} (GeV^{2});Events / %.3f GeV^{2}", Q2bins.default_width ), Q2bins );
  utils->bookHistos( h_enu, "h_enu", Form( "E (GeV);Reconstructed E_{#nu} (GeV);Events / %.3f GeV^{2}", Q2bins.default_width ), Enubins );
  utils->bookHistos( h_pmu, "h_pmu", "Reconstructed p_{#mu} (GeV)", UniformMuonPBins );
  utils->bookHistos( h_dthetaP, "h_dthetaP", "dThetaP",dthetaPerpbins );
  utils->bookHistos( h_dthetaR, "h_dthetaR", "dThetaR",dthetaReactbins );
  utils->bookHistos( h_pn, "h_pn", "pn", pnbins );


  cout<<"declaring h_q2qe_angles"<<endl;
  vector< vector<MnvH1D*> > h_q2qe_angles( 10, vector<MnvH1D*>(nHistos, NULL) );
  for( int i = 0; i< 10; i++ ) // h_q2qe_angle_01..10
    utils->bookHistos( &h_q2qe_angles[i][0], Form("h_q2qe_angle_%02d", i+1 ),Form( "Q^{2}_{QE};Reconstructed Q^{2}_{QE} (GeV^{2});Events / %.3f GeV^{2}", Q2bins.default_width ), Q2bins );

  MnvH2D *h_q2qe_ptmu[nHistos];
  utils->bookHistos( h_q2qe_ptmu, "h_q2qe_ptmu","q2qe:ptmu",Q2bins,muonPtbins);

  vector<string> region_names({"00","01","02","03","04","05","99","non99"});
  vector< vector<MnvH2D*> > h_q2qe_ptmu_regions(region_names.size(), vector<MnvH2D*>(nHistos,NULL) );
  vector< vector<MnvH2D*> > h_enu_q2qe_regions(region_names.size(), vector<MnvH2D*>(nHistos,NULL) );
  vector< vector<MnvH2D*> > h_muontheta_q2qe_regions(region_names.size(), vector<MnvH2D*>(nHistos,NULL) );

  for( unsigned int i = 0; i < region_names.size(); i++ )
  {
    utils->bookHistos( &h_q2qe_ptmu_regions[i][0], Form("h_q2qe_ptmu_region_%s", region_names[i].c_str()),"q2qe:ptmu",Q2bins,muonPtbins);
    utils->bookHistos( &h_enu_q2qe_regions[i][0], Form("h_enu_q2qe_region_%s", region_names[i].c_str()),"enu:q2qe",leptonmomentumbins, Q2bins);
    utils->bookHistos( &h_muontheta_q2qe_regions[i][0], Form("h_muontheta_q2qe_region_%s", region_names[i].c_str()),"muontheta:q2qe",MuonThetaBins,Q2bins );
  }



  MnvH2D *h_recoil_inc_q2qe[nHistos];
  MnvH2D *h_recoil_inc_q3[nHistos];
  utils->bookHistos( h_recoil_inc_q2qe, "h_recoil_inc_q2qe","recoil_inc:q2qe",UniformVisibleEBins, Q2bins);
  utils->bookHistos( h_recoil_inc_q3, "h_recoil_inc_q3","recoil_inc:q3",UniformVisibleEBins, UniformQ0Q3Bins);


  cout<<"declaring h_stki_xvars"<<endl;
  //STKI vs XVAR
  vector< string > stki_names({"dalphat","dphit", "pn","dpt","dptx","dpty","tp","ptheta"});
  vector< axis_binning > stki_bins({ dalphatbins, dphitbins, pnbins, dptbins, dptxbins, dptybins, protonKineticbins, protonThetabins});
  vector< vector< MnvH2D* > > h_stki_ptmu( stki_names.size(), vector< MnvH2D* >(nHistos, NULL) );
  vector< vector< MnvH2D* > > h_stki_q2qe( stki_names.size(), vector< MnvH2D* >(nHistos, NULL) );
  for( unsigned int i = 0; i < stki_names.size(); i++ )
  {
    string varname = stki_names[i];
    utils->bookHistos( &h_stki_ptmu[i][0], Form("h_%s_ptmu", varname.c_str()), Form("%s:ptmu", varname.c_str() ), stki_bins[i], muonPtbins );
    utils->bookHistos( &h_stki_q2qe[i][0], Form("h_%s_q2qe", varname.c_str()), Form("%s:q2qe", varname.c_str() ), stki_bins[i], Q2bins );
  }

  cout<<"declaring h_neut_xvars"<<endl;
  //NEUT vs XVAR
  vector< string > neut_names({"dthetaR","dthetaP","dPP","dPR","dPPi","dPRi"});
  vector< axis_binning > neut_bins({ dthetaReactbins, dthetaPerpbins,dptxbins,dptybins, dptxbins, dptybins });
  vector< vector< MnvH2D* > > h_neut_ptmu( neut_names.size(), vector< MnvH2D* >(nHistos, NULL) );
  vector< vector< MnvH2D* > > h_neut_q2qe( neut_names.size(), vector< MnvH2D* >(nHistos, NULL) );
  for( unsigned int i = 0; i < neut_names.size(); i++ )
  {
    string varname = neut_names[i];
    utils->bookHistos( &h_neut_ptmu[i][0], Form("h_%s_ptmu", varname.c_str()), Form("%s:ptmu", varname.c_str() ), neut_bins[i], muonPtbins );
    utils->bookHistos( &h_neut_q2qe[i][0], Form("h_%s_q2qe", varname.c_str()), Form("%s:q2qe", varname.c_str() ), neut_bins[i], Q2bins );
  }

  MnvH2D* h_dthetaP_center_q2qe[nHistos], *h_dthetaP_center_ptmu[nHistos];
  utils->bookHistos( h_dthetaP_center_ptmu, "h_dthetaP_center_ptmu", "dthetaPcenter:ptmu", dthetaPerpbins, muonPtbins );
  utils->bookHistos( h_dthetaP_center_q2qe, "h_dthetaP_center_q2qe", "dthetaPcenter:q2qe", dthetaPerpbins, Q2bins );


  //Misc
  cout<<"declaring MISC"<<endl;
  MnvH2D *h_dtheta2D_q2qe[nHistos],*h_dtheta2D_ptmu[nHistos];
  MnvH2D *h_dthetaP_dthetaR[nHistos], *h_dthetaP_dthetaR_Fine[nHistos];
  MnvH2D *h_dthetaPpos_center_q2qe[nHistos];

  utils->bookHistos( h_dtheta2D_q2qe, "h_dtheta2D_q2qe", "dtheta2D:q2qe", UniformDThetaBins, Q2bins );
  utils->bookHistos( h_dtheta2D_ptmu, "h_dtheta2D_ptmu", "dtheta2D:ptmu", UniformDThetaBins, muonPtbins );
  utils->bookHistos( h_dthetaP_dthetaR, "h_dthetaP_dthetaR", "dthetaP: dthetaR", dthetaPerpbins, dthetaReactbins );
  utils->bookHistos( h_dthetaP_dthetaR_Fine, "h_dthetaP_dthetaR_Fine", "dthetaP: dthetaR", dthetaFineBins, dthetaFineBins);
  utils->bookHistos( h_dthetaPpos_center_q2qe, "h_dthetaPpos_center_q2qe", "dthetaP: dthetaR", dthetaPerpbinsPositive, Q2bins );


  MnvH2D *h_node_nodeenergy[nHistos],*h_node_nodeenergy_long[nHistos],*h_node_nodeenergy_anchored[nHistos],*h_node_nodeenergy_vest[nHistos];
  MnvH2D *h_protonRes_nodeenergy_node01[nHistos],*h_protonRes_nodeenergy_node2[nHistos],*h_protonRes_nodeenergy_node3[nHistos],*h_protonRes_nodeenergy_node4[nHistos],*h_protonRes_nodeenergy_node5[nHistos];
  utils->bookHistos( h_node_nodeenergy, "h_node_nodeenergy","node:energy",50,0,50,500,0,50);
  utils->bookHistos( h_node_nodeenergy_long, "h_node_nodeenergy_long","node:energy",50,0,50,500,0,50);
  utils->bookHistos( h_node_nodeenergy_anchored, "h_node_nodeenergy_anchored","node:energy",50,0,50,500,0,50);
  utils->bookHistos( h_node_nodeenergy_vest, "h_node_nodeenergy_vest","node:energy",50,0,50,500,0,50);

  utils->bookHistos( h_protonRes_nodeenergy_node01, "h_protonRes_nodeenergy_node01","res:node",500,0,500,100,-1,1);
  utils->bookHistos( h_protonRes_nodeenergy_node2, "h_protonRes_nodeenergy_node2","res:node",500,0,500,100,-1,1);
  utils->bookHistos( h_protonRes_nodeenergy_node3, "h_protonRes_nodeenergy_node3","res:node",500,0,500,100,-1,1);
  utils->bookHistos( h_protonRes_nodeenergy_node4, "h_protonRes_nodeenergy_node4","res:node",500,0,500,100,-1,1);
  utils->bookHistos( h_protonRes_nodeenergy_node5, "h_protonRes_nodeenergy_node5","res:node",500,0,500,100,-1,1);

  MnvH2D *h_res_dthetaP_ptmu[nHistos],*h_res_dthetaP_q2qe[nHistos],*h_res_dthetaR_ptmu[nHistos],*h_res_dthetaR_q2qe[nHistos];
  MnvH2D *h_res_dtheta_ptmu[nHistos],*h_res_dtheta_q2qe[nHistos];
  MnvH2D *h_res_dthetaX_ptmu[nHistos],*h_res_dthetaX_q2qe[nHistos];
  MnvH2D *h_res_dthetaY_ptmu[nHistos],*h_res_dthetaY_q2qe[nHistos];
  MnvH2D *h_res_dthetaP[nHistos],*h_res_dthetaR[nHistos];
  utils->bookHistos( h_res_dthetaP_ptmu, "h_res_dthetaP_ptmu", "res:dthetaP_ptmu", 100,-50,50, 10,0,2.5);
  utils->bookHistos( h_res_dthetaR_ptmu, "h_res_dthetaR_ptmu", "res:dthetaR_ptmu", 100,-50,50, 10,0,2.5);
  utils->bookHistos( h_res_dtheta_ptmu, "h_res_dtheta_ptmu", "res:dtheta_ptmu", 100,-50,50, 10,0,2.5);
  utils->bookHistos( h_res_dthetaX_ptmu, "h_res_dthetaX_ptmu", "res:dthetaX_ptmu", 100,-50,50, 10,0,2.5);
  utils->bookHistos( h_res_dthetaY_ptmu, "h_res_dthetaY_ptmu", "res:dthetaY_ptmu", 100,-50,50, 10,0,2.5);

  utils->bookHistos( h_res_dthetaP_q2qe, "h_res_dthetaP_q2qe", "res:dthetaP_q2qe", 100,-50,50, 10,0,2.5);
  utils->bookHistos( h_res_dthetaR_q2qe, "h_res_dthetaR_q2qe", "res:dthetaR_q2qe", 100,-50,50, 10,0,2.5);
  utils->bookHistos( h_res_dtheta_q2qe, "h_res_dtheta_q2qe", "res:dtheta_q2qe", 100,-50,50, 10,0,2.5);
  utils->bookHistos( h_res_dthetaX_q2qe, "h_res_dthetaX_q2qe", "res:dthetaX_q2qe", 100,-50,50, 10,0,2.5);
  utils->bookHistos( h_res_dthetaY_q2qe, "h_res_dthetaY_q2qe", "res:dthetaY_q2qe", 100,-50,50, 10,0,2.5);

  utils->bookHistos( h_res_dthetaP, "h_res_dthetaP", "res:dthetaP", 100,-50,50, 72,-180,180);
  utils->bookHistos( h_res_dthetaR, "h_res_dthetaR", "res:dthetaR", 100,-50,50, 72,-180,180);


  //Some migration style plots

  MnvH2D *h_true_vs_reco_dthetaP[nHistos], *h_true_vs_reco_dthetaR[nHistos];
  utils->bookHistos( h_true_vs_reco_dthetaP, "h_true_vs_reco_dthetaP", "true:reco_dthetaP", 72,-180,180, 72,-180,180);
  utils->bookHistos( h_true_vs_reco_dthetaR, "h_true_vs_reco_dthetaR", "true:reco_dthetaR", 72,-180,180, 72,-180,180);



  cout<<"declaring HDL vs XVARS"<<endl;
  //dthetaRdthetaP vs XVAR
  axis_binning   dthetaPdthetaRBins;
  axis_binning   dthetaPdthetaRFineBins;

  MnvH2D   *h_dthetaPdthetaR_q2qe[nHistos];
  MnvH2D   *h_dthetaPdthetaR_q2qe_fine[nHistos];
  MnvH2D   *h_dthetaPdthetaR_q2qe_nc[nHistos];
  MnvH2D   *h_dthetaPdthetaR_q2qe_nc_fine[nHistos];

  cout<<"Declaring HDL"<<endl;
  HyperDimLinearizer *hdl_dthetaPdthetaR_q2qe = DeclareHDLHistos2D( utils, h_dthetaPdthetaR_q2qe, "h_dthetaPdthetaR_q2qe", "q2qe vs dthetaPdthetaR: dthetaP : dthetaR", dthetaPerpbins, Q2bins, dthetaReactbins, dthetaPdthetaRBins, RunCodeWithSystematics );
  HyperDimLinearizer *hdl_dthetaPdthetaR_q2qe_nc = DeclareHDLHistos2D( utils, h_dthetaPdthetaR_q2qe_nc, "h_dthetaPdthetaR_q2qe_nc", "q2qe vs dthetaPdthetaR with nc: dthetaP : dthetaR", dthetaPerpbins, Q2bins, dthetaReactbins, dthetaPdthetaRBins , RunCodeWithSystematics);
  HyperDimLinearizer *hdl_dthetaPdthetaR_q2qe_fine = DeclareHDLHistos2D( utils, h_dthetaPdthetaR_q2qe_fine, "h_dthetaPdthetaR_q2qe_fine", "q2qe vs dthetaPdthetaR: dthetaP : dthetaR", dthetaFineBins, Q2bins, dthetaFineBins, dthetaPdthetaRFineBins, RunCodeWithSystematics);
  HyperDimLinearizer *hdl_dthetaPdthetaR_q2qe_nc_fine = DeclareHDLHistos2D( utils, h_dthetaPdthetaR_q2qe_nc_fine, "h_dthetaPdthetaR_q2qe_nc_fine", "q2qe vs dthetaPdthetaR with nc: dthetaP : dthetaR", dthetaFineBins, Q2bins, dthetaFineBins, dthetaPdthetaRFineBins, RunCodeWithSystematics );
  cout<<"Declared HDL"<<endl;


  //Quadrant Based Histograms
  vector< string > stki_reduced_names({"dalphat","dphit", "pn","dpt","dptx","dpty"});
  vector< axis_binning > stki_reduced_bins({ dalphatbins, dphitbins, pnbins, dptbins, dptxbins, dptybins});

  vector<string> quad_names({"ll","lr","rl","rr"});
  vector< vector< vector<MnvH2D*> > > 
    quad_h_stki_q2qe( quad_names.size(), vector<vector<MnvH2D*>>( stki_reduced_names.size(), vector<MnvH2D*>(nHistos,NULL) ) ),
    quad_h_stki_ptmu( quad_names.size(), vector<vector<MnvH2D*>>( stki_reduced_names.size(), vector<MnvH2D*>(nHistos,NULL) ) );

  for( unsigned int i = 0; i< quad_names.size(); i++ )
  {
    //h_stki_xvar_ll
    string qn = quad_names[i];
    for( unsigned int j = 0; j < stki_reduced_names.size(); j++ )
    {
      string var = stki_reduced_names[j];

      utils->bookHistos( &quad_h_stki_q2qe[i][j][0],
          Form("h_%s_q2qe_%s", var.c_str(), qn.c_str() ),
          Form("%s_%s:q2qe", var.c_str(), qn.c_str() ),
          stki_reduced_bins[i],Q2bins );
      utils->bookHistos( &quad_h_stki_ptmu[i][j][0],
          Form("h_%s_ptmu_%s", var.c_str(), qn.c_str() ),
          Form("%s_%s:ptmu", var.c_str(), qn.c_str() ),
          stki_reduced_bins[i],muonPtbins );
    }
  }

  //Booking double transverse neutron plots
  MnvH2D *h_dphitt_q2qe[nHistos], *h_dphitt_ptmu[nHistos];
  MnvH2D *h_dthetaRt_q2qe[nHistos], *h_dthetaRt_ptmu[nHistos];
  MnvH2D *h_dthetaPt_q2qe[nHistos], *h_dthetaPt_ptmu[nHistos];
  MnvH2D *h_fitP_q2qe[nHistos], *h_fitP_ptmu[nHistos];
  MnvH2D *h_dFitP_q2qe[nHistos], *h_dFitP_ptmu[nHistos];
  MnvH2D *h_fitPRes_q2qe[nHistos], *h_fitPRes_ptmu[nHistos];
  MnvH2D *h_fitEnu_q2qe[nHistos], *h_fitEnu_ptmu[nHistos];
  MnvH2D *h_enuTrue_q2qe[nHistos], *h_enuTrue_ptmu[nHistos];
  MnvH2D *h_dFitEnu_q2qe[nHistos], *h_dFitEnu_ptmu[nHistos], *h_dEnu_q2qe[nHistos], *h_dEnu_ptmu[nHistos];
  utils->bookHistos( h_dphitt_q2qe,"h_dphitt_q2qe", "dphitt:q2", dphitbins, Q2bins );
  utils->bookHistos( h_dthetaRt_q2qe,"h_dthetaRt_q2qe","dthetaRt:q2",  dthetaReactbins, Q2bins );
  utils->bookHistos( h_dthetaPt_q2qe,"h_dthetaPt_q2qe","dthetapt:q2",  dthetaPerpbins, Q2bins );

  utils->bookHistos( h_dphitt_ptmu,"h_dphitt_ptmu", "dphitt:ptmu", dphitbins, muonPtbins );
  utils->bookHistos( h_dthetaRt_ptmu,"h_dthetaRt_ptmu", "dthetaR:ptmu", dthetaReactbins, muonPtbins );
  utils->bookHistos( h_dthetaPt_ptmu,"h_dthetaPt_ptmu", "dtheatP:ptmu", dthetaPerpbins, muonPtbins );


  utils->bookHistos( h_fitP_q2qe,"h_fitP_q2qe","fitP:q2",  kineticBins, Q2bins );
  utils->bookHistos( h_fitPRes_q2qe,"h_fitPRes_q2qe","fitPRes:q2",  kineticBins, Q2bins );
  utils->bookHistos( h_dFitP_q2qe,"h_dFitP_q2qe","dFitP:q2",  kineticResBins, Q2bins );
  utils->bookHistos( h_fitEnu_q2qe,"h_fitEnu_q2qe","fitEnu:q2",  leptonmomentumbins, Q2bins );
  utils->bookHistos( h_dFitEnu_q2qe,"h_dFitEnu_q2qe","dFitEnu:q2",  EnuResBins, Q2bins );
  utils->bookHistos( h_dEnu_q2qe,"h_dEnu_q2qe","dFitEnu:q2",  EnuResBins, Q2bins );
  utils->bookHistos( h_enuTrue_q2qe,"h_enuTrue_q2qe","enuTrue:q2",  leptonmomentumbins, Q2bins );

  utils->bookHistos( h_fitP_ptmu,"h_fitP_ptmu","fitP:ptmu",  kineticBins, muonPtbins );
  utils->bookHistos( h_fitPRes_ptmu,"h_fitPRes_ptmu","fitPRes:ptmu",  kineticBins, muonPtbins );
  utils->bookHistos( h_dFitP_ptmu,"h_dFitP_ptmu","dFitP:ptmu",  kineticResBins, muonPtbins );
  utils->bookHistos( h_fitEnu_ptmu,"h_fitEnu_ptmu","fitEnu:ptmu",  leptonmomentumbins, muonPtbins );
  utils->bookHistos( h_dFitEnu_ptmu,"h_dFitEnu_ptmu","dFitEnu:ptmu",  EnuResBins, muonPtbins );
  utils->bookHistos( h_dEnu_ptmu,"h_dEnu_ptmu","dFitEnu:ptmu",  EnuResBins, muonPtbins );
  utils->bookHistos( h_enuTrue_ptmu,"h_enuTrue_q2qe","enuTrue:q2",  leptonmomentumbins, muonPtbins );



  cout<<"booked quandrants"<<endl;


  //Vertex Energy Histograms
  MnvH2D *h_vtxEnergy_q2qe_lowRecoil[nHistos];
  utils->bookHistos( h_vtxEnergy_q2qe_lowRecoil, "h_vtxEnergy_q2qe_lowRecoil","vtxE : q2qe", VtxEnergyBins,Q2binsLow);


  MnvH2D *h_pRes_vs_node0[nHistos],*h_pRes_vs_node1[nHistos],*h_pRes_vs_node2[nHistos],*h_pRes_vs_node3[nHistos],*h_pRes_vs_nodeFirst[nHistos];
  utils->bookHistos( h_pRes_vs_node0,"h_pRes_vs_node0", "res:nodeE", 50,-1,1,10,0,30 );
  utils->bookHistos( h_pRes_vs_node1,"h_pRes_vs_node1", "res:nodeE", 50,-1,1,10,0,30 );
  utils->bookHistos( h_pRes_vs_node2,"h_pRes_vs_node2", "res:nodeE", 50,-1,1,10,0,30 );
  utils->bookHistos( h_pRes_vs_node3,"h_pRes_vs_node3", "res:nodeE", 50,-1,1,10,0,30 );
  utils->bookHistos( h_pRes_vs_nodeFirst,"h_pRes_vs_nodeFirst", "res:nodeE", 50,-1,1,10,0,30 );


  //Q2 tests: in angle_02
  MnvH2D *h_blobDist_q2qe[nHistos], *h_blobEnergy_q2qe[nHistos], *h_nBlobs_q2qe[nHistos],*h_n3DBlobs_q2qe[nHistos], *h_n2DBlobs_q2qe[nHistos], *h_muonTheta_q2qe[nHistos], *h_muonPhi_q2qe[nHistos], *h_muonE_q2qe[nHistos], *h_protonDEDX_q2qe[nHistos],*h_pionDEDX_q2qe[nHistos], *h_nonBlobEnergy_q2qe[nHistos], *h_hasTrack_q2qe[nHistos], *h_blobMaxE_q2qe[nHistos], *h_nClus_q2qe[nHistos], *h_blobEnergyLow_q2qe[nHistos];

  utils->bookHistos( h_blobDist_q2qe, "h_blobDist_q2qe", "dist:q2qe", BlobDistBins, Q2bins);
  utils->bookHistos( h_blobEnergy_q2qe, "h_blobEnergy_q2qe", "blobenergy:q2qe", BlobEnergyBins, Q2bins);
  utils->bookHistos( h_nonBlobEnergy_q2qe, "h_nonBlobEnergy_q2qe", "nonblobenergy:q2qe", BlobEnergyBins, Q2bins);
  utils->bookHistos( h_n3DBlobs_q2qe, "h_n3DBlobs_q2qe", "n3dblobs:q2qe", BlobNumBins, Q2bins );
  utils->bookHistos( h_n2DBlobs_q2qe, "h_n2DBlobs_q2qe", "n2dblobs:q2qe", BlobNumBins, Q2bins );
  utils->bookHistos( h_nBlobs_q2qe, "h_nBlobs_q2qe", "nblobs:q2qe", BlobNumBins, Q2bins );
  utils->bookHistos( h_muonTheta_q2qe, "h_muonTheta_q2qe", "muonTheta:q2qe", MuonThetaBins,Q2bins);
  utils->bookHistos( h_muonPhi_q2qe, "h_muonPhi_q2qe", "muonPhi:q2qe",MuonPhiBins, Q2bins );
  utils->bookHistos( h_muonE_q2qe, "h_muonE_q2qe", "muonE:q2qe",MuonEnergyBins, Q2bins );
  utils->bookHistos( h_protonDEDX_q2qe, "h_protonDEDX_q2qe", "protonDEDX:q2qe",10,-1,1,6,0.1,0.4);
  utils->bookHistos( h_pionDEDX_q2qe, "h_pionDEDX_q2qe", "pionDEDX:q2qe",10,-1,1,6,0.1,0.4);
  utils->bookHistos( h_hasTrack_q2qe, "h_hasTrack_q2qe", "hasTrack:q2qe",2,-.5,1.5,6,0.1,0.4);
  utils->bookHistos( h_blobMaxE_q2qe, "h_blobMaxE_q2qe", "blobMaxE:q2qe",BlobEnergyBins, Q2bins);
  utils->bookHistos( h_nClus_q2qe, "h_nClus_q2qe", "blobMaxE:q2qe",nClusBins, Q2bins);
  utils->bookHistos( h_blobEnergyLow_q2qe, "h_blobEnergyLow_q2qe", "blobenergy:q2qe", BlobEnergyBinsLow, Q2bins);



  MnvH2D *h_blobDist_q2qe_r1[nHistos], *h_blobEnergy_q2qe_r1[nHistos], *h_nBlobs_q2qe_r1[nHistos],*h_n3DBlobs_q2qe_r1[nHistos], *h_n2DBlobs_q2qe_r1[nHistos],  *h_nonBlobEnergy_q2qe_r1[nHistos], *h_blobMaxE_q2qe_r1[nHistos], *h_nClus_q2qe_r1[nHistos], *h_blobEnergyLow_q2qe_r1[nHistos];
  utils->bookHistos( h_blobDist_q2qe_r1, "h_blobDist_q2qe_r1", "dist:q2qe_r1", BlobDistBins, Q2bins);
  utils->bookHistos( h_blobEnergy_q2qe_r1, "h_blobEnergy_q2qe_r1", "blobenergy:q2qe_r1", BlobEnergyBins, Q2bins);
  utils->bookHistos( h_nonBlobEnergy_q2qe_r1, "h_nonBlobEnergy_q2qe_r1", "nonblobenergy:q2qe_r1", BlobEnergyBins, Q2bins);
  utils->bookHistos( h_n3DBlobs_q2qe_r1, "h_n3DBlobs_q2qe_r1", "n3dblobs:q2qe_r1", BlobNumBins, Q2bins );
  utils->bookHistos( h_n2DBlobs_q2qe_r1, "h_n2DBlobs_q2qe_r1", "n2dblobs:q2qe_r1", BlobNumBins, Q2bins );
  utils->bookHistos( h_nBlobs_q2qe_r1, "h_nBlobs_q2qe_r1", "nblobs:q2qe_r1", BlobNumBins, Q2bins );
  utils->bookHistos( h_blobMaxE_q2qe_r1, "h_blobMaxE_q2qe_r1", "blobMaxE:q2qe_r1",BlobEnergyBins, Q2bins);
  utils->bookHistos( h_nClus_q2qe_r1, "h_nClus_q2qe_r1", "blobMaxE:q2qe_r1",nClusBins, Q2bins);
  utils->bookHistos( h_blobEnergyLow_q2qe_r1, "h_blobEnergyLow_q2qe_r1", "blobenergy:q2qe", BlobEnergyBinsLow, Q2bins);

  MnvH2D *h_blobDist_q2qe_r5[nHistos], *h_blobEnergy_q2qe_r5[nHistos], *h_nBlobs_q2qe_r5[nHistos],*h_n3DBlobs_q2qe_r5[nHistos], *h_n2DBlobs_q2qe_r5[nHistos],  *h_nonBlobEnergy_q2qe_r5[nHistos], *h_blobMaxE_q2qe_r5[nHistos], *h_nClus_q2qe_r5[nHistos], *h_blobEnergyLow_q2qe_r5[nHistos];
  utils->bookHistos( h_blobDist_q2qe_r5, "h_blobDist_q2qe_r5", "dist:q2qe_r5", BlobDistBins, Q2bins);
  utils->bookHistos( h_blobEnergy_q2qe_r5, "h_blobEnergy_q2qe_r5", "blobenergy:q2qe_r5", BlobEnergyBins, Q2bins);
  utils->bookHistos( h_nonBlobEnergy_q2qe_r5, "h_nonBlobEnergy_q2qe_r5", "nonblobenergy:q2qe_r5", BlobEnergyBins, Q2bins);
  utils->bookHistos( h_n3DBlobs_q2qe_r5, "h_n3DBlobs_q2qe_r5", "n3dblobs:q2qe_r5", BlobNumBins, Q2bins );
  utils->bookHistos( h_n2DBlobs_q2qe_r5, "h_n2DBlobs_q2qe_r5", "n2dblobs:q2qe_r5", BlobNumBins, Q2bins );
  utils->bookHistos( h_nBlobs_q2qe_r5, "h_nBlobs_q2qe_r5", "nblobs:q2qe_r5", BlobNumBins, Q2bins );
  utils->bookHistos( h_blobMaxE_q2qe_r5, "h_blobMaxE_q2qe_r5", "blobMaxE:q2qe_r5",BlobEnergyBins, Q2bins);
  utils->bookHistos( h_nClus_q2qe_r5, "h_nClus_q2qe_r5", "blobMaxE:q2qe_r5",nClusBins, Q2bins);
  utils->bookHistos( h_blobEnergyLow_q2qe_r5, "h_blobEnergyLow_q2qe_r5", "blobenergy:q2qe_r5", BlobEnergyBinsLow, Q2bins);

  MnvH2D *h_blobDist_q2qe_r99[nHistos], *h_blobEnergy_q2qe_r99[nHistos], *h_nBlobs_q2qe_r99[nHistos],*h_n3DBlobs_q2qe_r99[nHistos], *h_n2DBlobs_q2qe_r99[nHistos],  *h_nonBlobEnergy_q2qe_r99[nHistos], *h_blobMaxE_q2qe_r99[nHistos], *h_nClus_q2qe_r99[nHistos], *h_blobEnergyLow_q2qe_r99[nHistos];
  utils->bookHistos( h_blobDist_q2qe_r99, "h_blobDist_q2qe_r99", "dist:q2qe_r99", BlobDistBins, Q2bins);
  utils->bookHistos( h_blobEnergy_q2qe_r99, "h_blobEnergy_q2qe_r99", "blobenergy:q2qe_r99", BlobEnergyBins, Q2bins);
  utils->bookHistos( h_nonBlobEnergy_q2qe_r99, "h_nonBlobEnergy_q2qe_r99", "nonblobenergy:q2qe_r99", BlobEnergyBins, Q2bins);
  utils->bookHistos( h_n3DBlobs_q2qe_r99, "h_n3DBlobs_q2qe_r99", "n3dblobs:q2qe_r99", BlobNumBins, Q2bins );
  utils->bookHistos( h_n2DBlobs_q2qe_r99, "h_n2DBlobs_q2qe_r99", "n2dblobs:q2qe_r99", BlobNumBins, Q2bins );
  utils->bookHistos( h_nBlobs_q2qe_r99, "h_nBlobs_q2qe_r99", "nblobs:q2qe_r99", BlobNumBins, Q2bins );
  utils->bookHistos( h_blobMaxE_q2qe_r99, "h_blobMaxE_q2qe_r99", "blobMaxE:q2qe_r99",BlobEnergyBins, Q2bins);
  utils->bookHistos( h_nClus_q2qe_r99, "h_nClus_q2qe_r99", "blobMaxE:q2qe_r99",nClusBins, Q2bins);
  utils->bookHistos( h_blobEnergyLow_q2qe_r99, "h_blobEnergyLow_q2qe_r99", "blobenergy:q2qe_r99", BlobEnergyBinsLow, Q2bins);



  //misc neutron diagnostics
  MnvH2D *h_blobDist_blobEnergy[nHistos]; utils->bookHistos( h_blobDist_blobEnergy, "h_blobDist_blobEnergy", "dist:energy", BlobDistBins, BlobEnergyBins);
  MnvH2D *h_blobDist_blobEnergyLow[nHistos]; utils->bookHistos( h_blobDist_blobEnergyLow, "h_blobDist_blobEnergyLow", "dist:energy", BlobDistBins, BlobEnergyBinsLow);
  MnvH2D *h_nHits_totalE[nHistos]; utils->bookHistos( h_nHits_totalE, "h_nHits_totalE", "nHits:histsE", nHitsBins, hitsEBins );



  MnvH2D *h_muontheta_q2qe[nHistos], *h_muontheta_ptmu[nHistos];
  MnvH2D *h_muonmomentum_q2qe[nHistos], *h_muonmomentum_ptmu[nHistos];
  MnvH2D *h_enu_q2qe[nHistos], *h_enu_ptmu[nHistos];
  utils->bookHistos( h_muontheta_q2qe, "h_muontheta_q2qe", "#theta_{#mu} (degree): Q^{2}_{QE} (GeV^{2})", thetabins, Q2bins );
  utils->bookHistos( h_muonmomentum_q2qe, "h_muonmomentum_q2qe", "p_{#mu} (GeV/c): Q^{2}_{QE} (GeV^{2})", leptonmomentumbins, Q2bins );
  utils->bookHistos( h_enu_q2qe, "h_enu_q2qe", "E_{#nu} (GeV/c): Q^{2}_{QE} (GeV^{2})", leptonmomentumbins, Q2bins );
  utils->bookHistos( h_muontheta_ptmu, "h_muontheta_ptmu", "#theta_{#mu} (degree): p_{T,#mu} (GeV/c))", thetabins, muonPtbins );
  utils->bookHistos( h_muonmomentum_ptmu, "h_muonmomentum_ptmu", "p_{#mu} (GeV/c): p_{T,#mu} (GeV/c))", leptonmomentumbins, muonPtbins );
  utils->bookHistos( h_enu_ptmu, "h_enu_ptmu", "E_{#nu} (GeV/c): p_{T,#mu} (GeV/c)", leptonmomentumbins, muonPtbins );


  MnvH2D* h_multiplicity_q2qe[nHistos];
  utils->bookHistos( h_multiplicity_q2qe, "h_multiplicity_q2qe", "mult : q2qe", BlobNumBins, Q2bins );

  //truth match
  MnvH2D* h_blobE_PartType[nHistos];
  MnvH2D* h_blobMaxE_PartType[nHistos];
  utils->bookHistos( h_blobE_PartType, "h_blobE_PartType", "energy:PartType", BlobEnergyBinsFine,PartTypeBins );
  utils->bookHistos( h_blobMaxE_PartType, "h_blobMaxE_PartType", "energy:PartType", BlobEnergyBinsFine,PartTypeBins );

  MnvH2D* h_blobE3D_PartType[nHistos];
  MnvH2D* h_blobMaxE3D_PartType[nHistos];
  utils->bookHistos( h_blobE3D_PartType, "h_blobE3D_PartType", "energy:PartType", BlobEnergyBinsFine,PartTypeBins );
  utils->bookHistos( h_blobMaxE3D_PartType, "h_blobMaxE3D_PartType", "energy:PartType", BlobEnergyBinsFine,PartTypeBins );

  MnvH2D* h_mainBlobE_PartType[nHistos];
  MnvH2D* h_mainBlobMaxE_PartType[nHistos];
  utils->bookHistos( h_mainBlobE_PartType, "h_mainBlobE_PartType", "energy:PartType", BlobEnergyBinsFine,PartTypeBins );
  utils->bookHistos( h_mainBlobMaxE_PartType, "h_mainBlobMaxE_PartType", "energy:PartType", BlobEnergyBinsFine,PartTypeBins );

  for( uint i = 0; i<nHistos; i++ )
  {
    for( uint j = 0; j< h_blobE_PartType[i]->GetNbinsY(); j++ )
    {
      h_blobE_PartType[i]->GetYaxis()->SetBinLabel( j+1, (returnPartName(j)).c_str() );
      h_blobMaxE_PartType[i]->GetYaxis()->SetBinLabel( j+1, (returnPartName(j)).c_str() );
      h_mainBlobE_PartType[i]->GetYaxis()->SetBinLabel( j+1, (returnPartName(j)).c_str() );
      h_mainBlobMaxE_PartType[i]->GetYaxis()->SetBinLabel( j+1, (returnPartName(j)).c_str() );
    }
  }

  //MnvH3D *h_dthetaP_dthetaR_q2qe[nHistos];
  //utils->bookHistos( h_dthetaP_dthetaR_q2qe, "h_dthetaP_dthetaR_q2qe"k


  //utils->SetHyperDim(my3d);  
  //recoil Study

 
  // Add Vertical and Lateral Error Bands
  // JO is removing error bands from all 1D because they are all 
  // present in 2D. Plz make projections directly from the 2D. 
  //--------------------------------------------------------------
  if(RunCodeWithSystematics){
      utils->addVertErrorBands( h_ptmu );
      utils->addLatErrorBands ( h_ptmu );
      
      utils->addVertErrorBands( h_q2qe );
      utils->addLatErrorBands ( h_q2qe );

      utils->addVertErrorBands( h_enu );
      utils->addLatErrorBands ( h_enu );

      utils->addVertErrorBands( h_dthetaP );
      utils->addVertErrorBands( h_dthetaR );
      utils->addVertErrorBands( h_pn );
      utils->addLatErrorBands( h_dthetaP );
      utils->addLatErrorBands( h_dthetaR );
      utils->addLatErrorBands( h_pn );

      utils->addVertErrorBands( h_q2qe_ptmu );
      utils->addLatErrorBands ( h_q2qe_ptmu );

      for( int i = 0; i< 10; i++ )
      {
        utils->addVertErrorBands( &h_q2qe_angles[i][0] );
        utils->addLatErrorBands( &h_q2qe_angles[i][0] );
      }


      for( unsigned int i = 0; i< region_names.size(); i++ )
      {
        utils->addVertErrorBands( &h_q2qe_ptmu_regions[i][0] );
        utils->addLatErrorBands ( &h_q2qe_ptmu_regions[i][0] );
        utils->addVertErrorBands( &h_enu_q2qe_regions[i][0] );
        utils->addLatErrorBands ( &h_enu_q2qe_regions[i][0] );
        utils->addVertErrorBands( &h_muontheta_q2qe_regions[i][0] );
        utils->addLatErrorBands ( &h_muontheta_q2qe_regions[i][0] );
      }



      //2-D
      //STKI vs XVAR
      for( unsigned int i = 0; i < stki_names.size(); i++ )
      {
        utils->addVertErrorBands( &h_stki_ptmu[i][0] );
        utils->addLatErrorBands( &h_stki_ptmu[i][0] );

        utils->addVertErrorBands( &h_stki_q2qe[i][0] );
        utils->addLatErrorBands( &h_stki_q2qe[i][0] );
      }
      //NeutVar vs XVAR
      for( unsigned int i = 0; i < neut_names.size(); i++ )
      {
        utils->addVertErrorBands( &h_neut_ptmu[i][0] );
        utils->addLatErrorBands( &h_neut_ptmu[i][0] );

        utils->addVertErrorBands( &h_neut_q2qe[i][0] );
        utils->addLatErrorBands( &h_neut_q2qe[i][0] );
      }


      utils->addVertErrorBands( h_dthetaP_center_ptmu );
      utils->addLatErrorBands( h_dthetaP_center_ptmu );

      utils->addVertErrorBands( h_dthetaP_center_q2qe );
      utils->addLatErrorBands( h_dthetaP_center_q2qe );


      utils->addVertErrorBands( h_dthetaP_dthetaR );
      utils->addVertErrorBands( h_dthetaP_dthetaR_Fine );
      utils->addVertErrorBands( h_dthetaPpos_center_q2qe );
      utils->addVertErrorBands( h_recoil_inc_q2qe );
      utils->addVertErrorBands( h_recoil_inc_q3 );

      utils->addLatErrorBands( h_dthetaP_dthetaR );
      utils->addLatErrorBands( h_dthetaP_dthetaR_Fine );
      utils->addLatErrorBands( h_dthetaPpos_center_q2qe );
      utils->addLatErrorBands( h_recoil_inc_q2qe );
      utils->addLatErrorBands( h_recoil_inc_q3 );


      utils->addVertErrorBands( h_muontheta_q2qe );
      utils->addVertErrorBands( h_muontheta_ptmu );
      utils->addVertErrorBands( h_muonmomentum_q2qe );
      utils->addVertErrorBands( h_muonmomentum_ptmu );
      utils->addVertErrorBands( h_enu_q2qe );
      utils->addVertErrorBands( h_enu_ptmu );

      utils->addLatErrorBands( h_muontheta_q2qe );
      utils->addLatErrorBands( h_muontheta_ptmu );
      utils->addLatErrorBands( h_muonmomentum_q2qe );
      utils->addLatErrorBands( h_muonmomentum_ptmu );
      utils->addLatErrorBands( h_enu_q2qe );
      utils->addLatErrorBands( h_enu_ptmu );

      utils->addVertErrorBands( h_multiplicity_q2qe );
      utils->addLatErrorBands( h_multiplicity_q2qe );

      
      utils->addVertErrorBands( h_dphitt_q2qe );
      utils->addVertErrorBands( h_dphitt_ptmu );
      utils->addVertErrorBands( h_dthetaRt_q2qe );
      utils->addVertErrorBands( h_dthetaRt_ptmu );
      utils->addVertErrorBands( h_dthetaPt_q2qe );
      utils->addVertErrorBands( h_dthetaPt_ptmu );
      utils->addLatErrorBands( h_dphitt_q2qe );
      utils->addLatErrorBands( h_dphitt_ptmu );
      utils->addLatErrorBands( h_dthetaRt_q2qe );
      utils->addLatErrorBands( h_dthetaRt_ptmu );
      utils->addLatErrorBands( h_dthetaPt_q2qe );
      utils->addLatErrorBands( h_dthetaPt_ptmu );



      utils->addVertErrorBands( h_fitP_q2qe );
      utils->addVertErrorBands( h_fitEnu_q2qe );
      utils->addVertErrorBands( h_fitPRes_q2qe );

      utils->addVertErrorBands( h_fitP_ptmu );
      utils->addVertErrorBands( h_fitEnu_ptmu );
      utils->addVertErrorBands( h_fitPRes_ptmu );

      utils->addLatErrorBands( h_fitP_q2qe );
      utils->addLatErrorBands( h_fitEnu_q2qe );
      utils->addLatErrorBands( h_fitPRes_q2qe );

      utils->addLatErrorBands( h_fitP_ptmu );
      utils->addLatErrorBands( h_fitEnu_ptmu );
      utils->addLatErrorBands( h_fitPRes_ptmu );


  }//end run with systematics




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
  
  //---------------------------------------------
  // Get Data Tree
  //---------------------------------------------
  if (n_datafiles==0) {
    cout<<"Data files won't be included"<<endl;
  } else {
    TChain* tree_data  = utils->getDataTree("CCQENu", n_datafiles );
    int entries_data   = tree_data->GetEntries();

    cout << "Data entries: " << entries_data << endl;
    bool isData = true;
    CCQENuEvent* data = new CCQENuEvent( tree_data, isData );
    isData = true;
    isTruth = false;
    data_events = AnalysisLooper( data, utils, cutter, ncutter, sample,  multiplicity, entries_data, data_counter, isData, isTruth, npMode );
  }

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
  cout<<" --------------------------------------------- " << endl;
  cout<<" DATA " << endl;
  cout<<" --------------------------------------------- " << endl;
  cout<<" ****************** " << endl;
  cout<<" DATA 1 Track Events " << endl;
  cout<<" ****************** " << endl;
  cout<< Form("nocuts = %i; cut1 = %i, cut2 = %i, cut3 = %i, cut4 = %i, cut5 = %i, cut6 = %i, cut6_sideband = %i",data_counter[10],data_counter[11], data_counter[12], data_counter[13],data_counter[14],data_counter[15],data_counter[16],data_counter[17]) << endl;
  cout<<" ****************** " << endl;
  cout<<" DATA 2 Track Events " << endl;
  cout<<" ****************** " << endl;
  cout<< Form("nocuts = %i; cut1 = %i, cut2 = %i, cut3 = %i, cut4 = %i, cut5 = %i, cut6 = %i, cut7 = %i, cut7_sideband = %i",data_counter[20],data_counter[21], data_counter[22], data_counter[23],data_counter[24],data_counter[25],data_counter[26],data_counter[27]) << endl;
  
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
  par.push_back("0");//multiplicity
  par.push_back("-1");//mc file
  par.push_back("-1");//ndata file
  par.push_back("0");//proton, neutron,proton/neutron mode,1,0,10
  par.push_back("0");//drop clusters?

  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }

  bool fluxConstraintHistoNeeded = ( par.at(4) == "1" ) ? true : false; 

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;

//int MuonSelectionHists(string filename, string sample, bool makeFluxConstraintHisto, int multiplicity, string playlist, int n_mcfiles = -1, int n_datafiles = -1, int npMode = 1)
  string output_path = par[1];
  string sample = par[3];
  int multiplicity = atoi(par[5].c_str() );
  string playlist = par[2];
  int n_mc = atoi(par[6].c_str() );
  int n_data = atoi(par[7].c_str() );
  int npmode = atoi(par[8].c_str() );
  dropClusters = bool(atoi(par[9].c_str() ) );

  GlobalParameters::Get().m_analysisType = (npmode == 10)? kPNTransverse : kNeutron;
  pass_neutronQsqCut = (npmode==10)? -1 : 1;
//
  //return MuonSelectionHists(output_path, sample, fluxConstraintHistoNeeded, multiplicity, par[2],atoi(par[6].c_str()), atoi(par[7].c_str()), atoi(par[8].c_str()) );
  return MuonSelectionHists(output_path, sample, fluxConstraintHistoNeeded, multiplicity, playlist, n_mc, n_data, npmode );
}
