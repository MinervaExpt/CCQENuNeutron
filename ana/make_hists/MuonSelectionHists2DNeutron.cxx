#include "include/GlobalParameters.h"
#include "include/CCQENuUtils.h"
#include "TParameter.h" 
#include "PlotUtils/HyperDimLinearizer.h"

#include "include/GeneralFunc.h" //in CCQENuNeutron/include/
#include "include/EventCuts.h"
#include "include/CommonBins.h"

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


double SolveNeutronP( XYZVector &neutronDir, XYZTVector &protonP )
{
  //XYZVector baseAxis(0,0,1);
  //AxisAngle R1( neutronDir.Cross( baseAxis ), ACos( baseAxis.Dot( neutronDir.Unit() ) ) );
  //XYZTVector protonR1 = (R1*protonP );
  double Mn = 0.9395654133, Mp=0.93827231;
  double Mn2 = Mn*Mn, Mp2 = Mp*Mp;
  double Ep = protonP.T(), Pz = protonP.Vect().Dot( neutronDir.Unit() );
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

template<class T>
double AnalysisLooperBase(T* evt, CCQENuUtils* utils, CCQENuCuts* cutter, NeutronBlobCuts* ncutter, string sample, int multiplicity, int n_entries, map<int,int> &cuts_counter, bool isData=false, bool isTruth=false, int npMode = 1)
{
  GeoUtils* geoUtils = new GeoUtils();
  PhysicsUtils* physicsUtils = new PhysicsUtils();

  double n_evts=0;

  //return -1;
  for( int i = 0; i< n_entries; ++i)
  {
    evt->GetEntry(i);
    if( IgnoreEvent( evt, isData )) continue;
    if (i%1000 == 0) cout<<"At Event "<<i<<" ("<<i*100./n_entries<<"%)"<<endl;
    //cout<<"at Event: "<<i<<endl;
    //==============Fill Neutron Blobs======================
    bool isMC = !isData;
    ncutter->UpdateCandidates( evt, isMC );
    if(isDebug) cout<<"Updated ncutter"<<endl;
    //==============Event Selection========================
    bool PassEventCut = (EventCut( evt, cutter, ncutter, cuts_counter, multiplicity, sample ,  isData, npMode));
    bool PassEventCutBase = EventCutBase( evt, cutter, ncutter, cuts_counter, isData, multiplicity, 1 , 1 , 1 , -1 , -1 , -1,-1,-1);
    if(isDebug) cout<<PassEventCutBase<<endl;
    if (!PassEventCutBase) continue;
    bool hasNC = false;
    bool hasProton = false;


    NeutronCandidates* ptrNeutronCandidates; 
    NeutronBlob* mainCandidate=NULL; 
    if (npMode != 1 ) 
    {
      ptrNeutronCandidates = ncutter->GetCandidates();
    }





    bool fill_common_histos = true;
    //===============Reconstruction ========================
    if(isDebug) cout<<"Getting CV Weight:"<<endl;
    double wgt = (isData)? 1: utils->GetCVWeight( evt, sample );
    n_evts+=wgt;
    if(isDebug) cout<<wgt<<endl;

    //-----------
    //Vertex
    XYZVector vertex( evt->vtx[0],evt->vtx[1],evt->vtx[2]);
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

    XYZVector xcoord(1,0,0);
    XYZVector ycoord = beam.Unit().Cross(xcoord);
    double muon_phi_detector = ACos( reco_muon_py/ reco_muon_px )*180./3.14159;


    //2D Neutron Variables:
    vector<NeutronBlob*> AllBlobs = ptrNeutronCandidates->AllBlobs;
    if(AllBlobs.size() == 0) continue;
    BlobSort::sort_blobVec(AllBlobs, false, "ClusMaxE" );

      NeutronBlob* blob = AllBlobs[0];
      if ( blob->Is3D ) continue;
      int view;
      XYZVector blobposT = ptrNeutronCandidates->GetClusterTPosClosestToVtx( blob, view );
      //XYZVector blobposT = blob->BlobPosT;
      blob->BlobPosT = blobposT;
      XYZVector vtx2D = ncutter->DetectorUtilsPtr()->ViewProjection2D( vertex, view ) ;
      XYZVector dPos = blobposT - vtx2D;

      //if( dPos.Z() < 100 ) continue;//downstream
      //if( abs(dPos.Z())<100 && abs(dPos.X())<100 ) continue; // not in vertex region
      XYZVector mu2D = (ncutter->DetectorUtilsPtr()->ViewProjection2D( reco_muon_3P, view )).Unit() ;
      if (ACos( mu2D.Dot( dPos.Unit() ) )*180/3.14159 < 15 ) continue;
      if ( mu2D.Cross( dPos ).R() < 100 ) continue;
      mainCandidate = blob;

    if (!mainCandidate) continue;

    


    //Proton Variables
    TVector3 protonVect;
    XYZVector protonVectXYZ;
    double protonAngle, protonMom, protonTn;

    if (npMode != 0 )
    {
      protonVect.SetXYZ(evt->CCQENu_proton_Px_fromdEdx,evt->CCQENu_proton_Py_fromdEdx,evt->CCQENu_proton_Pz_fromdEdx);//MeV
      protonVectXYZ.SetXYZ(evt->CCQENu_proton_Px_fromdEdx,evt->CCQENu_proton_Py_fromdEdx,evt->CCQENu_proton_Pz_fromdEdx);//MeV
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
    if (npMode !=1 || hasNC)
    {
      blobVtxDir = ncutter->MainNeutronCandidate()->BlobVtxDir;
      neutronDir.SetXYZ(blobVtxDir.X(),blobVtxDir.Y(), blobVtxDir.Z() ) ; 
    }

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
        fspart_res_angle = acos( true_fs_nucleon_p_comp.Unit().Dot( neutronDir.Unit() ) )*180/3.14159;
        fspart_res_angle_x = acos( true_fs_nucleon_p_comp.Unit().Dot( (neutronDir.Dot( hx )*hx).Unit() ) )*180/3.14159;
        fspart_res_angle_y = acos( true_fs_nucleon_p_comp.Unit().Dot( (neutronDir.Dot( hy )*hy).Unit() ) )*180/3.14159;

      }
    }


    //Expected nucleon kinematics

    // 2. binding and is/fs mass
    double binding_e = 0; 
    double is_part_mass = (neutrinoMode)? 0.9395654133: 0.93827231; //neutron/proton
    double fs_part_mass = (neutrinoMode)? 0.93827231: 0.9395654133;

    // 3. expected nucleon

    XYZTVector expectedNucleon4P = geoUtils->ComputeExpectedNucleon( beam, reco_muon_4P, is_part_mass, fs_part_mass, binding_e);
    XYZVector expectedNucleon3P = (XYZVector) expectedNucleon4P.Vect();


    // Set the primary nucleon momentum vec
    // proton/pn mode : proton. neutron mode: neutron direction
    TVector3 primaryNucleonVect = ( npMode == 0 )? neutronDir.Unit()*expectedNucleon3P.R() : protonVect ;
    //cout<<"Is proton primary? "<< (primaryNucleonVect==protonVect)<<endl;

    XYZVector primaryNucleonVect_XYZ( primaryNucleonVect.X(), primaryNucleonVect.Y(), primaryNucleonVect.Z() );

    //Calculating Q2qe:
    int charge = (neutrinoMode)? -1:1;
    double q2qe = physicsUtils->qsqCCQE( beam, reco_muon_4P, charge, binding_e ); //GeV^2
    double enu = physicsUtils->nuEnergyCCQE( beam, reco_muon_4P, charge, binding_e );



    //preparing for STKI and angular vars:
    double Mcarbon = 6*fs_part_mass+6*is_part_mass-92.163/1E3;
    double Enuc=sqrt(primaryNucleonVect_XYZ.Mag2()+fs_part_mass*fs_part_mass);
    //XYZTVector primaryNucleonVect_XYZT( protonVect.X(), protonVect.Y(), protonVect.Z(), Enuc );
    XYZTVector primaryNucleonVect_XYZT( primaryNucleonVect_XYZ.X(), primaryNucleonVect_XYZ.Y(), primaryNucleonVect_XYZ.Z(), Enuc );

    //calculating transverse variables
    //calculating angular variables
    std::vector<double> reco_neutron_vars = geoUtils->ComputeNeutronAngularVars( beam, expectedNucleon3P, primaryNucleonVect_XYZ );


    //Get 2D angulars
    int clusView = -1;
    XYZVector clusTPos = ptrNeutronCandidates->GetClusterTPosClosestToVtx( mainCandidate, clusView );
    mainCandidate->BlobPosT = clusTPos;
    //cout<<"ClusView: "<<clusView<<" --- "<<endl;
    //cout<<"TPos: "<<clusTPos.X()<<", "<< clusTPos.Y()<<", "<<clusTPos.Z() <<endl;
    //cout<<"MCPos: " << mainCandidate->MCTrackPos.X()<<", "<<mainCandidate->MCTrackPos.Y()<<", "<<mainCandidate->MCTrackPos.Z()<<", "<<endl;
    //cout<<"bTPos: "<<mainCandidate->BlobPosT.X()<<", "<<endl;
    XYZVector mc2DTpos =( ncutter->DetectorUtilsPtr()->ViewProjection2D(mainCandidate->MCTrackPos.Vect(), clusView ) ); 
    //cout<<"mTPos: "<<mc2DTpos.X()<<endl;
    //
    XYZVector blobDir2D = (mainCandidate->BlobPosT - vtx2D).Unit();
    //cout<<"blobDir2D: "<<blobDir2D.X()<<", " <<blobDir2D.Y()<<", " <<blobDir2D.Z()<<endl;
    XYZVector nu2D =( ncutter->DetectorUtilsPtr()->ViewProjection2D( beam, view ) ).Unit();
    XYZVector expected2D = ( ncutter->DetectorUtilsPtr()->ViewProjection2D( expectedNucleon3P, view ) ).Unit();

    XYZVector dtheta2DVec = blobDir2D.Cross( expected2D );
    double dtheta2D = ACos( blobDir2D.Dot( nu2D ) );
    double expTheta2D = ACos( expected2D.Dot( nu2D ) );
    dtheta2D-= expTheta2D;
    dtheta2D*=180/TMath::Pi();
    //int sign = dtheta2DVec.Y()>0? 1:-1;
    //int sign = dtheta2DVec.Y(>0? 1:-1;
    //double dtheta2D = sign*ASin(dtheta2DVec.R())*180/TMath::Pi();
    //if(abs(dtheta2D)>80) cout<<dtheta2D<<endl;
    //cout<<dtheta2D<<endl;
    //if( evt->mc_intType == 0 && evt->mc_targetA == 1 && mainCandidate->View == 1) 
    //{
    //  cout<<dtheta2D<<endl;
    //  cout<<mainCandidate->BlobPosT.X()<<", "<<mainCandidate->MCTrackPos.X()<<", "<<wgt<<endl;
    //}
    // test if the expected and muon cancel, they should

    double reco_dthetaP = reco_neutron_vars[vdthetaP];
    double reco_dthetaR = reco_neutron_vars[vdthetaR];
    double reco_dtheta = reco_neutron_vars[vdtheta];
    //if (isMC) cout<< evt->mc_intType<<endl;
    reco_dthetaP *= 180/3.14159;
    reco_dthetaR *= 180/3.14159;
    reco_dtheta *= 180/3.14159;
    //cout<<reco_dthetaP<<endl;
    double reco_dPp = reco_neutron_vars[vdPp];
    double reco_dPr = reco_neutron_vars[vdPr];
    double reco_dPp_infer = reco_neutron_vars[vdPpi];
    double reco_dPr_infer = reco_neutron_vars[vdPri];

    string quadrant="";
    if( hasNC && hasProton ) 
    {
      quadrant = GetQuadrantCode( beam, reco_muon_3P, protonVectXYZ, blobVtxDir );
      //cout<<quadrant<<endl;
    }

    if (isMC)
    {
      //get dthetaR/P resolution
      double true_fsp_p = true_fs_nucleon_p_comp.Mag(); 
      XYZVector true_fsp_3p( true_fs_nucleon_p_comp.X(), true_fs_nucleon_p_comp.Y(), true_fs_nucleon_p_comp.Z()  );
      //calculating angular variables
      std::vector<double> mc_neutron_vars = geoUtils->ComputeNeutronAngularVars( beam, expectedNucleon3P, true_fsp_3p );
      // test if the expected and muon cancel, they should

      double mc_dThetaP = mc_neutron_vars[vdthetaP];
      double mc_dThetaR = mc_neutron_vars[vdthetaR];
      //if (isMC) cout<< evt->mc_intType<<endl;
      mc_dThetaP *= 180/3.14159;
      mc_dThetaR *= 180/3.14159;
      fspart_res_dthetaR = reco_dthetaR - mc_dThetaR;
      fspart_res_dthetaP = reco_dthetaP - mc_dThetaP;


      if(fill_common_histos )
      {
        utils->fillHistosV3( utils->histos2D["h_true_vs_reco_dthetaP"], mc_dThetaP, reco_dthetaP, isData, evt, wgt);
        utils->fillHistosV3( utils->histos2D["h_true_vs_reco_dthetaR"], mc_dThetaR, reco_dthetaR, isData, evt, wgt);
      }

      
    }

    //Recoil Characteristic
    double recoil =ncutter->GetDefaultRecoilEnergy(evt);
    double vtxE  = evt->recoil_energy_nonmuon_vtx100mm;

    //Blob Characteristics
    double blobDist, blobEnergy, nBlobs, n2DBlobs, n3DBlobs;
    double nonBlobEnergy;
    double proton_score=-1, pion_score=-1;
    if( true )
    {
      blobDist = dPos.R();
      blobEnergy = mainCandidate->TotalE;
      nonBlobEnergy = recoil - blobEnergy;
      n3DBlobs = ptrNeutronCandidates->ThreeDBlobs.size();
      n2DBlobs = ptrNeutronCandidates->TwoDBlobs.size();
      nBlobs = ptrNeutronCandidates->AllBlobs.size();
      //cout<<blobDist<<", "<<blobEnergy<<endl;
      if( mainCandidate->hasTrack )
      {
        proton_score = mainCandidate->dEdX_proton;
        pion_score = mainCandidate->dEdX_pion;

      }
    }


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
      utils->fillHistosV3( utils->histos1D["h_pmu"], muon_p, isData, evt, wgt );
      //utils->fillHistosV3( utils->histos1D["h_enu"], reco_En,  isData, evt, wgt );
      
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
      //utils->fillHistosV3( utils->histos2D["h_sign_ptmu"],reco_sign,muon_pt_beam, isData, evt,wgt);
      //utils->fillHistosV3( utils->histos2D["h_signdalphat_ptmu"],reco_signed_dalphat,muon_pt_beam, isData, evt,wgt);
      //utils->fillHistosV3( utils->histos2D["h_signdphit_ptmu"],reco_signed_dphit,muon_pt_beam, isData, evt,wgt);

      utils->fillHistosV3( utils->histos2D["h_dthetaP_ptmu"],reco_dthetaP,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dthetaR_ptmu"],reco_dthetaR,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPP_ptmu"],reco_dPp,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPR_ptmu"],reco_dPr,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPPi_ptmu"],reco_dPp_infer,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPRi_ptmu"],reco_dPr_infer,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dtheta2D_ptmu"],   dtheta2D,muon_pt_beam, isData, evt,wgt);

      utils->fillHistosV3( utils->histos2D["h_dthetaP_q2qe"],reco_dthetaP,  q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dthetaR_q2qe"],reco_dthetaR,  q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPP_q2qe"],    reco_dPp,      q2qe, isData,evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPR_q2qe"],    reco_dPr,      q2qe, isData,evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPPi_q2qe"],   reco_dPp_infer,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dPRi_q2qe"],   reco_dPr_infer,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dtheta2D_q2qe"],   dtheta2D,q2qe, isData, evt,wgt);

      utils->fillHistosV3( utils->histos2D["h_dtheta2D_recoil"],   dtheta2D,recoil, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dtheta2D_vtxE"],   dtheta2D,vtxE, isData, evt,wgt);


      utils->fillHistosV3( utils->histos2D["h_dthetaP_dthetaR"],reco_dthetaP,  reco_dthetaR, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_dthetaP_dthetaR_Fine"],reco_dthetaP,  reco_dthetaR, isData, evt,wgt);

      utils->fillHistosV3( utils->histos2D["h_res_dthetaP_ptmu"],   fspart_res_dthetaP ,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dthetaR_ptmu"],   fspart_res_dthetaR ,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dtheta_ptmu"],   fspart_res_angle ,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dthetaX_ptmu"],   fspart_res_angle_x ,muon_pt_beam, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dthetaY_ptmu"],   fspart_res_angle_y ,muon_pt_beam, isData, evt,wgt);

      utils->fillHistosV3( utils->histos2D["h_res_dthetaP_q2qe"],   fspart_res_dthetaP ,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dthetaR_q2qe"],   fspart_res_dthetaR ,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dtheta_q2qe"],   fspart_res_angle ,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dthetaX_q2qe"],   fspart_res_angle_x ,q2qe, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dthetaY_q2qe"],   fspart_res_angle_y ,q2qe, isData, evt,wgt);

      utils->fillHistosV3( utils->histos2D["h_res_dthetaP"],   fspart_res_dthetaP ,reco_dthetaP, isData, evt,wgt);
      utils->fillHistosV3( utils->histos2D["h_res_dthetaR"],   fspart_res_dthetaR ,reco_dthetaR, isData, evt,wgt);


      if(hasNC && hasProton)
      {
        //utils->fillHistosV3( utils->histos2D[ Form("h_pn_q2qe_%s", quadrant.c_str())], reco_pn, q2qe, isData,evt,wgt);
        //utils->fillHistosV3( utils->histos2D[ Form("h_dphit_q2qe_%s", quadrant.c_str())], reco_dphit, q2qe, isData,evt,wgt);
        //utils->fillHistosV3( utils->histos2D[ Form("h_dalphat_q2qe_%s", quadrant.c_str())], reco_dalphat, q2qe, isData,evt,wgt);
        //utils->fillHistosV3( utils->histos2D[ Form("h_dpt_q2qe_%s", quadrant.c_str())], reco_dpt, q2qe, isData,evt,wgt);
        //utils->fillHistosV3( utils->histos2D[ Form("h_dptx_q2qe_%s", quadrant.c_str())], reco_dptx, q2qe, isData,evt,wgt);
        //utils->fillHistosV3( utils->histos2D[ Form("h_dpty_q2qe_%s", quadrant.c_str())], reco_dpty, q2qe, isData,evt,wgt);
      }

      //FillHDL2D( "h_dthetaRq2qe_ptmu", utils, reco_dthetaR, muon_pt_beam, q2qe, isData, evt, wgt, RunCodeWithSystematics );
      //FillHDL2D( "h_dthetaPq2qe_ptmu", utils, reco_dthetaP, muon_pt_beam, q2qe, isData, evt, wgt , RunCodeWithSystematics);
      //FillHDL2D( "h_tpq2qe_ptmu", utils, protonTn, muon_pt_beam, q2qe, isData, evt, wgt , RunCodeWithSystematics);

      FillHDL2D( "h_dthetaPdthetaR_q2qe", utils, reco_dthetaP, q2qe, reco_dthetaR, isData, evt, wgt, true);//no systematics
      if(hasNC) FillHDL2D( "h_dthetaPdthetaR_q2qe_nc", utils, reco_dthetaP, q2qe, reco_dthetaR, isData, evt, wgt, true);//no systematics
      if( abs(reco_dthetaR)<10 && abs(reco_dthetaP) < 10 )
      {
        FillHDL2D( "h_dthetaPdthetaR_q2qe_fine", utils, reco_dthetaP, q2qe, reco_dthetaR, isData, evt, wgt, true);//no systematics
        if(hasNC) FillHDL2D( "h_dthetaPdthetaR_q2qe_nc_fine", utils, reco_dthetaP, q2qe, reco_dthetaR, isData, evt, wgt, true);//no systematics
      }

      //if( abs(reco_dthetaR) < 10 ) FillHDL2D( "h_dthetaPq2qe_ptmu_center", utils, reco_dthetaP, muon_pt_beam, q2qe, isData, evt, wgt, RunCodeWithSystematics );
      //else FillHDL2D( "h_dthetaPq2qe_ptmu_side", utils, reco_dthetaP, muon_pt_beam, q2qe, isData, evt, wgt, RunCodeWithSystematics );
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
      if( abs(reco_dthetaR) < 10 ) 
      {
        utils->fillHistosV3( utils->histos2D["h_dthetaP_center_q2qe"], reco_dthetaP, q2qe, isData, evt, wgt);
        utils->fillHistosV3( utils->histos2D["h_dthetaPpos_center_q2qe"], abs(reco_dthetaP), q2qe, isData, evt, wgt);
        utils->fillHistosV3( utils->histos2D["h_dthetaP_center_ptmu"], reco_dthetaP, muon_pt_beam, isData, evt, wgt);

        int reco_dthetaP_pos = abs(reco_dthetaP);
        bool angle_01 = (abs(reco_dthetaP) < 1 && abs(reco_dthetaR) < 1 );
        bool angle_02 = (abs(reco_dthetaP) < 2 && abs(reco_dthetaR) < 2 );
        bool angle_03 = (abs(reco_dthetaP) < 3 && abs(reco_dthetaR) < 3 );
        bool angle_04 = (abs(reco_dthetaP) < 4 && abs(reco_dthetaR) < 4 );
        if( angle_01 ) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_01"], q2qe, isData, evt, wgt );
        if( angle_02 ) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_02"], q2qe, isData, evt, wgt );
        if( angle_03 ) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_03"], q2qe, isData, evt, wgt );
        if( angle_04 ) utils->fillHistosV3( utils->histos1D["h_q2qe_angle_04"], q2qe, isData, evt, wgt );
        if(RunCodeWithSystematics && !isData)
        {
          utils->fillVertErrorBands( utils->histos2D["h_dthetaP_center_ptmu"],reco_dthetaP,muon_pt_beam,evt);
          utils->fillVertErrorBands( utils->histos2D["h_dthetaP_center_q2qe"],reco_dthetaP,q2qe,evt);
          utils->fillVertErrorBands( utils->histos2D["h_dthetaPpos_center_q2qe"],abs(reco_dthetaP),q2qe,evt);

          if(angle_01) 
          {
            utils->fillVertErrorBands( utils->histos1D["h_q2qe_angle_01"], q2qe,  evt );
            utils->fillLatErrorBands(  utils->histos1D["h_q2qe_angle_01"],  "q2", q2qe, evt, true );
          }
          if(angle_02) 
          {
            utils->fillVertErrorBands( utils->histos1D["h_q2qe_angle_02"], q2qe,  evt );
            utils->fillLatErrorBands(  utils->histos1D["h_q2qe_angle_02"],  "q2", q2qe, evt, true );
          }
          if(angle_03) 
          {
            utils->fillVertErrorBands( utils->histos1D["h_q2qe_angle_03"], q2qe,  evt );
            utils->fillLatErrorBands(  utils->histos1D["h_q2qe_angle_03"],  "q2", q2qe, evt, true );
          }
          if(angle_04) 
          {
            utils->fillVertErrorBands( utils->histos1D["h_q2qe_angle_04"], q2qe,  evt );
            utils->fillLatErrorBands(  utils->histos1D["h_q2qe_angle_04"],  "q2", q2qe, evt, true );
          }
        }

        if(angle_02)
        {

        }
      }



      //------------------------------------------------------
      //Fill 2D Vertical Error Bands
      //Error Bands are not filled up for tmu_costheta 
      //and pmu_ptmu histos to avoid code slow down 
      //------------------------------------------------------
      if(RunCodeWithSystematics && !isData){

        utils->fillVertErrorBands( utils->histos1D["h_ptmu"],  muon_pt_beam, evt );
        utils->fillVertErrorBands( utils->histos1D["h_q2qe"],  q2qe, evt );
        utils->fillVertErrorBands( utils->histos1D["h_enu"],  enu, evt );
        utils->fillVertErrorBands( utils->histos2D["h_q2qe_ptmu"],  q2qe, muon_pt_beam, evt );
        
        //utils->fillVertErrorBands( utils->histos2D["h_dalphat_ptmu"],reco_dalphat,muon_pt_beam,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_dphit_ptmu"],reco_dphit,muon_pt_beam,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_pn_ptmu"],reco_pn,muon_pt_beam,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_dpt_ptmu"],reco_dpt,muon_pt_beam,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_dptx_ptmu"],reco_dptx,muon_pt_beam,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_dpty_ptmu"],reco_dpty,muon_pt_beam,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_tp_ptmu"],protonTn, muon_pt_beam,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_ptheta_ptmu"],protonAngle, muon_pt_beam,evt);

        //utils->fillVertErrorBands( utils->histos2D["h_dalphat_q2qe"],reco_dalphat,q2qe,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_dphit_q2qe"],reco_dphit,q2qe,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_pn_q2qe"],reco_pn,q2qe,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_dpt_q2qe"],reco_dpt,q2qe,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_dptx_q2qe"],reco_dptx,q2qe,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_dpty_q2qe"],reco_dpty,q2qe,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_tp_q2qe"],protonTn, q2qe,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_ptheta_q2qe"],protonAngle, q2qe,evt);


        //utils->fillVertErrorBands( utils->histos2D["h_sign_ptmu"],reco_sign,muon_pt_beam,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_signdalphat_ptmu"],reco_signed_dalphat,muon_pt_beam,evt);
        //utils->fillVertErrorBands( utils->histos2D["h_signdphit_ptmu"],reco_signed_dphit,muon_pt_beam,evt);

        utils->fillVertErrorBands( utils->histos2D["h_dthetaP_ptmu"],reco_dthetaP,muon_pt_beam,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dthetaR_ptmu"],reco_dthetaR,muon_pt_beam,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPP_ptmu"],reco_dPp,muon_pt_beam,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPR_ptmu"],reco_dPr,muon_pt_beam,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPPi_ptmu"],reco_dPp_infer,muon_pt_beam,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPRi_ptmu"],reco_dPr_infer,muon_pt_beam,evt);

        utils->fillVertErrorBands( utils->histos2D["h_dthetaP_q2qe"],reco_dthetaP,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dthetaR_q2qe"],reco_dthetaR,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPP_q2qe"],reco_dPp,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPR_q2qe"],reco_dPr,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPPi_q2qe"],reco_dPp_infer,q2qe,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dPRi_q2qe"],reco_dPr_infer,q2qe,evt);

        utils->fillVertErrorBands( utils->histos2D["h_dthetaP_dthetaR"],reco_dthetaP,reco_dthetaR,evt);
        utils->fillVertErrorBands( utils->histos2D["h_dthetaP_dthetaR_Fine"],reco_dthetaP,reco_dthetaR,evt);

        utils->fillVertErrorBands( utils->histos2D["h_dtheta2D_ptmu"], dtheta2D, muon_pt_beam, evt );
        utils->fillVertErrorBands( utils->histos2D["h_dtheta2D_q2qe"], dtheta2D, q2qe, evt );
        utils->fillVertErrorBands( utils->histos2D["h_dtheta2D_recoil"], dtheta2D, recoil, evt );
        utils->fillVertErrorBands( utils->histos2D["h_dtheta2D_vtxE"], dtheta2D, vtxE, evt );

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
  NeutronBlobCuts  *ncutter2D = new NeutronBlobCuts(false);


  //--------------------------------------------
  // Define histograms
  //--------------------------------------------

  GeoUtils *geoUtils = new GeoUtils();

  //---------------------------------------------
  // Book Histograms
  //---------------------------------------------
  //

  //1D
  MnvH1D *h_multiplicity[nHistos],*h_nodes[nHistos],*h_nodes_long[nHistos],*h_nodes_anchored[nHistos],*h_nodes_vest[nHistos],*h_ptmu[nHistos];
  MnvH1D *h_q2qe[nHistos], *h_enu[nHistos];
  MnvH1D* h_pmu[nHistos];

  utils->bookHistos( h_multiplicity, "h_multiplicity", "Multiplicity of Outgoing Tracks; Number of Outgoing Tracks; Number of Events", 9, 0., 9.);
  utils->bookHistos( h_nodes, "h_nodes", "nodes",50,0,50);
  utils->bookHistos( h_nodes_long, "h_nodes_long", "nodes",50,0,50);
  utils->bookHistos( h_nodes_anchored, "h_nodes_anchored", "nodes",50,0,50);
  utils->bookHistos( h_nodes_vest, "h_nodes_vest", "nodes",50,0,50);
  utils->bookHistos( h_ptmu, "h_ptmu", Form( "Muon Pt;Reconstructed Muon p_{T} (GeV);Events / %.3f GeV", muonPtbins.default_width ), muonPtbins );
  utils->bookHistos( h_q2qe, "h_q2qe", Form( "Q^{2}_{QE};Reconstructed Q^{2}_{QE} (GeV^{2});Events / %.3f GeV^{2}", Q2bins.default_width ), Q2bins );
  utils->bookHistos( h_enu, "h_enu", Form( "E (GeV);Reconstructed E_{#nu} (GeV);Events / %.3f GeV^{2}", Q2bins.default_width ), Enubins );
  utils->bookHistos( h_pmu, "h_pmu", "Reconstructed p_{#mu} (GeV)", UniformMuonPBins );

  MnvH1D *h_q2qe_angle_01[nHistos];
  MnvH1D *h_q2qe_angle_02[nHistos];
  MnvH1D *h_q2qe_angle_03[nHistos];
  MnvH1D *h_q2qe_angle_04[nHistos];
  utils->bookHistos( h_q2qe_angle_01, "h_q2qe_angle_01", Form( "Q^{2}_{QE};Reconstructed Q^{2}_{QE} (GeV^{2});Events / %.3f GeV^{2}", Q2bins.default_width ), Q2bins );
  utils->bookHistos( h_q2qe_angle_02, "h_q2qe_angle_02", Form( "Q^{2}_{QE};Reconstructed Q^{2}_{QE} (GeV^{2});Events / %.3f GeV^{2}", Q2bins.default_width ), Q2bins );
  utils->bookHistos( h_q2qe_angle_03, "h_q2qe_angle_03", Form( "Q^{2}_{QE};Reconstructed Q^{2}_{QE} (GeV^{2});Events / %.3f GeV^{2}", Q2bins.default_width ), Q2bins );
  utils->bookHistos( h_q2qe_angle_04, "h_q2qe_angle_04", Form( "Q^{2}_{QE};Reconstructed Q^{2}_{QE} (GeV^{2});Events / %.3f GeV^{2}", Q2bins.default_width ), Q2bins );





  //2-D

  MnvH2D *h_q2qe_ptmu[nHistos];
  utils->bookHistos( h_q2qe_ptmu, "h_q2qe_ptmu","q2qe:ptmu",Q2bins,muonPtbins);


  MnvH2D *h_dalphat_ptmu[nHistos],*h_dphit_ptmu[nHistos],*h_pn_ptmu[nHistos],*h_dpt_ptmu[nHistos],*h_dptx_ptmu[nHistos],*h_dpty_ptmu[nHistos],*h_tp_ptmu[nHistos],*h_ptheta_ptmu[nHistos],*h_sign_ptmu[nHistos],*h_signdalphat_ptmu[nHistos],*h_signdphit_ptmu[nHistos];

  MnvH2D *h_dalphat_q2qe[nHistos],*h_dphit_q2qe[nHistos],*h_pn_q2qe[nHistos],*h_dpt_q2qe[nHistos],*h_dptx_q2qe[nHistos],*h_dpty_q2qe[nHistos],*h_tp_q2qe[nHistos],*h_ptheta_q2qe[nHistos],*h_sign_q2qe[nHistos],*h_signdalphat_q2qe[nHistos],*h_signdphit_q2qe[nHistos];

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

  utils->bookHistos( h_dalphat_q2qe, "h_dalphat_q2qe","dat:q2qe",dalphatbins,Q2bins);
  utils->bookHistos( h_dphit_q2qe, "h_dphit_q2qe","dpt:q2qe",dphitbins,Q2bins);
  utils->bookHistos( h_pn_q2qe, "h_pn_q2qe","pn:q2qe",pnbins,Q2bins);
  utils->bookHistos( h_dpt_q2qe, "h_dpt_q2qe","dtp:q2qe",dptbins,Q2bins);
  utils->bookHistos( h_dptx_q2qe, "h_dptx_q2qe","dtpx:q2qe",dptxbins,Q2bins);
  utils->bookHistos( h_dpty_q2qe, "h_dpty_q2qe","dtpy:q2qe",dptybins,Q2bins);
  utils->bookHistos( h_tp_q2qe, "h_tp_q2qe","proton kin:q2qe",protonKineticbins,Q2bins);
  utils->bookHistos( h_ptheta_q2qe, "h_ptheta_q2qe", "protonang:q2qe",protonThetabins,Q2bins);

  //utils->bookHistos( h_sign_ptmu, "h_signed_ptmu","sign:ptmu",signbins,muonPtbins);
  //utils->bookHistos( h_signdalphat_ptmu,"h_signeddalphat_ptmu","sdalphat:ptm",signeddalphatbins,muonPtbins);
  //utils->bookHistos( h_signdphit_ptmu,"h_signeddphit_ptmu","sdphit:ptm",signeddphitbins,muonPtbins);

  MnvH2D *h_dthetaR_ptmu[nHistos], *h_dthetaP_ptmu[nHistos], *h_dPP_ptmu[nHistos], *h_dPR_ptmu[nHistos], *h_dPPi_ptmu[nHistos], *h_dPRi_ptmu[nHistos];
  MnvH2D *h_dthetaP_center_ptmu[nHistos];
  MnvH2D *h_dthetaR_q2qe[nHistos], *h_dthetaP_q2qe[nHistos], *h_dPP_q2qe[nHistos], *h_dPR_q2qe[nHistos], *h_dPPi_q2qe[nHistos], *h_dPRi_q2qe[nHistos];
  MnvH2D *h_dthetaP_center_q2qe[nHistos];
  MnvH2D *h_dthetaR[nHistos], *h_dthetaP[nHistos];

  MnvH2D *h_dthetaP_dthetaR[nHistos], *h_dthetaP_dthetaR_Fine[nHistos];

  MnvH2D *h_dthetaPpos_center_q2qe[nHistos];



  utils->bookHistos( h_dthetaR_ptmu, "h_dthetaR_ptmu", "dthetaR:ptmu", dthetaReactbins, muonPtbins);
  utils->bookHistos( h_dthetaP_ptmu, "h_dthetaP_ptmu", "dthetaP:ptmu", dthetaPerpbins, muonPtbins);
  utils->bookHistos( h_dthetaP_center_ptmu, "h_dthetaP_center_ptmu", "dthetaPcenter:ptmu", dthetaPerpbins, muonPtbins);
  utils->bookHistos( h_dPR_ptmu, "h_dPR_ptmu", "dPR:ptmu", dptybins, muonPtbins);
  utils->bookHistos( h_dPP_ptmu, "h_dPP_ptmu", "dPP:ptmu", dptxbins, muonPtbins);
  utils->bookHistos( h_dPRi_ptmu, "h_dPRi_ptmu", "dPRi:ptmu", dptybins, muonPtbins);
  utils->bookHistos( h_dPPi_ptmu, "h_dPPi_ptmu", "dPPi:ptmu", dptxbins, muonPtbins);

  utils->bookHistos( h_dthetaR_q2qe, "h_dthetaR_q2qe", "dthetaR:q2qe", dthetaReactbins, Q2bins);
  utils->bookHistos( h_dthetaP_q2qe, "h_dthetaP_q2qe", "dthetaP:q2qe", dthetaPerpbins, Q2bins);
  utils->bookHistos( h_dthetaP_center_q2qe, "h_dthetaP_center_q2qe", "dthetaPcenter:q2qe", dthetaPerpbins, Q2bins);
  utils->bookHistos( h_dPR_q2qe, "h_dPR_q2qe", "dPR:q2qe", dptybins, Q2bins);
  utils->bookHistos( h_dPP_q2qe, "h_dPP_q2qe", "dPP:q2qe", dptxbins, Q2bins);
  utils->bookHistos( h_dPRi_q2qe, "h_dPRi_q2qe", "dPRi:q2qe", dptybins, Q2bins);
  utils->bookHistos( h_dPPi_q2qe, "h_dPPi_q2qe", "dPPi:q2qe", dptxbins, Q2bins);

  utils->bookHistos( h_dthetaP_dthetaR, "h_dthetaP_dthetaR", "dthetaP: dthetaR", dthetaPerpbins, dthetaReactbins );
  utils->bookHistos( h_dthetaP_dthetaR_Fine, "h_dthetaP_dthetaR_Fine", "dthetaP: dthetaR", dthetaFineBins, dthetaFineBins);
  utils->bookHistos( h_dthetaPpos_center_q2qe, "h_dthetaPpos_center_q2qe", "dthetaP: dthetaR", dthetaPerpbinsPositive, Q2bins );

  utils->bookHistos( h_node_nodeenergy, "h_node_nodeenergy","node:energy",50,0,50,500,0,50);
  utils->bookHistos( h_node_nodeenergy_long, "h_node_nodeenergy_long","node:energy",50,0,50,500,0,50);
  utils->bookHistos( h_node_nodeenergy_anchored, "h_node_nodeenergy_anchored","node:energy",50,0,50,500,0,50);
  utils->bookHistos( h_node_nodeenergy_vest, "h_node_nodeenergy_vest","node:energy",50,0,50,500,0,50);

  utils->bookHistos( h_protonRes_nodeenergy_node01, "h_protonRes_nodeenergy_node01","res:node",500,0,500,100,-1,1);
  utils->bookHistos( h_protonRes_nodeenergy_node2, "h_protonRes_nodeenergy_node2","res:node",500,0,500,100,-1,1);
  utils->bookHistos( h_protonRes_nodeenergy_node3, "h_protonRes_nodeenergy_node3","res:node",500,0,500,100,-1,1);
  utils->bookHistos( h_protonRes_nodeenergy_node4, "h_protonRes_nodeenergy_node4","res:node",500,0,500,100,-1,1);
  utils->bookHistos( h_protonRes_nodeenergy_node5, "h_protonRes_nodeenergy_node5","res:node",500,0,500,100,-1,1);


  MnvH2D *h_dtheta2D_q2qe[nHistos];
  MnvH2D *h_dtheta2D_ptmu[nHistos];
  MnvH2D *h_dtheta2D_recoil[nHistos];
  MnvH2D *h_dtheta2D_vtxE[nHistos];
  utils->bookHistos( h_dtheta2D_q2qe, "h_dtheta2D_q2qe", "dtheta2D:q2qe", UniformDThetaBins, Q2bins );
  utils->bookHistos( h_dtheta2D_ptmu, "h_dtheta2D_ptmu", "dtheta2D:ptmu", UniformDThetaBins, muonPtbins );
  utils->bookHistos( h_dtheta2D_recoil, "h_dtheta2D_recoil", "dtheta2D:recoil", UniformDThetaBins, UniformVisibleEBins);
  utils->bookHistos( h_dtheta2D_vtxE, "h_dtheta2D_vtxE", "dtheta2D:vtxE", UniformDThetaBins, UniformVtxEBins );





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
  MnvH2D *h_true_vs_reco_neutronE[nHistos];
  utils->bookHistos( h_true_vs_reco_neutronE, "h_true_vs_reco_neutronE", "true:reco_neutE",200,0,1.2,200,0,1.2 );

  MnvH2D *h_true_vs_reco_dthetaP[nHistos], *h_true_vs_reco_dthetaR[nHistos];
  utils->bookHistos( h_true_vs_reco_dthetaP, "h_true_vs_reco_dthetaP", "true:reco_dthetaP", 72,-180,180, 72,-180,180);
  utils->bookHistos( h_true_vs_reco_dthetaR, "h_true_vs_reco_dthetaR", "true:reco_dthetaR", 72,-180,180, 72,-180,180);



  //axis_binning dthetaRq2qeBins, //Global X bins
  //             dthetaPq2qeBins, 
  //             tpq2qeBins, 
  axis_binning   dthetaPdthetaRBins;
  axis_binning   dthetaPdthetaRFineBins;
  //MnvH2D *h_dthetaRq2qe_ptmu[nHistos], 
  //       *h_dthetaPq2qe_ptmu[nHistos],
  //       *h_tpq2qe_ptmu[nHistos],
  MnvH2D   *h_dthetaPdthetaR_q2qe[nHistos];
  MnvH2D   *h_dthetaPdthetaR_q2qe_fine[nHistos];
  MnvH2D   *h_dthetaPdthetaR_q2qe_nc[nHistos];
  MnvH2D   *h_dthetaPdthetaR_q2qe_nc_fine[nHistos];



  //HyperDimLinearizer *hdl_dthetaRq2qe_ptmu = DeclareHDLHistos2D( utils, h_dthetaRq2qe_ptmu, "h_dthetaRq2qe_ptmu", "pt vs dthetaRq2qe: dthetaR : q2qe", dthetaReactbins, muonPtbins, Q2bins, dthetaRq2qeBins );
  //HyperDimLinearizer *hdl_dthetaPq2qe_ptmu = DeclareHDLHistos2D( utils, h_dthetaPq2qe_ptmu, "h_dthetaPq2qe_ptmu", "pt vs dthetaPq2qe: dthetaP : q2qe", dthetaPerpbins, muonPtbins, Q2bins, dthetaPq2qeBins, true );
  ////HyperDimLinearizer *hdl_tpq2qe_ptmu = DeclareHDLHistos2D( utils, h_tpq2qe_ptmu, "h_tpq2qe_ptmu", "pt vs tpq2qe: tp : q2qe", protonKineticbins, muonPtbins, Q2bins, tpq2qeBins );
  HyperDimLinearizer *hdl_dthetaPdthetaR_q2qe = DeclareHDLHistos2D( utils, h_dthetaPdthetaR_q2qe, "h_dthetaPdthetaR_q2qe", "q2qe vs dthetaPdthetaR: dthetaP : dthetaR", dthetaPerpbins, Q2bins, dthetaReactbins, dthetaPdthetaRBins );
  HyperDimLinearizer *hdl_dthetaPdthetaR_q2qe_nc = DeclareHDLHistos2D( utils, h_dthetaPdthetaR_q2qe_nc, "h_dthetaPdthetaR_q2qe_nc", "q2qe vs dthetaPdthetaR with nc: dthetaP : dthetaR", dthetaPerpbins, Q2bins, dthetaReactbins, dthetaPdthetaRBins );
  HyperDimLinearizer *hdl_dthetaPdthetaR_q2qe_fine = DeclareHDLHistos2D( utils, h_dthetaPdthetaR_q2qe_fine, "h_dthetaPdthetaR_q2qe_fine", "q2qe vs dthetaPdthetaR: dthetaP : dthetaR", dthetaFineBins, Q2bins, dthetaFineBins, dthetaPdthetaRFineBins);
  HyperDimLinearizer *hdl_dthetaPdthetaR_q2qe_nc_fine = DeclareHDLHistos2D( utils, h_dthetaPdthetaR_q2qe_nc_fine, "h_dthetaPdthetaR_q2qe_nc_fine", "q2qe vs dthetaPdthetaR with nc: dthetaP : dthetaR", dthetaFineBins, Q2bins, dthetaFineBins, dthetaPdthetaRFineBins );




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

  //MnvH2D *h_dthetaPq2qe_ptmu_center[nHistos],*h_dthetaPq2qe_ptmu_side[nHistos];//test dthetaR cut
  //DeclareHDLHistos2D( utils, h_dthetaPq2qe_ptmu_center, "h_dthetaPq2qe_ptmu_center", "pt vs dthetaPq2qe: dthetaP : q2qe", dthetaPerpbins, muonPtbins, Q2bins, dthetaPq2qeBins, true );
  //DeclareHDLHistos2D( utils, h_dthetaPq2qe_ptmu_side, "h_dthetaPq2qe_ptmu_side", "pt vs dthetaPq2qe: dthetaP : q2qe", dthetaPerpbins, muonPtbins, Q2bins, dthetaPq2qeBins, true );

  MnvH2D *h_pRes_vs_node0[nHistos],*h_pRes_vs_node1[nHistos],*h_pRes_vs_node2[nHistos],*h_pRes_vs_node3[nHistos],*h_pRes_vs_nodeFirst[nHistos];
  utils->bookHistos( h_pRes_vs_node0,"h_pRes_vs_node0", "res:nodeE", 50,-1,1,10,0,30 );
  utils->bookHistos( h_pRes_vs_node1,"h_pRes_vs_node1", "res:nodeE", 50,-1,1,10,0,30 );
  utils->bookHistos( h_pRes_vs_node2,"h_pRes_vs_node2", "res:nodeE", 50,-1,1,10,0,30 );
  utils->bookHistos( h_pRes_vs_node3,"h_pRes_vs_node3", "res:nodeE", 50,-1,1,10,0,30 );
  utils->bookHistos( h_pRes_vs_nodeFirst,"h_pRes_vs_nodeFirst", "res:nodeE", 50,-1,1,10,0,30 );


  //Q2 tests: in angle_02
  MnvH2D *h_blobDist_q2qe[nHistos], *h_blobEnergy_q2qe[nHistos], *h_nBlobs_q2qe[nHistos],*h_n3DBlobs_q2qe[nHistos], *h_n2DBlobs_q2qe[nHistos], *h_muonTheta_q2qe[nHistos], *h_muonPhi_q2qe[nHistos], *h_muonE_q2qe[nHistos], *h_protonDEDX_q2qe[nHistos],*h_pionDEDX_q2qe[nHistos], *h_nonBlobEnergy_q2qe[nHistos], *h_hasTrack_q2qe[nHistos];

  MnvH2D *h_blobMaxE_q2qe[nHistos];
  MnvH2D *h_nClus_q2qe[nHistos];
  utils->bookHistos( h_blobDist_q2qe, "h_blobDist_q2qe", "dist:q2qe", BlobDistBins, Q2bins);
  utils->bookHistos( h_blobEnergy_q2qe, "h_blobEnergy_q2qe", "blobenergy:q2qe", BlobEnergyBins, Q2bins);
  utils->bookHistos( h_nonBlobEnergy_q2qe, "h_nonBlobEnergy_q2qe", "nonblobenergy:q2qe", BlobEnergyBins, Q2bins);
  utils->bookHistos( h_n3DBlobs_q2qe, "h_n3DBlobs_q2qe", "n3dblobs:q2qe", BlobNumBins, Q2bins );
  utils->bookHistos( h_n2DBlobs_q2qe, "h_n2DBlobs_q2qe", "n2dblobs:q2qe", BlobNumBins, Q2bins );
  utils->bookHistos( h_nBlobs_q2qe, "h_nBlobs_q2qe", "n2dblobs:q2qe", BlobNumBins, Q2bins );
  utils->bookHistos( h_muonTheta_q2qe, "h_muonTheta_q2qe", "muonTheta:q2qe", MuonThetaBins,Q2bins);
  utils->bookHistos( h_muonPhi_q2qe, "h_muonPhi_q2qe", "muonPhi:q2qe",MuonPhiBins, Q2bins );
  utils->bookHistos( h_muonE_q2qe, "h_muonE_q2qe", "muonE:q2qe",MuonEnergyBins, Q2bins );
  utils->bookHistos( h_protonDEDX_q2qe, "h_protonDEDX_q2qe", "protonDEDX:q2qe",10,-1,1,6,0.1,0.4);
  utils->bookHistos( h_pionDEDX_q2qe, "h_pionDEDX_q2qe", "pionDEDX:q2qe",10,-1,1,6,0.1,0.4);
  utils->bookHistos( h_hasTrack_q2qe, "h_hasTrack_q2qe", "hasTrack:q2qe",2,-.5,1.5,6,0.1,0.4);
  utils->bookHistos( h_blobMaxE_q2qe, "h_blobMaxE_q2qe", "blobMaxE:q2qe",BlobEnergyBins, Q2bins);
  utils->bookHistos( h_nClus_q2qe, "h_nClus_q2qe", "NClusters:q2qe",nClusBins, Q2bins);


  MnvH2D *h_nHits_totalE[nHistos];
  utils->bookHistos( h_nHits_totalE, "h_nHits_totalE", "nHits:histsE", nHitsBins, hitsEBins );


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

      utils->addVertErrorBands( h_q2qe_angle_01 );
      utils->addLatErrorBands ( h_q2qe_angle_01 );
      utils->addVertErrorBands( h_q2qe_angle_02 );
      utils->addLatErrorBands ( h_q2qe_angle_02 );
      utils->addVertErrorBands( h_q2qe_angle_03 );
      utils->addLatErrorBands ( h_q2qe_angle_03 );
      utils->addVertErrorBands( h_q2qe_angle_04 );
      utils->addLatErrorBands ( h_q2qe_angle_04 );


      utils->addVertErrorBands( h_q2qe_ptmu );
      utils->addLatErrorBands ( h_q2qe_ptmu );
      //2-D
      utils->addVertErrorBands( h_dalphat_ptmu );
      utils->addVertErrorBands( h_dphit_ptmu );
      utils->addVertErrorBands( h_pn_ptmu );
      utils->addVertErrorBands( h_dpt_ptmu );
      utils->addVertErrorBands( h_dptx_ptmu );
      utils->addVertErrorBands( h_dpty_ptmu );
      utils->addVertErrorBands( h_tp_ptmu );
      utils->addVertErrorBands( h_ptheta_ptmu );

      utils->addVertErrorBands( h_dalphat_q2qe );
      utils->addVertErrorBands( h_dphit_q2qe );
      utils->addVertErrorBands( h_pn_q2qe );
      utils->addVertErrorBands( h_dpt_q2qe );
      utils->addVertErrorBands( h_dptx_q2qe );
      utils->addVertErrorBands( h_dpty_q2qe );
      utils->addVertErrorBands( h_tp_q2qe );
      utils->addVertErrorBands( h_ptheta_q2qe );


      ////hutils->addVertErrorBands( h_sign_ptmu );
      ////hutils->addVertErrorBands( h_signdalphat_ptmu );
      ////hutils->addVertErrorBands( h_signdphit_ptmu );
      utils->addVertErrorBands( h_dthetaR_ptmu );
      utils->addVertErrorBands( h_dthetaP_ptmu );
      utils->addVertErrorBands( h_dthetaP_center_ptmu );
      utils->addVertErrorBands( h_dPP_ptmu );
      utils->addVertErrorBands( h_dPR_ptmu );
      utils->addVertErrorBands( h_dPPi_ptmu );
      utils->addVertErrorBands( h_dPRi_ptmu );
      utils->addVertErrorBands( h_dtheta2D_ptmu );

      utils->addVertErrorBands( h_dthetaR_q2qe );
      utils->addVertErrorBands( h_dthetaP_q2qe );
      utils->addVertErrorBands( h_dthetaP_center_q2qe );
      utils->addVertErrorBands( h_dPP_q2qe );
      utils->addVertErrorBands( h_dPR_q2qe );
      utils->addVertErrorBands( h_dPPi_q2qe );
      utils->addVertErrorBands( h_dPRi_q2qe );
      utils->addVertErrorBands( h_dtheta2D_q2qe );

      utils->addVertErrorBands( h_dthetaP_dthetaR );
      utils->addVertErrorBands( h_dthetaP_dthetaR_Fine );
      utils->addVertErrorBands( h_dthetaPpos_center_q2qe );


      utils->addVertErrorBands( h_dtheta2D_recoil );
      utils->addVertErrorBands( h_dtheta2D_vtxE );
      
      utils->addLatErrorBands( h_dalphat_ptmu );
      utils->addLatErrorBands( h_dphit_ptmu );
      utils->addLatErrorBands( h_pn_ptmu );
      utils->addLatErrorBands( h_dpt_ptmu );
      utils->addLatErrorBands( h_dptx_ptmu );
      utils->addLatErrorBands( h_dpty_ptmu );
      utils->addLatErrorBands( h_tp_ptmu );
      utils->addLatErrorBands( h_ptheta_ptmu );

      utils->addLatErrorBands( h_dalphat_q2qe );
      utils->addLatErrorBands( h_dphit_q2qe );
      utils->addLatErrorBands( h_pn_q2qe );
      utils->addLatErrorBands( h_dpt_q2qe );
      utils->addLatErrorBands( h_dptx_q2qe );
      utils->addLatErrorBands( h_dpty_q2qe );
      utils->addLatErrorBands( h_tp_q2qe );
      utils->addLatErrorBands( h_ptheta_q2qe );



      ////hutils->addLatErrorBands( h_sign_ptmu );
      ////hutils->addLatErrorBands( h_signdalphat_ptmu );
      ////hutils->addLatErrorBands( h_signdphit_ptmu );
      utils->addLatErrorBands( h_dthetaR_ptmu );
      utils->addLatErrorBands( h_dthetaP_ptmu );
      utils->addLatErrorBands( h_dthetaP_center_ptmu );
      utils->addLatErrorBands( h_dPP_ptmu );
      utils->addLatErrorBands( h_dPR_ptmu );
      utils->addLatErrorBands( h_dPPi_ptmu );
      utils->addLatErrorBands( h_dPRi_ptmu );
      utils->addLatErrorBands( h_dtheta2D_ptmu );

      utils->addLatErrorBands( h_dthetaR_q2qe );
      utils->addLatErrorBands( h_dthetaP_q2qe );
      utils->addLatErrorBands( h_dthetaP_center_q2qe );
      utils->addLatErrorBands( h_dPP_q2qe );
      utils->addLatErrorBands( h_dPR_q2qe );
      utils->addLatErrorBands( h_dPPi_q2qe );
      utils->addLatErrorBands( h_dPRi_q2qe );
      utils->addLatErrorBands( h_dtheta2D_q2qe );

      utils->addLatErrorBands( h_dthetaP_dthetaR );
      utils->addLatErrorBands( h_dthetaP_dthetaR_Fine );
      utils->addLatErrorBands( h_dthetaPpos_center_q2qe );


      utils->addLatErrorBands( h_dtheta2D_recoil );
      utils->addLatErrorBands( h_dtheta2D_vtxE );
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
  mc_events=AnalysisLooper( mc, utils, cutter, ncutter2D, sample,  multiplicity, entries_mc, mc_counter, isData, isTruth, npMode );
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
    data_events = AnalysisLooper( data, utils, cutter, ncutter2D, sample,  multiplicity, entries_data, data_counter, isData, isTruth, npMode );
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
  par.push_back("1");//proton, neutron,proton/neutron mode

  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }

  bool fluxConstraintHistoNeeded = ( par.at(4) == "1" ) ? true : false; 

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;

//int MuonSelectionHists(string filename, string sample, bool makeFluxConstraintHisto, int multiplicity, string playlist, int n_mcfiles = -1, int n_datafiles = -1, int npMode = 1)
  return MuonSelectionHists(par[1], par[3], fluxConstraintHistoNeeded, atoi(par[5].c_str()), par[2],atoi(par[6].c_str()), atoi(par[7].c_str()), atoi(par[8].c_str()) );
}
