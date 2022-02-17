#include "include/CCQENuUtils.h"
#include "include/ComputeUtils.h"
#include "TParameter.h" 
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"

using namespace CCQENU_ANA;
using namespace Neutron_ANA;
using namespace ROOT::Math;
using namespace ROOT;

bool DEBUG=false;
double MuonCylinderCut = 100;//mm
double NeutronZCut = 100;//mm
double MuonConeAngleDeg = 15;//degree

DetectorUtils* detUtils = new DetectorUtils();

bool CommonCut( CCQENuEvent *evt, CCQENuCuts *cutter, int multiplicity = 0, bool isSignal = false )
{
  if( !cutter->passInteractionVertex( evt ) )       return false;
  if( !cutter->passDeadTimeCuts( evt ) )            return false;
  if( !cutter->passNuHelicityCut( evt ) )           return false;  
  //if( !cutter->passExtraTracksProtonsCut ( evt ) )  return false;
  if( !cutter->passTrackAngleCut( evt ) )           return false;
  //Fiducial Cuts
  double x = evt->vtx[0];
  double y = evt->vtx[1];
  double z = evt->vtx[2];
  if( z < 5980 || z > 8422 )                        return false;
  if( !detUtils->IsInHexagon( x, y, 850 ) )         return false;

  if (evt->multiplicity != multiplicity)            return false;
  return true;
}

bool CommonNeutronCut( CCQENuEvent* evt, NeutronBlobCuts* ncutter )
{
  
  bool distToMuonCut = ncutter->DistanceToMuonCut( evt, MuonCylinderCut );
  bool cutForwardNeutron = true;
  bool forwardCandidateCut = ncutter->ForwardNeutronCut( evt, ncutter->MainNeutronCandidate(), cutForwardNeutron,NeutronZCut );
  bool muonConeCut = ncutter->MuonConeCut( evt, ncutter->MainNeutronCandidate(), MuonConeAngleDeg, true, 0 );
  bool energyExist = ncutter->MainNeutronCandidate()->TotalE > 0;
  return distToMuonCut*forwardCandidateCut*muonConeCut*energyExist;
}

bool EventSelectionDataMC( CCQENuEvent *evt, CCQENuCuts *cutter, NeutronBlobCuts* ncutter, bool isSignal = true, bool ISMC=true, int multiplicity = 0 , bool pass_michel = true, bool pass_nblobs = true, bool pass_neutronMainCandidate = true, bool pass_MainCandidateCut = true) //all true --> signal
{
  //Repopulate the neutron candidates class first
  ncutter->UpdateCandidates( evt, ISMC );
  //mic side band: michel=true, nblobs = false
  if (!CommonCut(evt, cutter, multiplicity, isSignal ) ) return false; //no longer has signal function
  // neutron blobbing and cuts
  if(DEBUG) cout<<"Start main candidate cut"<<endl;
  bool hasMainCandidate = ncutter->GetMainCandidate( evt );
  if (!hasMainCandidate) return false;
  if (!CommonNeutronCut( evt, ncutter ) ) return false;
  return true;
}

double AnalysisLooper( CCQENuEvent* evt, CCQENuUtils* utils, CCQENuCuts *cutter, NeutronBlobCuts* ncutter, string sample, int multiplicity, int n_entries, bool isMC = true )
{
  cout<<"Start AnalysisLooper"<<endl;
  GeoUtils * geoUtils = ncutter->GeoUtilsPtr();
  PhysicsUtils * physicsUtils = new PhysicsUtils();
  cout<<"Got GeoUtils*"<<endl;
  double n_outputevent = 0;
  for( int i = 0; i<n_entries; ++ i )
  {
    if( (i)%100000 == 0 ) 
    {
      double percent = 100*(i+1)/n_entries;
      if( isMC )  cout << "Reading MC entry : " ;
      else cout<< "Reading Data entry : " ;
      cout<< i+1 << " - " << 100*(i+1)/n_entries << " % " << endl; 
    }
    evt->GetEntry(i);

    // MC1 continue if outside mc weight range
    if (isMC) 
    {
      
      if(evt->truth_genie_wgt_MFP_N[1]>50. || evt->truth_genie_wgt_MFP_N[2]>50. || evt->truth_genie_wgt_MFP_N[3]>50. || 
      evt->truth_genie_wgt_MFP_N[4]>50. || evt->truth_genie_wgt_MFP_N[5]>50. || evt->truth_genie_wgt_MFP_N[6]>50. || evt->truth_genie_wgt_MFP_N[7]>50.) 
      {
        cout << "Event ignored - truth_genie_wgt_MFP_N[1] == " << evt->truth_genie_wgt_MFP_N[1]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_N[2] == " << evt->truth_genie_wgt_MFP_N[2]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_N[3] == " << evt->truth_genie_wgt_MFP_N[3]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_N[4] == " << evt->truth_genie_wgt_MFP_N[4]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_N[5] == " << evt->truth_genie_wgt_MFP_N[5]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_N[6] == " << evt->truth_genie_wgt_MFP_N[6]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_N[7] == " << evt->truth_genie_wgt_MFP_N[7]<< endl;
        continue;
      }
      if(evt->truth_genie_wgt_MFP_pi[1]>50. || evt->truth_genie_wgt_MFP_pi[2]>50. || evt->truth_genie_wgt_MFP_pi[3]>50. || evt->truth_genie_wgt_MFP_pi[4]>50. 
      || evt->truth_genie_wgt_MFP_pi[5]>50. || evt->truth_genie_wgt_MFP_pi[6]>50. || evt->truth_genie_wgt_MFP_pi[7]>50.) 
      {
        cout << "Event ignored - truth_genie_wgt_MFP_pi[1] == " << evt->truth_genie_wgt_MFP_pi[1]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_pi[2] == " << evt->truth_genie_wgt_MFP_pi[2]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_pi[3] == " << evt->truth_genie_wgt_MFP_pi[3]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_pi[4] == " << evt->truth_genie_wgt_MFP_pi[4]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_pi[5] == " << evt->truth_genie_wgt_MFP_pi[5]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_pi[6] == " << evt->truth_genie_wgt_MFP_pi[6]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_pi[7] == " << evt->truth_genie_wgt_MFP_pi[7]<< endl;
        continue;
      }
    } //end MC1
    
    // event selection
    if (sample == "Signal" )
    {
      if( ! EventSelectionDataMC( evt, cutter, ncutter, true,isMC, multiplicity ) ) continue;
    }
    else if (sample == "MichelSideBand")
    {
      if( ! EventSelectionDataMC( evt, cutter, ncutter, false,isMC, multiplicity, false,true, true,true ) ) continue;
    }
    else if (sample == "BlobSideBand" )
    {
      if( ! EventSelectionDataMC( evt, cutter, ncutter, false,isMC, multiplicity, true,false, true,true ) ) continue;
    }
    else if (sample == "MicBlobSideBand")
    {
      if( ! EventSelectionDataMC( evt, cutter, ncutter, false,isMC, multiplicity, false,false, true,true ) ) continue;
    }
    else if( sample == "SideBand" )
    {
      if (multiplicity==0)
      {
        if( !cutter->passCCQESelection(evt,"Sideband") )    continue;
      }
      else
      {
        if(!EventSelectionDataMC(evt, cutter, ncutter, true, isMC, multiplicity) &&!cutter->passSidebandFunction( evt, 0, 0 )) continue;
      }
    }
    else {
      cout<<"Wrong Sample Selection."<<endl;
      cout<<"Select: Signal, MichelSideBand, SideBand"<<endl;
      return 0;
    }
    // end event selection
    //bool EventSelectionDataMC( CCQENuEvent *evt, CCQENuCuts *cutter, NeutronBlobCuts* ncutter,bool IsSignal, bool ISMC = true, int multiplicity = 0 , bool pass_michel = true, bool pass_nblobs = true, bool pass_neutronMainCandidate = true, bool pass_MainCandidateCut = true) //all true --> signal
    
    double wgt = (isMC)? utils->GetCVWeight( evt , sample): 1;    
    n_outputevent+=wgt;
    //-----------
    // This is the old method. Because of angular bias introduced from hadronic energy we need to pick out the right angle from a downstream node
    // New method until general reco is fixed is to get the reco P, a downstream theta (already corrected to the beam direction), and calculate
    // the pt and pz from these two quantities
    //-----------
 //double muon_theta     =  evt->CCQENu_muon_theta;
    double muon_theta     = evt->muon_theta_allNodes[20]; // get theta from 20th node
    double cos_muon_theta =  cos(muon_theta);
    double muon_T         =  evt->CCQENu_muon_T / pow(10,3); //GeV
    
    double reco_muon_px   =  evt->CCQENu_leptonE[0];
    double reco_muon_py   =  evt->CCQENu_leptonE[1];
    double reco_muon_pz   =  evt->CCQENu_leptonE[2];
      /*  
    double muon_pz_beam   = utils->GetLongitudinalMomentumWRTBeam( reco_muon_px, reco_muon_py, reco_muon_pz ) / pow(10,3); //GeV
    */
    double muon_p         = utils->GetTotalMomentum( reco_muon_px, reco_muon_py, reco_muon_pz ) / pow(10,3); //GeV

    //New method (maybe temporary)
    double muon_pt_beam = sin(muon_theta)*muon_p;
    //double muon_pt_beam   = utils->GetTransverseMomentumWRTBeam( reco_muon_px, reco_muon_py, reco_muon_pz ) / pow(10,3); //GeV
    double muon_pz_beam = cos(muon_theta)*muon_p;
    //double muon_pz_beam = reco_muon_pz;
    
    //Convert to degrees
    muon_theta*= 180. / 3.14159;
    

    // Uncomment out if you need to make trkChi2 plot   
    // double mnv_trk_chi2PerDoF  = evt->muon_minerva_trk_chi2PerDoF; 
    
    //Change muon momentum to GeV
    reco_muon_px /= pow(10,3);
    reco_muon_py /= pow(10,3);
    reco_muon_pz /= pow(10,3);


    //reconstructed invariant mass
    double q2_gen = 2.0*(evt->CCQENu_leptonE[3]+evt->recoil_energy_nonmuon_nonvtx0mm-muon_p*1000.0*cos(muon_theta))-(105*105);
    double reco_w = (938.*938.)+(2*938.*evt->recoil_energy_nonmuon_nonvtx0mm)-q2_gen;
    reco_w = (reco_w<0) ? 0.0 : sqrt(reco_w);

    //Cone Neutron energy
    XYZVector muonDirWRTDetector = XYZVector( reco_muon_px, reco_muon_py, reco_muon_pz ).Unit();
    XYZVector vtx(evt->vtx[0], evt->vtx[1], evt->vtx[2]);
    XYZVector blobVtxDirTemp = ncutter->MainNeutronCandidate()->BlobVtxDir.Unit();
    Cone muonCone( vtx, muonDirWRTDetector, 15./180*3.14);


    Cone neutronCone05( vtx, blobVtxDirTemp, 5./180*3.14 );
    Cone neutronCone10( vtx, blobVtxDirTemp, 10./180*3.14 );
    Cone neutronCone15( vtx, blobVtxDirTemp, 15./180*3.14 );
    Cone neutronCone20( vtx, blobVtxDirTemp, 20./180*3.14 );
    Cone neutronCone30( vtx, blobVtxDirTemp, 30./180*3.14 );
    Cone neutronCone45( vtx, blobVtxDirTemp, 45./180*3.14 );
    double neutron_cone_E_15 = 0, neutron_cone_E_20 = 0, neutron_cone_E_05 = 0, neutron_cone_E_10 = 0, neutron_cone_E_30 = 0, neutron_cone_E_45=0;
    double neutron_cone2D_E_15 = 0, neutron_cone2D_E_20 = 0, neutron_cone2D_E_05 = 0, neutron_cone2D_E_10 = 0, neutron_cone2D_E_30 = 0, neutron_cone2D_E_45=0;
    double neutron_cone3D_E_15 = 0, neutron_cone3D_E_20 = 0, neutron_cone3D_E_05 = 0, neutron_cone3D_E_10 = 0, neutron_cone3D_E_30 = 0, neutron_cone3D_E_45=0;
    double muon_E=0, muon_E2D = 0, muon_E3D=0;

    vector<NeutronBlob*> all_blobs = ncutter->GetCandidates()->AllBlobs;
    for ( vector<NeutronBlob*>::iterator itBlob = all_blobs.begin(); itBlob != all_blobs.end(); ++itBlob )
    {
      XYZVector blobPos( (*itBlob)->BlobBegPos.X(),(*itBlob)->BlobBegPos.Y(),(*itBlob)->BlobBegPos.Z());
      double E = (*itBlob)->TotalE;
      if (E==0) continue;
      int view = (*itBlob)->View;
      if (view >3 || view < 1) continue;
      if ((*itBlob)->Is3D)
      { 
        if (muonCone.InsideConeAbsPos( blobPos )) { muon_E+=E; muon_E3D+=E;continue; }
        if (neutronCone05.InsideConeAbsPos( blobPos )) {neutron_cone3D_E_05+=E; neutron_cone_E_05+=E;}
        if (neutronCone10.InsideConeAbsPos( blobPos )) {neutron_cone3D_E_10+=E; neutron_cone_E_10+=E;}
        if (neutronCone15.InsideConeAbsPos( blobPos )) {neutron_cone3D_E_15+=E; neutron_cone_E_15+=E;}
        if (neutronCone20.InsideConeAbsPos( blobPos )) {neutron_cone3D_E_20+=E; neutron_cone_E_20+=E;}
        if (neutronCone30.InsideConeAbsPos( blobPos )) {neutron_cone3D_E_30+=E; neutron_cone_E_30+=E;}
        if (neutronCone45.InsideConeAbsPos( blobPos )) {neutron_cone3D_E_45+=E; neutron_cone_E_45+=E;}
      }
      else
      {
        double tpos = (*itBlob)->BlobPosT.X(), zpos = (*itBlob)->BlobPosT.Z() ;
        if (muonCone.InsideConeAbsPos( view, tpos, zpos )) { muon_E+=E; muon_E2D+=E;continue; }
        if (neutronCone05.InsideConeAbsPos( view, tpos, zpos  )) {neutron_cone2D_E_05+=E; neutron_cone_E_05+=E;}
        if (neutronCone10.InsideConeAbsPos( view, tpos, zpos  )) {neutron_cone2D_E_10+=E; neutron_cone_E_10+=E;}
        if (neutronCone15.InsideConeAbsPos( view, tpos, zpos  )) {neutron_cone2D_E_15+=E; neutron_cone_E_15+=E;}
        if (neutronCone20.InsideConeAbsPos( view, tpos, zpos  )) {neutron_cone2D_E_20+=E; neutron_cone_E_20+=E;}
        if (neutronCone30.InsideConeAbsPos( view, tpos, zpos  )) {neutron_cone2D_E_30+=E; neutron_cone_E_30+=E;}
        if (neutronCone45.InsideConeAbsPos( view, tpos, zpos  )) {neutron_cone2D_E_45+=E; neutron_cone_E_45+=E;}
      }
    }

    //recoil variants MeV
    double recoil =ncutter->GetDefaultRecoilEnergy(evt);
    double recoil_inc = evt->recoil_summed_energy[0];
    //double recoil_NoMC = recoil -  ncutter->MainNeutronCandidate()->TotalE;
    double recoil_NoMC = recoil -  neutron_cone_E_30;
    double recoil_inc_NoMC = recoil_inc - ncutter->MainNeutronCandidate()->TotalE;
    double recoil_neutron_3D = ncutter->GetNCRecoilEnergy(3);
    double recoil_neutron_2D = ncutter->GetNCRecoilEnergy(2);
    double recoil_neutron_All = ncutter->GetNCRecoilEnergy(1);

    double recoil_all_nonMuon = evt->recoil_energy_nonmuon_nonvtx100mm;

    

    //vtx energy
    double vtx_energy = evt->vtx_blobs_energy;

    //

    //vtx energy fraction
    double vtx_energy_fraction = vtx_energy/(vtx_energy+recoil);
    if(vtx_energy+recoil <=0)vtx_energy_fraction=-0.005; 
    //When there is absolutely no reconstructed recoil we need to avoid NAN.
    if(vtx_energy/(vtx_energy+recoil) > 1) cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
    //Neutrino Energies
    double nu_cal_E = (evt->CCQENu_muon_E+recoil_inc)/1e3;
    double nu_qe_E = evt->CCQENu_enu_muon/1e3;

    

    //------------------------------------------------------------------
    // Since this event passed main neutron candidate selection,
    // it has a main candidate
    //------------------------------------------------------------------
    // muon 4P
    XYZTVector muon( reco_muon_px, reco_muon_py, reco_muon_pz, evt->CCQENu_muon_E/1e3 );
    // neutrino direction
    XYZVector beamAxis = geoUtils->BeamAxis( 0 );
    // blob direction wrt vertex
    XYZVector blobVtxDir = ncutter->MainNeutronCandidate()->BlobVtxDir;
    // calculated theoretical neutron direction assuming Hydrogen, Eb = 0
    XYZTVector expectedFSP = geoUtils->ComputeExpectedNucleon( beamAxis, muon );
    XYZVector expectedFSPp = (XYZVector) expectedFSP.Vect();
    
    // Neutron Speed
    double beta_nc = expectedFSP.Beta();
    double speed_nc = beta_nc*TMath::C();

    // NC distance
    XYZTVector blobVtxVec= ncutter->MainNeutronCandidate()->BlobVtxVec;
    double blobVtxZ = ncutter->MainNeutronCandidate()->BlobVtxVec.Z();
    double blobvtxdR=blobVtxVec.R();
    double blobvtxdT=blobVtxVec.T();

    double expectedR= expectedFSPp.R()*blobVtxZ/expectedFSPp.Z();
    double expecteddT = expectedR/1000/speed_nc*10e9;
    double dT = blobvtxdT-expecteddT;
    //cout<<blobvtxdR<<endl;
    //cout<<blobvtxdT<<" "<<expecteddT<<" "<<blobvtxdT<<endl;

    // Perpendicular and reaction plane angles
    std::vector<double> neutAngularVars = geoUtils->ComputeNeutronAngularVars( beamAxis, expectedFSPp,blobVtxDir );
    // Inferred neutron momentum
    XYZVector inferredNeutP = blobVtxDir.Unit() * expectedFSP.R();
    XYZTVector inferredNeut4P( inferredNeutP.X(), inferredNeutP.Y(), inferredNeutP.Z(), expectedFSP.E() );

    double PerpAngle = neutAngularVars[0], ReactAngle = neutAngularVars[1], dTheta = neutAngularVars[2];
    double Q2 = evt->CCQENu_Q2/1e6 ;
    int charge = (neutrinoMode)? -1 : +1;
    double bindingE = 0;
    double Q2_p = physicsUtils->qsqCCQE( beamAxis, muon, charge, bindingE) ;

    //---------------------------------
    // Transverse Variables //only phiT for neutrons
    //---------------------------------
    double transverseEb = geoUtils->EbAntinuGeV() + geoUtils->MnGeV() - geoUtils->MpGeV(); //make Remnant Nucleon Mass compliant with anti neutrino
    std::vector<double> neutTransverseVars = geoUtils->GetTransverseVariables( beamAxis, muon, inferredNeut4P, geoUtils->C12GeV(), transverseEb );
    double Enu=      neutTransverseVars[0];
    double Pn=       neutTransverseVars[1];
    double dpT=      neutTransverseVars[2];
    double dpTx=     neutTransverseVars[3];
    double dpTy=     neutTransverseVars[4];
    double dalphaT=  neutTransverseVars[5];
    double phiTDeg=     neutTransverseVars[6];
    double sign = neutTransverseVars[7];
    double phiTRadian= phiTDeg/180*3.141592654*sign;


    //--------------------------------------
    // recoil and vtx e
    // -------------------------------------
    double recoil_total = ncutter->GetNCRecoilEnergy(1);
    double recoil_total_3D = ncutter->GetNCRecoilEnergy(3);
    double recoil_total_2D = ncutter->GetNCRecoilEnergy(2);

    double recoil_ratio = ncutter->MainNeutronCandidate()->TotalE/recoil_total;
    //cout<<ncutter->MainNeutronCandidate()->TotalE<<" "<<recoil_total<<" "<<recoil_ratio<<endl;
    recoil_total-= ncutter->MainNeutronCandidate()->TotalE; 
    recoil_total_3D -= ncutter->MainNeutronCandidate()->TotalE;

    double vtxE=evt->vtx_blobs_energy;
    recoil_ratio = vtxE/recoil_total;
    //--------------------------------------
    // test cut
    //--------------------------------------
    //test cut
    //if( abs(PerpAngle) > 0.2 ) continue;
    //if( vtxE > 10 ) continue;


    //-------------------------------------------------------------------
    // Check that the CV of events pass the proton and recoil cuts 
    // This is because we have systematics due to these cuts 
    // and a change in their values can affect the composition 
    // of the selected sample. 
    //-------------------------------------------------------------------
    // only for MC
    //bool doesCVpass_SingleProtonCut = (evt->multiplicity==1)? true : cutter->passSingleProtonCut( evt, 0, 0 );
    //bool doesCVpass_ExtraProtonsCut = (evt->multiplicity==1)? true : cutter->passExtraProtonsCut( evt, NULL, 0 );
    //bool doesCVpass_RecoilCut = cutter->passSignalFunction( evt, 0, 0 );   

    //bool fillHisto = (isMC)? ( doesCVpass_SingleProtonCut && doesCVpass_ExtraProtonsCut && doesCVpass_RecoilCut ) : true;

    bool fillHisto = true;
    if (DEBUG) cout<<"pre-fill"<<endl;
    if( fillHisto )
    {
            //2D-neutron histograms
      //-------------------------------------------------------
      //Fill 2-D Plots for Neutrons
      //-------------------------------------------------------
      if(DEBUG) cout<<"fillHisto Start Neutron"<<endl;
      utils->fillHistos( utils->histos2D["h_neutron_dtheta_qsq"], dTheta, Q2, isMC, evt, wgt);
      utils->fillHistos( utils->histos2D["h_neutron_dreactplane_qsq"], ReactAngle, Q2, isMC, evt, wgt);
      utils->fillHistos( utils->histos2D["h_neutron_dperp_qsq"], PerpAngle, Q2, isMC, evt, wgt);
      utils->fillHistos( utils->histos2D["h_neutron_dperp_dreactplane"], PerpAngle, ReactAngle, isMC, evt, wgt);

      utils->fillHistos( utils->histos2D["h_neutron_phit_dtheta"], phiTRadian, dTheta, isMC, evt, wgt);
      utils->fillHistos( utils->histos2D["h_neutron_phit_dperp"], phiTRadian, PerpAngle, isMC, evt, wgt);
      utils->fillHistos( utils->histos2D["h_neutron_phit_dreactplane"], phiTRadian, ReactAngle, isMC, evt, wgt);
      utils->fillHistos( utils->histos2D["h_neutron_phit_qsq"], phiTRadian, Q2 , isMC, evt, wgt);

      utils->fillHistos( utils->histos2D["h_neutron_dreactplane_dt"], ReactAngle, dT, isMC, evt, wgt);
      utils->fillHistos( utils->histos2D["h_neutron_dperp_dt"], PerpAngle, dT, isMC, evt, wgt);


      utils->fillHistos( utils->histos2D["h_neutron_dreactplane_recoilE"], ReactAngle, recoil_total, isMC, evt, wgt);
      utils->fillHistos( utils->histos2D["h_neutron_dreactplane_recoilE_ratio"], ReactAngle, recoil_ratio, isMC, evt, wgt);
      utils->fillHistos( utils->histos2D["h_neutron_dreactplane_vtxE"], ReactAngle, evt->vtx_blobs_energy, isMC, evt, wgt);

      utils->fillHistos( utils->histos2D["h_study_recoilE_qsq"], recoil_NoMC, Q2, isMC, evt, wgt );
      utils->fillHistos( utils->histos2D["h_study_recoilEDefault_qsq"], recoil, Q2, isMC, evt, wgt );
      utils->fillHistos( utils->histos2D["h_study_recoilE_pz"], recoil_NoMC, muon_pz_beam, isMC, evt, wgt );
      utils->fillHistos( utils->histos2D["h_study_recoilE_pt"], recoil_NoMC, muon_pt_beam, isMC, evt, wgt );

      //utils->fillHistos( utils->histos2D["h_study_recoilE05_neutE05"], recoil_all_nonMuon-muon_E-neutron_cone_E_05, neutron_cone_E_05, isMC, evt, wgt );
      //utils->fillHistos( utils->histos2D["h_study_recoilE10_neutE10"], recoil_all_nonMuon-muon_E-neutron_cone_E_10, neutron_cone_E_10, isMC, evt, wgt );
      //utils->fillHistos( utils->histos2D["h_study_recoilE15_neutE15"], recoil_all_nonMuon-muon_E-neutron_cone_E_15, neutron_cone_E_15, isMC, evt, wgt );
      //utils->fillHistos( utils->histos2D["h_study_recoilE20_neutE20"], recoil_all_nonMuon-muon_E-neutron_cone_E_20, neutron_cone_E_20, isMC, evt, wgt );
      //utils->fillHistos( utils->histos2D["h_study_recoilE30_neutE30"], recoil_all_nonMuon-muon_E-neutron_cone_E_30, neutron_cone_E_30, isMC, evt, wgt );
      //utils->fillHistos( utils->histos2D["h_study_recoilE45_neutE45"], recoil_all_nonMuon-muon_E-neutron_cone_E_45, neutron_cone_E_45, isMC, evt, wgt );

      utils->fillHistos( utils->histos2D["h_study_recoilE05_neutE05"], recoil-muon_E-neutron_cone_E_05, neutron_cone_E_05, isMC, evt, wgt );
      utils->fillHistos( utils->histos2D["h_study_recoilE10_neutE10"], recoil-muon_E-neutron_cone_E_10, neutron_cone_E_10, isMC, evt, wgt );
      utils->fillHistos( utils->histos2D["h_study_recoilE15_neutE15"], recoil-muon_E-neutron_cone_E_15, neutron_cone_E_15, isMC, evt, wgt );
      utils->fillHistos( utils->histos2D["h_study_recoilE20_neutE20"], recoil-muon_E-neutron_cone_E_20, neutron_cone_E_20, isMC, evt, wgt );
      utils->fillHistos( utils->histos2D["h_study_recoilE30_neutE30"], recoil-muon_E-neutron_cone_E_30, neutron_cone_E_30, isMC, evt, wgt );
      utils->fillHistos( utils->histos2D["h_study_recoilE45_neutE45"], recoil-muon_E-neutron_cone_E_45, neutron_cone_E_45, isMC, evt, wgt );



      if(DEBUG) cout<<utils->histos2D["h_neutron_phit_qsq"][kMC]->Integral()<<endl;
      if(DEBUG) cout<<"fillHisto end Neutron"<<endl;
      //3-D
      utils->fillHistos( utils->histos3D["h_enucal_enulep_pt"], nu_cal_E, nu_qe_E, muon_pt_beam, isMC, evt, wgt);
      
      if(DEBUG) cout<<"fillHisto part 1"<<endl;
      // Now filling MC only histos and error bands
      if (isMC)
      {
        if(DEBUG) cout<<"fillHisto part MC"<<endl;
        double true_muon_px   = evt->mc_primFSLepton[0]/pow(10,3);//GeV
        double true_muon_py   = evt->mc_primFSLepton[1]/pow(10,3);
        double true_muon_pz   = evt->mc_primFSLepton[2]/pow(10,3);
        double true_muon_p = utils->GetTotalMomentum( true_muon_px, true_muon_py, true_muon_pz ); //GeV
        double true_theta     = utils->GetTheta( true_muon_px, true_muon_py, true_muon_pz );
        double true_muon_pt_beam   = utils->GetTransverseMomentumWRTBeam( true_muon_px, true_muon_py, true_muon_pz ); //GeV 
        double true_muon_pz_beam   = utils->GetLongitudinalMomentumWRTBeam( true_muon_px, true_muon_py, true_muon_pz ); //GeV
        double true_q2_tmk    = utils->GetTrueCCQEQ2( true_muon_p*1000.0, true_theta, bindingE );//MeV NEED INPUTS IN MEV
        double true_enu       = utils->GetTrueCCQENeutrinoEnergy( true_muon_p*1000, true_theta, bindingE );//MeV NEED INPUTS IN MEV

        //------------------------------------------------------
        //Fill 2D Vertical Error Bands
        //Error Bands are not filled up for tmu_costheta 
        //and pmu_ptmu histos to avoid code slow down 
        //------------------------------------------------------
        if(DEBUG) cout<<"fillHisto part MC 2"<<endl;
        if(RunCodeWithSystematics)
        {
          if(DEBUG) cout<<"fillHisto Sys  part MC 2"<<endl;
          bool doVertErrBand = true, doLatErrBand = true;
          if ( doVertErrBand )
          {
            /*
            utils->fillVertErrorBands( utils->histos2D["h_pzmu_ptmu"], muon_pz_beam, muon_pt_beam, evt );
            utils->fillVertErrorBands( utils->histos2D["h_enu_ptmu"], evt->CCQENu_enu_muon/1e3, muon_pt_beam, evt );
            utils->fillVertErrorBands( utils->histos2D["h_q2_ptmu"], evt->CCQENu_Q2/1e6, muon_pt_beam, evt );
            utils->fillVertErrorBands( utils->histos2D["h_recoil_inc_ptmu"], recoil_inc, muon_pt_beam, evt);
            utils->fillVertErrorBands( utils->histos2D["h_recoil_ptmu"], recoil, muon_pt_beam, evt);
            if(evt->CCQENu_Q2/1e6 < 0.2){
              utils->fillVertErrorBands( utils->histos2D["h_vtx_old_low_q2_300mm"], evt->recoil_energy_nonmuon_vtx300mm, muon_pt_beam, evt);
              utils->fillVertErrorBands( utils->histos2D["h_vtx_old_low_q2_250mm"], evt->recoil_energy_nonmuon_vtx250mm, muon_pt_beam, evt);
              utils->fillVertErrorBands( utils->histos2D["h_vtx_old_low_q2_200mm"], evt->recoil_energy_nonmuon_vtx200mm, muon_pt_beam, evt);
              utils->fillVertErrorBands( utils->histos2D["h_vtx_old_low_q2_150mm"], evt->recoil_energy_nonmuon_vtx150mm, muon_pt_beam, evt);
              utils->fillVertErrorBands( utils->histos2D["h_vtx_new_low_q2"], vtx_energy, muon_pt_beam, evt);
            }
            else{
              utils->fillVertErrorBands( utils->histos2D["h_vtx_old_high_q2_300mm"], evt->recoil_energy_nonmuon_vtx300mm, muon_pt_beam, evt);
              utils->fillVertErrorBands( utils->histos2D["h_vtx_old_high_q2_250mm"], evt->recoil_energy_nonmuon_vtx250mm, muon_pt_beam, evt);
              utils->fillVertErrorBands( utils->histos2D["h_vtx_old_high_q2_200mm"], evt->recoil_energy_nonmuon_vtx200mm, muon_pt_beam, evt);
              utils->fillVertErrorBands( utils->histos2D["h_vtx_old_high_q2_150mm"], evt->recoil_energy_nonmuon_vtx150mm, muon_pt_beam, evt);
              utils->fillVertErrorBands( utils->histos2D["h_vtx_new_high_q2"], vtx_energy, muon_pt_beam, evt);
            }
            utils->fillVertErrorBands( utils->histos2D["h_vtx_old_300mm"], evt->recoil_energy_nonmuon_vtx300mm, muon_pt_beam, evt);
            utils->fillVertErrorBands( utils->histos2D["h_vtx_old_250mm"], evt->recoil_energy_nonmuon_vtx250mm, muon_pt_beam, evt);
            utils->fillVertErrorBands( utils->histos2D["h_vtx_old_200mm"], evt->recoil_energy_nonmuon_vtx200mm, muon_pt_beam, evt);
            utils->fillVertErrorBands( utils->histos2D["h_vtx_old_150mm"], evt->recoil_energy_nonmuon_vtx150mm, muon_pt_beam, evt);
            */
            //neutron
            utils->fillVertErrorBands( utils->histos2D["h_neutron_dtheta_qsq"], dTheta, Q2, evt );
            utils->fillVertErrorBands( utils->histos2D["h_neutron_dreactplane_qsq"], ReactAngle, Q2, evt );
            utils->fillVertErrorBands( utils->histos2D["h_neutron_dperp_qsq"], PerpAngle, Q2, evt );
            utils->fillVertErrorBands( utils->histos2D["h_neutron_dperp_dreactplane"], PerpAngle, ReactAngle, evt );

            utils->fillVertErrorBands( utils->histos2D["h_neutron_phit_dtheta"], phiTRadian, dTheta, evt );
            utils->fillVertErrorBands( utils->histos2D["h_neutron_phit_dperp"], phiTRadian, PerpAngle, evt );
            utils->fillVertErrorBands( utils->histos2D["h_neutron_phit_dreactplane"], phiTRadian, ReactAngle, evt );
            utils->fillVertErrorBands( utils->histos2D["h_neutron_phit_qsq"], phiTRadian, Q2 , evt );

            utils->fillVertErrorBands( utils->histos2D["h_vtx_new"], vtx_energy, muon_pt_beam, evt);
            //utils->fillVertErrorBands( h_tmu_theta, muon_T, muon_theta, evt );
            //utils->fillVertErrorBands( h_tmu_costheta, muon_T, cos_muon_theta, evt );
            //utils->fillVertErrorBands( h_pmu_ptmu, muon_p, muon_pt, evt );
            if(DEBUG) cout<<"fillHisto Sys  part MC 2"<<endl;
          }
         
          if ( doLatErrBand )
          {
            /*
         
            if(DEBUG) cout<<"fillHisto lat Sys  part MC 2"<<endl;
            //utils->fillLatErrorBands( h_theta, "theta", muon_theta, evt, true );
            //utils->fillLatErrorBands( h_tmu,   "tmu",   muon_T, evt, true );
            utils->fillLatErrorBands( utils->histos1D["h_ptmu"],  "pT",    muon_pt_beam, evt, true );
            utils->fillLatErrorBands( utils->histos1D["h_pzmu"],  "pZ",    muon_pz_beam, evt, true );
            utils->fillLatErrorBands( utils->histos1D["h_vtx_energy"],  "vtx_energy",    vtx_energy, evt, true );
            
            //----------------------------------------------------
            //Fill 2D Lateral Error Bands
            //Error Bands are not filled up for tmu_costheta 
            //and pmu_ptmu histos to avoid code slow down 
            //----------------------------------------------------
            //  utils->fillLatErrorBands( h_tmu_theta, "tmu", "theta", muon_T, muon_theta, evt, true );
            //utils->fillLatErrorBands( h_tmu_costheta, "tmu", "theta", muon_T, cos_muon_theta, evt, true );
            //utils->fillLatErrorBands( h_pmu_ptmu, "p", "pT", muon_p, muon_pt, evt, true );
            utils->fillLatErrorBands( utils->histos2D["h_pzmu_ptmu"], "pZ", "pT", muon_pz_beam, muon_pt_beam, evt, true );
            utils->fillLatErrorBands( utils->histos2D["h_enu_ptmu"], "enu", "pT", evt->CCQENu_enu_muon/1e3, muon_pt_beam, evt, true );
            utils->fillLatErrorBands( utils->histos2D["h_q2_ptmu"], "q2", "pT", evt->CCQENu_Q2/1e6, muon_pt_beam, evt, true );
            utils->fillLatErrorBands( utils->histos2D["h_recoil_inc_ptmu"],"recoil","pT", recoil_inc, muon_pt_beam, evt, true);
            utils->fillLatErrorBands( utils->histos2D["h_recoil_ptmu"],"recoil","pT", recoil, muon_pt_beam, evt, true);
            utils->fillLatErrorBands( utils->histos2D["h_recoil_ptmu_before"],"recoil","pT", recoil, muon_pt_beam, evt, true);
            if(evt->CCQENu_Q2/1e6 < 0.2){
              utils->fillLatErrorBands( utils->histos2D["h_vtx_old_low_q2_300mm"],"vtx","pT",evt->recoil_energy_nonmuon_vtx300mm, muon_pt_beam, evt, true);
              utils->fillLatErrorBands( utils->histos2D["h_vtx_old_low_q2_250mm"],"vtx","pT",evt->recoil_energy_nonmuon_vtx250mm, muon_pt_beam, evt, true);
              utils->fillLatErrorBands( utils->histos2D["h_vtx_old_low_q2_200mm"],"vtx","pT",evt->recoil_energy_nonmuon_vtx200mm, muon_pt_beam, evt, true);
              utils->fillLatErrorBands( utils->histos2D["h_vtx_old_low_q2_150mm"],"vtx","pT",evt->recoil_energy_nonmuon_vtx150mm, muon_pt_beam, evt, true);
              utils->fillLatErrorBands( utils->histos2D["h_vtx_new_low_q2"],"vtx","pT",vtx_energy, muon_pt_beam, evt, true);
            }
            else{
              utils->fillLatErrorBands( utils->histos2D["h_vtx_old_high_q2_300mm"],"vtx","pT",evt->recoil_energy_nonmuon_vtx300mm, muon_pt_beam, evt, true);
              utils->fillLatErrorBands( utils->histos2D["h_vtx_old_high_q2_250mm"],"vtx","pT",evt->recoil_energy_nonmuon_vtx250mm, muon_pt_beam, evt, true);
              utils->fillLatErrorBands( utils->histos2D["h_vtx_old_high_q2_200mm"],"vtx","pT",evt->recoil_energy_nonmuon_vtx200mm, muon_pt_beam, evt, true);
              utils->fillLatErrorBands( utils->histos2D["h_vtx_old_high_q2_150mm"],"vtx","pT",evt->recoil_energy_nonmuon_vtx150mm, muon_pt_beam, evt, true);
              utils->fillLatErrorBands( utils->histos2D["h_vtx_new_high_q2"],"vtx","pT",vtx_energy, muon_pt_beam, evt, true);
            }
            utils->fillLatErrorBands( utils->histos2D["h_vtx_old_300mm"],"vtx","pT",evt->recoil_energy_nonmuon_vtx300mm, muon_pt_beam, evt, true);
            utils->fillLatErrorBands( utils->histos2D["h_vtx_old_250mm"],"vtx","pT",evt->recoil_energy_nonmuon_vtx250mm, muon_pt_beam, evt, true);
            utils->fillLatErrorBands( utils->histos2D["h_vtx_old_200mm"],"vtx","pT",evt->recoil_energy_nonmuon_vtx200mm, muon_pt_beam, evt, true);
            utils->fillLatErrorBands( utils->histos2D["h_vtx_old_150mm"],"vtx","pT",evt->recoil_energy_nonmuon_vtx150mm, muon_pt_beam, evt, true);
            utils->fillLatErrorBands( utils->histos2D["h_vtx_new"],"vtx","pT",vtx_energy, muon_pt_beam, evt, true);
            */
          }
        }
      }
    } 
    if( isMC)
    {
        if(DEBUG) cout<<"fillHisto more mc"<<endl;
        /*
        utils->fillHistos( utils->histos2D["h_recoil_ptmu_before"], recoil, muon_pt_beam, evt, wgt);
        utils->fillVertErrorBands( utils->histos2D["h_recoil_ptmu_before"], recoil, muon_pt_beam, evt);
        */
    }
  }//end looping
  if (isMC) cout << "Finished looping over all MC entries. Number of MC events satisfying all selection criteria = " << n_outputevent << endl; 
  else cout << "Finished looping over all DATA entries. Number of DATA events satisfying all selection criteria = " << n_outputevent << endl; 
  return n_outputevent;
}


int NeutronSelectionHists(string filename, string sample, bool makeFluxConstraintHisto, int multiplicity, string playlist, int n_mcfiles = -1, int n_datafiles = -1)
{
  //---------------------------------------------
  // create file to store histograms
  //---------------------------------------------
  //TFile *f = new TFile( filename.c_str(), "RECREATE" )
  cout<<"Begin NeutronSelectionHists"<<endl;

  CCQENuUtils  *utils  = new CCQENuUtils( false, makeFluxConstraintHisto );
  utils->setPlaylist(playlist);
  utils->setFluxReweighterPlaylist();
  CCQENuCuts   *cutter = new CCQENuCuts();
  AnaBinning *binner = new AnaBinning();
  CCQENuBinning *minmodbinner = new CCQENuBinning();

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
  axis_binning fineEnu;
  axis_binning fineEnuDiff;
  axis_binning vtxEnergy;
  axis_binning vtxEnergyFraction;





  vector<double> bins;
  for(float i = 0.0; i<=20.0; i+=0.2) bins.push_back(i);
  fineEnu.bin_edges = bins;
  fineEnu.nbins = bins.size()-1;
  fineEnu.min = bins.front();
  fineEnu.max = bins.back();
  bins.clear();
  for(float i = -2.0; i<=2.0; i+=0.1) bins.push_back(i);
  fineEnuDiff.bin_edges = bins;
  fineEnuDiff.nbins = bins.size()-1;
  fineEnuDiff.min = bins.front();
  fineEnuDiff.max = bins.back();
  bins.clear();
  //New binning
  if(multiplicity==1){
    bins.push_back(0);
    bins.push_back(10);
    bins.push_back(20);
    bins.push_back(30);
    bins.push_back(40);
    bins.push_back(50);
    bins.push_back(60);
    bins.push_back(70);
    bins.push_back(80);
    bins.push_back(90);
    bins.push_back(100);
    bins.push_back(120);
    bins.push_back(140);
    bins.push_back(160);
    bins.push_back(180);
    bins.push_back(200);
    bins.push_back(250);
    bins.push_back(300);
    bins.push_back(400);
    bins.push_back(500);
  }
  else{
    bins.push_back(0);
    bins.push_back(20);
    bins.push_back(40);
    bins.push_back(60);
    bins.push_back(80);
    bins.push_back(100);
    bins.push_back(200);
    bins.push_back(400);
    bins.push_back(500);  
  }



  vtxEnergy.bin_edges = bins;
  vtxEnergy.nbins = bins.size()-1;
  vtxEnergy.min = bins.front();
  vtxEnergy.max = bins.back();

  //now for recoil fraction
  bins.clear();
  for(float i=-0.01;i<=1.01;i+=0.01) bins.push_back(i);
  vtxEnergyFraction.bin_edges = bins;
  vtxEnergyFraction.nbins = bins.size()-1;
  vtxEnergyFraction.min = bins.front();
  vtxEnergyFraction.max = bins.back();



  axis_binning uniformbins_0_2GeV;
  bins.clear();
  for(float i=0;i<=2;i+=0.05) bins.push_back(i);
  uniformbins_0_2GeV.bin_edges = bins;
  uniformbins_0_2GeV.nbins = bins.size()-1;
  uniformbins_0_2GeV.min = bins.front();
  uniformbins_0_2GeV.max = bins.back();

   axis_binning uniformbins_0_8GeV;
  bins.clear();
  for(float i=0;i<=8;i+=0.05) bins.push_back(i);
  uniformbins_0_8GeV.bin_edges = bins;
  uniformbins_0_8GeV.nbins = bins.size()-1;
  uniformbins_0_8GeV.min = bins.front();
  uniformbins_0_8GeV.max = bins.back();

  

  axis_binning uniformbins_0_500MeV;
  bins.clear();
  for(float i=0;i<=500;i+=5) bins.push_back(i);
  uniformbins_0_500MeV.bin_edges = bins;
  uniformbins_0_500MeV.nbins = bins.size()-1;
  uniformbins_0_500MeV.min = bins.front();
  uniformbins_0_500MeV.max = bins.back();
  

 


 

  cout<<"BOOK HISTOGRAMS"<<endl;
   //Neutrons

  NeutronBlobCuts *ncutter = new NeutronBlobCuts();
  ncutter->SetDebug(0);
  NeutronBlobBinning *neutblobbinner = new NeutronBlobBinning();

  axis_binning neutreactplanebins = neutblobbinner->GetThetaReactionPlaneBinsRadian();
  axis_binning neutperpendicularbins = neutblobbinner->GetThetaPerpBinsRadian();
  axis_binning neutthetabins = neutblobbinner->GetThetaBinsRadian();
  axis_binning neutphitbins = neutblobbinner->GetPhiTBinsRadian();
  axis_binning neutdTbins = neutblobbinner->GetNeutronTimingNs();

  MnvH2D *h_neutron_dperp_dreactplane[nHistos], *h_neutron_dtheta_qsq[nHistos], *h_neutron_dperp_qsq[nHistos], *h_neutron_dreactplane_qsq[nHistos];
  
  MnvH2D *h_neutron_phit_qsq[nHistos], *h_neutron_phit_dtheta[nHistos], *h_neutron_phit_dperp[nHistos], *h_neutron_phit_dreactplane[nHistos];

  MnvH2D *h_neutron_dperp_dt[nHistos], *h_neutron_dreactplane_dt[nHistos];



  utils->bookHistos( h_neutron_dperp_qsq, "h_neutron_dperp_qsq","Neutron Candidate #theta_{Perp} vs Q^{2}_{QE}: #theta_{Perp} (radian) : Q^{2}_{QE} (GeV^{2})", neutperpendicularbins, Q2bins );
  utils->bookHistos( h_neutron_dreactplane_qsq, "h_neutron_dreactplane_qsq","Neutron Candidate #theta_{React} vs Q^{2}_{QE}: #theta_{React} (radian) : Q^{2}_{QE} (GeV^{2})", neutreactplanebins, Q2bins );
  utils->bookHistos( h_neutron_dtheta_qsq, "h_neutron_dtheta_qsq","Neutron Candidate #theta vs Q^{2}_{QE}: #theta (radian) : Q^{2}_{QE} (GeV^{2})", neutthetabins, Q2bins );
  utils->bookHistos( h_neutron_dperp_dreactplane, "h_neutron_dperp_dreactplane","Neutron Candidate #theta_{Perp} vs #theta_{React}: #theta_{Perp} (radian) :  #theta_{Perp} (radian)", neutperpendicularbins, neutthetabins );
  utils->bookHistos( h_neutron_phit_qsq, "h_neutron_phit_qsq","Neutron Candidate #phi_{T} vs Q^{2}_{QE}: #phi_{T} (radian): Q^{2}_{QE} (GeV^{2})", neutphitbins, Q2bins );
  utils->bookHistos( h_neutron_phit_dtheta, "h_neutron_phit_dtheta","Neutron Candidate #phi_{T} vs #theta: #phi_{T} (radian): #theta (radian)", neutphitbins, neutthetabins );
  utils->bookHistos( h_neutron_phit_dperp, "h_neutron_phit_dperp","Neutron Candidate #phi_{T} vs #theta_{Perp}: #phi_{T} (radian): #theta_{Perp} (radian)", neutphitbins,neutperpendicularbins );
  utils->bookHistos( h_neutron_phit_dreactplane, "h_neutron_phit_dreactplane","Neutron Candidate #phi_{T} vs #theta_{React}: #phi_{T} (radian): #theta_{React} (radian)", neutphitbins, neutreactplanebins );
  cout<<"Stage 1"<<endl;

  utils->bookHistos( h_neutron_dperp_dt, "h_neutron_dperp_dt", "NC #theta_{P} vs dT: #theta_{P} : dT (ns)", neutperpendicularbins, neutdTbins );
  utils->bookHistos( h_neutron_dreactplane_dt, "h_neutron_dreactplane_dt", "NC #theta_{R} vs dT: #theta_{R} : dT (ns)", neutreactplanebins, neutdTbins );


  MnvH2D *h_neutron_dreactplane_vtxE[nHistos], *h_neutron_dreactplane_recoilE[nHistos],*h_neutron_dreactplane_recoilE_ratio[nHistos];
  utils->bookHistos( h_neutron_dreactplane_vtxE, "h_neutron_dreactplane_vtxE", "NC #theta_{R} vs vtxE: #theta_{R} : vtxE (GeV)", neutreactplanebins, vtxEnergy);
  utils->bookHistos( h_neutron_dreactplane_recoilE, "h_neutron_dreactplane_recoilE", "NC #theta_{R} vs Recoil: #theta_{R} : recoil (GeV)", neutreactplanebins, blobEnergyBinsMeV);
  utils->bookHistos( h_neutron_dreactplane_recoilE_ratio, "h_neutron_dreactplane_recoilE_ratio", "NC #theta_{R} vs RecoilRatio: #theta_{R} : recoil ratio", neutreactplanebins, vtxEnergyFraction);

  cout<<"Stage 2"<<endl;
  MnvH2D *h_study_recoilE_pz[nHistos];
  MnvH2D *h_study_recoilE_pt[nHistos];
  MnvH2D *h_study_recoilE_qsq[nHistos];
  MnvH2D *h_study_recoilEDefault_qsq[nHistos];

  cout<<"Stage 3"<<endl;
  utils->bookHistos( h_study_recoilE_qsq, "h_study_recoilE_qsq", "recoil,qsq: E (GeV) : Q^{2}_{CCQE} GeV^{2}", uniformbins_0_500MeV, uniformbins_0_2GeV );
  utils->bookHistos( h_study_recoilEDefault_qsq, "h_study_recoilEDefault_qsq", "recoil,qsq: E (GeV) : Q^{2}_{CCQE} GeV^{2}", uniformbins_0_500MeV, uniformbins_0_2GeV );
  utils->bookHistos( h_study_recoilE_pz, "h_study_recoilE_pz", "recoil,pz : E (GeV) : pz (GeV)", uniformbins_0_500MeV, uniformbins_0_8GeV);
  utils->bookHistos( h_study_recoilE_pt, "h_study_recoilE_pt", "recoil,pt : E (GeV) : pt (GeV)", uniformbins_0_500MeV, uniformbins_0_2GeV);


  MnvH2D *h_study_recoilE05_neutE05[nHistos];
  MnvH2D *h_study_recoilE10_neutE10[nHistos];
  MnvH2D *h_study_recoilE15_neutE15[nHistos];
  MnvH2D *h_study_recoilE20_neutE20[nHistos];
  MnvH2D *h_study_recoilE30_neutE30[nHistos];
  MnvH2D *h_study_recoilE45_neutE45[nHistos];
  utils->bookHistos( h_study_recoilE05_neutE05, "h_study_recoilE05_neutE05", "recoil,neut: E (MeV) : E (MeV)", uniformbins_0_500MeV, uniformbins_0_500MeV );
  utils->bookHistos( h_study_recoilE15_neutE15, "h_study_recoilE15_neutE15", "recoil,neut: E (MeV) : E (MeV)", uniformbins_0_500MeV, uniformbins_0_500MeV );
  utils->bookHistos( h_study_recoilE10_neutE10, "h_study_recoilE10_neutE10", "recoil,neut: E (MeV) : E (MeV)", uniformbins_0_500MeV, uniformbins_0_500MeV );
  utils->bookHistos( h_study_recoilE20_neutE20, "h_study_recoilE20_neutE20", "recoil,neut: E (MeV) : E (MeV)", uniformbins_0_500MeV, uniformbins_0_500MeV );
  utils->bookHistos( h_study_recoilE30_neutE30, "h_study_recoilE30_neutE30", "recoil,neut: E (MeV) : E (MeV)", uniformbins_0_500MeV, uniformbins_0_500MeV );
  utils->bookHistos( h_study_recoilE45_neutE45, "h_study_recoilE45_neutE45", "recoil,neut: E (MeV) : E (MeV)", uniformbins_0_500MeV, uniformbins_0_500MeV );


  cout<<"Stage 4"<<endl;
  //3-D
  MnvH3D *h_enucal_enulep_pt[nHistos];

  utils->bookHistos( h_enucal_enulep_pt, "h_enucal_enulep_pt", "enucal enulep pt; enucal;enulep;pt",fineEnu,fineEnu,muonPtbins);

  cout<<"END BOOK HISTOGRAMS"<<endl;
  //--------------------------------------------------------------
  // Add Vertical and Lateral Error Bands
  // JO is removing error bands from all 1D because they are all 
  // present in 2D. Plz make projections directly from the 2D. 
  //--------------------------------------------------------------
  if(RunCodeWithSystematics){
    cout<<"ADD SYS. HISTO"<<endl;
    //Neutron 
    utils->addVertErrorBands( h_neutron_dperp_dreactplane );
    utils->addVertErrorBands( h_neutron_dperp_qsq );
    utils->addVertErrorBands( h_neutron_dreactplane_qsq );
    utils->addVertErrorBands( h_neutron_dtheta_qsq );
    utils->addVertErrorBands( h_neutron_phit_qsq );
    utils->addVertErrorBands( h_neutron_phit_dreactplane );
    utils->addVertErrorBands( h_neutron_phit_dperp );
    utils->addVertErrorBands( h_neutron_phit_dtheta );
    utils->addVertErrorBands( h_neutron_dperp_dt );
    utils->addVertErrorBands( h_neutron_dreactplane_dt );

    utils->addLatErrorBands( h_neutron_dperp_dreactplane );
    utils->addLatErrorBands( h_neutron_dperp_qsq );
    utils->addLatErrorBands( h_neutron_dreactplane_qsq );
    utils->addLatErrorBands( h_neutron_dtheta_qsq );
    utils->addLatErrorBands( h_neutron_phit_qsq );
    utils->addLatErrorBands( h_neutron_phit_dreactplane );
    utils->addLatErrorBands( h_neutron_phit_dperp );
    utils->addLatErrorBands( h_neutron_phit_dtheta );
    utils->addLatErrorBands( h_neutron_dperp_dt );
    utils->addLatErrorBands( h_neutron_dreactplane_dt );





    cout<<"ADDED SYS. HISTO"<<endl;

  }//end run with systematics
  //  return 0;

  double data_events, mc_events;
  //---------------------------------------------
  // Get MC Tree
  //---------------------------------------------
  cout<<"Getting MC Trees"<<endl;
  TChain *truth_mc = utils->getMCTree("Truth", n_mcfiles );
  TChain* tree_mc   =  utils->getMCTree("CCQENu", n_mcfiles );
  int entries_mc   = tree_mc->GetEntries();
  cout << "MC entries: " << entries_mc << endl;
  utils->setmnvHadronReweightTruthTree(truth_mc);
  utils->setmnvHadronReweightDataTree(tree_mc);
  CCQENuEvent* mc = new CCQENuEvent( tree_mc );
  //int AnalysisLooper( CCQENuEvent* evt, CCQENuUtils* utils, CCQENuCuts *cutter, NeutronBlobCuts* ncutter, string sample, int multiplicity, int n_entries, bool isMC = true )
  
  cout<<"Running AnalysisLooper"<<endl;
  mc_events = AnalysisLooper( mc, utils, cutter, ncutter, sample, multiplicity, entries_mc, true );
  
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
    data_events = AnalysisLooper( data, utils, cutter, ncutter, sample, multiplicity,entries_data, false );
    
    //cout << "Finished looping over all DATA entries. Number of DATA events satisfying all selection criteria = " << data_events << endl; 
    delete data; 
    delete tree_data; 
  }

  //Print out some information
  /*
  cout<<" --------------------------------------------- " << endl;
  cout<<" MC " << endl;
  cout<<" --------------------------------------------- " << endl;
  cout<<" ****************** " << endl;
  cout<<" MC 1 Track Events " << endl;
  cout<<" ****************** " << endl;
  cout<< Form("nocuts = %i; cut1 = %i, cut2 = %i, cut3 = %i, cut4 = %i, cut5 = %i, cut6 = %i, cut6_sideband = %i", mc_nocuts_1track, mc_cut1_1track, mc_cut2_1track, mc_cut3_1track, mc_cut4_1track, mc_cut5_1track, mc_cut6_1track, mc_cut6_1track_sideband) << endl;
  cout<<" ****************** " << endl;
  cout<<" MC 2 Track Events " << endl;
  cout<<" ****************** " << endl;
  cout<< Form("nocuts = %i; cut1 = %i, cut2 = %i, cut3 = %i, cut4 = %i, cut5 = %i, cut6 = %i, cut7 = %i, cut7_sideband = %i", mc_nocuts_2track, mc_cut1_2track, mc_cut2_2track, mc_cut3_2track, mc_cut4_2track, mc_cut5_2track, mc_cut6_2track, mc_cut7_2track, mc_cut7_2track_sideband) << endl;
  cout<<" --------------------------------------------- " << endl;
  cout<<" DATA " << endl;
  cout<<" --------------------------------------------- " << endl;
  cout<<" ****************** " << endl;
  cout<<" DATA 1 Track Events " << endl;
  cout<<" ****************** " << endl;
  cout<< Form("nocuts = %i; cut1 = %i, cut2 = %i, cut3 = %i, cut4 = %i, cut5 = %i, cut6 = %i, cut6_sideband = %i", data_nocuts_1track, data_cut1_1track, data_cut2_1track, data_cut3_1track, data_cut4_1track, data_cut5_1track, data_cut6_1track, data_cut6_1track_sideband) << endl;
  cout<<" ****************** " << endl;
  cout<<" DATA 2 Track Events " << endl;
  cout<<" ****************** " << endl;
  cout<< Form("nocuts = %i; cut1 = %i, cut2 = %i, cut3 = %i, cut4 = %i, cut5 = %i, cut6 = %i, cut7 = %i, cut7_sideband = %i", data_nocuts_2track, data_cut1_2track, data_cut2_2track, data_cut3_2track, data_cut4_2track, data_cut5_2track, data_cut6_2track, data_cut7_2track, data_cut7_2track_sideband) << endl;
  */
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
  
  f->cd();
  for( unsigned int i=kData; i<nHistos; ++i ){
    /*
    h_theta[i]->Write();
    h_tmu[i]->Write();
    h_multiplicity[i]->Write();
    h_ptmu[i]->Write();
    h_pzmu[i]->Write();
    h_tmu_theta[i]->Write();
    h_tmu_costheta[i]->Write();
    h_pmu_ptmu[i]->Write();
    h_pzmu_ptmu[i]->Write();
    h_enu_ptmu[i]->Write();
    h_q2_ptmu[i]->Write();
    h_true_w[i]->Write();
    h_w[i]->Write();
    h_recoil_inc_ptmu[i]->Write();
    h_recoil_ptmu[i]->Write();
    h_recoil_ptmu_before[i]->Write();
    h_true_w_ptmu[i]->Write();
    h_vtx_energy[i]->Write();
    h_enucal_enulep_pt[i]->Write();
    h_Ediff_ptmu[i]->Write();
    h_vtxrecoilratio_pt[i]->Write();
    h_vtxrecoilratio[i]->Write();
    h_vtxrecoilratio_lowq2[i]->Write();
    h_vtxrecoilratio_highq2[i]->Write();
    h_vtx_old_low_q2_300mm[i]->Write();
    h_vtx_old_low_q2_250mm[i]->Write();
    h_vtx_old_low_q2_200mm[i]->Write();
    h_vtx_old_low_q2_150mm[i]->Write();
    h_vtx_new_low_q2[i]->Write();
    h_vtx_old_high_q2_300mm[i]->Write();
    h_vtx_old_high_q2_250mm[i]->Write();
    h_vtx_old_high_q2_200mm[i]->Write();
    h_vtx_old_high_q2_150mm[i]->Write();
    h_vtx_new_high_q2[i]->Write();
    h_vtx_old_300mm[i]->Write();
    h_vtx_old_250mm[i]->Write();
    h_vtx_old_200mm[i]->Write();
    h_vtx_old_150mm[i]->Write();
    h_vtx_new[i]->Write();
    */
        
    h_neutron_dperp_dreactplane[i]->Write(); 
    h_neutron_dtheta_qsq[i]->Write(); 
    h_neutron_dperp_qsq[i]->Write(); 
    h_neutron_dreactplane_qsq[i]->Write();

    h_neutron_phit_qsq[i]->Write(); 
    h_neutron_phit_dtheta[i]->Write(); 
    h_neutron_phit_dperp[i]->Write(); 
    h_neutron_phit_dreactplane[i]->Write();

    h_neutron_dperp_dt[i]->Write();
    h_neutron_dreactplane_dt[i]->Write();
    
    h_neutron_dreactplane_vtxE[i]->Write();
    h_neutron_dreactplane_recoilE[i]->Write();
    h_neutron_dreactplane_recoilE_ratio[i]->Write();
    
    //h_study_recoil3D_recoil2D_recoilAll[i]->Write();
    h_study_recoilE_pz[i]->Write();
    h_study_recoilE_pt[i]->Write();
    h_study_recoilE_qsq[i]->Write();
    h_study_recoilEDefault_qsq[i]->Write();

    h_study_recoilE05_neutE05[i]->Write();
    h_study_recoilE10_neutE10[i]->Write();
    h_study_recoilE15_neutE15[i]->Write();
    h_study_recoilE20_neutE20[i]->Write();
    h_study_recoilE30_neutE30[i]->Write();
    h_study_recoilE45_neutE45[i]->Write();


    // Uncomment out if you need to make this plot   
    // h_trkChi2PerDoF[i]->Write(); 
    // h_trkChi2PerDoF_inclusive[i]->Write(); 
    //
    //h_neutron_phit_qsq[i]->Write(); 
    //h_theta[i]->Write();
  }
   
  f->Close();
  delete f;
  delete utils;
  delete cutter; 
  delete ncutter;
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
      "\t-./NeutronSelectionHists  Name_and_path_to_Output_file Playlist Selection_Signal  Make_Flux_Constraint_Histo  Multiplicity  Number_MC_files  Number_DATA_files \n\n"<<
      "\t-Name_and_path_to_Output_file\t =\t Name of and path to the Output ROOT file that will be created \n"<<
      "\t-Playlist\t =\t Name of the playlist you want to run over e.g. minerva1 \n"<<
      "\t-Selection_Signal\t =\t Can be: Signal, MichelSideBand, SideBand (this is:non-vtx vs Q2 sideband) \n"<<
      "\t-Make_Flux_Constraint_Histo\t =\t If TRUE Enter 1; If FALSE Enter 0 \n"<<
      "\t-Multiplicity\t =\t Enter 0, 1 or 2. Here 0: 1 or 2 tracks; 1: 1-track events only; 2: 2-track events only \n"<<
      "\t-Number_MC_files\t =\t Number of MonteCarlo files. To use all files, set this to: -1 \n"<<
      "\t-Number_DATA_files\t =\t Number of Data files. To use all files, set this to: -1" << std::endl;
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    return 0; 
  }
  
  //! Default parameters
  std::vector<std::string> par;
  par.push_back("NeutronSelectionHists");
  par.push_back( Form("%s/ana/rootfiles/NeutronSelectionHists.root",getenv("CCQENUROOT") ) );
  par.push_back("minerva1");
  par.push_back("Signal");
  par.push_back("0");
  par.push_back("0");
  par.push_back("-1");
  par.push_back("-1");

  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }

  bool fluxConstraintHistoNeeded = ( par.at(4) == "1" ) ? true : false; 

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;

  return NeutronSelectionHists(par[1], par[3], fluxConstraintHistoNeeded, atoi(par[5].c_str()), par[2],atoi(par[6].c_str()), atoi(par[7].c_str()) );
}
