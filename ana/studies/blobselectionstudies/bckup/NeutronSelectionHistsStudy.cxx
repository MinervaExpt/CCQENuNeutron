#include "include/CCQENuUtils.h"
#include "include/CommonBins.h"
#include "include/EventCuts.h"
#include "TParameter.h" 
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"

using namespace CCQENU_ANA;
using namespace Neutron_ANA;
using namespace ROOT::Math;

bool DEBUG=false;
double MuonCylinderCut = 100;//mm
double NeutronZCut = 100;//mm
double MuonConeAngleDeg = 15;//degree

DetectorUtils* detUtils = new DetectorUtils();
double AnalysisLooper( CCQENuEvent* evt, CCQENuUtils* utils, CCQENuCuts *cutter, NeutronBlobCuts* ncutter, string sample, int multiplicity, int n_entries, bool isMC = true )
{
  cout<<"Start AnalysisLooper"<<endl;
  GeoUtils * geoUtils = ncutter->GeoUtilsPtr();
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
    
    map<int,int> cuts_counter;
    bool pass_all = EventCutBase( evt, cutter, ncutter, cuts_counter, !isMC, 1);
    bool pass_nc = EventCutBase( evt, cutter, ncutter, cuts_counter, !isMC, 1, -1,-1,-1,-1, -1, 1,1,-1);
    bool pass_nc_forward = EventCutBase( evt, cutter, ncutter, cuts_counter, !isMC, 1, -1,-1,-1,-1, -1, 1,1,1);
    bool pass_nc_backward = pass_nc && !pass_nc_forward;

    // event selection
    //if (sample == "Signal" )
    //{
    //  if( ! EventSelectionDataMC( evt, cutter, ncutter, true,isMC, multiplicity ) ) continue;
    //}
    //else if (sample == "MichelSideBand")
    //{
    //  if( ! EventSelectionDataMC( evt, cutter, ncutter, false,isMC, multiplicity, false,true, true,true ) ) continue;
    //}
    //else if (sample == "BlobSideBand" )
    //{
    //  if( ! EventSelectionDataMC( evt, cutter, ncutter, false,isMC, multiplicity, true,false, true,true ) ) continue;
    //}
    //else if (sample == "MicBlobSideBand")
    //{
    //  if( ! EventSelectionDataMC( evt, cutter, ncutter, false,isMC, multiplicity, false,false, true,true ) ) continue;
    //}
    //else if( sample == "SideBand" )
    //{
    //  if (multiplicity==0)
    //  {
    //    if( !cutter->passCCQESelection(evt,"Sideband") )    continue;
    //  }
    //  else
    //  {
    //    if(!EventSelectionDataMC(evt, cutter, ncutter, true, isMC, multiplicity) &&!cutter->passSidebandFunction( evt, 0, 0 )) continue;
    //  }
    //}
    //else {
    //  cout<<"Wrong Sample Selection."<<endl;
    //  cout<<"Select: Signal, MichelSideBand, SideBand"<<endl;
    //  return 0;
    //}
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

    //recoil variants
    double recoil = evt->nonvtx_iso_blobs_energy + evt->dis_id_energy;
    double recoil_inc = evt->recoil_summed_energy[0];

    //vtx energy
    double vtx_energy = evt->vtx_blobs_energy;

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
    // Perpendicular and reaction plane angles
    std::vector<double> neutAngularVars = geoUtils->ComputeNeutronAngularVars( beamAxis, expectedFSPp,blobVtxDir );
    // Inferred neutron momentum
    XYZVector inferredNeutP = blobVtxDir.Unit() * expectedFSP.R();
    XYZTVector inferredNeut4P( inferredNeutP.X(), inferredNeutP.Y(), inferredNeutP.Z(), expectedFSP.E() );

    double PerpAngle = neutAngularVars[0], ReactAngle = neutAngularVars[1], dTheta = neutAngularVars[2];
    double Q2 = evt->CCQENu_Q2/1e6 ;


    //---------------------------------
    // Transverse Variables //only phiT for neutrons
    //---------------------------------
    double transverseEb = geoUtils->EbAntinuGeV() + geoUtils->MnGeV() - geoUtils->MpGeV(); //make Remnant Nucleon Mass compliant with anti neutrino
    std::vector<double> neutTransverseVars = geoUtils->GetTransverseVariables( beamAxis, muon, inferredNeut4P, geoUtils->C12GeV(), transverseEb );
    double Enu=      neutTransverseVars[vEnu];
    double Pn=       neutTransverseVars[vPn];
    double dpT=      neutTransverseVars[vdpT];
    double dpTx=     neutTransverseVars[vdpTx];
    double dpTy=     neutTransverseVars[vdpTy];
    double dalphaT=  neutTransverseVars[vdalphaT];
    double phiTDeg=     neutTransverseVars[vdphiT];
    double phiTRadian= phiTDeg/180*3.141592654;

    

    //-------------------------------------------------------------------
    // Check that the CV of events pass the proton and recoil cuts 
    // This is because we have systematics due to these cuts 
    // and a change in their values can affect the composition 
    // of the selected sample. 
    //-------------------------------------------------------------------
    // only for MC
    bool doesCVpass_SingleProtonCut = (evt->multiplicity==1)? true : cutter->passSingleProtonCut( evt, 0, 0 );
    bool doesCVpass_ExtraProtonsCut = (evt->multiplicity==1)? true : cutter->passExtraProtonsCut( evt, NULL, 0 );
    bool doesCVpass_RecoilCut = cutter->passSignalFunction( evt, 0, 0 );   

    bool fillHisto = (isMC)? ( doesCVpass_SingleProtonCut && doesCVpass_ExtraProtonsCut && doesCVpass_RecoilCut ) : true;

    if (DEBUG) cout<<"pre-fill"<<endl;
    if( fillHisto )
    {
      /*
      if(DEBUG) cout<<"fillHisto"<<endl;
      utils->fillHistos( utils->histos1D["h_theta"], muon_theta, isMC, evt, wgt );
      if(DEBUG) cout<<"fillHisto Done 1D 0"<<endl;
      utils->fillHistos( utils->histos1D["h_tmu"], muon_T, isMC, evt, wgt );
      utils->fillHistos( utils->histos1D["h_ptmu"], muon_pt_beam, isMC, evt, wgt );
      utils->fillHistos( utils->histos1D["h_pzmu"], muon_pz_beam, isMC, evt, wgt );
      if(DEBUG) cout<<"fillHisto Done 1D 2"<<endl;
      utils->fillHistos( utils->histos1D["h_multiplicity"], evt->multiplicity, isMC, evt, wgt );
      utils->fillHistos( utils->histos1D["h_w"],reco_w, isMC, evt, wgt );
      if(DEBUG) cout<<"fillHisto Done 1D 3"<<endl;
      utils->fillHistos( utils->histos1D["h_vtx_energy"], vtx_energy, isMC, evt, wgt );
      utils->fillHistos( utils->histos1D["h_vtxrecoilratio"], vtx_energy_fraction, isMC, evt, wgt);
      if(evt->CCQENu_Q2/1e6 < 0.2) utils->fillHistos( utils->histos1D["h_vtxrecoilratio_lowq2"], vtx_energy_fraction, isMC, evt, wgt);
      else utils->fillHistos( utils->histos1D["h_vtxrecoilratio_highq2"], vtx_energy_fraction, isMC, evt, wgt);
      // Uncomment out if you need to make this plot   
      // utils->fillHistos( h_trkChi2PerDoF, mnv_trk_chi2PerDoF ); 
      
      if(DEBUG) cout<<"fillHisto Done 1D"<<endl;
      //2-D
      utils->fillHistos( utils->histos2D["h_tmu_theta"], muon_T, muon_theta, isMC, evt, wgt );
      utils->fillHistos( utils->histos2D["h_tmu_costheta"], muon_T, cos_muon_theta, isMC, evt, wgt );
      utils->fillHistos( utils->histos2D["h_pmu_ptmu"], muon_p, muon_pt_beam, isMC, evt, wgt );
      utils->fillHistos( utils->histos2D["h_pzmu_ptmu"], muon_pz_beam, muon_pt_beam, isMC, evt, wgt );
      utils->fillHistos( utils->histos2D["h_enu_ptmu"], evt->CCQENu_enu_muon/1e3, muon_pt_beam, isMC, evt, wgt );
      utils->fillHistos( utils->histos2D["h_q2_ptmu"], evt->CCQENu_Q2/1e6, muon_pt_beam, isMC, evt, wgt );
      utils->fillHistos( utils->histos2D["h_recoil_inc_ptmu"], recoil_inc, muon_pt_beam, isMC, evt, wgt );
      utils->fillHistos( utils->histos2D["h_recoil_ptmu"], recoil, muon_pt_beam , isMC, evt, wgt);
      utils->fillHistos( utils->histos2D["h_recoil_ptmu_before"], recoil, muon_pt_beam , isMC, evt, wgt);
      utils->fillHistos( utils->histos2D["h_Ediff_ptmu"], nu_qe_E-nu_cal_E, muon_pt_beam , isMC, evt, wgt);
      utils->fillHistos( utils->histos2D["h_vtxrecoilratio_pt"],vtx_energy_fraction,muon_pt_beam, isMC, evt, wgt);      

      if(evt->CCQENu_Q2/1e6 < 0.2){
        utils->fillHistos( utils->histos2D["h_vtx_old_low_q2_300mm"], evt->recoil_energy_nonmuon_vtx300mm, muon_pt_beam, isMC, evt, wgt);
        utils->fillHistos( utils->histos2D["h_vtx_old_low_q2_250mm"], evt->recoil_energy_nonmuon_vtx250mm, muon_pt_beam, isMC, evt, wgt);
        utils->fillHistos( utils->histos2D["h_vtx_old_low_q2_200mm"], evt->recoil_energy_nonmuon_vtx200mm, muon_pt_beam, isMC, evt, wgt);
        utils->fillHistos( utils->histos2D["h_vtx_old_low_q2_150mm"], evt->recoil_energy_nonmuon_vtx150mm, muon_pt_beam, isMC, evt, wgt);
        utils->fillHistos( utils->histos2D["h_vtx_new_low_q2"], vtx_energy, muon_pt_beam, isMC, evt, wgt);
      }
      else{
        utils->fillHistos( utils->histos2D["h_vtx_old_high_q2_300mm"], evt->recoil_energy_nonmuon_vtx300mm, muon_pt_beam, isMC, evt, wgt);
        utils->fillHistos( utils->histos2D["h_vtx_old_high_q2_250mm"], evt->recoil_energy_nonmuon_vtx250mm, muon_pt_beam, isMC, evt, wgt);
        utils->fillHistos( utils->histos2D["h_vtx_old_high_q2_200mm"], evt->recoil_energy_nonmuon_vtx200mm, muon_pt_beam, isMC, evt, wgt);
        utils->fillHistos( utils->histos2D["h_vtx_old_high_q2_150mm"], evt->recoil_energy_nonmuon_vtx150mm, muon_pt_beam, isMC, evt, wgt);
        utils->fillHistos( utils->histos2D["h_vtx_new_high_q2"], vtx_energy, muon_pt_beam, isMC, evt, wgt);
      }
      utils->fillHistos( utils->histos2D["h_vtx_old_300mm"], evt->recoil_energy_nonmuon_vtx300mm, muon_pt_beam, isMC, evt, wgt);
      utils->fillHistos( utils->histos2D["h_vtx_old_250mm"], evt->recoil_energy_nonmuon_vtx250mm, muon_pt_beam, isMC, evt, wgt);
      utils->fillHistos( utils->histos2D["h_vtx_old_200mm"], evt->recoil_energy_nonmuon_vtx200mm, muon_pt_beam, isMC, evt, wgt);
      utils->fillHistos( utils->histos2D["h_vtx_old_150mm"], evt->recoil_energy_nonmuon_vtx150mm, muon_pt_beam, isMC, evt, wgt);
      utils->fillHistos( utils->histos2D["h_vtx_new"], vtx_energy, muon_pt_beam, isMC, evt, wgt);
      */
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
        // Neutron Blob PDG
        //------------------------------------------------------
        double pdg_maincandidate =  ncutter->MainNeutronCandidate()->MCPID;
        double nblobs = ncutter->GetCandidates()->AllBlobs.size();

        utils->fillHistos( utils->histos2D["h_neutron_truth_mcpid_nblobs"], pdg_maincandidate, nblobs, isMC, evt, wgt);
        utils->fillHistos( utils->histos2D["h_neutron_truth_dreactplane_nblobs"], ReactAngle, nblobs, isMC, evt, wgt);
        utils->fillHistos( utils->histos2D["h_neutron_truth_dreactplane_mcpid"], ReactAngle, pdg_maincandidate, isMC, evt, wgt);

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
    if( isMC && doesCVpass_SingleProtonCut && doesCVpass_ExtraProtonsCut )
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
  NeutronBlobCuts *ncutter = new NeutronBlobCuts();

  cout<<"BOOK HISTOGRAMS"<<endl;
  //---------------------------------------------
  // Book Histograms
  //---------------------------------------------
  //Neutrons

  //Histograms we need: in bins of recoil energy, ptmu, ptz
  //recoil energy, number of blobs or number of tracks
  //recoil energy-MC energy, number of blobs/tracks
  //proton score cut. 
  //Define recoil = total blob energies
  MnvH3D* h_recoil_ptmu_ptz_0[nHistos],h_recoil_ptmu_ptz_1[nHistos];
  MnvH2D* h_recoil_nBlobs_0[nHistos],h_recoil_nTracks_0[nHistos];
  MnvH2D* h_recoil_nBlobs_1[nHistos],h_recoil_nTracks_1[nHistos];
  TH2D* h_pid_score = new TH2D("h_score_pid","h_score_pid", 4,0,4,100,-10,10);   //nucleon: 1, pion+-: 2, pion0:3, em:4

 utils->bookHistos( h_recoil_ptmu_ptz_0,"h_recoil_ptmu_ptz_0","recoil_ptmu_ptz0", recoilbins3D, ptbins3D, pzbins3D );
 utils->bookHistos( h_recoil_ptmu_ptz_1,"h_recoil_ptmu_ptz_1","recoil_ptmu_ptz1", recoilbins3D, ptbins3D, pzbins3D );

 utils->bookHistos( h_recoil_nBlobs_0, "h_recoil_nBlobs_0", "recoil_nBlobs0", recoilbins3D, nDiscreteBins );
 utils->bookHistos( h_recoil_nBlobs_1, "h_recoil_nBlobs_1", "recoil_nBlobs0", recoilbins3D, nDiscreteBins );
 utils->bookHistos( h_recoil_nTracks_0, "h_recoil_nTracks_0", "recoil_nTracks0", recoilbins3D, nDiscreteBins );
 utils->bookHistos( h_recoil_nTracks_1, "h_recoil_nTracks_1", "recoil_nTracks0", recoilbins3D, nDiscreteBins );

   // Neutron Truth Plots:

  //MnvH2D *h_neutron_truth_mcpid_nblobs[nHistos];
  //utils->bookHistos( h_neutron_truth_mcpid_nblobs, "h_neutron_truth_mcpid_nblobs","Neutron Candidate MCPID vs NBlob: PID : NBlobs", neut_mcpdg, neut_nblobs );

  //MnvH2D *h_neutron_truth_dreactplane_nblobs[nHistos];
  //utils->bookHistos( h_neutron_truth_dreactplane_nblobs, "h_neutron_truth_dreactplane_nblobs","Neutron Candidate #theta_{React} vs NBlob: #theta_{React} : NBlobs", neutthetabins, neut_nblobs );

  //MnvH2D *h_neutron_truth_dreactplane_mcpid[nHistos];
  //utils->bookHistos( h_neutron_truth_dreactplane_mcpid, "h_neutron_truth_dreactplane_mcpid","Neutron Candidate #theta_{React} vs MCPID: #theta_{React} : NMCPID", neutthetabins, neut_mcpdg );
  //
  ////Other diagnostic plots
  //MnvH1D *h_neutron_truth_mcpid_allblobs[nHistos],*h_neutron_truth_mcpid_2Dblobs[nHistos],*h_neutron_truth_mcpid_3Dblobs[nHistos];
  //utils->bookHistos(h_neutron_truth_mcpid_allblobs,"h_neutron_truth_mcpid_allblobs","",neut_mcpdg);
  //utils->bookHistos(h_neutron_truth_mcpid_2Dblobs,"h_neutron_truth_mcpid_2Dblobs","",neut_mcpdg);
  //utils->bookHistos(h_neutron_truth_mcpid_3Dblobs,"h_neutron_truth_mcpid_3Dblobs","",neut_mcpdg);

  //MnvH1D *h_neutron_truth_2DblobE_total[nHistos],*h_neutron_truth_2DblobE_max[nHistos],*h_neutron_truth_2DblobE_avg[nHistos];
  //utils->bookHistos(h_neutron_truth_2DblobE_total,"h_neutron_truth_2DblobE_total","",neut_mcpdg);

  //MnvH1D *h_neutron_truth_3DblobE_total[nHistos],*h_neutron_truth_3DblobE_max[nHistos],*h_neutron_truth_3DblobE_avg[nHistos];
  //MnvH2D *h_neutron_truth_angleMuon_mcpid_3Dblobs[nHistos];
  //MnvH2D *h_neutron_truth_angleMuon_mcpid_mainCandidate[nHistos];

  //3-D
  //MnvH3D *h_enucal_enulep_pt[nHistos];

  //utils->bookHistos( h_enucal_enulep_pt, "h_enucal_enulep_pt", "enucal enulep pt; enucal;enulep;pt",fineEnu,fineEnu,muonPtbins);

  cout<<"END BOOK HISTOGRAMS"<<endl;
  //--------------------------------------------------------------
  // Add Vertical and Lateral Error Bands
  // JO is removing error bands from all 1D because they are all 
  // present in 2D. Plz make projections directly from the 2D. 
  //--------------------------------------------------------------
  if(RunCodeWithSystematics){
    cout<<"ADD SYS. HISTO"<<endl;
    
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
        
    /*
    h_neutron_dperp_dreactplane[i]->Write(); 
    h_neutron_dtheta_qsq[i]->Write(); 
    h_neutron_dperp_qsq[i]->Write(); 
    h_neutron_dreactplane_qsq[i]->Write();

    h_neutron_phit_qsq[i]->Write(); 
    h_neutron_phit_dtheta[i]->Write(); 
    h_neutron_phit_dperp[i]->Write(); 
    h_neutron_phit_dreactplane[i]->Write();
    */
    h_neutron_truth_dreactplane_mcpid[i]->Write();
    h_neutron_truth_dreactplane_nblobs[i]->Write();
    h_neutron_truth_mcpid_nblobs[i]->Write();
    
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
