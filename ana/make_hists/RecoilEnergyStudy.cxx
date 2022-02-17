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

vector<double> neutronConeAngles({5,10,15,20,25,30,35,40,45});

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

  for( int i = 0; i< n_entries; ++i)
  {
    evt->GetEntry(i);
    if( IgnoreEvent( evt, isData )) continue;
    if (i%5000 == 0) cout<<"At Event "<<i<<" -- "<<float(i)/n_entries*100<<"%s"<<endl;
    //==============Fill Neutron Blobs======================
    bool isMC = !isData;
    ncutter->UpdateCandidates( evt, isMC );
    if(isDebug) cout<<"Updated ncutter"<<endl;
    //==============Event Selection========================
    bool useRecoil = false;
    bool PassEventCut = (EventCut( evt, cutter, ncutter, cuts_counter, multiplicity, sample ,  isData, npMode,useRecoil ));
    if(isDebug) cout<<PassEventCut<<endl;


    //now add in neutron-->
    std::vector<bool> commonNeutronCuts = CommonNeutronCuts(evt,ncutter);
    bool hasNC = commonNeutronCuts[0] && commonNeutronCuts[1];
    bool isForward = commonNeutronCuts[2] && hasNC;
    bool isBackward = commonNeutronCuts[3] && hasNC;
    bool hasProton = (npMode!=0);


    if (!PassEventCut) continue;
    NeutronCandidates* ptrNeutronCandidates; 
    if (npMode != 1 || hasNC) ptrNeutronCandidates = ncutter->GetCandidates();
    bool fill_common_histos = true;
    bool passRecoil = cutter->passSignalFunction( evt, 0, 0 ) ;

    //===============Reconstruction ========================
    if(isDebug) cout<<"Getting CV Weight:"<<endl;
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

    //Neutron Variables
    TVector3 neutronDir(-1,-1,-1);
    XYZVector blobVtxDir(-1,-1,-1);
    if (npMode !=1 || hasNC)
    {
      blobVtxDir = ncutter->MainNeutronCandidate()->BlobVtxDir;
      neutronDir.SetXYZ(blobVtxDir.X(),blobVtxDir.Y(), blobVtxDir.Z() ) ; 
    }
    //Expected nucleon kinematics

    // 2. binding and is/fs mass
    double binding_e = 0; 
    int charge = 1;
    double is_part_mass = (neutrinoMode)? 0.9395654133: 0.93827231; //neutron/proton
    double fs_part_mass = (neutrinoMode)? 0.93827231: 0.9395654133;

    // 3. expected nucleon
    XYZTVector expectedNucleon4P = geoUtils->ComputeExpectedNucleon( beam, reco_muon_4P, is_part_mass, fs_part_mass, binding_e);
    XYZVector expectedNucleon3P = (XYZVector) expectedNucleon4P.Vect();

    TVector3 primaryNucleonVect =  neutronDir.Unit()*expectedNucleon3P.R();

    XYZVector primaryNucleonVect_XYZ( primaryNucleonVect.X(), primaryNucleonVect.Y(), primaryNucleonVect.Z() );

    //Calculating Q2qe:
    double q2qe = physicsUtils->qsqCCQE( beam, reco_muon_4P, charge, binding_e ); //GeV^2
    double enu = physicsUtils->nuEnergyCCQE( beam, reco_muon_4P, charge, binding_e );

    //Calculate recoil
      //Cone Neutron energy
    XYZVector muonDirWRTDetector = (XYZVector) reco_muon_3P.Unit();
    XYZVector vtx(evt->vtx[0], evt->vtx[1], evt->vtx[2]);
    XYZVector blobVtxDirTemp = ncutter->MainNeutronCandidate()->BlobVtxDir.Unit();

    Cone muonCone( vtx, muonDirWRTDetector, 15./180*3.14);
    double muon_E = 0, muon_E_2D=0, muon_E_3D=0;

    //define neutron cones
    //vector<double> neutronConeAngles({5,10,15,20,25,40,40,45});
    int nAngles = neutronConeAngles.size();
    vector<Cone> neutronCones;
    for(auto angle : neutronConeAngles) neutronCones.push_back( Cone( vtx, blobVtxDirTemp, angle/180*3.1415 ) );
    vector<double> neutronConeEnergy(nAngles,0);
    vector<double> neutronConeEnergy2D(nAngles,0);
    vector<double> neutronConeEnergy3D(nAngles,0);

    vector<NeutronBlob*> all_blobs = ncutter->GetCandidates()->AllBlobs;

    NeutronBlob* mainCandidate = ncutter->MainNeutronCandidate();
    double mcDist = (mainCandidate->BlobBegPos.Vect() - vtx ).R();
    NeutronBlob* newMainCandidate = NULL;
    double newE = 0;
    for ( auto itBlob = all_blobs.begin(); itBlob != all_blobs.end(); ++itBlob )
    {
      XYZVector blobPos=(XYZVector) (*itBlob)->BlobBegPos.Vect();
      double E = (*itBlob)->TotalE;
      if (E==0) continue;
      int view = (*itBlob)->View;
      if (view >3 || view < 1) continue;

      if ((*itBlob)->Is3D)
      { 
        if (muonCone.InsideCone( blobPos )) { muon_E+=E; muon_E_3D+=E;continue; }
        bool inCone = false;
        int i = 0;
        while( !inCone && i <nAngles)
        {
            inCone = neutronCones[i].InsideCone( blobPos );
            if(!inCone) i++;
        }
        int I = i;
        for( ; i<nAngles; i++ ){
          neutronConeEnergy[i]+=E;
          neutronConeEnergy3D[i]+=E;
        }


        int dist = (blobPos - vtx).R();
        if(
            I<= 5 && 
            dist< mcDist &&
            dist > 100 &&
            E > newE &&
            E != mainCandidate->TotalE
          ) {newMainCandidate = *itBlob;newE=E;}
      }
      else
      {
        XYZVector tpos = (*itBlob)->BlobPosT;
        if (muonCone.InsideCone( view, tpos)) { muon_E+=E; muon_E_2D+=E;continue; }
        bool inCone = false;
        int i = 0;
        while( !inCone && i <nAngles)
        {
            inCone = neutronCones[i].InsideCone( view, tpos );
            if(!inCone) i++;
        }
        for( ; i<nAngles; i++ ){
            neutronConeEnergy[i]+=E;
            neutronConeEnergy2D[i]+=E;
        }
      }
    }


    if (newMainCandidate == NULL) newMainCandidate = mainCandidate;

    double mainNC_E = newMainCandidate->TotalE;

    //recoil variants MeV
    double recoil =ncutter->GetDefaultRecoilEnergy(evt);
    double recoil_inc = evt->recoil_summed_energy[0];
    vector<double> recoil_neutronCones(nAngles, recoil );
    for( int i = 0; i< nAngles; i++)  recoil_neutronCones[i]-=neutronConeEnergy[i];
    double recoil_inc_NoMC = recoil_inc - ncutter->MainNeutronCandidate()->TotalE;
    double recoil_neutron_3D = ncutter->GetNCRecoilEnergy(3);
    double recoil_neutron_2D = ncutter->GetNCRecoilEnergy(2);
    double recoil_neutron_All = ncutter->GetNCRecoilEnergy(1);

    double recoil_all_nonMuon = evt->recoil_energy_nonmuon_nonvtx100mm;
    //cout<<recoil<<endl;
    //cout<<neutronConeEnergy[3]<<endl;


    //calculating angular variables
    std::vector<double> reco_neutron_vars = geoUtils->ComputeNeutronAngularVars( beam, expectedNucleon3P, primaryNucleonVect_XYZ );
    double reco_dthetaP = reco_neutron_vars[vdthetaP];
    double reco_dthetaR = reco_neutron_vars[vdthetaR];
    double reco_dtheta = reco_neutron_vars[vdtheta];
    //if (isMC) cout<< evt->mc_intType<<endl;
    reco_dthetaP *= 180/3.14159;
    reco_dthetaR *= 180/3.14159;
    reco_dtheta *= 180/3.14159;

    XYZVector updatedNucleonVect = (XYZVector) (newMainCandidate->BlobBegPos.Vect() - vtx);
    std::vector<double> reco_updated_neutron_vars = geoUtils->ComputeNeutronAngularVars( beam, expectedNucleon3P, updatedNucleonVect );
    double reco_updated_dthetaP = reco_updated_neutron_vars[vdthetaP];
    double reco_updated_dthetaR = reco_updated_neutron_vars[vdthetaR];
    double reco_updated_dtheta = reco_updated_neutron_vars[vdtheta];
    //if (isMC) cout<< evt->mc_intType<<endl;
    reco_updated_dthetaP *= 180/3.14159;
    reco_updated_dthetaR *= 180/3.14159;
    reco_updated_dtheta *= 180/3.14159;


    //cout<<"dpTx - dPp = "<<abs(reco_dptx) - abs(reco_dPp)<<endl;

    //I feel there is also a need  to include the left-right uncertainty of the beam, since LR-direction is pretty important....

    if (fill_common_histos)
    {
      utils->fillHistosV3(utils->histos2D["h_dthetaP_dthetaR"], reco_dthetaP, reco_dthetaR, isData, evt, wgt );
      utils->fillHistosV3(utils->histos2D["h_dthetaP_dthetaR_fine"], reco_dthetaP, reco_dthetaR, isData, evt, wgt );
      utils->fillHistosV3(utils->histos2D["h_dthetaP_dthetaR_update"], reco_updated_dthetaP, reco_updated_dthetaR, isData, evt, wgt );
      utils->fillHistosV3(utils->histos2D["h_dthetaP_dthetaR_fine_update"], reco_updated_dthetaP, reco_updated_dthetaR, isData, evt, wgt );


      utils->fillHistosV3(utils->histos2D["h_q2qe_recoil"], q2qe, recoil, isData, evt, wgt );
      if(passRecoil) utils->fillHistosV3(utils->histos2D["h_q2qe_recoil_cut"], q2qe, recoil, isData, evt, wgt );
      for( int i = 0; i<nAngles;i++) 
      {
        utils->fillHistosV3(utils->histos2D[Form("h_q2qe_recoil_cone_%02d", 5+i*5)], q2qe, recoil_neutronCones[i], isData, evt, wgt );
        if(passRecoil) utils->fillHistosV3(utils->histos2D[Form("h_q2qe_recoil_cone_%02d_cut", 5+i*5)], q2qe, recoil_neutronCones[i], isData, evt, wgt );
      }




    } // End of fill_common_histos
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
  vector<double> bins;
  axis_binning uniformbins;
  bins.clear();
  for(float i=0;i<=1000;i+=5.) bins.push_back(i);
  uniformbins.bin_edges = bins;
  uniformbins.nbins = bins.size()-1;
  uniformbins.min = bins.front();
  uniformbins.max = bins.back();
 
  //2-D
  MnvH2D* h_q2qe_recoil[nHistos],* h_q2qe_recoil_cut[nHistos];
  MnvH2D* h_q2qe_recoil_cone_05[nHistos], *h_q2qe_recoil_cone_05_cut[nHistos];
  MnvH2D* h_q2qe_recoil_cone_10[nHistos], *h_q2qe_recoil_cone_10_cut[nHistos];
  MnvH2D* h_q2qe_recoil_cone_15[nHistos], *h_q2qe_recoil_cone_15_cut[nHistos];
  MnvH2D* h_q2qe_recoil_cone_20[nHistos], *h_q2qe_recoil_cone_20_cut[nHistos];
  MnvH2D* h_q2qe_recoil_cone_25[nHistos], *h_q2qe_recoil_cone_25_cut[nHistos];
  MnvH2D* h_q2qe_recoil_cone_30[nHistos], *h_q2qe_recoil_cone_30_cut[nHistos];
  MnvH2D* h_q2qe_recoil_cone_35[nHistos], *h_q2qe_recoil_cone_35_cut[nHistos];
  MnvH2D* h_q2qe_recoil_cone_40[nHistos], *h_q2qe_recoil_cone_40_cut[nHistos];
  MnvH2D* h_q2qe_recoil_cone_45[nHistos], *h_q2qe_recoil_cone_45_cut[nHistos];

  utils->bookHistos( h_q2qe_recoil, "h_q2qe_recoil","q2qe:recoil", Q2bins, uniformbins ); utils->bookHistos( h_q2qe_recoil_cut, "h_q2qe_recoil_cut","q2qe:recoil", Q2bins, uniformbins );
  utils->bookHistos( h_q2qe_recoil_cone_05, "h_q2qe_recoil_cone_05","q2qe:recoil", Q2bins, uniformbins ); utils->bookHistos( h_q2qe_recoil_cone_05_cut, "h_q2qe_recoil_cone_05_cut","q2qe:recoil", Q2bins, uniformbins );
  utils->bookHistos( h_q2qe_recoil_cone_10, "h_q2qe_recoil_cone_10","q2qe:recoil", Q2bins, uniformbins ); utils->bookHistos( h_q2qe_recoil_cone_10_cut, "h_q2qe_recoil_cone_10_cut","q2qe:recoil", Q2bins, uniformbins );
  utils->bookHistos( h_q2qe_recoil_cone_15, "h_q2qe_recoil_cone_15","q2qe:recoil", Q2bins, uniformbins ); utils->bookHistos( h_q2qe_recoil_cone_15_cut, "h_q2qe_recoil_cone_15_cut","q2qe:recoil", Q2bins, uniformbins );
  utils->bookHistos( h_q2qe_recoil_cone_20, "h_q2qe_recoil_cone_20","q2qe:recoil", Q2bins, uniformbins ); utils->bookHistos( h_q2qe_recoil_cone_20_cut, "h_q2qe_recoil_cone_20_cut","q2qe:recoil", Q2bins, uniformbins );
  utils->bookHistos( h_q2qe_recoil_cone_25, "h_q2qe_recoil_cone_25","q2qe:recoil", Q2bins, uniformbins ); utils->bookHistos( h_q2qe_recoil_cone_25_cut, "h_q2qe_recoil_cone_25_cut","q2qe:recoil", Q2bins, uniformbins );
  utils->bookHistos( h_q2qe_recoil_cone_30, "h_q2qe_recoil_cone_30","q2qe:recoil", Q2bins, uniformbins ); utils->bookHistos( h_q2qe_recoil_cone_30_cut, "h_q2qe_recoil_cone_30_cut","q2qe:recoil", Q2bins, uniformbins );
  utils->bookHistos( h_q2qe_recoil_cone_35, "h_q2qe_recoil_cone_35","q2qe:recoil", Q2bins, uniformbins ); utils->bookHistos( h_q2qe_recoil_cone_35_cut, "h_q2qe_recoil_cone_35_cut","q2qe:recoil", Q2bins, uniformbins );
  utils->bookHistos( h_q2qe_recoil_cone_40, "h_q2qe_recoil_cone_40","q2qe:recoil", Q2bins, uniformbins ); utils->bookHistos( h_q2qe_recoil_cone_40_cut, "h_q2qe_recoil_cone_40_cut","q2qe:recoil", Q2bins, uniformbins );
  utils->bookHistos( h_q2qe_recoil_cone_45, "h_q2qe_recoil_cone_45","q2qe:recoil", Q2bins, uniformbins ); utils->bookHistos( h_q2qe_recoil_cone_45_cut, "h_q2qe_recoil_cone_45_cut","q2qe:recoil", Q2bins, uniformbins );
  //2-D
  axis_binning   dthetaPdthetaRBins;
  axis_binning   dthetaPdthetaRFineBins;
  //MnvH2D *h_dthetaRq2qe_ptmu[nHistos], 
  //       *h_dthetaPq2qe_ptmu[nHistos],
  //       *h_tpq2qe_ptmu[nHistos],
  MnvH2D   *h_dthetaP_dthetaR[nHistos];
  MnvH2D   *h_dthetaP_dthetaR_fine[nHistos];
  MnvH2D   *h_dthetaP_dthetaR_update[nHistos];
  MnvH2D   *h_dthetaP_dthetaR_fine_update[nHistos];
  utils->bookHistos( h_dthetaP_dthetaR, "h_dthetaP_dthetaR", "dthetaP:dthetaR",dthetaPerpbins, dthetaReactbins );
  utils->bookHistos( h_dthetaP_dthetaR_fine, "h_dthetaP_dthetaR_fine","dthetaP:dthetaR", dthetaFineBins, dthetaFineBins );
  utils->bookHistos( h_dthetaP_dthetaR_update, "h_dthetaP_dthetaR_update","dthetaP:dthetaR", dthetaPerpbins, dthetaReactbins );
  utils->bookHistos( h_dthetaP_dthetaR_fine_update, "h_dthetaP_dthetaR_fine_update","dthetaP:dthetaR", dthetaFineBins, dthetaFineBins );



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
    //data_events = AnalysisLooper( data, utils, cutter, ncutter, sample,  multiplicity, entries_data, data_counter, isData, isTruth, npMode );
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
