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


int testEvent = -1;
bool testMode=false;

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
  if(isData && testMode) return 0;
  GeoUtils* geoUtils = new GeoUtils();
  PhysicsUtils* physicsUtils = new PhysicsUtils();
  KinematicsModifier* kineMod = new KinematicsModifier();

  double n_evts=0;
  //return 0;

  for( int i = 0; i< n_entries; ++i)
  {
    evt->GetEntry(i);
    if( IgnoreEvent( evt, isData )) continue;
    if (i%1000 == 0 && !testMode) cout<<"At Event "<<i<<" ("<<(i+1)*100./float(n_entries)<<"%)"<<endl;
    //==============Fill Neutron Blobs======================
    bool isMC = !isData;
    if(testMode)
    {
      if( !cutter->passTrueCCQE( evt ) ) continue;
      if( evt->mc_targetZ != 1 ) continue;
    }
    
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

    //bool PassEventCut = (EventCut( evt, cutter, ncutter, cuts_counter, multiplicity, sample ,  isData, npMode, useRecoilCut, orderPar));
    int pass_signalFunc = -1;
    int pass_michel = -1;
    int pass_nblobs = -1;
    int pass_blobHitEnergy = -1;
    int pass_singleProton = -1;
    int pass_extraProtons = -1;
    int pass_neutronCut = -1;
    int pass_neutronMuonCut = -1;
    int pass_neutronVtxCut = -1;
    int pass_neutronClusMaxE = -1;
    
    bool PassEventCut = EventCutBase( evt, cutter, ncutter, cuts_counter, isData, multiplicity, pass_signalFunc , pass_michel , pass_nblobs , pass_blobHitEnergy, pass_singleProton, pass_extraProtons , pass_neutronCut, pass_neutronMuonCut, pass_neutronVtxCut, pass_neutronClusMaxE, orderPar );

    if(isDebug) cout<<PassEventCut<<endl;

    if (!PassEventCut) continue;

    pass_signalFunc = 1;
    bool pass_side = EventCutBase( evt, cutter, ncutter, cuts_counter, isData, multiplicity, pass_signalFunc , pass_michel , pass_nblobs , pass_blobHitEnergy, pass_singleProton, pass_extraProtons , pass_neutronCut, pass_neutronMuonCut, pass_neutronVtxCut, pass_neutronClusMaxE, orderPar );

    pass_nblobs = 1;
    pass_blobHitEnergy = 1;
    bool pass_side_blob = EventCutBase( evt, cutter, ncutter, cuts_counter, isData, multiplicity, pass_signalFunc , pass_michel , pass_nblobs , pass_blobHitEnergy, pass_singleProton, pass_extraProtons , pass_neutronCut, pass_neutronMuonCut, pass_neutronVtxCut, pass_neutronClusMaxE, orderPar );

    pass_neutronCut = 1;
    pass_neutronMuonCut = 1;
    bool pass_side_blob_neutron = EventCutBase( evt, cutter, ncutter, cuts_counter, isData, multiplicity, pass_signalFunc , pass_michel , pass_nblobs , pass_blobHitEnergy, pass_singleProton, pass_extraProtons , pass_neutronCut, pass_neutronMuonCut, pass_neutronVtxCut, pass_neutronClusMaxE, orderPar );



    pass_nblobs = -1;
    pass_blobHitEnergy = -1;
    bool pass_side_neutron = EventCutBase( evt, cutter, ncutter, cuts_counter, isData, multiplicity, pass_signalFunc , pass_michel , pass_nblobs , pass_blobHitEnergy, pass_singleProton, pass_extraProtons , pass_neutronCut, pass_neutronMuonCut, pass_neutronVtxCut, pass_neutronClusMaxE, orderPar );


    pass_signalFunc = -1;
    pass_nblobs = 1;
    pass_blobHitEnergy = 1;
    pass_nblobs = 1;
    pass_blobHitEnergy = 1;
    bool pass_blob_neutron = EventCutBase( evt, cutter, ncutter, cuts_counter, isData, multiplicity, pass_signalFunc , pass_michel , pass_nblobs , pass_blobHitEnergy, pass_singleProton, pass_extraProtons , pass_neutronCut, pass_neutronMuonCut, pass_neutronVtxCut, pass_neutronClusMaxE, orderPar );
    
    bool fill_common_histos = true;
    //===============Reconstruction ========================
    //-----------
    //Vertex


   XYZTVector reco_muon_XYZT(evt->CCQENu_leptonE[0]/1000.,
                     evt->CCQENu_leptonE[1]/1000.,
                     evt->CCQENu_leptonE[2]/1000.,
                     evt->CCQENu_leptonE[3]/1000.);

   
    
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


    //Expected nucleon kinematics

    // 2. binding and is/fs mass
    int charge = (neutrinoMode)? -1:1;
    double binding_e = (npMode==0)? 0 : bindingE;
    double is_part_mass = (neutrinoMode)? 0.9395654133: 0.93827231; //neutron/proton
    double fs_part_mass = (neutrinoMode)? 0.93827231: 0.9395654133;

    //Calculating Q2qe:
    double q2qe = physicsUtils->qsqCCQE( beam, reco_muon_XYZT, charge, binding_e ); //GeV^2
    double ccqenu_q2 = evt->CCQENu_Q2/1000/1000;
    double enu = physicsUtils->nuEnergyCCQE( beam, reco_muon_XYZT, charge, binding_e );
    double q0 = enu - reco_muon_XYZT.E();
    double q3 = pow(q0*q0+q2qe,.5);

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

    double wgt = (isData)? 1: utils->GetCVWeight( evt, sample );
      //cout<<"Fill start"<<endl;
    recoil/=1000;
    if (fill_common_histos)
    {
      q2qe = ccqenu_q2;
      utils->fillHistosV3( utils->histos2D["h_q2qe_recoil"], q2qe, recoil, isData, evt, wgt ); 
      if(pass_side) utils->fillHistosV3( utils->histos2D["h_q2qe_recoil_SF"], q2qe,recoil,  isData, evt, wgt ); 
      if(pass_side_blob) utils->fillHistosV3( utils->histos2D["h_q2qe_recoil_SF_blob"], q2qe,recoil,  isData, evt, wgt ); 
      if(pass_side_blob_neutron) utils->fillHistosV3( utils->histos2D["h_q2qe_recoil_SF_blob_neutron"], q2qe,recoil,  isData, evt, wgt ); 
      if(pass_side_neutron) utils->fillHistosV3( utils->histos2D["h_q2qe_recoil_SF_neutron"], q2qe,recoil,  isData, evt, wgt ); 
      if(pass_blob_neutron) utils->fillHistosV3( utils->histos2D["h_q2qe_recoil_blob_neutron"], q2qe,recoil,  isData, evt, wgt ); 

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


  axis_binning UniformVisibleEBinsGeV = GetBins(0,1.2,120); //GeV
  axis_binning UniformQ2BinsGeV = GetBins(0,1.6,1280); //GeV

  MnvH2D *h_q2qe_recoil[nHistos];
  MnvH2D *h_q2qe_recoil_SF[nHistos];
  MnvH2D *h_q2qe_recoil_SF_neutron[nHistos];
  MnvH2D *h_q2qe_recoil_SF_blob[nHistos];
  MnvH2D *h_q2qe_recoil_SF_blob_neutron[nHistos];
  MnvH2D *h_q2qe_recoil_blob_neutron[nHistos];
  utils->bookHistos( h_q2qe_recoil, "h_q2qe_recoil","recoil:q2qe",UniformQ2BinsGeV, UniformVisibleEBinsGeV);
  utils->bookHistos( h_q2qe_recoil_SF, "h_q2qe_recoil_SF","recoil:q2qe",UniformQ2BinsGeV, UniformVisibleEBinsGeV);
  utils->bookHistos( h_q2qe_recoil_SF_blob, "h_q2qe_recoil_SF_blob","recoil:q2qe",UniformQ2BinsGeV, UniformVisibleEBinsGeV);
  utils->bookHistos( h_q2qe_recoil_SF_blob_neutron, "h_q2qe_recoil_SF_blob_neutron","recoil:q2qe",UniformQ2BinsGeV, UniformVisibleEBinsGeV);
  utils->bookHistos( h_q2qe_recoil_SF_neutron, "h_q2qe_recoil_SF_neutron","recoil:q2qe",UniformQ2BinsGeV, UniformVisibleEBinsGeV);
  utils->bookHistos( h_q2qe_recoil_blob_neutron, "h_q2qe_recoil_blob_neutron","recoil:q2qe",UniformQ2BinsGeV, UniformVisibleEBinsGeV);



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
