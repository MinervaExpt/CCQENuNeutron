#include "include/GlobalParameters.h"
#include "include/GlobalIncludes.h"
#include "include/CCQENuUtils.h"
#include "TParameter.h" 
#include "PlotUtils/HyperDimLinearizer.h"

#include "include/GeneralFunc.h" //in CCQENuNeutron/include/
#include "include/EventCuts.h"
#include "include/CommonBins.h"

#include "TF1.h"
#include "TH1D.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
using namespace CCQENU_ANA;

vector<double> dropFrac({0, 0.1,0.2,0.3});
vector<double> dropEnergy({0,5, 10, 15,20});
vector< pair<double, double> > dropParameters;
double modMin = 0;
double modMax = 5000;
double angle = 5;
bool UseAbfalterer = false;
bool SetResToInelastic = false;

int rwType = 0;


string neut_par = "ClusTotalE";
//string neut_par = "ClusMaxE";

//forward declaration
//HyperDimLinearizer* DeclareHDLHistos2D( CCQENuUtils* utils, MnvH2D **h, std::string name, std::string title, axis_binning &xbins, axis_binning &ybins, axis_binning &zbins, axis_binning &global_x );
//void FillHDL2D( std::string name, CCQENuUtils *utils, double x, double y, double z, bool isData, CCQENuEvent* evt, double wgt );
//================================================

//int npMode = 1;//0:neutron, 1:proton, 10:neutron+proton

template<class T>
double AnalysisLooperBase(T* evt, CCQENuUtils* &utils, CCQENuCuts* cutter, NeutronBlobCuts* ncutter, string sample, int multiplicity, int n_entries, map<int,int> &cuts_counter, bool isData=false, bool isTruth=false, int npMode = 10)
{
  GeoUtils* geoUtils = new GeoUtils();
  PhysicsUtils* physicsUtils = new PhysicsUtils();

  map< double, TH1D* > weightMaps;


  double n_evts=0;
  //return 0;

  for( int i = 0; i< n_entries; ++i)
  {
    evt->GetEntry(i);
    bool isMC = !isData;
    if( IgnoreEvent( evt, isData )) continue;
    if (i%1000 == 0) cout<<"At Event "<<i<<" ("<<(i+1)*100./float(n_entries)<<"%)"<<endl;
    //==============Apply Precut to speed things up========
    double wgt = (isData)? 1: utils->GetCVWeight( evt, sample );
    n_evts+=wgt;

    ncutter->UpdateCandidates( evt, isMC, -1, -1 );
    bool useRecoilCut = true;
    bool PassEventCut = (EventCut( evt, cutter, ncutter, cuts_counter, multiplicity, sample ,  isData, npMode, useRecoilCut, neut_par));
    if (!PassEventCut) continue;
    //cout<<"Passed Event: "<<evt->mc_incoming<<endl;
    //==============Fill Neutron Blobs======================

    double beam_bias = 0;
    XYZVector beam = geoUtils->BeamAxis( beam_bias ); 

    double reco_muon_px   =  evt->CCQENu_leptonE[0]/1000.; //GeV
    double reco_muon_py   =  evt->CCQENu_leptonE[1]/1000.; //GeV
    double reco_muon_pz   =  evt->CCQENu_leptonE[2]/1000.; //GeV
    XYZTVector reco_muon_4P( reco_muon_px, reco_muon_py, reco_muon_pz ,evt->CCQENu_muon_E/1e3);
    XYZVector reco_muon_3P( reco_muon_px, reco_muon_py, reco_muon_pz );

    int charge = (neutrinoMode)? -1:1;

    // 2. binding and is/fs mass
    double binding_e = 27.13/1E3;
    double is_part_mass = (neutrinoMode)? 0.9395654133: 0.93827231; //neutron/proton
    double fs_part_mass = (neutrinoMode)? 0.93827231: 0.9395654133;


    double q2qe = physicsUtils->qsqCCQE( beam, reco_muon_4P, charge, binding_e ); //GeV^2


    


    std::vector<bool> commonNeutronCuts = CommonNeutronCuts(evt,ncutter, neut_par);
    bool hasNC = commonNeutronCuts[0] && commonNeutronCuts[1];
    if(! hasNC ) continue;
    //cout<<"--has NC: "<<endl;

    NeutronCandidates* ptrNeutronCandidates; 
    NeutronBlob* mainCandidate; 
    if (npMode != 1 || hasNC) 
    {
      ptrNeutronCandidates = ncutter->GetCandidates();
      mainCandidate = ncutter->MainNeutronCandidate();
    }
    bool fill_common_histos = true;
    double blobDist, blobEnergy, nBlobs, n2DBlobs, n3DBlobs;
    double nonBlobEnergy;
    double proton_score=-1, pion_score=-1;
    int blobType = -1;
    bool inVtx = false;
    if( hasNC )
    {
      blobDist = mainCandidate->BlobVtxVec.R();
      blobEnergy = mainCandidate->TotalE;
      //nonBlobEnergy = recoil - blobEnergy;
      n3DBlobs = ptrNeutronCandidates->ThreeDBlobs.size();
      n2DBlobs = ptrNeutronCandidates->TwoDBlobs.size();
      nBlobs = ptrNeutronCandidates->AllBlobs.size();
      inVtx = blobDist<200;
      if( mainCandidate->hasTrack )
      {
        proton_score = mainCandidate->dEdX_proton;
        pion_score = mainCandidate->dEdX_pion;

      }
      if( isMC )
      {
        blobType = parsePartType( mainCandidate->MCPID );
      }
    }




    // 3. expected nucleon

    XYZTVector expectedNucleon4P = geoUtils->ComputeExpectedNucleon( beam, reco_muon_4P, is_part_mass, fs_part_mass, binding_e);
    XYZVector expectedNucleon3P = (XYZVector) expectedNucleon4P.Vect();



    //Neutron
    XYZVector neutronDirXYZ = ncutter->MainNeutronCandidate()->BlobVtxDir;
    TVector3  neutronDir( neutronDirXYZ.X(), neutronDirXYZ.Y(), neutronDirXYZ.Z() );

    //Proton
    TVector3 protonVect(evt->CCQENu_proton_Px_fromdEdx,evt->CCQENu_proton_Py_fromdEdx,evt->CCQENu_proton_Pz_fromdEdx);
    protonVect*=1E-3;
    XYZVector protonXYZ(evt->CCQENu_proton_Px_fromdEdx,evt->CCQENu_proton_Py_fromdEdx,evt->CCQENu_proton_Pz_fromdEdx);
    protonXYZ*=1E-3;
    double protonMass = 0.93827231;
    double Eproton = sqrt(protonXYZ.Mag2()+protonMass*protonMass);
    TLorentzVector proton4P( protonVect.X(), protonVect.Y(), protonVect.Z(), Eproton );
    XYZTVector     protonXYZT( protonVect.X(), protonVect.Y(), protonVect.Z(), Eproton );

    //Fit Neutron
    vector<double> fitMomentum = geoUtils->FitMomentum( beam, reco_muon_3P, protonXYZ, (XYZVector) neutronDirXYZ.Unit(), geoUtils->MnGeV() );
    double pneutron = fitMomentum[0];
    double Mneutron = geoUtils->MnGeV();
    double Eneutron = pow(pneutron*pneutron+Mneutron*Mneutron,0.5);

    XYZVector neutronXYZ = neutronDirXYZ*pneutron;

    XYZTVector neutronFitXYZT( neutronXYZ.X(),neutronXYZ.Y(),neutronXYZ.Z(),Eneutron);
    XYZVector neutronFitXYZ = neutronXYZ;

    double Mcarbon = 6*fs_part_mass+6*is_part_mass-92.163/1E3;
    std::vector<double> ddtransverseVars = geoUtils->GetTransverseVariables( beam, reco_muon_4P+protonXYZT, neutronFitXYZT, Mcarbon, 27.13/1E3);
    double dd_En      = ddtransverseVars[vEnu];
    double dd_pn      = ddtransverseVars[vPn];
    double dd_dpt     = ddtransverseVars[vdpT];
    double dd_dptx    = ddtransverseVars[vdpTx];
    double dd_dpty    = ddtransverseVars[vdpTy];
    double dd_dalphat = ddtransverseVars[vdalphaT];
    double dd_dphit   = ddtransverseVars[vdphiT];
    double dd_sign    = ddtransverseVars[vsign];
    double dd_signed_dalphat = dd_sign*dd_dalphat ;
    double dd_signed_dphit = dd_dphit *dd_sign;

    dd_dphit*=180./TMath::Pi();

    //Neutrino Energy
    double Enu_True = evt->mc_incomingE/1000.;
    double Enu_CCQE = physicsUtils->nuEnergyCCQE( beam, reco_muon_4P, charge, binding_e );
    double Enu_Fit = beam.Dot( reco_muon_3P+protonXYZ+neutronFitXYZ );
    double Enu_NoFit = beam.Dot( reco_muon_3P+protonXYZ );

    double dEnu_CCQE = Enu_CCQE - Enu_True;
    double dEnu_Fit = Enu_Fit - Enu_True;
    double dEnu_NoFit = Enu_NoFit - Enu_True;

    double dEnuFrac_CCQE = dEnu_CCQE/Enu_True;
    double dEnuFrac_Fit = dEnu_CCQE/Enu_True;
    double dEnuFrac_NoFit = dEnu_CCQE/Enu_True;

    utils->fillHistosV3( utils->histos2D["h_enuTrue_dphitt"], Enu_True, dd_dphit, isData, evt, wgt );
    utils->fillHistosV3( utils->histos2D["h_enuCCQE_dphitt"], Enu_CCQE, dd_dphit, isData, evt, wgt );
    utils->fillHistosV3( utils->histos2D["h_enuFit_dphitt"], Enu_Fit, dd_dphit, isData, evt, wgt );
    utils->fillHistosV3( utils->histos2D["h_enuNoFit_dphitt"], Enu_NoFit, dd_dphit, isData, evt, wgt );

    utils->fillHistosV3( utils->histos2D["h_denuCCQE_dphitt"], dEnu_CCQE, dd_dphit, isData, evt, wgt );
    utils->fillHistosV3( utils->histos2D["h_denuFit_dphitt"], dEnu_Fit, dd_dphit, isData, evt, wgt );
    utils->fillHistosV3( utils->histos2D["h_denuNoFit_dphitt"], dEnu_NoFit, dd_dphit, isData, evt, wgt );

    utils->fillHistosV3( utils->histos2D["h_denuFracCCQE_dphitt"], dEnuFrac_CCQE, dd_dphit, isData, evt, wgt );
    utils->fillHistosV3( utils->histos2D["h_denuFracFit_dphitt"], dEnuFrac_Fit, dd_dphit, isData, evt, wgt );
    utils->fillHistosV3( utils->histos2D["h_denuFracNoFit_dphitt"], dEnuFrac_NoFit, dd_dphit, isData, evt, wgt );

    utils->fillHistosV3( utils->histos2D["h_enuTrue_q2qe"], Enu_True, q2qe, isData, evt, wgt );
    utils->fillHistosV3( utils->histos2D["h_enuCCQE_q2qe"], Enu_CCQE, q2qe, isData, evt, wgt );
    utils->fillHistosV3( utils->histos2D["h_enuFit_q2qe"], Enu_Fit, q2qe, isData, evt, wgt );
    utils->fillHistosV3( utils->histos2D["h_enuNoFit_q2qe"], Enu_NoFit, q2qe, isData, evt, wgt );

    utils->fillHistosV3( utils->histos2D["h_denuCCQE_q2qe"], dEnu_CCQE, q2qe, isData, evt, wgt );
    utils->fillHistosV3( utils->histos2D["h_denuFit_q2qe"], dEnu_Fit, q2qe, isData, evt, wgt );
    utils->fillHistosV3( utils->histos2D["h_denuNoFit_q2qe"], dEnu_NoFit, q2qe, isData, evt, wgt );

    utils->fillHistosV3( utils->histos2D["h_denuFracCCQE_q2qe"], dEnuFrac_CCQE, q2qe, isData, evt, wgt );
    utils->fillHistosV3( utils->histos2D["h_denuFracFit_q2qe"], dEnuFrac_Fit, q2qe, isData, evt, wgt );
    utils->fillHistosV3( utils->histos2D["h_denuFracNoFit_q2qe"], dEnuFrac_NoFit, q2qe, isData, evt, wgt );


  }
  return n_evts;
}



double AnalysisLooper(CCQENuEvent* evt, CCQENuUtils* utils, CCQENuCuts* cutter, NeutronBlobCuts* ncutter, string sample, int multiplicity, int n_entries, map<int,int>&cutsC, bool isData=false, bool isTruth=false, int npMode = 1 )
{ 
  return AnalysisLooperBase(evt, utils, cutter, ncutter, sample, multiplicity,n_entries, cutsC, isData, isTruth, npMode );
}

int MuonSelectionHists(string filename, string sample, bool makeFluxConstraintHisto, int multiplicity, string playlist, int n_mcfiles = -1, int n_datafiles = -1, int npMode = 10)
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
  //utils->SetNeutronReweight();
  CCQENuCuts   *cutter = new CCQENuCuts();
  NeutronBlobCuts  *ncutter = new NeutronBlobCuts();
  GeoUtils *geoUtils = new GeoUtils();


  //--------------------------------------------
  // Define histograms
  //--------------------------------------------
  MnvH2D* h_enuTrue_dphitt[nHistos];
  MnvH2D* h_enuCCQE_dphitt[nHistos];
  MnvH2D* h_enuFit_dphitt[nHistos];
  MnvH2D* h_enuNoFit_dphitt[nHistos];
  MnvH2D* h_denuCCQE_dphitt[nHistos];
  MnvH2D* h_denuFit_dphitt[nHistos];
  MnvH2D* h_denuNoFit_dphitt[nHistos];
  MnvH2D* h_denuFracCCQE_dphitt[nHistos];
  MnvH2D* h_denuFracFit_dphitt[nHistos];
  MnvH2D* h_denuFracNoFit_dphitt[nHistos];

  utils->bookHistos( h_enuTrue_dphitt, "h_enuTrue_dphitt", "EnuTrue(GeV): dphitt(degree)", leptonmomentumbins,dphitbins );
  utils->bookHistos( h_enuCCQE_dphitt, "h_enuCCQE_dphitt", "EnuCCQE(GeV): dphitt(degree)", leptonmomentumbins,dphitbins );
  utils->bookHistos( h_enuFit_dphitt, "h_enuFit_dphitt", "EnuFit(GeV): dphitt(degree)", leptonmomentumbins,dphitbins );
  utils->bookHistos( h_enuNoFit_dphitt, "h_enuNoFit_dphitt", "EnuNoFit(GeV): dphitt(degree)", leptonmomentumbins,dphitbins );

  axis_binning dEnuBins=GetBins(-2,2,100);
  utils->bookHistos( h_denuCCQE_dphitt, "h_denuCCQE_dphitt", "dEnuCCQE(GeV): dphitt(degree)", dEnuBins,dphitbins );
  utils->bookHistos( h_denuFit_dphitt, "h_denuFit_dphitt", "dEnuFit(GeV): dphitt(degree)", dEnuBins,dphitbins );
  utils->bookHistos( h_denuNoFit_dphitt, "h_denuNoFit_dphitt", "dEnuNoFit(GeV): dphitt(degree)", dEnuBins,dphitbins );
  utils->bookHistos( h_denuFracCCQE_dphitt, "h_denuFracCCQE_dphitt", "dEnuCCQE(GeV): dphitt(degree)", dEnuBins,dphitbins );
  utils->bookHistos( h_denuFracFit_dphitt, "h_denuFracFit_dphitt", "dEnuFit(GeV): dphitt(degree)", dEnuBins,dphitbins );
  utils->bookHistos( h_denuFracNoFit_dphitt, "h_denuFracNoFit_dphitt", "dEnuNoFit(GeV): dphitt(degree)", dEnuBins,dphitbins );
  
  MnvH2D* h_enuTrue_q2qe[nHistos];
  MnvH2D* h_enuCCQE_q2qe[nHistos];
  MnvH2D* h_enuFit_q2qe[nHistos];
  MnvH2D* h_enuNoFit_q2qe[nHistos];
  MnvH2D* h_denuCCQE_q2qe[nHistos];
  MnvH2D* h_denuFit_q2qe[nHistos];
  MnvH2D* h_denuNoFit_q2qe[nHistos];
  MnvH2D* h_denuFracCCQE_q2qe[nHistos];
  MnvH2D* h_denuFracFit_q2qe[nHistos];
  MnvH2D* h_denuFracNoFit_q2qe[nHistos];

  utils->bookHistos( h_enuTrue_q2qe, "h_enuTrue_q2qe", "EnuTrue(GeV): q2qe(GeV^{2})", leptonmomentumbins,Q2bins );
  utils->bookHistos( h_enuCCQE_q2qe, "h_enuCCQE_q2qe", "EnuCCQE(GeV): q2qe(GeV^{2})", leptonmomentumbins,Q2bins );
  utils->bookHistos( h_enuFit_q2qe, "h_enuFit_q2qe", "EnuFit(GeV): q2qe(GeV^{2})", leptonmomentumbins,Q2bins );
  utils->bookHistos( h_enuNoFit_q2qe, "h_enuNoFit_q2qe", "EnuNoFit(GeV): q2qe(GeV^{2})", leptonmomentumbins,Q2bins );

  utils->bookHistos( h_denuCCQE_q2qe, "h_denuCCQE_q2qe", "dEnuCCQE(GeV): q2qe(GeV^{2})", dEnuBins,Q2bins );
  utils->bookHistos( h_denuFit_q2qe, "h_denuFit_q2qe", "dEnuFit(GeV): q2qe(GeV^{2})", dEnuBins,Q2bins );
  utils->bookHistos( h_denuNoFit_q2qe, "h_denuNoFit_q2qe", "dEnuNoFit(GeV): q2qe(GeV^{2})", dEnuBins,Q2bins );
  utils->bookHistos( h_denuFracCCQE_q2qe, "h_denuFracCCQE_q2qe", "dEnuCCQE(GeV): q2qe(GeV^{2})", dEnuBins,Q2bins );
  utils->bookHistos( h_denuFracFit_q2qe, "h_denuFracFit_q2qe", "dEnuFit(GeV): q2qe(GeV^{2})", dEnuBins,Q2bins );
  utils->bookHistos( h_denuFracNoFit_q2qe, "h_denuFracNoFit_q2qe", "dEnuNoFit(GeV): q2qe(GeV^{2})", dEnuBins,Q2bins );


  TVector2 nevents;

  //--------------------------------------------
  // Define GeoUtils
  //--------------------------------------------

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
  mc_events=AnalysisLooper( mc, utils, cutter, ncutter, sample,  multiplicity, entries_mc, mc_counter, isData, isTruth, npMode );
  delete mc; 
  delete tree_mc; 
  
  isData = true;
  isTruth= false;
  if( n_datafiles != 0 )
  {
    TChain* tree_data  = utils->getDataTree("CCQENu", n_datafiles );
    int entries_data   = tree_data->GetEntries();
    CCQENuEvent* data = new CCQENuEvent( tree_data, isData );
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
  cout<<"MC Events: "<<mc_events<<endl;
  
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
  par.push_back("1");//proton, neutron,proton/neutron mode
  par.push_back("0"); //modMin
  par.push_back("0"); //modMax
  par.push_back("0"); //UseAbfalterer
  par.push_back("0"); //SetResToInelastic
  par.push_back("1"); //RWType 1-> Drop Cluster, other-> reweight

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
  modMin = stod(par[9].c_str() );
  modMax = stod(par[10].c_str() );
  UseAbfalterer = stoi( par[11].c_str() );
  SetResToInelastic = stoi(par[12].c_str() );
  rwType = stoi(par[13].c_str() );
  

  pass_neutronQsqCut = (npmode==10)?-1:1;

  
//
  //return MuonSelectionHists(output_path, sample, fluxConstraintHistoNeeded, multiplicity, par[2],atoi(par[6].c_str()), atoi(par[7].c_str()), atoi(par[8].c_str()) );
  return MuonSelectionHists(output_path, sample, fluxConstraintHistoNeeded, multiplicity, playlist, n_mc, n_data, npmode );
}
