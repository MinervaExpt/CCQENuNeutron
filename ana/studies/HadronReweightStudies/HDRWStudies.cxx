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
    if( i == testEvent ) cout<<"At testEvent "<<testEvent<<endl;
    if( IgnoreEvent( evt, isData )) continue;
    //if( !cutter->passTrueCCQELike( evt ) ) continue;
    //if( !cutter->passTrueCCQE( evt ) ) continue;
    //if( evt->mc_targetZ != 1 ) continue;
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

    vector<XYZTVector> true_particles_XYZT;
    bool passTrue = GetTrueParticles(  evt, utils, kineMod, physicsUtils, geoUtils, isMC, npMode, true_particles_XYZT );
    XYZTVector neutrino = true_particles_XYZT[0];
    XYZTVector muon = true_particles_XYZT[1];   XYZVector muon_XYZ = (XYZVector) muon.Vect();
    XYZTVector proton = true_particles_XYZT[2]; XYZVector proton_XYZ = (XYZVector) proton.Vect();
    XYZTVector neutron = true_particles_XYZT[3]; XYZVector neutron_XYZ = (XYZVector) neutron.Vect();

    TVector3 proton_V3 = Convert( proton_XYZ );
    TVector3 neutron_V3 = Convert( neutron_XYZ );

    bool has_proton = proton.E()>0;
    bool has_neutron = neutron.E()>0;


    XYZTVector nu( evt->mc_incomingPartVec[0]/1000.,evt->mc_incomingPartVec[1]/1000.,evt->mc_incomingPartVec[2]/1000.,evt->mc_incomingPartVec[3]/1000.);

    double cos_muon_theta =  beam.Unit().Dot( muon.Vect().Unit() );
    double muon_theta= ACos( cos_muon_theta );
    double muon_T    = (muon.E()-muon.M());
    double muon_px   = muon.Px();
    double muon_py   = muon.Py();
    double muon_pz   = muon.Pz();
    double muon_p    = muon.P();

    muon_theta*=180/TMath::Pi();


    XYZTVector primary_nucleon_XYZT = (npMode == 0)? neutron : proton;
    XYZVector primary_nucleon_XYZ = primary_nucleon_XYZT.Vect();

    XYZTVector secondary_nucleon = (npMode == 0)? proton : neutron;

    bool has_primary = primary_nucleon_XYZT.E()>0;
    bool has_secondary = secondary_nucleon.E()>0;

    if( i == testEvent )
    {
      cout<<"Muon momentum: "<<muon.X()<<" "<<muon.Y()<<" "<<muon.Z()<<endl;
      cout<<"Beam direction: "<<beam.X()<<" "<<beam.Y()<<" "<<beam.Z()<<endl;
      cout<<"Muon theta: "<<muon_theta<<endl;

    }


    //=========Truth Kinematic Cuts==========================
    //if( !evt->truth_is_fiducial ) continue;
    //if( !DoFullKinematicRegion )
    //{
    //  if( muon_theta > 20 ) continue; //in beam cooridnate
    //  if( ( muon_p < MuonPMin || muon_p > MuonPMax ) && useMuonPCut ) continue; // a muon momentum cut
    //}

    //if(npMode == 0 && has_proton) continue; //for neutron only mode
    //if ( npMode !=0 && !NucleonKinematicsCut( beam, (XYZVector) proton.Vect() ) ) continue;
    //if(!has_primary) continue;

    //cout<<i<<", "<<endl;
    //Calculations

    //New method (maybe temporary)
    double muon_pt_beam = sin(muon_theta)*muon_p;
    double muon_pz_beam = cos(muon_theta)*muon_p;

    //Convert to degrees
    muon_theta*= 180. / TMath::Pi();

    // Binding and is/fs mass
    int charge = (neutrinoMode)? -1:1;
    double binding_e = (npMode==0)? 0 : bindingE;
    double is_part_mass = (neutrinoMode)? 0.9395654133: 0.93827231; //neutron/proton
    double fs_part_mass = (neutrinoMode)? 0.93827231: 0.9395654133;


    // 3. expected nucleon
    XYZTVector expected_nucleon_XYZT = geoUtils->ComputeExpectedNucleon( beam, muon, is_part_mass, fs_part_mass, binding_e);
    XYZVector expected_nucleon_XYZ = (XYZVector) expected_nucleon_XYZT.Vect();

    //Calculating Q2qe:
    double q2qe = physicsUtils->qsqCCQE( beam, muon, charge, binding_e ); //GeV^2
    double enu = physicsUtils->nuEnergyCCQE( beam, muon, charge, binding_e ); //GeV^2

    //preparing for STKI and angular vars:
    double Mcarbon = 6*fs_part_mass+6*is_part_mass-92.163/1E3;
    double Enuc=primary_nucleon_XYZT.E();


    //calculating transverse variables
    std::vector<double> true_transverse_vars = geoUtils->GetTransverseVariables( beam, muon, primary_nucleon_XYZT, Mcarbon, 27.13/1E3);
    //calculating angular variables
    std::vector<double> true_neutron_vars = geoUtils->ComputeNeutronAngularVars( beam, expected_nucleon_XYZ, primary_nucleon_XYZ );
    // test if the expected and muon cancel, they should

    double En      = true_transverse_vars[vEnu];
    double pn      = true_transverse_vars[vPn];
    double dpt     = true_transverse_vars[vdpT];
    double dptx    = true_transverse_vars[vdpTx];
    double dpty    = true_transverse_vars[vdpTy];
    double dalphat = true_transverse_vars[vdalphaT];
    double dphit   = true_transverse_vars[vdphiT];
    double sign    = true_transverse_vars[vsign];
    double signed_dalphat = sign*dalphat ;
    double signed_dphit = dphit *sign;

    double dthetaP = true_neutron_vars[vdthetaP];
    double dthetaR = true_neutron_vars[vdthetaR];
    double dtheta = true_neutron_vars[vdtheta];
    //if (isMC) cout<< evt->mc_intType<<endl;
    dthetaP *= 180/TMath::Pi();
    dthetaR *= 180/TMath::Pi();
    dtheta *= 180/TMath::Pi();
    //cout<<dthetaP<<endl;
    double dPp = true_neutron_vars[vdPp];
    double dPr = true_neutron_vars[vdPr];
    double dPp_infer = true_neutron_vars[vdPpi];
    double dPr_infer = true_neutron_vars[vdPri];

    double protonTn = ( proton.E()- proton.M());
    double protonAngle = ACos( proton.Vect().Unit().Dot( beam ) ) * 180/TMath::Pi();

    bool q2Cut = passQ2( q2qe );

    double neutronWeight = utils->GetNeutronCVRW( evt );

    if (fill_common_histos)
    {
      //---------------------------------------------
      //Fill 1-D Plots
      //---------------------------------------------
      //cout<<"Fill start"<<endl;

      utils->fillHistosV3( utils->histos1D["h_q2qe"],q2qe, isData, evt, wgt );
      utils->fillHistosV3( utils->histos1D["h_q2qe_fine"],q2qe, isData, evt, wgt );
      utils->fillHistosV3( utils->histos1D["h_q2qe_neutroncvrw"],q2qe, isData, evt, wgt*neutronWeight );
      utils->fillHistosV3( utils->histos1D["h_q2qe_fine_neutroncvrw"],q2qe, isData, evt, wgt*neutronWeight );

      utils->fillVertErrorBands(utils->histos1D["h_q2qe"],q2qe, evt);
      utils->fillVertErrorBands(utils->histos1D["h_q2qe_fine"],q2qe, evt);
      utils->fillVertErrorBands(utils->histos1D["h_q2qe_neutroncvrw"],q2qe, evt);
      utils->fillVertErrorBands(utils->histos1D["h_q2qe_fine_neutroncvrw"],q2qe, evt);
      
    } // End of fillHistosV3 and fillVertErrorBands where CV of events pass the proton and recoil cuts 
  }
  return n_evts;
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


  axis_binning UniformQ2BinsGeV = GetBins(0,10,2000); //GeV

  MnvH1D *h_q2qe[nHistos], *h_q2qe_fine[nHistos];
  utils->bookHistos( h_q2qe, "h_q2qe","q2qe",Q2bins );
  utils->bookHistos( h_q2qe_fine, "h_q2qe_fine","q2qe_fine",UniformQ2BinsGeV );

  MnvH1D *h_q2qe_neutroncvrw[nHistos];
  MnvH1D *h_q2qe_fine_neutroncvrw[nHistos];
  utils->bookHistos( h_q2qe_neutroncvrw, "h_q2qe_neutroncvrw","q2qe_ncvrw",Q2bins );
  utils->bookHistos( h_q2qe_fine_neutroncvrw, "h_q2qe_fine_neutroncvrw","q2qe_fine_ncvrw",UniformQ2BinsGeV );

  utils->addVertErrorBands( h_q2qe );
  utils->addVertErrorBands( h_q2qe_fine );
  utils->addVertErrorBands( h_q2qe_neutroncvrw );
  utils->addVertErrorBands( h_q2qe_fine_neutroncvrw );



  TChain *truth_mc = utils->getMCTree("Truth", n_mcfiles );
  
  utils->setmnvHadronReweightTruthTree(truth_mc);
  utils->SetNeutronReweight();

  int entries_truth = truth_mc->GetEntries();
  cout << "Truth entries: " << entries_truth << endl;
  CCQENuTruth* truth = new CCQENuTruth( truth_mc );
  map<int,int> cuts_counter;
  double truth_events = EfficiencyLooperTruth(truth, utils,cutter, ncutter, "Signal", multiplicity, entries_truth, cuts_counter, npMode);
 
 
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
