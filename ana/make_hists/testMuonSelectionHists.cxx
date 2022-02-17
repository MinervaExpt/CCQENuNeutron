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
double AnalysisLooperBase(T* evt, CCQENuUtils* &utils, CCQENuCuts* cutter, NeutronBlobCuts* ncutter, string sample, int multiplicity, int n_entries, map<int,int> &cuts_counter, bool isData=false, bool isTruth=false, int npMode = 1)
{
  GeoUtils* geoUtils = new GeoUtils();
  PhysicsUtils* physicsUtils = new PhysicsUtils();

  double n_evts=0;
  //return 0;

  for( int i = 0; i< n_entries; ++i)
  {
    evt->GetEntry(i);
    if( IgnoreEvent( evt, isData )) continue;
    if (i%1000 == 0) cout<<"At Event "<<i<<" ("<<(i+1)*100./float(n_entries)<<"%)"<<endl;
    //==============Fill Neutron Blobs======================
    bool isMC = !isData;
    ncutter->UpdateCandidates( evt, isMC );
    if(isDebug) cout<<"Updated ncutter"<<endl;
    //==============Event Selection========================
    bool PassEventCut = (EventCut( evt, cutter, ncutter, cuts_counter, multiplicity, sample ,  isData, npMode));
    if(isDebug) cout<<PassEventCut<<endl;


    //now add in neutron-->
    std::vector<bool> commonNeutronCuts = CommonNeutronCuts(evt,ncutter);
    bool hasNC = commonNeutronCuts[0] && commonNeutronCuts[1];
    bool isForward = commonNeutronCuts[2] && hasNC;
    bool isBackward = commonNeutronCuts[3] && hasNC;
    bool hasProton = (npMode!=0);

    //double neutronwgt = utils->GetNeutronCVRW( evt );
    //cout<<"NeutronWgt: "<<neutronwgt<<endl;


    if (!PassEventCut) continue;
    NeutronCandidates* ptrNeutronCandidates; 
    NeutronBlob* mainCandidate; 
    if (npMode != 1 || hasNC) 
    {
      ptrNeutronCandidates = ncutter->GetCandidates();
      mainCandidate = ncutter->MainNeutronCandidate();
    }
    if(!isData)
    {
      double w = utils->GetNeutronCVRW( evt );
      double ke = 0;
      double mass = 939.56542052 ;
      for( uint i = 0; i< evt->mc_nFSPart; i++ )
      {
        if ( evt->mc_FSPartPDG[i] != 2112 ) continue;
        if( evt->mc_FSPartE[i]-mass > ke ) ke = evt->mc_FSPartE[i]-mass;
      }

      //cout<<ke<<", "<<w<<endl;
      if(ke>0)
      {
        utils->fillHistosV3( utils->histos2D["h_ke_weight"], ke, w, isData, evt, 1 );
      }
    }


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
  MnvH2D* h_ke_weight[nHistos];
  utils->bookHistos( h_ke_weight,"h_ke_weight","KE:weight",100,0,1000,100,0,10);


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
  string output_path = par[1];
  string sample = par[3];
  int multiplicity = atoi(par[5].c_str() );
  string playlist = par[2];
  int n_mc = atoi(par[6].c_str() );
  int n_data = atoi(par[7].c_str() );
  int npmode = atoi(par[8].c_str() );
//
  //return MuonSelectionHists(output_path, sample, fluxConstraintHistoNeeded, multiplicity, par[2],atoi(par[6].c_str()), atoi(par[7].c_str()), atoi(par[8].c_str()) );
  return MuonSelectionHists(output_path, sample, fluxConstraintHistoNeeded, multiplicity, playlist, n_mc, n_data, npmode );
}
