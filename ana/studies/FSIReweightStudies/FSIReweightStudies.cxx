#include "include/GlobalIncludes.h"
#include "include/CCQENuUtils.h"
#include "include/CCQENuCuts.h"
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
CCQENuCuts *cutter = new CCQENuCuts(); 

weight_fsi* wfsi = new weight_fsi("", fsi_pdg, true, 2);
double AnalysisLooper( CCQENuTruth* evt, CCQENuUtils* utils, int n_entries)
{
  bool isMC = true;
  cout<<"Start AnalysisLooper"<<endl;
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
   
    // event selection
    // end event selection
    //bool EventSelectionDataMC( CCQENuEvent *evt, CCQENuCuts *cutter, NeutronBlobCuts* ncutter,bool IsSignal, bool ISMC = true, int multiplicity = 0 , bool pass_michel = true, bool pass_nblobs = true, bool pass_neutronMainCandidate = true, bool pass_MainCandidateCut = true) //all true --> signal
      
    double wgt = 1;    
    n_outputevent+=wgt;
    //Get True KE
    double KE=-1000;
    int mode=wfsi->getQEFSIMode( evt->mc_intType, evt->mc_targetA, evt->mc_targetZ, evt->mc_er_nPart, evt->mc_er_ID, evt->mc_er_status, evt->mc_er_FD, evt->mc_er_LD, evt->mc_er_mother, evt->mc_er_Px, evt->mc_er_Py, evt->mc_er_Pz, evt->mc_er_E, KE);
    double weight=wfsi->getWeightPoly( evt->mc_intType, evt->mc_targetA, evt->mc_targetZ, evt->mc_er_nPart, evt->mc_er_ID, evt->mc_er_status, evt->mc_er_FD, evt->mc_er_LD, evt->mc_er_mother, evt->mc_er_Px, evt->mc_er_Py, evt->mc_er_Pz, evt->mc_er_E);
    KE/=1000;//converts to GeV
    if ( cutter->passTrueCCQE(evt) && evt->mc_targetA==12)
    {
      cout<<"isEvent CCQE? "<<cutter->passTrueCCQE(evt)<<endl;
      cout<<"Event A: "<<evt->mc_targetA<<endl;
      cout<<i<<" "<<mode<<" "<<KE<<endl;
      utils->fillHistos( utils->histos1D["h_neutron_true_KE"], KE, isMC, evt, wgt);
      utils->fillHistos( utils->histos2D["h_neutron_true_WeightVsKE"], KE, weight, isMC, evt, wgt);
    }


  }
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

  axis_binning muonPbins      = binner->GetMuonEnergyUniformBinsGeV();
  axis_binning trueKEBins;
  axis_binning trueWeightBins;

  vector<double> bins;
  for( float i = 0.0; i< 2; i+=0.05 ) bins.push_back(i);
  trueKEBins.bin_edges=bins;
  trueKEBins.nbins = bins.size()-1;
  trueKEBins.min = bins.front();
  trueKEBins.max = bins.back();
  bins.clear();
  for( float i = 0.0; i< 4; i+=0.01 ) bins.push_back(i);
  trueWeightBins.bin_edges=bins;
  trueWeightBins.nbins = bins.size()-1;
  trueWeightBins.min = bins.front();
  trueWeightBins.max = bins.back();


  

  cout<<"BOOK HISTOGRAMS"<<endl;
  //---------------------------------------------
  // Book Histograms
  //---------------------------------------------
  MnvH1D *h_neutron_true_KE[nHistos];
  MnvH2D *h_neutron_true_WeightVsKE[nHistos];

  utils->bookHistos( h_neutron_true_KE, "h_neutron_true_KE","Neutron True KE (GeV)", trueKEBins);
  utils->bookHistos( h_neutron_true_WeightVsKE, "h_neutron_true_WeightVsKE","WeightVsKE", trueKEBins, trueWeightBins);
  cout<<"END BOOK HISTOGRAMS"<<endl;
  //--------------------------------------------------------------
  // Add Vertical and Lateral Error Bands
  // JO is removing error bands from all 1D because they are all 
  // present in 2D. Plz make projections directly from the 2D. 
  //--------------------------------------------------------------
    double data_events, mc_events;
  //---------------------------------------------
  // Get MC Tree
  //---------------------------------------------
  cout<<"Getting MC Trees"<<endl;
  TChain *truth_mc = utils->getMCTree("Truth", n_mcfiles );
  int entries_mc   = truth_mc->GetEntries();
  cout << "MC entries: " << entries_mc << endl;
  utils->setmnvHadronReweightTruthTree(truth_mc);
  CCQENuTruth* mc = new CCQENuTruth( truth_mc );
  //int AnalysisLooper( CCQENuEvent* evt, CCQENuUtils* utils, CCQENuCuts *cutter, NeutronBlobCuts* ncutter, string sample, int multiplicity, int n_entries, bool isMC = true )
  
  cout<<"Running AnalysisLooper"<<endl;
  mc_events = AnalysisLooper( mc, utils, entries_mc );
  data_events = 0;
  
  delete mc; 
  delete truth_mc;
  
  //---------------------------------------------
  // Get Data Tree
  //---------------------------------------------
   
  //Write MC and DATA POT to file
  
  TFile *f = new TFile( filename.c_str(), "RECREATE" );
   
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
    h_neutron_true_KE[i]->Write();
    h_neutron_true_WeightVsKE[i]->Write();
  }
   
  f->Close();
  delete f;
  delete utils;
  delete cutter; 
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
