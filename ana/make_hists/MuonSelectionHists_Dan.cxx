#include "include/CCQENuUtils.h"
#include "TParameter.h" 
#include "PlotUtils/HyperDimLinearizer.h"
using namespace CCQENU_ANA;

int MuonSelectionHists(string filename, string sample, bool makeFluxConstraintHisto, int multiplicity, string playlist, int n_mcfiles = -1, int n_datafiles = -1)
{

  //---------------------------------------------
  // create file to store histograms
  //---------------------------------------------
  //TFile *f = new TFile( filename.c_str(), "RECREATE" );
  makeFluxConstraintHisto = false;//turn off nue
  CCQENuUtils  *utils  = new CCQENuUtils( false, makeFluxConstraintHisto );
  utils->setPlaylist(playlist);
  utils->setFluxReweighterPlaylist();
  utils->setAnalysisType(kIMD);
  CCQENuCuts   *cutter = new CCQENuCuts();
  AnaBinning *binner = new AnaBinning();
  CCQENuBinning *ccqenubinner = new CCQENuBinning();

  //Original 2D LE binning
  axis_binning muonPtbins     = ccqenubinner->GetInclusiveMuonPtBinsGeV();
  axis_binning muonPzbins     = ccqenubinner->GetInclusiveMuonPzBinsGeV();

  

  vector<double> et_bins;
  vector<double> full_recoil_bins;
  vector<double> vtx_recoil_bins;

  for(int i=0;i<21;i++) et_bins.push_back(0.00025*i);
  for(int i=0;i<11;i++) full_recoil_bins.push_back(20*i);
  for(int i=0;i<11;i++) vtx_recoil_bins.push_back(10*i);

  //Temp 3D cv only
  //E*theta*theta 
  vector<vector<double> > imd;
  imd.push_back(vtx_recoil_bins);
  imd.push_back(et_bins);
  imd.push_back(full_recoil_bins);
  HyperDimLinearizer *imd3d = new HyperDimLinearizer(imd,0);
  int n_bins_long_imd = (vtx_recoil_bins.size()+1)*(full_recoil_bins.size()+1);
  vector<double> bins_x_global_imd;
  for(int i=0;i<n_bins_long_imd;i++) bins_x_global_imd.push_back(i);

  axis_binning glob_x_long_imd;
  glob_x_long_imd.bin_edges = bins_x_global_imd;
  glob_x_long_imd.nbins = bins_x_global_imd.size()-1;
  glob_x_long_imd.min = bins_x_global_imd.front();
  glob_x_long_imd.max = bins_x_global_imd.back();

  axis_binning etbin;
  etbin.bin_edges = et_bins;
  etbin.nbins = et_bins.size()-1;
  etbin.min = et_bins.front();
  etbin.max = et_bins.back();


  MnvH2D *h_vtxetfull[nHistosInc];
  utils->bookHistos( h_vtxetfull, "h_vtxetfull","et vs vtxfull",glob_x_long_imd,etbin);
  


  //---------------------------------------------
  // Book Histograms
  //---------------------------------------------
  //
  //2-D
  MnvH2D *h_E_etheta[nHistosInc],*h_nuParent_nuParentP[nHistosInc],*h_nuParent_nuParentPt[nHistosInc],*h_kaon_parPz_parPt[nHistosInc],*h_pion_parPz_parPt[nHistosInc];

  utils->bookHistos( h_E_etheta, "h_E_etheta","E vs eth2",45,5,50,40,0,0.01);
  utils->bookHistos( h_nuParent_nuParentP, "h_nuParent_nuParentP","parent vs P",800,-400,400,480,0,120);
  utils->bookHistos( h_nuParent_nuParentPt, "h_nuParent_nuParentPt","parent vs Pt",800,-400,400,200,0,1);
  utils->bookHistos( h_kaon_parPz_parPt, "h_kaon_parPz_parPt","Kaon pt vs pz",400,0,120,200,0,1);
  utils->bookHistos( h_pion_parPz_parPt, "h_pion_parPz_parPt","Pion pt vs pz",400,0,120,200,0,1);

  //--------------------------------------------------------------
  // Add Vertical and Lateral Error Bands
  // JO is removing error bands from all 1D because they are all 
  // present in 2D. Plz make projections directly from the 2D. 
  //--------------------------------------------------------------
  if(RunCodeWithSystematics){
    cout << "Add systematics" << endl;
    utils->addVertErrorBands( h_E_etheta );
    utils->addLatErrorBands( h_E_etheta );
    
  }//end run with systematics
  //  return 0;
  TChain *truth_mc = utils->getMCTree("Truth", n_mcfiles );
    //---------------------------------------------
  // Get MC Tree
  //---------------------------------------------
  TChain* tree_mc   =  utils->getMCTree("CCQENu", n_mcfiles );
  int entries_mc   = tree_mc->GetEntries();
  cout << "MC entries: " << entries_mc << endl;
  CCQENuEvent* mc = new CCQENuEvent( tree_mc );
  cout<<"Set MC"<<endl;
  CCQENuTruth* truth = new CCQENuTruth( truth_mc );
  utils->setmnvHadronReweightDataTree(tree_mc);
  utils->setmnvHadronReweightTruthTree(truth_mc);
  cout<<"Set Truth"<<endl;

  //-----------------------------------------------------------
  // Running over event selection criteria for each MC event 
  //-----------------------------------------------------------
  double mc_events = 0., data_events = 0.;
  
  int count_high_vals = 0;

  for (int i=0; i<entries_mc ; ++i) {
    if( (i+1)%100000 == 0 ) 
      cout << "Reading MC entry : " << i+1 << " - " << 100*(i+1)/entries_mc << " % " << endl; 
    mc->GetEntry(i);
    //    if(mc->mc_intType!=6) continue;
    if ( sample == "Signal") {
      if( multiplicity==0) { //1+2tracks
	if( !cutter->passInteractionVertex( mc ) )     continue;	
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )   	  continue;  
	if( !cutter->passTrackAngleCut( mc ) )          continue;
      }
    } else {
      cout<<"Wrong Sample Selection."<<endl;
      cout<<"Select: Signal "<<endl;
      return 1;
    }

    //Get CVweight
    double wgt = utils->GetCVWeight( mc , sample);    
    mc_events+=wgt;
    //-----------
    // This is the old method. Because of angular bias introduced from hadronic energy we need to pick out the right angle from a downstream node
    // New method until general reco is fixed is to get the reco P, a downstream theta (already corrected to the beam direction), and calculate
    // the pt and pz from these two quantities
    //-----------
    
    //double muon_theta     =  mc->CCQENu_muon_theta;
    double muon_theta     = mc->muon_theta;
    double cos_muon_theta =  cos(muon_theta);
    double muon_T         =  mc->CCQENu_muon_T / pow(10,3); //GeV
    
    double reco_muon_px   =  mc->CCQENu_leptonE[0];
    double reco_muon_py   =  mc->CCQENu_leptonE[1];
    double reco_muon_pz   =  mc->CCQENu_leptonE[2];
        
    double muon_pt_beam   = utils->GetTransverseMomentumWRTBeam( reco_muon_px, reco_muon_py, reco_muon_pz ) / pow(10,3); //GeV
    double muon_pz_beam   = utils->GetLongitudinalMomentumWRTBeam( reco_muon_px, reco_muon_py, reco_muon_pz ) / pow(10,3); //GeV
    
    double muon_p         = utils->GetTotalMomentum( reco_muon_px, reco_muon_py, reco_muon_pz ) / pow(10,3); //GeV

    //Change muon momentum to GeV
    reco_muon_px /= pow(10,3);
    reco_muon_py /= pow(10,3);
    reco_muon_pz /= pow(10,3);


    //reconstructed invariant mass
    double q2_gen = 2.0*(mc->CCQENu_leptonE[3]+mc->recoil_energy_nonmuon_nonvtx0mm-muon_p*1000.0*cos(muon_theta))-(105*105);
    double reco_w = (938.*938.)+(2*938.*mc->recoil_energy_nonmuon_nonvtx0mm)-q2_gen;
    reco_w = (reco_w<0) ? 0.0 : sqrt(reco_w);

    //recoil variants
    double recoil = mc->nonvtx_iso_blobs_energy + mc->dis_id_energy;
    double recoil_inc = mc->recoil_summed_energy[0];
    double muon_em_blobs = mc->muon_iso_blobs_energy > 0 ? mc->muon_iso_blobs_energy:0.0;

    //vtx energy
    double vtx_energy = mc->vtx_blobs_energy;

    //vtx energy fraction
    double vtx_energy_fraction = vtx_energy/(vtx_energy+recoil);
    if(vtx_energy+recoil <=0)vtx_energy_fraction=-0.005; //When there is absolutely no reconstructed recoil we need to avoid NAN.

    double etheta2 =  muon_theta*muon_theta*mc->CCQENu_leptonE[3]/1000;
    //    if(vtx_energy_fraction <0) cout << vtx_energy << "\t" << recoil << "\t" << vtx_energy_fraction << endl;
    /*
    cout << "-------------------------------------------------------------------------" << endl;
    cout << vtx_energy_fraction << "\t" << mc->recoil_energy_nonmuon_vtx300mm << "\t" << mc->recoil_energy_nonmuon_nonvtx300mm << endl;
    cout << "\t" << vtx_energy/(vtx_energy+recoil)<<"\t" << vtx_energy<<"\t" << vtx_energy+recoil << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    */
    if(vtx_energy/(vtx_energy+recoil) > 1) cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
    //Neutrino Energies
    double nu_cal_E = (mc->CCQENu_muon_E+recoil_inc)/1e3;
    double nu_qe_E = mc->CCQENu_enu_muon/1e3;
    

    
    vector<double> val_vect_imd;
    val_vect_imd.push_back(vtx_energy);
    val_vect_imd.push_back(etheta2);
    val_vect_imd.push_back(recoil_inc);
    float globx_long_imd = imd3d->GetBin(val_vect_imd).first;

    //-------------------------------------------------------------------
    // Check that the CV of events pass the proton and recoil cuts 
    // This is because we have systematics due to these cuts 
    // and a change in their values can affect the composition 
    // of the selected sample. 
    //-------------------------------------------------------------------
    /*
    if(muon_pz_beam > 10.0){
      utils->fillHistos( h_vtxetfull,globx_long_imd+0.0001,etheta2,mc,wgt);
    }
    */
    //vtx and full recoil cuts... (x-axis has E from 5 to 50 so I can cut in histograms and avoid the blasted later cut issue
    if(recoil_inc<80 && vtx_energy <10){
      utils->fillHistos( h_E_etheta,mc->CCQENu_leptonE[3]/1000,etheta2,mc,wgt);
      //These are for the neutrino ancestor plots which we want to have passed the reco cuts
      if(mc->CCQENu_leptonE[3]/1000>10 && etheta2<0.00075){

	//	cout << mc->multiplicity << "\t" << mc->mc_run << "\t" << mc->mc_subrun << "\t" << mc->mc_nthEvtInFile << endl;
	double pdpz = mc->mc_fr_nuParentProdP[2];
	double pdpx = mc->mc_fr_nuParentProdP[0]*pdpz;//Apparently these are stored as the fraction of Pz... undo that
	double pdpy = mc->mc_fr_nuParentProdP[1]*pdpz;

	double parentP = TMath::Sqrt(pdpx*pdpx+pdpy*pdpy+pdpz*pdpz)/1000.;
	double parentPt = TMath::Sqrt(pdpx*pdpx+pdpy*pdpy)/1000.;

	int parentPdg = mc->mc_fr_nuParentID;

	utils->fillHistos( h_nuParent_nuParentP,parentPdg,parentP,mc,wgt);
	utils->fillHistos( h_nuParent_nuParentPt,parentPdg,parentPt,mc,wgt);

	if(parentPdg==211) utils->fillHistos( h_pion_parPz_parPt,pdpz/1000.,parentPt,mc,wgt);
	else if(parentPdg==321) utils->fillHistos( h_kaon_parPz_parPt,pdpz/1000.,parentPt,mc,wgt);
	else cout << "Found a parent which is not a pion or kaon with pdg="<<parentPdg << endl;

      }
      //------------------------------------------------------
      //Fill 2D Vertical Error Bands
      //Error Bands are not filled up for tmu_costheta 
      //and pmu_ptmu histos to avoid code slow down 
      //------------------------------------------------------
      if(RunCodeWithSystematics){
	utils->fillVertErrorBands( h_E_etheta,mc->CCQENu_leptonE[3]/1000,etheta2, mc );	    
      }//end run with systematics  
      // End of fillHistos and fillVertErrorBands where CV of events pass the proton and recoil cuts 
    }//End CV cuts
    //---------------------------------------------------------------------------------------
    //Now fillLatErrorBands for events whose CV did not pass the proton and recoil cuts  
    //What if the shifted proton score and recoil pass the selection cuts ? 
    //Those events should also be included in the LatErrorBands
    //fillLatErrorBands have their own checks for shifted proton scores and recoils  
    //---------------------------------------------------------------------------------------

    //----------------------------------------------------
    //Fill Lateral Error Bands
    //----------------------------------------------------
    if(RunCodeWithSystematics){
      utils->fillLatErrorBands( h_E_etheta, "emu", "etheta", mc->CCQENu_leptonE[3]/1000,etheta2 , mc, true );
    }//end run with systematics

  } //End of for loop over MC entries 
  
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

  //---------------------------------------------
  // Fill DATA Histograms
  //---------------------------------------------
    for (int i=0; i<entries_data ; ++i) {
      if( (i+1)%10000 == 0 ) 
	cout << "Reading DATA entry : " << i+1 << " - " << 100*(i+1)/entries_data << " % " << endl; 

      data->GetEntry(i);

      if( sample == "Signal") {
	if( multiplicity==0) { //1+2tracks
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) ) 	    continue;
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	}
      }
      else{
	cout << "You gave an invalid sample" << endl;
	return 1;
      }
      
      data_events++; 
      //-----------
      // This is the old method. Because of angular bias introduced from hadronic energy we need to pick out the right angle from a downstream node
      // New method until general reco is fixed is to get the reco P, a downstream theta (already corrected to the beam direction), and calculate
      // the pt and pz from these two quantities
      //-----------
      
      //double muon_theta     =  data->CCQENu_muon_theta;
      double muon_theta = data->muon_theta;
      double cos_muon_theta =  cos(muon_theta);

      double muon_T = data->CCQENu_muon_T / pow(10,3); //GeV
      
      double reco_muon_px   =  data->CCQENu_leptonE[0];
      double reco_muon_py   =  data->CCQENu_leptonE[1];
      double reco_muon_pz   =  data->CCQENu_leptonE[2];
      
      double muon_pt_beam   = utils->GetTransverseMomentumWRTBeam( reco_muon_px, reco_muon_py, reco_muon_pz ) / pow(10,3); //GeV 
      double muon_pz_beam   = utils->GetLongitudinalMomentumWRTBeam( reco_muon_px, reco_muon_py, reco_muon_pz ) / pow(10,3); //GeV 
      
      double muon_p    = utils->GetTotalMomentum( reco_muon_px, reco_muon_py, reco_muon_pz) / pow(10,3); //GeV

      //Change muon momentum to GeV
      reco_muon_px /= pow(10,3);
      reco_muon_py /= pow(10,3);
      reco_muon_pz /= pow(10,3);
      //reconstructed invariant mass
      double q2_gen = 2.0*(data->CCQENu_leptonE[3]+data->recoil_energy_nonmuon_nonvtx0mm-muon_p*1000.0*cos(muon_theta))-(105*105);
      double reco_w = (938.*938.)+(2*938.*data->recoil_energy_nonmuon_nonvtx0mm)-q2_gen;
      reco_w = (reco_w<0) ? 0.0 : sqrt(reco_w);

      //recoil variants
      double recoil = data->nonvtx_iso_blobs_energy + data->dis_id_energy;
      double recoil_inc = data->recoil_summed_energy[0];

      //vtx energy
      double vtx_energy = data->vtx_blobs_energy;

      //vtx energy fration
      double vtx_energy_fraction = vtx_energy/(vtx_energy+recoil);
      if(vtx_energy+recoil <=0) vtx_energy_fraction=-0.005; //When there is absolutely no reconstructed recoil we need to avoid NAN.
      //Neutrino Energies
      double nu_cal_E = (data->CCQENu_muon_E+recoil_inc)/1e3;
      double nu_qe_E = data->CCQENu_enu_muon/1e3;

      double etheta2 =  muon_theta*muon_theta*data->CCQENu_leptonE[3]/1000;

      vector<double> val_vect_imd;
      val_vect_imd.push_back(vtx_energy);
      val_vect_imd.push_back(etheta2);
      val_vect_imd.push_back(recoil_inc);
      float globx_long_imd = imd3d->GetBin(val_vect_imd).first;

      if(recoil_inc<80 && vtx_energy <10) utils->fillHistos( h_E_etheta,data->CCQENu_leptonE[3]/1000,etheta2);
      //utils->fillHistos( h_vtxetfull,globx_long_imd+0.0001,etheta2);
      

    } //End of for loop over DATA entries 
    
    cout << "Finished looping over all DATA entries. Number of DATA events satisfying all selection criteria = " << data_events << endl; 
    delete data; 
    delete tree_data; 
  }

  
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

  for( unsigned int i=kData; i<nHistosInc; ++i ){
    h_vtxetfull[i]->Write();
    h_E_etheta[i]->Write();
    h_nuParent_nuParentP[i]->Write();
    h_nuParent_nuParentPt[i]->Write();
    h_pion_parPz_parPt[i]->Write();
    h_kaon_parPz_parPt[i]->Write();
  }
   
  f->Close();
  delete f;
  delete utils;
  delete cutter; 
  delete binner;
  delete ccqenubinner; 

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
      "\t-Selection_Signal\t =\t Can be: Signal \n"<<
      "\t-Make_Flux_Constraint_Histo\t =\t If TRUE Enter 1; If FALSE Enter 0 \n"<<
      "\t-Multiplicity\t =\t Enter 0, 1 or 2. Here 0: 1 or 2 tracks; 1: 1-track events only; 2: 2-track events only \n"<<
      "\t-Number_MC_files\t =\t Number of MonteCarlo files. To use all files, set this to: -1 \n"<<
      "\t-Number_DATA_files\t =\t Number of Data files. To use all files, set this to: -1" << std::endl;
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

  return MuonSelectionHists(par[1], par[3], fluxConstraintHistoNeeded, atoi(par[5].c_str()), par[2],atoi(par[6].c_str()), atoi(par[7].c_str()) );
}
