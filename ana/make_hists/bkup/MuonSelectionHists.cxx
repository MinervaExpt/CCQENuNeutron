#include "include/CCQENuUtils.h"
#include "TParameter.h" 
#include "PlotUtils/HyperDimLinearizer.h"
using namespace CCQENU_ANA;
using namespace Neutron_ANA;

int MuonSelectionHists(string filename, string sample, bool makeFluxConstraintHisto, int multiplicity, string playlist, int n_mcfiles = -1, int n_datafiles = -1, int xtype=-1, int ztype=-1)
{
  //---------------------------------------------
  // create file to store histograms
  //---------------------------------------------
  //TFile *f = new TFile( filename.c_str(), "RECREATE" );
  bool typeCheck = xtype==-1 && ztype==-1;//Do normal
  CCQENuUtils  *utils  = new CCQENuUtils( false, makeFluxConstraintHisto );
  utils->setPlaylist(playlist);
  utils->setFluxReweighterPlaylist();
  CCQENuCuts   *cutter = new CCQENuCuts();
  AnaBinning *binner = new AnaBinning();
  CCQENuBinning *minmodbinner = new CCQENuBinning();

  NeutronBlobBinning *neutbinner = new NeutronBlobBinning();
  GeoUtils *geoUtils = new GeoUtils();

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
  axis_binning dalphatbins    = minmodbinner->GetDalphaT();
  axis_binning dphitbins    = minmodbinner->GetDphiT();
  axis_binning pnbins    = minmodbinner->GetPn();
  axis_binning dptbins    = minmodbinner->GetDeltaPt();
  axis_binning dptxbins    = minmodbinner->GetDeltaPtx();
  axis_binning dptybins    = minmodbinner->GetDeltaPty();
  axis_binning protonThetabins = minmodbinner->GetProtonThetaBins();
  axis_binning protonKineticbins = minmodbinner->GetProtonEnergyBinsGeV();
  axis_binning signbins          = minmodbinner->GetSignedBins();
  axis_binning signeddalphatbins = minmodbinner->GetSignedDeltaAlphaT();
  axis_binning signeddphitbins = minmodbinner->GetSignedDeltaPhiT();


  axis_binning dthetaPerpbins = neutbinner->GetThetaPerpBinsDegree();
  axis_binning dthetaReactbins = neutbinner->GetThetaReactionPlaneBinsDegree();
  
  axis_binning dalphatbins_TKI    = minmodbinner->GetDalphaT_TKI();
  axis_binning dphitbins_TKI    = minmodbinner->GetDphiT_TKI();
  axis_binning pnbins_TKI    = minmodbinner->GetPn_TKI();
  axis_binning dptbins_TKI    = minmodbinner->GetDeltaPt_TKI();
  axis_binning dptxbins_TKI    = minmodbinner->GetDeltaPtx_TKI();
  axis_binning dptybins_TKI    = minmodbinner->GetDeltaPty_TKI();
  axis_binning protonThetabins_TKI = minmodbinner->GetProtonThetaBins_TKI();
  axis_binning protonKineticbins_TKI = minmodbinner->GetProtonEnergyBinsGeV_TKI();
  axis_binning signeddalphatbins_TKI = minmodbinner->GetSignedDeltaAlphaT_TKI(); 
  axis_binning signeddphitbins_TKI = minmodbinner->GetSignedDeltaPhiT_TKI();




  //---------------------------------------------
  // Book Histograms
  //---------------------------------------------
  //

  //1D
  MnvH1D *h_multiplicity[nHistos],*h_nodes[nHistos],*h_nodes_long[nHistos],*h_nodes_anchored[nHistos],*h_nodes_vest[nHistos],*h_ptmu[nHistos];

  utils->bookHistos( h_multiplicity, "h_multiplicity", "Multiplicity of Outgoing Tracks; Number of Outgoing Tracks; Number of Events", 9, 0., 9.);
  utils->bookHistos( h_nodes, "h_nodes", "nodes",50,0,50);
  utils->bookHistos( h_nodes_long, "h_nodes_long", "nodes",50,0,50);
  utils->bookHistos( h_nodes_anchored, "h_nodes_anchored", "nodes",50,0,50);
  utils->bookHistos( h_nodes_vest, "h_nodes_vest", "nodes",50,0,50);
  utils->bookHistos( h_ptmu, "h_ptmu", Form( "Muon Pt;Reconstructed Muon p_{T} (GeV);Events / %.3f GeV", muonPtbins.default_width ), muonPtbins );
  //2-D
  MnvH2D *h_dalphat_ptmu[nHistos],*h_dphit_ptmu[nHistos],*h_pn_ptmu[nHistos],*h_dpt_ptmu[nHistos],*h_dptx_ptmu[nHistos],*h_dpty_ptmu[nHistos],*h_tp_ptmu[nHistos],*h_ptheta_ptmu[nHistos],*h_sign_ptmu[nHistos],*h_signdalphat_ptmu[nHistos],*h_signdphit_ptmu[nHistos];
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
  utils->bookHistos( h_sign_ptmu, "h_signed_ptmu","sign:ptmu",signbins,muonPtbins);
  utils->bookHistos( h_signdalphat_ptmu,"h_signeddalphat_ptmu","sdalphat:ptm",signeddalphatbins,muonPtbins);
  utils->bookHistos( h_signdphit_ptmu,"h_signeddphit_ptmu","sdphit:ptm",signeddphitbins,muonPtbins);

  MnvH2D *h_dthetaR_ptmu[nHistos], *h_dthetaP_ptmu[nHistos], *h_dPP_ptmu[nHistos], *h_dPR_ptmu[nHistos], *h_dPPi_ptmu[nHistos], *h_dPRi_ptmu[nHistos];

  utils->bookHistos( h_dthetaR_ptmu, "h_dthetaR_ptmu", "dthetaR:ptmu", dthetaReactbins, muonPtbins);
  utils->bookHistos( h_dthetaP_ptmu, "h_dthetaP_ptmu", "dthetaP:ptmu", dthetaPerpbins, muonPtbins);
  utils->bookHistos( h_dPR_ptmu, "h_dPR_ptmu", "dPR:ptmu", dptybins, muonPtbins);
  utils->bookHistos( h_dPP_ptmu, "h_dPP_ptmu", "dPP:ptmu", dptxbins, muonPtbins);
  utils->bookHistos( h_dPRi_ptmu, "h_dPRi_ptmu", "dPRi:ptmu", dptybins, muonPtbins);
  utils->bookHistos( h_dPPi_ptmu, "h_dPPi_ptmu", "dPPi:ptmu", dptxbins, muonPtbins);

  utils->bookHistos( h_node_nodeenergy, "h_node_nodeenergy","node:energy",50,0,50,500,0,50);
  utils->bookHistos( h_node_nodeenergy_long, "h_node_nodeenergy_long","node:energy",50,0,50,500,0,50);
  utils->bookHistos( h_node_nodeenergy_anchored, "h_node_nodeenergy_anchored","node:energy",50,0,50,500,0,50);
  utils->bookHistos( h_node_nodeenergy_vest, "h_node_nodeenergy_vest","node:energy",50,0,50,500,0,50);

  utils->bookHistos( h_protonRes_nodeenergy_node01, "h_protonRes_nodeenergy_node01","res:node",500,0,500,100,-1,1);
  utils->bookHistos( h_protonRes_nodeenergy_node2, "h_protonRes_nodeenergy_node2","res:node",500,0,500,100,-1,1);
  utils->bookHistos( h_protonRes_nodeenergy_node3, "h_protonRes_nodeenergy_node3","res:node",500,0,500,100,-1,1);
  utils->bookHistos( h_protonRes_nodeenergy_node4, "h_protonRes_nodeenergy_node4","res:node",500,0,500,100,-1,1);
  utils->bookHistos( h_protonRes_nodeenergy_node5, "h_protonRes_nodeenergy_node5","res:node",500,0,500,100,-1,1);

  //Special 3D for TKI vs TKI
  string spec_title = "";
  vector<vector<double> > full3D;
  //Now to pick the variables
  map<int,vector<double> > axis_bin_map;
  map<int,string > axis_name;
  axis_bin_map[0] = muonPzbins.bin_edges;
  axis_bin_map[1] = dalphatbins_TKI.bin_edges;
  axis_bin_map[2] = dphitbins_TKI.bin_edges;
  axis_bin_map[3] = pnbins_TKI.bin_edges;
  axis_bin_map[4] = dptbins_TKI.bin_edges;
  axis_bin_map[5] = dptxbins_TKI.bin_edges;
  axis_bin_map[6] = dptybins_TKI.bin_edges;
  axis_bin_map[7] = protonKineticbins_TKI.bin_edges;
  axis_bin_map[8] = protonThetabins_TKI.bin_edges;
  axis_bin_map[9] = signbins.bin_edges;
  axis_bin_map[10] = signeddalphatbins_TKI.bin_edges;
  axis_bin_map[11] = signeddphitbins_TKI.bin_edges;

  axis_name[0] = "pzmu";
  axis_name[1] = "dalphat";
  axis_name[2] = "dphit";
  axis_name[3] = "pn";
  axis_name[4] = "dpt";
  axis_name[5] = "dptx";
  axis_name[6] = "dpty";
  axis_name[7] = "tp";
  axis_name[8] = "prottheta";
  axis_name[9] = "sign";
  axis_name[10] = "signeddalphat";
  axis_name[11] = "signeddphit";

  vector<double> dummy;
  dummy.push_back(0);
  dummy.push_back(1);
  if(!typeCheck){
    spec_title = "h_"+axis_name[xtype]+"_"+axis_name[ztype];
    cout << "I will be running with " << spec_title << endl;
    full3D.push_back(axis_bin_map[xtype]);
    full3D.push_back(muonPtbins.bin_edges);
    full3D.push_back(axis_bin_map[ztype]);
  }
  else{
    full3D.push_back(dummy);
    full3D.push_back(dummy);
    full3D.push_back(dummy);
  }

  HyperDimLinearizer *my3d = new HyperDimLinearizer(full3D,0);
  utils->SetHyperDim(my3d);  
  vector<double> bins_x_global_long;
  int n_bins_long =2;
  if(!typeCheck) n_bins_long = (axis_bin_map[xtype].size()+1)*(axis_bin_map[ztype].size()+1);
  for(int i=0;i<n_bins_long;i++) bins_x_global_long.push_back(i);
  axis_binning glob_x_long;
  glob_x_long.bin_edges = bins_x_global_long;
  glob_x_long.nbins = bins_x_global_long.size()-1;
  glob_x_long.min = bins_x_global_long.front();
  glob_x_long.max = bins_x_global_long.back();
  
  MnvH2D *h_ptmu_tkitki[nHistos];
  utils->bookHistos( h_ptmu_tkitki,spec_title.c_str(),"pt vs globx long: GlobLong:Pt",glob_x_long,muonPtbins);


  
  

  //--------------------------------------------------------------
  // Add Vertical and Lateral Error Bands
  // JO is removing error bands from all 1D because they are all 
  // present in 2D. Plz make projections directly from the 2D. 
  //--------------------------------------------------------------
  if(RunCodeWithSystematics){
    if(typeCheck){
      utils->addVertErrorBands( h_ptmu );
      utils->addLatErrorBands ( h_ptmu );
      
      //2-D
      utils->addVertErrorBands( h_dalphat_ptmu );
      utils->addVertErrorBands( h_dphit_ptmu );
      utils->addVertErrorBands( h_pn_ptmu );
      utils->addVertErrorBands( h_dpt_ptmu );
      utils->addVertErrorBands( h_dptx_ptmu );
      utils->addVertErrorBands( h_dpty_ptmu );
      utils->addVertErrorBands( h_tp_ptmu );
      utils->addVertErrorBands( h_ptheta_ptmu );
      utils->addVertErrorBands( h_sign_ptmu );
      utils->addVertErrorBands( h_signdalphat_ptmu );
      utils->addVertErrorBands( h_signdphit_ptmu );
      utils->addVertErrorBands( h_dthetaR_ptmu );
      utils->addVertErrorBands( h_dthetaP_ptmu );
      utils->addVertErrorBands( h_dPP_ptmu );
      utils->addVertErrorBands( h_dPR_ptmu );
      utils->addVertErrorBands( h_dPPi_ptmu );
      utils->addVertErrorBands( h_dPRi_ptmu );
      
      utils->addLatErrorBands( h_dalphat_ptmu );
      utils->addLatErrorBands( h_dphit_ptmu );
      utils->addLatErrorBands( h_pn_ptmu );
      utils->addLatErrorBands( h_dpt_ptmu );
      utils->addLatErrorBands( h_dptx_ptmu );
      utils->addLatErrorBands( h_dpty_ptmu );
      utils->addLatErrorBands( h_tp_ptmu );
      utils->addLatErrorBands( h_ptheta_ptmu );
      utils->addLatErrorBands( h_sign_ptmu );
      utils->addLatErrorBands( h_signdalphat_ptmu );
      utils->addLatErrorBands( h_signdphit_ptmu );
      utils->addLatErrorBands( h_dthetaR_ptmu );
      utils->addLatErrorBands( h_dthetaP_ptmu );
      utils->addLatErrorBands( h_dPP_ptmu );
      utils->addLatErrorBands( h_dPR_ptmu );
      utils->addLatErrorBands( h_dPPi_ptmu );
      utils->addLatErrorBands( h_dPRi_ptmu );
 
    }
    else{
      utils->addVertErrorBands( h_ptmu_tkitki );
      utils->addLatErrorBands( h_ptmu_tkitki );
    }
  }//end run with systematics
  //  return 0;
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
  double mc_events = 0., data_events = 0.;
  int mc_nocuts_1track = 0, mc_cut1_1track = 0, mc_cut2_1track =0, mc_cut3_1track = 0, mc_cut4_1track = 0, mc_cut5_1track = 0, mc_cut6_1track = 0, mc_cut6_1track_sideband = 0;
  int mc_nocuts_2track = 0, mc_cut1_2track = 0, mc_cut2_2track =0, mc_cut3_2track = 0, mc_cut4_2track = 0, mc_cut5_2track = 0, mc_cut6_2track = 0, mc_cut7_2track = 0, mc_cut7_2track_sideband = 0;
  int data_nocuts_1track = 0, data_cut1_1track = 0, data_cut2_1track =0, data_cut3_1track = 0, data_cut4_1track = 0, data_cut5_1track = 0, data_cut6_1track = 0, data_cut6_1track_sideband = 0;
  int data_nocuts_2track = 0, data_cut1_2track = 0, data_cut2_2track =0, data_cut3_2track = 0, data_cut4_2track = 0, data_cut5_2track = 0, data_cut6_2track = 0, data_cut7_2track = 0, data_cut7_2track_sideband = 0;
  
  for (int i=0; i<entries_mc ; ++i) {
    if( (i+1)%100000 == 0 ) 
      cout << "Reading MC entry : " << i+1 << " - " << 100*(i+1)/entries_mc << " % " << endl; 
    mc->GetEntry(i);    
    if ( sample == "Signal") {
      if( multiplicity==0) { //1+2tracks
        if( !cutter->passInteractionVertex( mc ) )       continue;	
        if( !cutter->passCCQESelection(mc) )             continue;
      }
      else if( multiplicity==1) { //1 track only 
        if( !cutter->passInteractionVertex( mc ) )       continue;	mc_nocuts_1track++; 
        if( !cutter->passDeadTimeCuts( mc ) )            continue;	mc_cut1_1track++; 
        if( !cutter->passNuHelicityCut( mc ) )   	       continue;  
        if( !(mc->multiplicity==1) )                     continue;	mc_cut3_1track++; 
        if( !cutter->passImprovedMichelCut( mc ) )       continue;	mc_cut4_1track++; 
        if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( !cutter->passNBlobs( mc ) )                  continue;	mc_cut5_1track++; 
        if ( !cutter->passTrackAngleCut( mc ) )          continue;	mc_cut6_1track++; 
      }
      else if( multiplicity==2) { //2 tracks only 
        if( !cutter->passInteractionVertex( mc ) )       continue; 	mc_nocuts_2track++;  
        if( !cutter->passDeadTimeCuts( mc ) )            continue;	mc_cut1_2track++; 
        if( !cutter->passNuHelicityCut( mc ) )           continue;	mc_cut2_2track++; 
        if( !(mc->multiplicity==2) )                     continue;	mc_cut3_2track++; 
        if( !cutter->passImprovedMichelCut( mc ) )       continue;	mc_cut5_2track++; 
        if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( !cutter->passNBlobs( mc ) )                  continue;	mc_cut6_2track++; 
        if( !cutter->passTrackAngleCut( mc ) )           continue;	mc_cut7_2track++; 
      }
    }
    else if( sample == "MichelSideBand") {
      if( multiplicity==0) {//1+2tracks
        if( !cutter->passInteractionVertex( mc ) )       continue;	
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !cutter->passMultiplicity( mc ) )            continue;
        if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( !cutter->passNBlobs( mc ) )                  continue;
        if( cutter->passImprovedMichelCut( mc ) )        continue;
        if( !cutter->passTrackAngleCut( mc ) )           continue;
      }
      else if( multiplicity==1) { //1 track only 
	if( !cutter->passInteractionVertex( mc ) )       continue;	
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !(mc->multiplicity==1) )                     continue;
        if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( !cutter->passNBlobs( mc ) )                  continue;
        if( cutter->passImprovedMichelCut( mc ) )        continue;
        if( !cutter->passTrackAngleCut( mc ) )           continue;
      }
      else if( multiplicity==2) { //2 tracks only 
        if( !cutter->passInteractionVertex( mc ) )       continue;	
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !(mc->multiplicity==2) )                     continue;
        if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( !cutter->passNBlobs( mc ) )                  continue;
        if( cutter->passImprovedMichelCut( mc ) )        continue;
        if( !cutter->passTrackAngleCut( mc ) )           continue;
      }
    }
    else if( sample == "BlobSideBand") {
      if( multiplicity==0) {//1+2tracks
        if( !cutter->passInteractionVertex( mc ) )       continue;	
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !cutter->passMultiplicity( mc ) )            continue;
        if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( cutter->passNBlobs( mc ) )                   continue;
        if( !cutter->passImprovedMichelCut( mc ) )       continue;
        if( !cutter->passTrackAngleCut( mc ) )           continue;
      }
      else if( multiplicity==1) { //1 track only 
        if( !cutter->passInteractionVertex( mc ) )       continue;	
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !(mc->multiplicity==1) )                     continue;
        if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( cutter->passNBlobs( mc ) )                   continue;
        if( !cutter->passImprovedMichelCut( mc ) )       continue;
        if( !cutter->passTrackAngleCut( mc ) )           continue;
      }
      else if( multiplicity==2) { //2 tracks only 
        if( !cutter->passInteractionVertex( mc ) )       continue;	
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !(mc->multiplicity==2) )                     continue;
        if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( cutter->passNBlobs( mc ) )                   continue;
        if( !cutter->passImprovedMichelCut( mc ) )       continue;
        if( !cutter->passTrackAngleCut( mc ) )           continue;
      }
    }
    else if( sample == "MicBlobSideBand") {
      if( multiplicity==0) {//1+2tracks
        if( !cutter->passInteractionVertex( mc ) )       continue;
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !cutter->passMultiplicity( mc ) )            continue;
        if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( cutter->passNBlobs( mc ) )                   continue;
        if( cutter->passImprovedMichelCut( mc ) )        continue;
        if( !cutter->passTrackAngleCut( mc ) )           continue;
      }
      else if( multiplicity==1) { //1 track only 
        if( !cutter->passInteractionVertex( mc ) )       continue;
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !(mc->multiplicity==1) )                     continue;
        if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( cutter->passNBlobs( mc ) )                   continue;
        if( cutter->passImprovedMichelCut( mc ) )        continue;
        if( !cutter->passTrackAngleCut( mc ) )          continue;
      }
      else if( multiplicity==2) { //2 tracks only      
        if( !cutter->passInteractionVertex( mc ) )       continue;
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !(mc->multiplicity==2) )                     continue;
        if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( cutter->passNBlobs( mc ) )                   continue;
        if( cutter->passImprovedMichelCut( mc ) )        continue;
        if( !cutter->passTrackAngleCut( mc ) )           continue;
      }
    } else {
      cout<<"Wrong Sample Selection."<<endl;
      cout<<"Select: Signal, MichelSideBand, BlobSideBand, MicBlobSideBand"<<endl;
      continue;
    }

    //Get CVweight
    double wgt = utils->GetCVWeight( mc , sample);    
    mc_events+=wgt;
    //-----------
    //Muon Variables
    double muon_theta     = mc->muon_theta_allNodes[20]; // get theta from 20th node, note NX kludge fixes this...
    double cos_muon_theta =  cos(muon_theta);
    double muon_T         =  mc->CCQENu_muon_T / pow(10,3); //GeV
    double reco_muon_px   =  mc->CCQENu_leptonE[0];
    double reco_muon_py   =  mc->CCQENu_leptonE[1];
    double reco_muon_pz   =  mc->CCQENu_leptonE[2];
    double muon_p         = utils->GetTotalMomentum( reco_muon_px, reco_muon_py, reco_muon_pz ) / pow(10,3); //GeV
    //New method (maybe temporary)
    double muon_pt_beam = sin(muon_theta)*muon_p;
    double muon_pz_beam = cos(muon_theta)*muon_p;

    //Proton Variables
    TVector3 protonVect(mc->CCQENu_proton_Px_fromdEdx,mc->CCQENu_proton_Py_fromdEdx,mc->CCQENu_proton_Pz_fromdEdx);//MeV
    double protonAngle = mc->CCQENu_proton_theta*180./3.14159;
    double protonMom = protonVect.Mag();//MeV
    double protonTn = TMath::Sqrt(protonMom*protonMom+utils->M_p*utils->M_p)-utils->M_p;//Done in MeV

    //Truth Proton Variables
    //TVector3 true_proton_p_comp = utils->GetHighestTrueProtonMomentumComp(mc);
    TVector3 true_proton_p_comp(mc->proton_prong_4p[0],mc->proton_prong_4p[1],mc->proton_prong_4p[2]);
    double true_protonMom = true_proton_p_comp.Mag();
    double proton_res = 1-protonMom/true_protonMom;
    
    //Convert to degrees, GeV
    muon_theta*= 180. / 3.14159;
    reco_muon_px /= pow(10,3);
    reco_muon_py /= pow(10,3);
    reco_muon_pz /= pow(10,3);
    protonVect *= pow(10,-3);
    protonTn /= pow(10,3);


    double reco_dalphat = utils->GetDeltaAlphaT(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
    double reco_dphit = utils->GetDeltaPhiT(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
    double reco_pn = utils->GetPn(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
    double reco_dpt = utils->GetDeltaPt(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
    double reco_dptx = utils->GetDeltaPtx(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
    double reco_dpty = utils->GetDeltaPty(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
    double reco_sign = utils->GetMuonProtonSign(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
    double reco_signed_dalphat = utils->GetSignedDeltaAlphaT(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
    double reco_signed_dphit = utils->GetSignedDeltaPhiT(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);


    //Directional "transverse" variables, Tejin
    // beam_bias default to 0, using this method allows calculation of angle uncertainties. 
    // 1. Define the beam direction, and muon 4P
    double beam_bias = 0;
    XYZVector beam = geoUtils->BeamAxis( beam_bias ); 
    XYZTVector reco_muon_4P( reco_muon_px, reco_muon_py, reco_muon_pz ,mc->CCQENu_muon_E/1e3);
    
    // 2. Compute expected stationary nucleon scattering momentum
    // XYZTVector ComputeExpectedNucleon ( XYZVector &nu, XYZTVector &muon, double ISMass = 0.93827231,double FSMass =   0.93956563, double BindingE = 0.00 ); // vbar p --> mu n, Initial Mass, Final state mass
    double binding_e = 0; // theoretically we should include binding E for carbon.. but I'm comparing to hydrogen in anti-nu, so probably set it to 0.
    double is_neutron_mass = 0.9395654133;
    double fs_proton_mass = 0.93827231;

    XYZTVector expectedProton4P = geoUtils->ComputeExpectedNucleon( beam, reco_muon_4P, is_neutron_mass, fs_proton_mass, binding_e);
    XYZVector expectedProton3P = (XYZVector) expectedProton4P.Vect();

    // 3. 
    // ComputeNeutronAngularVars( XYZVector &nu, XYZVector &expVec, XYZVector &targetVec ), variables are : dThetaP, dThetaR, dThetaTotal,  dPp, dPr, inferred dPp and dPr
    XYZVector protonVect_XYZ( protonVect.X(), protonVect.Y(), protonVect.Z() );
    vector<double> reco_vars = geoUtils->ComputeNeutronAngularVars( beam, expectedProton3P, protonVect_XYZ );
    double reco_dThetaP = reco_vars[0]*189/TMath::Pi();
    double reco_dThetaR = reco_vars[1]*189/TMath::Pi();
    double reco_dTheta = reco_vars[2]*189/TMath::Pi();
    double reco_dPp = reco_vars[3];
    double reco_dPr = reco_vars[4];
    double reco_dPp_infer = reco_vars[5];
    double reco_dPr_infer = reco_vars[6];

    //I feel there is also a need  to include the left-right uncertainty of the beam, since LR-direction is pretty important....




    map<int,double> varmap;
    varmap[0] = muon_pz_beam;//muonPzbins;
    varmap[1] = reco_dalphat;//dalphatbins;
    varmap[2] = reco_dphit;//dphitbins;
    varmap[3] = reco_pn;//pnbins;
    varmap[4] = reco_dpt;//dptbins;
    varmap[5] = reco_dptx;//dptxbins;
    varmap[6] = reco_dpty;//dptybins;
    varmap[7] = protonTn;//protonKineticbins;
    varmap[8] = protonAngle;//protonThetabins;
    varmap[9] = reco_sign;//signbins;
    varmap[10] = reco_signed_dalphat;//signeddalphatbins;
    varmap[11] = reco_signed_dphit;//signeddphitbins;
    varmap[12] = reco_dThetaP;  
    varmap[13] = reco_dThetaR;  
    varmap[14] = reco_dTheta;
    varmap[15] = reco_dPp;      
    varmap[16] = reco_dPr;     
    varmap[17] = reco_dPp_infer;
    varmap[18] = reco_dPr_infer;

 

    vector<double> val_vect;
    if(!typeCheck){
      val_vect.push_back(varmap[xtype]);
      val_vect.push_back(muon_pt_beam);
      val_vect.push_back(varmap[ztype]);
    }
    double globx_long = 0;
    if(!typeCheck)globx_long = my3d->GetBin(val_vect).first;
    
    

    //-------------------------------------------------------------------
    // Check that the CV of events pass the proton and recoil cuts 
    // This is because we have systematics due to these cuts 
    // and a change in their values can affect the composition 
    // of the selected sample. 
    //-------------------------------------------------------------------
    bool doesCVpass_SingleProtonCut = (mc->multiplicity==1)? true : cutter->passSingleProtonCut( mc, 0, 0 );
    bool doesCVpass_ExtraProtonsCut = (mc->multiplicity==1)? true : cutter->passExtraProtonsCut( mc, NULL, 0 );
    bool doesCVpass_RecoilCut = cutter->passSignalFunction( mc, 0, 0 ); 

    if( doesCVpass_SingleProtonCut && doesCVpass_ExtraProtonsCut && doesCVpass_RecoilCut ) {     
      //---------------------------------------------
      //Fill 1-D Plots
      //---------------------------------------------
      if(typeCheck){
	utils->fillHistos( h_multiplicity, mc->multiplicity, mc, wgt ); 
	utils->fillHistos( h_ptmu, muon_pt_beam, mc, wgt );
	
	int n_nodes = mc->CCQENu_proton_nodes_nodesNormE_sz;
	int pattrec = mc->CCQENu_proton_patternRec;
	utils->fillHistos(h_nodes,n_nodes,mc,wgt);
	if(pattrec==1 || pattrec==2)utils->fillHistos(h_nodes_long,n_nodes,mc,wgt);
	else if(pattrec==3)         utils->fillHistos(h_nodes_anchored,n_nodes,mc,wgt);
	else if (pattrec==6)        utils->fillHistos(h_nodes_vest,n_nodes,mc,wgt);
	
	double summed_01 = 0.0;
	for(int node=0;node<n_nodes;node++){
	  utils->fillHistos(h_node_nodeenergy,node,mc->CCQENu_proton_nodes_nodesNormE[node],mc,wgt);
	  if(pattrec==1||pattrec==2)utils->fillHistos(h_node_nodeenergy_long,node,mc->CCQENu_proton_nodes_nodesNormE[node],mc,wgt);
	  else if(pattrec==3)utils->fillHistos(h_node_nodeenergy_anchored,node,mc->CCQENu_proton_nodes_nodesNormE[node],mc,wgt);
	  else if(pattrec==6) utils->fillHistos(h_node_nodeenergy_vest,node,mc->CCQENu_proton_nodes_nodesNormE[node],mc,wgt);
	  
	  if(node==0) summed_01+=mc->CCQENu_proton_nodes_nodesNormE[node];
	  if(node==1){
	    summed_01+=mc->CCQENu_proton_nodes_nodesNormE[node];
	    utils->fillHistos(h_protonRes_nodeenergy_node01,summed_01,proton_res,mc,wgt);
	    summed_01=0;
	  }
	  if(node==2) utils->fillHistos(h_protonRes_nodeenergy_node2,mc->CCQENu_proton_nodes_nodesNormE[node],proton_res,mc,wgt);
	  if(node==3) utils->fillHistos(h_protonRes_nodeenergy_node3,mc->CCQENu_proton_nodes_nodesNormE[node],proton_res,mc,wgt);
	  if(node==4) utils->fillHistos(h_protonRes_nodeenergy_node4,mc->CCQENu_proton_nodes_nodesNormE[node],proton_res,mc,wgt);
	  if(node==5) utils->fillHistos(h_protonRes_nodeenergy_node5,mc->CCQENu_proton_nodes_nodesNormE[node],proton_res,mc,wgt);
	  
	  
	}
	
	//-------------------------------------------------------
	//Fill 2-D Plots
	//-------------------------------------------------------
	utils->fillHistos( h_dalphat_ptmu,reco_dalphat,muon_pt_beam,mc,wgt);
	utils->fillHistos( h_dphit_ptmu,reco_dphit,muon_pt_beam,mc,wgt);
	utils->fillHistos( h_pn_ptmu,reco_pn,muon_pt_beam,mc,wgt);
	utils->fillHistos( h_dpt_ptmu,reco_dpt,muon_pt_beam,mc,wgt);
	utils->fillHistos( h_dptx_ptmu,reco_dptx,muon_pt_beam,mc,wgt);
	utils->fillHistos( h_dpty_ptmu,reco_dpty,muon_pt_beam,mc,wgt);
	utils->fillHistos( h_tp_ptmu,protonTn, muon_pt_beam,mc,wgt);
	utils->fillHistos( h_ptheta_ptmu,protonAngle, muon_pt_beam,mc,wgt);
	utils->fillHistos( h_sign_ptmu,reco_sign,muon_pt_beam,mc,wgt);
	utils->fillHistos( h_signdalphat_ptmu,reco_signed_dalphat,muon_pt_beam,mc,wgt);
	utils->fillHistos( h_signdphit_ptmu,reco_signed_dphit,muon_pt_beam,mc,wgt);

	utils->fillHistos( h_dthetaP_ptmu,reco_dThetaP,muon_pt_beam,mc,wgt);
	utils->fillHistos( h_dthetaR_ptmu,reco_dThetaR,muon_pt_beam,mc,wgt);
	utils->fillHistos( h_dPP_ptmu,reco_dPp,muon_pt_beam,mc,wgt);
	utils->fillHistos( h_dPR_ptmu,reco_dPr,muon_pt_beam,mc,wgt);
	utils->fillHistos( h_dPPi_ptmu,reco_dPp_infer,muon_pt_beam,mc,wgt);
	utils->fillHistos( h_dPRi_ptmu,reco_dPr_infer,muon_pt_beam,mc,wgt);
      }
      else{
	utils->fillHistos( h_ptmu_tkitki,globx_long+0.0001, muon_pt_beam,mc,wgt);
      }
      //------------------------------------------------------
      //Fill 2D Vertical Error Bands
      //Error Bands are not filled up for tmu_costheta 
      //and pmu_ptmu histos to avoid code slow down 
      //------------------------------------------------------
      if(RunCodeWithSystematics){
	if(typeCheck){
	  utils->fillVertErrorBands( h_ptmu,  muon_pt_beam, mc );
	  
	  utils->fillVertErrorBands( h_dalphat_ptmu,reco_dalphat,muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_dphit_ptmu,reco_dphit,muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_pn_ptmu,reco_pn,muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_dpt_ptmu,reco_dpt,muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_dptx_ptmu,reco_dptx,muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_dpty_ptmu,reco_dpty,muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_tp_ptmu,protonTn, muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_ptheta_ptmu,protonAngle, muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_sign_ptmu,reco_sign,muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_signdalphat_ptmu,reco_signed_dalphat,muon_pt_beam,mc);
	  utils->fillVertErrorBands( h_signdphit_ptmu,reco_signed_dphit,muon_pt_beam,mc);

    utils->fillVertErrorBands( h_dthetaP_ptmu,reco_dThetaP,muon_pt_beam,mc);
    utils->fillVertErrorBands( h_dthetaR_ptmu,reco_dThetaR,muon_pt_beam,mc);
    utils->fillVertErrorBands( h_dPP_ptmu,reco_dPp,muon_pt_beam,mc);
    utils->fillVertErrorBands( h_dPR_ptmu,reco_dPr,muon_pt_beam,mc);
    utils->fillVertErrorBands( h_dPPi_ptmu,reco_dPp_infer,muon_pt_beam,mc);
    utils->fillVertErrorBands( h_dPRi_ptmu,reco_dPr_infer,muon_pt_beam,mc);

    
	}
	else{
	  utils->fillVertErrorBands( h_ptmu_tkitki,globx_long+0.0001, muon_pt_beam,mc);
	}
      }//end run with systematics        
    } // End of fillHistos and fillVertErrorBands where CV of events pass the proton and recoil cuts 
  
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
    
    if(RunCodeWithSystematics){
      if(typeCheck){
	utils->fillLatErrorBands( h_ptmu,  "pT",    muon_pt_beam, mc, true );
      }

      //----------------------------------------------------
      //Fill 2D Lateral Error Bands
      //Error Bands are not filled up for tmu_costheta 
      //and pmu_ptmu histos to avoid code slow down 
      //----------------------------------------------------
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
	  if( !cutter->passInteractionVertex( data ) )     continue;	
	  if( !cutter->passCCQESelection( data ) )         continue;
	}
	else if( multiplicity==1) {
	  if( !cutter->passInteractionVertex( data ) )      continue;	  data_nocuts_1track++; 
	  if( !cutter->passDeadTimeCuts( data ) )           continue;	  data_cut1_1track++; 
	  if( !cutter->passNuHelicityCut( data ) ) 	    continue;     data_cut2_1track++; 
	  if( !(data->multiplicity==1) )                    continue;	  data_cut3_1track++; 
	  if( !cutter->passImprovedMichelCut( data ) )      continue; 	  data_cut4_1track++; 
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( !cutter->passNBlobs( data ) )                 continue;	  data_cut5_1track++; 
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;	  data_cut6_1track++; 
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	}
	else if( multiplicity==2) {
	  if( !cutter->passInteractionVertex( data ) )      continue;	  data_nocuts_2track++; 
	  if( !cutter->passDeadTimeCuts( data ) )           continue;	  data_cut1_2track++; 
	  if( !cutter->passNuHelicityCut( data ) )          continue;	  data_cut2_2track++; 
	  if( !(data->multiplicity==2) )                    continue;	  data_cut3_2track++; 
	  if( !cutter->passSingleProtonCut( data ) )        continue;
	  if( !cutter->passExtraProtonsCut( data ) )        continue;	  data_cut4_2track++; 
	  if( !cutter->passImprovedMichelCut( data ) )      continue;	  data_cut5_2track++; 
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( !cutter->passNBlobs( data ) )                 continue;	  data_cut6_2track++; 
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;	  data_cut7_2track++; 
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	}
      }
      else if( sample == "MichelSideBand") {
	if( multiplicity==0) {//1+2tracks
	  bool passProtonCut = (data->multiplicity==1)? true : cutter->passSingleProtonCut( data );
	  bool pass_ExtraProtonsCut = (data->multiplicity==1)? true : cutter->passExtraProtonsCut( data );
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !cutter->passMultiplicity( data ) )           continue;
	  if( !passProtonCut )                              continue;
	  if( !pass_ExtraProtonsCut )                       continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( !cutter->passNBlobs( data ) )                 continue;
	  if( cutter->passImprovedMichelCut( data ) )       continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	}
	else if( multiplicity==1) {
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !(data->multiplicity==1) )                    continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( !cutter->passNBlobs( data ) )                 continue;
	  if( cutter->passImprovedMichelCut( data ) )       continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	}
	else if( multiplicity==2) {
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !(data->multiplicity==2) )                    continue;
	  if( !cutter->passSingleProtonCut( data ) )        continue;
	  if( !cutter->passExtraProtonsCut( data ) )        continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( !cutter->passNBlobs( data ) )                 continue;
	  if( cutter->passImprovedMichelCut( data ) )       continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
	  if ( !cutter->passTrackAngleCut( data ) )         continue;
	}
      }

      else if( sample == "BlobSideBand") {
	if( multiplicity==0) {//1+2tracks
	  bool passProtonCut = (data->multiplicity==1)? true : cutter->passSingleProtonCut( data );
	  bool pass_ExtraProtonsCut = (data->multiplicity==1)? true : cutter->passExtraProtonsCut( data );
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !cutter->passMultiplicity( data ) )           continue;
	  if( !passProtonCut )                              continue;
	  if( !pass_ExtraProtonsCut )                       continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( cutter->passNBlobs( data ) )                  continue;
	  if( !cutter->passImprovedMichelCut( data ) )      continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	}
	else if( multiplicity==1) {
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !(data->multiplicity==1) )                    continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( cutter->passNBlobs( data ) )                  continue;
	  if( !cutter->passImprovedMichelCut( data ) )      continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
	  if ( !cutter->passTrackAngleCut( data ) )         continue;
	}
	else if( multiplicity==2) {
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !(data->multiplicity==2) )                    continue;
	  if( !cutter->passSingleProtonCut( data ) )        continue;
	  if( !cutter->passExtraProtonsCut( data ) )        continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( cutter->passNBlobs( data ) )                  continue;
	  if( !cutter->passImprovedMichelCut( data ) )      continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	}
      }
      else if( sample == "MicBlobSideBand") {
	if( multiplicity==0) {//1+2tracks
	  bool passProtonCut = (data->multiplicity==1)? true : cutter->passSingleProtonCut( data );
	  bool pass_ExtraProtonsCut = (data->multiplicity==1)? true : cutter->passExtraProtonsCut( data );
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !cutter->passMultiplicity( data ) )           continue;
	  if( !passProtonCut )                              continue;
	  if( !pass_ExtraProtonsCut )                       continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( cutter->passNBlobs( data ) )                  continue;
	  if( cutter->passImprovedMichelCut( data ) )       continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	}
	else if( multiplicity==1) {
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !(data->multiplicity==1) )                    continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( cutter->passNBlobs( data ) )                  continue;
	  if( cutter->passImprovedMichelCut( data ) )       continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	}
	else if( multiplicity==2) {
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !(data->multiplicity==2) )                    continue;
	  if( !cutter->passSingleProtonCut( data ) )        continue;
	  if( !cutter->passExtraProtonsCut( data ) )        continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( cutter->passNBlobs( data ) )                  continue;
	  if( cutter->passImprovedMichelCut( data ) )       continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
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
      
      //Muon Variables
      double muon_theta = data->muon_theta_allNodes[20];// get theta from 20th node
      double cos_muon_theta =  cos(muon_theta);
      double muon_T = data->CCQENu_muon_T / pow(10,3); //GeV
      double reco_muon_px   =  data->CCQENu_leptonE[0];
      double reco_muon_py   =  data->CCQENu_leptonE[1];
      double reco_muon_pz   =  data->CCQENu_leptonE[2];
      double muon_p    = utils->GetTotalMomentum( reco_muon_px, reco_muon_py, reco_muon_pz) / pow(10,3); //GeV
      double muon_pt_beam = sin(muon_theta)*muon_p;
      double muon_pz_beam = cos(muon_theta)*muon_p;

      //Proton Variables
      TVector3 protonVect(data->CCQENu_proton_Px_fromdEdx,data->CCQENu_proton_Py_fromdEdx,data->CCQENu_proton_Pz_fromdEdx);//MeV
      double protonAngle = data->CCQENu_proton_theta*180./3.14159;
      double protonMom = protonVect.Mag();//MeV
      double protonTn = TMath::Sqrt(protonMom*protonMom+utils->M_p*utils->M_p)-utils->M_p;//MeV

      //Convert to degrees, GeV
      muon_theta*= 180. / 3.14159;
      reco_muon_px /= pow(10,3);
      reco_muon_py /= pow(10,3);
      reco_muon_pz /= pow(10,3);
      protonVect *= pow(10,-3);
      protonTn /= pow(10,3);

      //Transverse variables all inputs are GeV
      double reco_dalphat = utils->GetDeltaAlphaT(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
      double reco_dphit = utils->GetDeltaPhiT(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
      double reco_pn = utils->GetPn(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
      double reco_dpt = utils->GetDeltaPt(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
      double reco_dptx = utils->GetDeltaPtx(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
      double reco_dpty = utils->GetDeltaPty(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
      double reco_sign = utils->GetMuonProtonSign(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
      double reco_signed_dalphat = utils->GetSignedDeltaAlphaT(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
      double reco_signed_dphit = utils->GetSignedDeltaPhiT(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);

      //Directional "transverse" variables, Tejin
      // beam_bias default to 0, using this method allows calculation of angle uncertainties. 
      // 1. Define the beam direction, and muon 4P
      double beam_bias = 0;
      XYZVector beam = geoUtils->BeamAxis( beam_bias ); 
      XYZTVector reco_muon_4P( reco_muon_px, reco_muon_py, reco_muon_pz ,data->CCQENu_muon_E/1e3);
      
      // 2. Compute expected stationary nucleon scattering momentum
      // XYZTVector ComputeExpectedNucleon ( XYZVector &nu, XYZTVector &muon, double ISMass = 0.93827231,double FSMass =   0.93956563, double BindingE = 0.00 ); // vbar p --> mu n, Initial Mass, Final state mass
      double binding_e = 0; // theoretically we should include binding E for carbon.. but I'm comparing to hydrogen in anti-nu, so probably set it to 0.
      double is_neutron_mass = 0.9395654133;
      double fs_proton_mass = 0.93827231;

      XYZTVector expectedProton4P = geoUtils->ComputeExpectedNucleon( beam, reco_muon_4P, is_neutron_mass, fs_proton_mass, binding_e);
      XYZVector expectedProton3P = (XYZVector) expectedProton4P.Vect();

      // 3. 
      // ComputeNeutronAngularVars( XYZVector &nu, XYZVector &expVec, XYZVector &targetVec ), variables are : dThetaP, dThetaR, dThetaTotal,  dPp, dPr, inferred dPp and dPr
      XYZVector protonVect_XYZ( protonVect.X(), protonVect.Y(), protonVect.Z() );
      vector<double> reco_vars = geoUtils->ComputeNeutronAngularVars( beam, expectedProton3P, protonVect_XYZ );
      double reco_dThetaP = reco_vars[0]*180/TMath::Pi();
      double reco_dThetaR = reco_vars[1]*180/TMath::Pi();
      double reco_dTheta = reco_vars[2]*180/TMath::Pi();
      double reco_dPp = reco_vars[3];
      double reco_dPr = reco_vars[4];
      double reco_dPp_infer = reco_vars[5];
      double reco_dPr_infer = reco_vars[6];

      //I feel there is also a need  to include the left-right uncertainty of the beam, since LR-direction is pretty important....




      map<int,double> varmap;
      varmap[0] = muon_pz_beam;//muonPzbins;
      varmap[1] = reco_dalphat;//dalphatbins;
      varmap[2] = reco_dphit;//dphitbins;
      varmap[3] = reco_pn;//pnbins;
      varmap[4] = reco_dpt;//dptbins;
      varmap[5] = reco_dptx;//dptxbins;
      varmap[6] = reco_dpty;//dptybins;
      varmap[7] = protonTn;//protonKineticbins;
      varmap[8] = protonAngle;//protonThetabins;
      varmap[9] = reco_sign;//signbins;
      varmap[10] = reco_signed_dalphat;//signeddalphatbins;
      varmap[11] = reco_signed_dphit;//signeddphitbins;
      varmap[12] = reco_dThetaP;  
      varmap[13] = reco_dThetaR;  
      varmap[14] = reco_dTheta;
      varmap[15] = reco_dPp;      
      varmap[16] = reco_dPr;     
      varmap[17] = reco_dPp_infer;
      varmap[18] = reco_dPr_infer;


      
      vector<double> val_vect;
      double globx_long = 0;
      if(!typeCheck){
	val_vect.push_back(varmap[xtype]);
	val_vect.push_back(muon_pt_beam);
	val_vect.push_back(varmap[ztype]);
	globx_long = my3d->GetBin(val_vect).first;
      }



      if(typeCheck){
	
	//FillHistograms
	utils->fillHistos( h_multiplicity, data->multiplicity );
	//2-D
	int n_nodes = data->CCQENu_proton_nodes_nodesNormE_sz;
	int pattrec = data->CCQENu_proton_patternRec;
	utils->fillHistos(h_nodes,n_nodes);
	if(pattrec==1 || pattrec==2)utils->fillHistos(h_nodes_long,n_nodes);
	else if(pattrec==3)         utils->fillHistos(h_nodes_anchored,n_nodes);
	else if (pattrec==6)        utils->fillHistos(h_nodes_vest,n_nodes);
	
	for(int node=0;node<n_nodes;node++){
	  utils->fillHistos(h_node_nodeenergy,node,data->CCQENu_proton_nodes_nodesNormE[node]);
	  if(pattrec==1||pattrec==2)utils->fillHistos(h_node_nodeenergy_long,node,data->CCQENu_proton_nodes_nodesNormE[node]);
	  else if(pattrec==3)utils->fillHistos(h_node_nodeenergy_anchored,node,data->CCQENu_proton_nodes_nodesNormE[node]);
	  else if(pattrec==6) utils->fillHistos(h_node_nodeenergy_vest,node,data->CCQENu_proton_nodes_nodesNormE[node]);
	}
	utils->fillHistos( h_ptmu, muon_pt_beam );
	
	
	utils->fillHistos( h_dalphat_ptmu,reco_dalphat,muon_pt_beam);
	utils->fillHistos( h_dphit_ptmu,reco_dphit,muon_pt_beam);
	utils->fillHistos( h_pn_ptmu,reco_pn,muon_pt_beam);
	utils->fillHistos( h_dpt_ptmu,reco_dpt,muon_pt_beam);
	utils->fillHistos( h_dptx_ptmu,reco_dptx,muon_pt_beam);
	utils->fillHistos( h_dpty_ptmu,reco_dpty,muon_pt_beam);
	utils->fillHistos( h_tp_ptmu,protonTn, muon_pt_beam);
	utils->fillHistos( h_ptheta_ptmu,protonAngle, muon_pt_beam);
	utils->fillHistos( h_sign_ptmu,reco_sign,muon_pt_beam);
	utils->fillHistos( h_signdalphat_ptmu,reco_signed_dalphat,muon_pt_beam);
	utils->fillHistos( h_signdphit_ptmu,reco_signed_dphit,muon_pt_beam);

	utils->fillHistos( h_dthetaP_ptmu,reco_dThetaP,muon_pt_beam);
	utils->fillHistos( h_dthetaR_ptmu,reco_dThetaR,muon_pt_beam);
	utils->fillHistos( h_dPP_ptmu,reco_dPp,muon_pt_beam);
	utils->fillHistos( h_dPR_ptmu,reco_dPr,muon_pt_beam);
	utils->fillHistos( h_dPPi_ptmu,reco_dPp_infer,muon_pt_beam);
	utils->fillHistos( h_dPRi_ptmu,reco_dPr_infer,muon_pt_beam);

      }
      else{
	utils->fillHistos( h_ptmu_tkitki,globx_long+0.0001, muon_pt_beam);
      }
    } //End of for loop over DATA entries 
    
    cout << "Finished looping over all DATA entries. Number of DATA events satisfying all selection criteria = " << data_events << endl; 
    delete data; 
    delete tree_data; 
  }

  //Print out some information
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
  
  //==================================================================
  // Create ROOT file to store histograms
  // Write to file all the created histograms 
  //==================================================================
  if(!typeCheck) filename = filename+"_"+axis_name[xtype]+"_"+axis_name[ztype]+".root";
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
    if(typeCheck){
      h_multiplicity[i]->Write();
      h_node_nodeenergy[i]->Write();
      h_nodes[i]->Write();
      h_ptmu[i]->Write();
      
      h_node_nodeenergy_long[i]->Write();
      h_nodes_long[i]->Write();
      h_node_nodeenergy_anchored[i]->Write();
      h_nodes_anchored[i]->Write();
      h_node_nodeenergy_vest[i]->Write();
      h_nodes_vest[i]->Write();
      
      h_protonRes_nodeenergy_node01[i]->Write();
      h_protonRes_nodeenergy_node2[i]->Write();
      h_protonRes_nodeenergy_node3[i]->Write();
      h_protonRes_nodeenergy_node4[i]->Write();
      h_protonRes_nodeenergy_node5[i]->Write();
      
      
      h_dalphat_ptmu[i]->Write();
      h_dphit_ptmu[i]->Write();
      h_pn_ptmu[i]->Write();
      h_dpt_ptmu[i]->Write();
      h_dptx_ptmu[i]->Write();
      h_dpty_ptmu[i]->Write();
      h_tp_ptmu[i]->Write();
      h_ptheta_ptmu[i]->Write();
      h_sign_ptmu[i]->Write();
      h_signdalphat_ptmu[i]->Write();
      h_signdphit_ptmu[i]->Write();

      h_dthetaR_ptmu[i]->Write();
      h_dthetaP_ptmu[i]->Write();
      h_dPP_ptmu[i]->Write();
      h_dPR_ptmu[i]->Write();
      h_dPPi_ptmu[i]->Write();
      h_dPRi_ptmu[i]->Write();
    }
    else{
      h_ptmu_tkitki[i]->Write();
    }
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
      "\t-./MuonSelectionHists  Name_and_path_to_Output_file  Selection_Signal  Make_Flux_Constraint_Histo  Multiplicity  Number_MC_files  Number_DATA_files \n\n"<<
      "\t-Name_and_path_to_Output_file\t =\t Name of and path to the Output ROOT file that will be created \n"<<
      "\t-Playlist\t =\t Name of the playlist you want to run over e.g. minerva1 \n"<<
      "\t-Selection_Signal\t =\t Can be: Signal, MichelSideBand, SideBand (this is:non-vtx vs Q2 sideband) \n"<<
      "\t-Make_Flux_Constraint_Histo\t =\t If TRUE Enter 1; If FALSE Enter 0 \n"<<
      "\t-Multiplicity\t =\t Enter 0, 1 or 2. Here 0: 1 or 2 tracks; 1: 1-track events only; 2: 2-track events only \n"<<
      "\t-Number_MC_files\t =\t Number of MonteCarlo files. To use all files, set this to: -1 \n"<<
      "\t-Number_DATA_files\t =\t Number of Data files. To use all files, set this to: -1\n" <<
      "\t-Variable combo\t = \tChoose TKI vs TKI combinations, default is off (-1)" << std::endl;
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
  par.push_back("-1");//no special variable combinations
  par.push_back("-1");

  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }

  bool fluxConstraintHistoNeeded = ( par.at(4) == "1" ) ? true : false; 

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;

  return MuonSelectionHists(par[1], par[3], fluxConstraintHistoNeeded, atoi(par[5].c_str()), par[2],atoi(par[6].c_str()), atoi(par[7].c_str()), atoi(par[8].c_str()), atoi(par[9].c_str()) );
}
