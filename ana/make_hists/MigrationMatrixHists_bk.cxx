#include "include/CCQENuUtils.h"
#include "PlotUtils/HyperDimLinearizer.h"
#include "TParameter.h"
#include "include/GeneralFunc.h" //in CCQENuNeutron/include/
#include "include/CommonBins.h" //in CCQENuNeutron/include/
#include "include/EventCuts.h"

using namespace CCQENU_ANA;

int MigrationMatrixHists( string filename, string playlist, bool makeFluxConstraintHisto, bool applyFluxConstraint, int multiplicity = 0, string sample = "pzmu", int n_mcfiles = -1 , int xtype = -1, int ytype = -1, int npMode = 1, string xvarname = "PT") 
{
  //---------------------------------------------
  // create file to store histograms
  //---------------------------------------------
  bool typeCheck = xtype==-1 && ytype==-1;//Do normal
  CCQENuUtils  *utils  = new CCQENuUtils( false, makeFluxConstraintHisto );
  utils->setPlaylist(playlist);
  utils->setFluxReweighterPlaylist();
  CCQENuCuts   *cutter = new CCQENuCuts();
  AnaBinning *binner = new AnaBinning();

  PhysicsUtils* physicsUtils = new PhysicsUtils();

  GeoUtils *geoUtils = new GeoUtils();

  //--------------------------------------------
  // Xvars
  //--------------------------------------------
  axis_binning xbins = muonPtbins;
  if (xvarname == "PT") xbins = muonPtbins;
  if (xvarname == "Q2") xbins = Q2bins;
  string xvartitle = "ptmu";
  if (xvarname == "Q2") xvartitle = "q2qe";

  //---------------------------------------------
  // Book Histograms
  //---------------------------------------------
  //
  // 2-D Smearing histograms for 1-D unfolding
  // This is only used to check the smearing in Reco vs Truth from 1-D variables.


  //2-D Smearing histograms for 2-D unfolding
  MinervaUnfold::MnvResponse *response_dalphat = NULL;
  MnvH2D *h_dalphat_xvar_migration = NULL;
  MnvH2D *h_dalphat_xvar = NULL;
  MnvH2D *h_dalphat_xvar_generated = NULL;
  //Dummy object: Just not to have data as a 0x0 object
  // utils->bookDataHisto( h_dalphat_xvar_migration, "h_dalphat_xvar_migration", "Dummy Migration", DALPHATbins, xbins );
  // utils->bookDataHisto( h_dalphat_xvar, "h_dalphat_xvar", "Dummy Migration", DALPHATbins, xbins );
  // utils->bookDataHisto( h_dalphat_xvar_generated, "h_dalphat_xvar_generated", "Dummy Migration", DALPHATbins, xbins );
  if(sample=="dalphat") utils->setupResponse(response_dalphat, "h_dalphat_"+xvartitle, "Delta AlphaT vs "+xvarname,dalphatbins,xbins,dalphatbins,xbins );

  MinervaUnfold::MnvResponse *response_dphit = NULL;
  MnvH2D *h_dphit_xvar_migration = NULL;
  MnvH2D *h_dphit_xvar = NULL;
  MnvH2D *h_dphit_xvar_generated = NULL;
  //Dummy object: Just not to have data as a 0x0 object
  // utils->bookDataHisto( h_dphit_xvar_migration, "h_dphit_xvar_migration", "Dummy Migration", DPHITbins, xbins );
  // utils->bookDataHisto( h_dphit_xvar, "h_dphit_xvar", "Dummy Migration", DPHITbins, xbins );
  // utils->bookDataHisto( h_dphit_xvar_generated, "h_dphit_xvar_generated", "Dummy Migration", DPHITbins, xbins );
  if(sample=="dphit") utils->setupResponse(response_dphit, "h_dphit_"+xvartitle, "Delta PhiT vs "+xvarname,dphitbins,xbins,dphitbins,xbins );

  MinervaUnfold::MnvResponse *response_pn = NULL;
  MnvH2D *h_pn_xvar_migration = NULL;
  MnvH2D *h_pn_xvar = NULL;
  MnvH2D *h_pn_xvar_generated = NULL;
  //Dummy object: Just not to have data as a 0x0 object
  // utils->bookDataHisto( h_pn_xvar_migration, "h_pn_xvar_migration", "Dummy Migration", PNbins, xbins );
  // utils->bookDataHisto( h_pn_xvar, "h_pn_xvar", "Dummy Migration", PNbins, xbins );
  // utils->bookDataHisto( h_pn_xvar_generated, "h_pn_xvar_generated", "Dummy Migration", PNbins, xbins );
  if(sample=="pn") utils->setupResponse(response_pn, "h_pn_"+xvartitle, "pn vs "+xvarname,pnbins,xbins,pnbins,xbins);


  MinervaUnfold::MnvResponse *response_dpt = NULL;
  MnvH2D *h_dpt_xvar_migration = NULL;
  MnvH2D *h_dpt_xvar = NULL;
  MnvH2D *h_dpt_xvar_generated = NULL;
  //Dummy object: Just not to have data as a 0x0 object
  // utils->bookDataHisto( h_dpt_xvar_migration, "h_dpt_xvar_migration", "Dummy Migration", DPTbins, xbins );
  // utils->bookDataHisto( h_dpt_xvar, "h_dpt_xvar", "Dummy Migration", DPTbins, xbins );
  // utils->bookDataHisto( h_dpt_xvar_generated, "h_dpt_xvar_generated", "Dummy Migration", DPTbins, xbins );
  if(sample=="dpt") utils->setupResponse(response_dpt, "h_dpt_"+xvartitle, "dpt vs "+xvarname,dptbins,xbins,dptbins,xbins);


  MinervaUnfold::MnvResponse *response_dptx = NULL;
  MnvH2D *h_dptx_xvar_migration = NULL;
  MnvH2D *h_dptx_xvar = NULL;
  MnvH2D *h_dptx_xvar_generated = NULL;
  //Dummy object: Just not to have data as a 0x0 object
  // utils->bookDataHisto( h_dptx_xvar_migration, "h_dptx_xvar_migration", "Dummy Migration", DPTXbins, xbins );
  // utils->bookDataHisto( h_dptx_xvar, "h_dptx_xvar", "Dummy Migration", DPTXbins, xbins );
  // utils->bookDataHisto( h_dptx_xvar_generated, "h_dptx_xvar_generated", "Dummy Migration", DPTXbins, xbins );
  if(sample=="dptx") utils->setupResponse(response_dptx, "h_dptx_"+xvartitle, "dptx vs "+xvarname,dptxbins,xbins,dptxbins,xbins);


  MinervaUnfold::MnvResponse *response_dpty = NULL;
  MnvH2D *h_dpty_xvar_migration = NULL;
  MnvH2D *h_dpty_xvar = NULL;
  MnvH2D *h_dpty_xvar_generated = NULL;
  //Dummy object: Just not to have data as a 0x0 object
  // utils->bookDataHisto( h_dpty_xvar_migration, "h_dpty_xvar_migration", "Dummy Migration", DPTYbins, xbins );
  // utils->bookDataHisto( h_dpty_xvar, "h_dpty_xvar", "Dummy Migration", DPTYbins, xbins );
  // utils->bookDataHisto( h_dpty_xvar_generated, "h_dpty_xvar_generated", "Dummy Migration", DPTYbins, xbins );
  if(sample=="dpty") utils->setupResponse(response_dpty, "h_dpty_"+xvartitle, "dpty vs "+xvarname,dptybins,xbins,dptybins,xbins);


  MinervaUnfold::MnvResponse *response_tp = NULL;
  MnvH2D *h_tp_xvar_migration = NULL;
  MnvH2D *h_tp_xvar = NULL;
  MnvH2D *h_tp_xvar_generated = NULL;
  //Dummy object: Just not to have data as a 0x0 object
  // utils->bookDataHisto( h_tp_xvar_migration, "h_tp_xvar_migration", "Dummy Migration", TPbins, xbins );
  // utils->bookDataHisto( h_tp_xvar, "h_tp_xvar", "Dummy Migration", TPbins, xbins );
  // utils->bookDataHisto( h_tp_xvar_generated, "h_tp_xvar_generated", "Dummy Migration", TPbins, xbins );
  if(sample=="tp") utils->setupResponse(response_tp, "h_tp_"+xvartitle, "tp vs "+xvarname,protonKineticbins,xbins,protonKineticbins,xbins);


  MinervaUnfold::MnvResponse *response_ptheta = NULL;
  MnvH2D *h_ptheta_xvar_migration = NULL;
  MnvH2D *h_ptheta_xvar = NULL;
  MnvH2D *h_ptheta_xvar_generated = NULL;
  //Dummy object: Just not to have data as a 0x0 object
  // utils->bookDataHisto( h_ptheta_xvar_migration, "h_ptheta_xvar_migration", "Dummy Migration", PTHETAbins, xbins );
  // utils->bookDataHisto( h_ptheta_xvar, "h_ptheta_xvar", "Dummy Migration", PTHETAbins, xbins );
  // utils->bookDataHisto( h_ptheta_xvar_generated, "h_ptheta_xvar_generated", "Dummy Migration", PTHETAbins, xbins );
  if(sample=="ptheta") utils->setupResponse(response_ptheta, "h_ptheta_"+xvartitle, "proton theta vs "+xvarname,protonThetabins,xbins,protonThetabins,xbins);


  MinervaUnfold::MnvResponse *response_signed = NULL;
  MnvH2D *h_signed_xvar_migration = NULL;
  MnvH2D *h_signed_xvar = NULL;
  MnvH2D *h_signed_xvar_generated = NULL;
  //Dummy object: Just not to have data as a 0x0 object
  // utils->bookDataHisto( h_signed_xvar_migration, "h_signed_xvar_migration", "Dummy Migration", Signbins, xbins );
  // utils->bookDataHisto( h_signed_xvar, "h_signed_xvar", "Dummy Migration", Signbins, xbins );
  // utils->bookDataHisto( h_signed_xvar_generated, "h_signed_xvar_generated", "Dummy Migration", Signbins, xbins );
  if(sample=="signed") utils->setupResponse(response_signed, "h_signed_"+xvartitle, "proton theta vs "+xvarname,signbins,xbins,signbins,xbins);


  MinervaUnfold::MnvResponse *response_signeddalphat = NULL;
  MnvH2D *h_signeddalphat_xvar_migration = NULL;
  MnvH2D *h_signeddalphat_xvar = NULL;
  MnvH2D *h_signeddalphat_xvar_generated = NULL;
  //Dummy object: Just not to have data as a 0x0 object
  // utils->bookDataHisto( h_signeddalphat_xvar_migration, "h_signeddalphat_xvar_migration", "Dummy Migration", SIGNEDDALPHATbins, xbins );
  // utils->bookDataHisto( h_signeddalphat_xvar, "h_signeddalphat_xvar", "Dummy Migration", SIGNEDDALPHATbins, xbins );
  // utils->bookDataHisto( h_signeddalphat_xvar_generated, "h_signeddalphat_xvar_generated", "Dummy Migration", SIGNEDDALPHATbins, xbins );
  if(sample=="signeddalphat") utils->setupResponse(response_signeddalphat, "h_signeddalphat_"+xvartitle, "proton theta vs "+xvarname,signeddalphatbins,xbins,signeddalphatbins,xbins);


  MinervaUnfold::MnvResponse *response_signeddphit = NULL;
  MnvH2D *h_signeddphit_xvar_migration = NULL;
  MnvH2D *h_signeddphit_xvar = NULL;
  MnvH2D *h_signeddphit_xvar_generated = NULL;
  //Dummy object: Just not to have data as a 0x0 object
  // utils->bookDataHisto( h_signeddphit_xvar_migration, "h_signeddphit_xvar_migration", "Dummy Migration", SIGNEDDPHITbins, xbins );
  // utils->bookDataHisto( h_signeddphit_xvar, "h_signeddphit_xvar", "Dummy Migration", SIGNEDDPHITbins, xbins );
  // utils->bookDataHisto( h_signeddphit_xvar_generated, "h_signeddphit_xvar_generated", "Dummy Migration", SIGNEDDPHITbins, xbins );
  if(sample=="signeddphit") utils->setupResponse(response_signeddphit, "h_signeddphit_"+xvartitle, "proton theta vs "+xvarname,signeddphitbins,xbins,signeddphitbins,xbins);

  //Tejin's thetaP/R variables
  MinervaUnfold::MnvResponse *response_dthetaR = NULL;
  MnvH2D *h_dthetaR_xvar_migration = NULL;
  MnvH2D *h_dthetaR_xvar = NULL;
  MnvH2D *h_dthetaR_xvar_generated = NULL;
  if(sample=="dthetaR") utils->setupResponse(response_dthetaR, "h_dthetaR_"+xvartitle, "proton thetaR vs "+xvarname,dthetaReactbins,xbins,dthetaReactbins,xbins);

  MinervaUnfold::MnvResponse *response_dthetaP = NULL;
  MnvH2D *h_dthetaP_xvar_migration = NULL;
  MnvH2D *h_dthetaP_xvar = NULL;
  MnvH2D *h_dthetaP_xvar_generated = NULL;
  if(sample=="dthetaP") utils->setupResponse(response_dthetaP, "h_dthetaP_"+xvartitle, "proton thetaP vs "+xvarname,dthetaPerpbins,xbins,dthetaPerpbins,xbins);

  MinervaUnfold::MnvResponse *response_dPR = NULL;
  MnvH2D *h_dPR_xvar_migration = NULL;
  MnvH2D *h_dPR_xvar = NULL;
  MnvH2D *h_dPR_xvar_generated = NULL;
  if(sample=="dPR") utils->setupResponse(response_dPR, "h_dPR_"+xvartitle, "proton PR vs "+xvarname,dptxbins,xbins,dptxbins,xbins);

  MinervaUnfold::MnvResponse *response_dPP = NULL;
  MnvH2D *h_dPP_xvar_migration = NULL;
  MnvH2D *h_dPP_xvar = NULL;
  MnvH2D *h_dPP_xvar_generated = NULL;
  if(sample=="dPP") utils->setupResponse(response_dPP, "h_dPP_"+xvartitle, "proton dPP vs "+xvarname,dptxbins,xbins,dptxbins,xbins);


  MinervaUnfold::MnvResponse *response_dPRi = NULL;
  MnvH2D *h_dPRi_xvar_migration = NULL;
  MnvH2D *h_dPRi_xvar = NULL;
  MnvH2D *h_dPRi_xvar_generated = NULL;
  if(sample=="dPRi") utils->setupResponse(response_dPRi, "h_dPRi_"+xvartitle, "proton PR vs "+xvarname,dptxbins,xbins,dptxbins,xbins);

  MinervaUnfold::MnvResponse *response_dPPi = NULL;
  MnvH2D *h_dPPi_xvar_migration = NULL;
  MnvH2D *h_dPPi_xvar = NULL;
  MnvH2D *h_dPPi_xvar_generated = NULL;
  if(sample=="dPPi") utils->setupResponse(response_dPPi, "h_dPPi_"+xvartitle, "proton dPPi vs "+xvarname,dptxbins,xbins,dptxbins,xbins);




  //special case
 //Special 3D for TKI vs TKI
  string spec_title = "";
  vector<vector<double> > full3D;
  //Now to pick the variables
  map<int,axis_binning > axis_bin_map;
  map<int,string > axis_name;
  axis_bin_map[0] = muonPzbins;
  axis_bin_map[1] = dalphatbins_TKI;
  axis_bin_map[2] = dphitbins_TKI;
  axis_bin_map[3] = pnbins_TKI;
  axis_bin_map[4] = dptbins_TKI;
  axis_bin_map[5] = dptxbins_TKI;
  axis_bin_map[6] = dptybins_TKI;
  axis_bin_map[7] = protonKineticbins_TKI;
  axis_bin_map[8] = protonThetabins_TKI;
  axis_bin_map[9] = signbins;
  axis_bin_map[10] = signeddalphatbins_TKI;
  axis_bin_map[11] = signeddphitbins_TKI;

  axis_name[0] = "muonPz";
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

  if(!typeCheck) spec_title = "h_"+axis_name[xtype]+"_"+axis_name[ytype];
  if(sample=="tkitki") cout << "I will be running " << spec_title << endl;

  MinervaUnfold::MnvResponse *response_tkitki = NULL;
  MnvH2D *h_tkitki_xvar_migration = NULL;
  MnvH2D *h_tkitki_xvar = NULL;
  MnvH2D *h_tkitki_xvar_generated = NULL;
  //Dummy object: Just not to have data as a 0x0 object
  // utils->bookDataHisto( h_tkitki_xvar_migration, "h_tkitki_xvar_migration", "Dummy Migration", TKITKIbins, xbins );
  // utils->bookDataHisto( h_tkitki_xvar, "h_tkitki_xvar", "Dummy Migration", TKITKIbins, xbins );
  // utils->bookDataHisto( h_tkitki_xvar_generated, "h_tkitki_xvar_generated", "Dummy Migration", TKITKIbins, xbins );
  if(sample=="tkitki") utils->setupResponse(response_tkitki, spec_title.c_str(), "tkitki",axis_bin_map[xtype],axis_bin_map[ytype],axis_bin_map[xtype],axis_bin_map[ytype]);








  //  cout << response_q2 << endl;
  //---------------------------------------------
  // Add Error Bands
  //---------------------------------------------
  if(RunCodeWithSystematics){
    if(sample=="dalphat"){
    //Response dalphat/pt
    utils->addVertErrorBands( response_dalphat );
    utils->addLatErrorBands(  response_dalphat );
    }

    if(sample=="dphit"){
    //Response dphit/pt
    utils->addVertErrorBands( response_dphit );
    utils->addLatErrorBands(  response_dphit );
    }

    if(sample=="pn"){
    //Response pn/pt
    utils->addVertErrorBands( response_pn );
    utils->addLatErrorBands(  response_pn );
    }


    if(sample=="dpt"){
    //Response dpt/pt
    utils->addVertErrorBands( response_dpt );
    utils->addLatErrorBands(  response_dpt );
    }

    if(sample=="dptx"){
    //Response dptx/pt
    utils->addVertErrorBands( response_dptx );
    utils->addLatErrorBands(  response_dptx );
    }

    if(sample=="dpty"){
    //Response dpty/pt
    utils->addVertErrorBands( response_dpty );
    utils->addLatErrorBands(  response_dpty );
    }

    if(sample=="tp"){
    //Response tp/pt
    utils->addVertErrorBands( response_tp );
    utils->addLatErrorBands(  response_tp );
    }

    if(sample=="ptheta"){
    //Response ptheta/pt
    utils->addVertErrorBands( response_ptheta );
    utils->addLatErrorBands(  response_ptheta );
    }


    if(sample=="signed"){
    //Response signed/pt
    utils->addVertErrorBands( response_signed );
    utils->addLatErrorBands(  response_signed );
    }


    if(sample=="signeddalphat"){
    //Response signeddalphat/pt
    utils->addVertErrorBands( response_signeddalphat );
    utils->addLatErrorBands(  response_signeddalphat );
    }


    if(sample=="signeddphit"){
    //Response signeddphit/pt
    utils->addVertErrorBands( response_signeddphit );
    utils->addLatErrorBands(  response_signeddphit );
    }

    if(sample=="dthetaR"){
    //Response signeddphit/pt
    utils->addVertErrorBands( response_dthetaR );
    utils->addLatErrorBands(  response_dthetaR );
    }

    if(sample=="dthetaP"){
    //Response signeddphit/pt
    utils->addVertErrorBands( response_dthetaP );
    utils->addLatErrorBands(  response_dthetaP );
    }

    if(sample=="dPR"){
    //Response signeddphit/pt
    utils->addVertErrorBands( response_dPR );
    utils->addLatErrorBands(  response_dPR );
    }

    if(sample=="dPP"){
    //Response signeddphit/pt
    utils->addVertErrorBands( response_dPP );
    utils->addLatErrorBands(  response_dPP );
    }

    if(sample=="dPRi"){
    //Response signeddphit/pt
    utils->addVertErrorBands( response_dPRi );
    utils->addLatErrorBands(  response_dPRi );
    }

    if(sample=="dPPi"){
    //Response signeddphit/pt
    utils->addVertErrorBands( response_dPPi );
    utils->addLatErrorBands(  response_dPPi );
    }




    if(sample=="tkitki"){
    //Response tkitki/pt
    utils->addVertErrorBands( response_tkitki );
    utils->addLatErrorBands(  response_tkitki );
    }

  }//end run with sysetmatics



  TChain *truth_mc = utils->getMCTree("Truth", n_mcfiles );

  //---------------------------------------------
  // Get MC Tree
  //---------------------------------------------
  TChain* tree_mc  = utils->getMCTree("CCQENu", n_mcfiles );
  int entries_mc   = tree_mc->GetEntries();
  utils->setmnvHadronReweightTruthTree(truth_mc);
  utils->setmnvHadronReweightDataTree(tree_mc);
  cout << "MC entries: " << entries_mc << endl;
  CCQENuEvent* mc = new CCQENuEvent( tree_mc );

  //---------------------------------------------
  // Fill MC Histograms
  //---------------------------------------------
  double mc_events = 0., data_events = 0.;

  for (int i=0; i<entries_mc ; ++i) {
    if( (i+1)%100000 == 0 ) 
      cout << "Reading MC entry : " << i+1 << " - " << 100*(i+1)/entries_mc << " % " << endl; 
    
    mc->GetEntry(i);

    if(!cutter->passTrueCCQELike(mc)) continue;

    if (multiplicity==0){ //1+2tracks
      if( !cutter->passInteractionVertex( mc ) )     continue;	
      if( !cutter->passCCQESelection(mc) )           continue;
    }
    else if( multiplicity==1){ //1 track only 
      if( !cutter->passInteractionVertex( mc ) )     continue;	
      if( !cutter->passDeadTimeCuts( mc ) )          continue;
      if( !cutter->passNuHelicityCut( mc ) )         continue;
      if( !(mc->multiplicity==1) )                   continue;
      if( !cutter->passImprovedMichelCut( mc ) )     continue;
      if( !cutter->passExtraTracksProtonsCut ( mc ) )continue;
      if( !cutter->passNBlobs( mc ) )                continue;
      if( !cutter->passTrackAngleCut( mc ) )         continue;
    }
    else if( multiplicity==2){ //2 tracks only 
      if( !cutter->passInteractionVertex( mc ) )     continue;	
      if( !cutter->passDeadTimeCuts( mc ) )          continue;
      if( !cutter->passNuHelicityCut( mc ) )         continue;
      if( !(mc->multiplicity==2) )                   continue;
      if( !cutter->passImprovedMichelCut( mc ) )     continue;
      if( !cutter->passExtraTracksProtonsCut ( mc ) )continue;
      if( !cutter->passNBlobs( mc ) )                continue;
      if( !cutter->passTrackAngleCut( mc ) )         continue;
    }

    //--------
    // THis is the old method. Beacause of angular bias introduced from hadronic energy we need to pick out the right angle from the downstream node
    // New method until general reco is fixed is to get the reco P, a downstream theta (already corrected to the beam direction), and calculate
    // the pt and pz from these two quantities
    //-------


    //Muon Variables
    double reco_muon_theta     = mc->muon_theta_allNodes[20]; // get theta from 20th node, note NX kludge fixes this...
    double reco_cos_muon_theta = cos(reco_muon_theta);
    double reco_muon_T         = mc->CCQENu_muon_T;
    double reco_muon_px   = mc->CCQENu_leptonE[0];
    double reco_muon_py        = mc->CCQENu_leptonE[1];
    double reco_muon_pz        = mc->CCQENu_leptonE[2];
    double muon_p              = utils->GetTotalMomentum( reco_muon_px, reco_muon_py, reco_muon_pz );
    double reco_muon_pt_beam = sin(reco_muon_theta)*muon_p;
    double reco_muon_pz_beam = cos(reco_muon_theta)*muon_p;

    //Proton Variables
    TVector3 protonVect(mc->CCQENu_proton_Px_fromdEdx,mc->CCQENu_proton_Py_fromdEdx,mc->CCQENu_proton_Pz_fromdEdx);//MeV
    double reco_protonAngle = mc->CCQENu_proton_theta*180./3.14159;
    double reco_protonMom = protonVect.Mag();//MeV
    double reco_protonTn = TMath::Sqrt(reco_protonMom*reco_protonMom+utils->M_p*utils->M_p)-utils->M_p;//Done in MeV


    


    //Truth Proton Variables
    
    //Convert to degrees, GeV
    reco_muon_theta*= 180. / 3.14159;
    reco_muon_px /= pow(10,3);
    reco_muon_py /= pow(10,3);
    reco_muon_pz /= pow(10,3);
    protonVect *= pow(10,-3);
    reco_protonTn /= pow(10,3);
    reco_muon_T /= pow(10,3);
    reco_muon_pt_beam /= pow(10,3);
    reco_muon_pz_beam /= pow(10,3);

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
    //
    double binding_e = 0; // theoretically we should include binding E for carbon.. but I'm comparing to hydrogen in anti-nu, so probably set it to 0.
    double is_part_mass = (neutrinoMode)? 0.9395654133: 0.93827231; //neutron/proton
    double fs_part_mass = (neutrinoMode)? 0.93827231: 0.9395654133;
    int charge = (neutrinoMode)? -1 : 1;



    double beam_bias = 0;
    XYZVector beam = geoUtils->BeamAxis( beam_bias ); 
    XYZTVector reco_muon_4P( reco_muon_px, reco_muon_py, reco_muon_pz ,mc->CCQENu_muon_E/1e3);

    double reco_q2qe = physicsUtils->qsqCCQE( beam, reco_muon_4P, charge, binding_e ); //GeV^2
    
    XYZTVector expectedPRoton4P = geoUtils->ComputeExpectedNucleon( beam, reco_muon_4P, is_part_mass, fs_part_mass, binding_e);
    XYZVector expectedPRoton3P = (XYZVector) expectedPRoton4P.Vect();

    XYZVector protonVect_XYZ( protonVect.X(), protonVect.Y(), protonVect.Z() );
    vector<double> reco_vars = geoUtils->ComputeNeutronAngularVars( beam, expectedPRoton3P, protonVect_XYZ );
    double reco_dthetaP = reco_vars[0]*180/TMath::Pi();
    double reco_dthetaR = reco_vars[1]*180/TMath::Pi();
    double reco_dtheta = reco_vars[2]*180/TMath::Pi();
    double reco_dPP = reco_vars[3];
    double reco_dPR = reco_vars[4];
    double reco_dPPi = reco_vars[5];
    double reco_dPRi = reco_vars[6];



    //Get True values
    const double bindingE = 34.0;

    TVector3 true_proton_p_comp = utils->GetHighestTrueProtonMomentumComp(mc);
    double true_proton_px = true_proton_p_comp.X();
    double true_proton_py = true_proton_p_comp.Y();
    double true_proton_pz = true_proton_p_comp.Z();

    double true_proton_P = TMath::Sqrt(true_proton_px*true_proton_px+
				       true_proton_py*true_proton_py+
				       true_proton_pz*true_proton_pz);

    double true_protonTn = TMath::Sqrt(true_proton_P*true_proton_P+(utils->M_p)*(utils->M_p))-utils->M_p;
    double true_protonAngle = TMath::ACos( TMath::Sqrt(true_proton_px*true_proton_px+true_proton_py*true_proton_py)/(true_proton_pz))*180.0/3.14159;

    double true_muon_px   = mc->mc_primFSLepton[0];
    double true_muon_py   = mc->mc_primFSLepton[1];
    double true_muon_pz   = mc->mc_primFSLepton[2];
    double true_muon_E    = mc->mc_primFSLepton[3];
    double true_muon_p    = utils->GetTotalMomentum( true_muon_px, true_muon_py, true_muon_pz ); //GeV
    double true_theta     = utils->GetTheta( true_muon_px, true_muon_py, true_muon_pz );
    double true_muon_pt_beam   = utils->GetTransverseMomentumWRTBeam( true_muon_px, true_muon_py, true_muon_pz ); //GeV 
    double true_muon_pz_beam   = utils->GetLongitudinalMomentumWRTBeam( true_muon_px, true_muon_py, true_muon_pz ); //GeV


    //Add true theta cut
    if(true_theta*180./3.14159 >20) continue;
    //Scales
    true_proton_px /= pow(10,3);
    true_proton_py /= pow(10,3);
    true_proton_pz /= pow(10,3);
    true_proton_P  /= pow(10,3);
    true_protonTn        /= pow(10,3);
    true_muon_px   /= pow(10,3);
    true_muon_py   /= pow(10,3);
    true_muon_pz   /= pow(10,3);
    true_muon_E   /= pow(10,3);
    true_muon_p    /= pow(10,3);
    true_muon_pt_beam /= pow(10,3);
    true_muon_pz_beam /= pow(10,3);
  
    double true_dalphat   = utils->GetDeltaAlphaT(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_dphit     = utils->GetDeltaPhiT(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_pn        = utils->GetPn(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_dpt       = utils->GetDeltaPt(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_dptx      = utils->GetDeltaPtx(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_dpty      = utils->GetDeltaPty(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_sign      = utils->GetMuonProtonSign(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_signed_dalphat = utils->GetSignedDeltaAlphaT(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_signed_dphit   = utils->GetSignedDeltaPhiT(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);

    //Tejin recalculate using truth

    XYZTVector true_muon_4P( true_muon_px, true_muon_py, true_muon_pz ,true_muon_E);
    double true_q2qe = physicsUtils->qsqCCQE( beam, true_muon_4P, charge, binding_e ); //GeV^2
    expectedPRoton4P = geoUtils->ComputeExpectedNucleon( beam, true_muon_4P, is_part_mass, fs_part_mass, binding_e);
    expectedPRoton3P = (XYZVector) expectedPRoton4P.Vect();
    protonVect_XYZ.SetXYZ( true_proton_px, true_proton_py, true_proton_pz );

    vector<double> true_vars = geoUtils->ComputeNeutronAngularVars( beam, expectedPRoton3P, protonVect_XYZ );
    double true_dthetaP   = true_vars[0]*180/TMath::Pi();
    double true_dthetaR   = true_vars[1]*180/TMath::Pi();
    double true_dtheta    = true_vars[2]*180/TMath::Pi();
    double true_dPP       = true_vars[3];
    double true_dPR       = true_vars[4];
    double true_dPPi = true_vars[5];
    double true_dPRi = true_vars[6];




    map<int,double> reco_varmap;
    reco_varmap[0] = reco_muon_pz_beam;//muonPzbins;
    reco_varmap[1] = reco_dalphat;//dalphatbins;
    reco_varmap[2] = reco_dphit;//dphitbins;
    reco_varmap[3] = reco_pn;//pnbins;
    reco_varmap[4] = reco_dpt;//dptbins;
    reco_varmap[5] = reco_dptx;//dptxbins;
    reco_varmap[6] = reco_dpty;//dptybins;
    reco_varmap[7] = reco_protonTn;//protonKineticbins;
    reco_varmap[8] = reco_protonAngle;//protonThetabins;
    reco_varmap[9] = reco_sign;//signbins;
    reco_varmap[10] = reco_signed_dalphat;//signeddalphatbins;
    reco_varmap[11] = reco_signed_dphit;//signeddphitbins;

    reco_varmap[12] = reco_dthetaP;  
    reco_varmap[13] = reco_dthetaR;  
    reco_varmap[14] = reco_dtheta;
    reco_varmap[15] = reco_dPP;      
    reco_varmap[16] = reco_dPR;     
    reco_varmap[17] = reco_dPPi;
    reco_varmap[18] = reco_dPRi;


 

    map<int,double> true_varmap;
    true_varmap[0] = true_muon_pz_beam;//muonPzbins;
    true_varmap[1] = true_dalphat;//dalphatbins;
    true_varmap[2] = true_dphit;//dphitbins;
    true_varmap[3] = true_pn;//pnbins;
    true_varmap[4] = true_dpt;//dptbins;
    true_varmap[5] = true_dptx;//dptxbins;
    true_varmap[6] = true_dpty;//dptybins;
    true_varmap[7] = true_protonTn;//protonKineticbins;
    true_varmap[8] = true_protonAngle;//protonThetabins;
    true_varmap[9] = true_sign;//signbins;
    true_varmap[10] = true_signed_dalphat;//signeddalphatbins;
    true_varmap[11] = true_signed_dphit;//signeddphitbins;

    true_varmap[12] = true_dthetaP;  
    true_varmap[13] = true_dthetaR;  
    true_varmap[14] = true_dtheta;
    true_varmap[15] = true_dPP;      
    true_varmap[16] = true_dPR;     
    true_varmap[17] = true_dPPi;
    true_varmap[18] = true_dPRi;

    //Set X var
    double reco_xvar = -99999;
    double true_xvar = -99999;

    if( xvarname == "PT" )
    {
      reco_xvar = reco_muon_pt_beam;
      true_xvar = true_muon_pt_beam;
    }
    if( xvarname == "Q2" )
    {
      reco_xvar = reco_q2qe;
      true_xvar = true_q2qe;
    }




    //Get cvweight
    double wgt = utils->GetCVWeight( mc );
    mc_events+=wgt;
    
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
      //      cout << mc->CCQENu_nuHelicity << endl;
      //---------------------------------------------
      // Fill Histograms
      //---------------------------------------------
      //Start with the 2-D migration histograms for 1-D unfolding

      //Then, fill response object for the 2-D migration histograms for 2-D unfolding
      if(sample=="dalphat")   utils->fillResponse( response_dalphat, reco_dalphat, reco_xvar, true_dalphat, true_xvar, mc, wgt );
      if(sample=="dphit")     utils->fillResponse( response_dphit, reco_dphit, reco_xvar, true_dphit, true_xvar, mc, wgt );
      if(sample=="pn")        utils->fillResponse( response_pn, reco_pn, reco_xvar, true_pn, true_xvar, mc, wgt );
      if(sample=="dpt")       utils->fillResponse( response_dpt, reco_dpt, reco_xvar, true_dpt, true_xvar, mc, wgt );
      if(sample=="dptx")      utils->fillResponse( response_dptx, reco_dptx, reco_xvar, true_dptx, true_xvar, mc, wgt );
      if(sample=="dpty")      utils->fillResponse( response_dpty, reco_dpty, reco_xvar, true_dpty, true_xvar, mc, wgt );
      if(sample=="tp")        utils->fillResponse( response_tp, reco_protonTn, reco_xvar, true_protonTn, true_xvar, mc, wgt );
      if(sample=="ptheta")    utils->fillResponse( response_ptheta, reco_protonAngle, reco_xvar, true_protonAngle, true_xvar, mc, wgt );
      if(sample=="signed")    utils->fillResponse( response_signed, reco_sign, reco_xvar, true_sign, true_xvar, mc, wgt );
      if(sample=="signeddphit")    utils->fillResponse( response_signeddphit, reco_signed_dphit, reco_xvar, true_signed_dphit, true_xvar, mc, wgt );
      if(sample=="signeddalphat")    utils->fillResponse( response_signeddalphat, reco_signed_dalphat, reco_xvar, true_signed_dalphat, true_xvar, mc, wgt );

      //Tejinc
      if(sample=="dthetaR")    utils->fillResponse( response_dthetaR, reco_dthetaR, reco_xvar, true_dthetaR, true_xvar, mc, wgt );
      if(sample=="dthetaP")    utils->fillResponse( response_dthetaP, reco_dthetaP, reco_xvar, true_dthetaP, true_xvar, mc, wgt );
      if(sample=="dPR")    utils->fillResponse( response_dPR, reco_dPR, reco_xvar, true_dPR, true_xvar, mc, wgt );
      if(sample=="dPP")    utils->fillResponse( response_dPP, reco_dPP, reco_xvar, true_dPP, true_xvar, mc, wgt );
      if(sample=="dPRi")    utils->fillResponse( response_dPRi, reco_dPRi, reco_xvar, true_dPRi, true_xvar, mc, wgt );
      if(sample=="dPPi")    utils->fillResponse( response_dPPi, reco_dPPi, reco_xvar, true_dPPi, true_xvar, mc, wgt );

      if(sample=="tkitki")    utils->fillResponse( response_tkitki, reco_varmap[xtype], reco_varmap[ytype], true_varmap[xtype], true_varmap[ytype], mc, wgt );

      //---------------------------
      //Fill Vertical Errors
      //---------------------------
      //Fill Response object (the one used for 2-D unfolding) only
      if(RunCodeWithSystematics){
	if(sample=="dalphat")  utils->fillVertErrorBands( response_dalphat, reco_dalphat, reco_xvar, true_dalphat, true_xvar,  mc );
	if(sample=="dphit")    utils->fillVertErrorBands( response_dphit, reco_dphit, reco_xvar, true_dphit, true_xvar,  mc );
	if(sample=="pn")       utils->fillVertErrorBands( response_pn, reco_pn, reco_xvar, true_pn, true_xvar,  mc );
	if(sample=="dpt")      utils->fillVertErrorBands( response_dpt, reco_dpt, reco_xvar, true_dpt, true_xvar,  mc );
	if(sample=="dptx")     utils->fillVertErrorBands( response_dptx, reco_dptx, reco_xvar, true_dptx, true_xvar,  mc );
	if(sample=="dpty")     utils->fillVertErrorBands( response_dpty, reco_dpty, reco_xvar, true_dpty, true_xvar,  mc );
	if(sample=="tp")       utils->fillVertErrorBands( response_tp, reco_protonTn, reco_xvar, true_protonTn, true_xvar,  mc );
	if(sample=="ptheta")   utils->fillVertErrorBands( response_ptheta, reco_protonAngle, reco_xvar, true_protonAngle, true_xvar,  mc );
	if(sample=="signed")    utils->fillVertErrorBands( response_signed, reco_sign, reco_xvar, true_sign, true_xvar, mc );
	if(sample=="signeddphit")    utils->fillVertErrorBands( response_signeddphit, reco_signed_dphit, reco_xvar, true_signed_dphit, true_xvar, mc );
	if(sample=="signeddalphat")    utils->fillVertErrorBands( response_signeddalphat, reco_signed_dalphat, reco_xvar, true_signed_dalphat, true_xvar, mc );

  if(sample=="dthetaR")    utils->fillVertErrorBands( response_dthetaR, reco_dthetaR, reco_xvar, true_dthetaR, true_xvar, mc );
  if(sample=="dthetaP")    utils->fillVertErrorBands( response_dthetaP, reco_dthetaP, reco_xvar, true_dthetaP, true_xvar, mc );
  if(sample=="dPR")    utils->fillVertErrorBands( response_dPR, reco_dPR, reco_xvar, true_dPR, true_xvar, mc );
  if(sample=="dPP")    utils->fillVertErrorBands( response_dPP, reco_dPP, reco_xvar, true_dPP, true_xvar, mc );
  if(sample=="dPRi")    utils->fillVertErrorBands( response_dPRi, reco_dPRi, reco_xvar, true_dPRi, true_xvar, mc );
  if(sample=="dPPi")    utils->fillVertErrorBands( response_dPPi, reco_dPPi, reco_xvar, true_dPPi, true_xvar, mc );



	if(sample=="tkitki")    utils->fillVertErrorBands( response_tkitki, reco_varmap[xtype], reco_varmap[ytype], true_varmap[xtype], true_varmap[ytype], mc);
      }//end run with systematics
    } // End of fillHistos/fillResponse and fillVertErrorBands where CV of events pass the proton and recoil cuts 

    //---------------------------------------------------------------------------------------
    //Now fillLatErrorBands for events whose CV did not pass the proton and recoil cuts  
    //What if the shifted proton score and recoil pass the selection cuts ? 
    //Those events should also be included in the LatErrorBands
    //fillLatErrorBands have their own checks for shifted proton scores and recoils  
    //---------------------------------------------------------------------------------------
    if(RunCodeWithSystematics){


      //Will Add in transverse


    }//end run with sysetmatics
  } //End of for loop over MC entries 
 
  cout << "Finished looping over all MC entries " << endl; 
  delete mc; 
  delete tree_mc; 

  //-----------------------------------------------
  // Getting h_migration, h_reco, h_truth
  //-----------------------------------------------
  if(sample=="dalphat"){
    utils->getResponseObjects(response_dalphat, h_dalphat_xvar_migration, h_dalphat_xvar, h_dalphat_xvar_generated);
    delete response_dalphat;
  }

  if(sample=="dphit"){
    utils->getResponseObjects(response_dphit, h_dphit_xvar_migration, h_dphit_xvar, h_dphit_xvar_generated);
    delete response_dphit;
  }

  if(sample=="pn"){
    utils->getResponseObjects(response_pn, h_pn_xvar_migration, h_pn_xvar, h_pn_xvar_generated);
    delete response_pn;
  }

  if(sample=="dpt"){
    utils->getResponseObjects(response_dpt, h_dpt_xvar_migration, h_dpt_xvar, h_dpt_xvar_generated);
    delete response_dpt;
  }

  if(sample=="dptx"){
    utils->getResponseObjects(response_dptx, h_dptx_xvar_migration, h_dptx_xvar, h_dptx_xvar_generated);
    delete response_dptx;
  }

  if(sample=="dpty"){
    utils->getResponseObjects(response_dpty, h_dpty_xvar_migration, h_dpty_xvar, h_dpty_xvar_generated);
    delete response_dpty;
  }

  if(sample=="tp"){
    utils->getResponseObjects(response_tp, h_tp_xvar_migration, h_tp_xvar, h_tp_xvar_generated);
    delete response_tp;
  }

  if(sample=="ptheta"){
    utils->getResponseObjects(response_ptheta, h_ptheta_xvar_migration, h_ptheta_xvar, h_ptheta_xvar_generated);
    delete response_ptheta;
  }

  if(sample=="signed"){
    utils->getResponseObjects(response_signed, h_signed_xvar_migration, h_signed_xvar, h_signed_xvar_generated);
    delete response_signed;
  }

  if(sample=="signeddalphat"){
    utils->getResponseObjects(response_signeddalphat, h_signeddalphat_xvar_migration, h_signeddalphat_xvar, h_signeddalphat_xvar_generated);
    delete response_signeddalphat;
  }

  if(sample=="signeddphit"){
    utils->getResponseObjects(response_signeddphit, h_signeddphit_xvar_migration, h_signeddphit_xvar, h_signeddphit_xvar_generated);
    delete response_signeddphit;
  }

  if(sample=="dthetaR"){
    utils->getResponseObjects(response_dthetaR, h_dthetaR_xvar_migration, h_dthetaR_xvar, h_dthetaR_xvar_generated);
    delete response_dthetaR;
  }
  if(sample=="dthetaP"){
    utils->getResponseObjects(response_dthetaP, h_dthetaP_xvar_migration, h_dthetaP_xvar, h_dthetaP_xvar_generated);
    delete response_dthetaP;
  }
  if(sample=="dPR"){
    utils->getResponseObjects(response_dPR, h_dPR_xvar_migration, h_dPR_xvar, h_dPR_xvar_generated);
    delete response_dPR;
  }
  if(sample=="dPP"){
    utils->getResponseObjects(response_dPP, h_dPP_xvar_migration, h_dPP_xvar, h_dPP_xvar_generated);
    delete response_dPP;
  }
  if(sample=="dPRi"){
    utils->getResponseObjects(response_dPRi, h_dPRi_xvar_migration, h_dPRi_xvar, h_dPRi_xvar_generated);
    delete response_dPRi;
  }
  if(sample=="dPPi"){
    utils->getResponseObjects(response_dPPi, h_dPPi_xvar_migration, h_dPPi_xvar, h_dPPi_xvar_generated);
    delete response_dPPi;
  }


  if(sample=="tkitki"){
    utils->getResponseObjects(response_tkitki, h_tkitki_xvar_migration, h_tkitki_xvar, h_tkitki_xvar_generated);
    delete response_tkitki;
  }

  
  MnvH2D *h_X_xvar_migration[200] = {NULL};

  //Need to split response to lat/vert to beable to save to disk...
  if(sample=="dalphat"){
    utils->splitObject(h_dalphat_xvar_migration, h_X_xvar_migration);
    delete h_dalphat_xvar_migration;
  }
  if(sample=="dphit"){
    utils->splitObject(h_dphit_xvar_migration, h_X_xvar_migration);
    delete h_dphit_xvar_migration;
  }
  if(sample=="pn"){
    utils->splitObject(h_pn_xvar_migration, h_X_xvar_migration);
    delete h_pn_xvar_migration;
  }

  if(sample=="dpt"){
    utils->splitObject(h_dpt_xvar_migration, h_X_xvar_migration);
    delete h_dpt_xvar_migration;
  }

  if(sample=="dptx"){
    utils->splitObject(h_dptx_xvar_migration, h_X_xvar_migration);
    delete h_dptx_xvar_migration;
  }

  if(sample=="dpty"){
    utils->splitObject(h_dpty_xvar_migration, h_X_xvar_migration);
    delete h_dpty_xvar_migration;
  }
 
  if(sample=="tp"){
    utils->splitObject(h_tp_xvar_migration, h_X_xvar_migration);
    delete h_tp_xvar_migration;
  }

  if(sample=="ptheta"){
    utils->splitObject(h_ptheta_xvar_migration, h_X_xvar_migration);
    delete h_ptheta_xvar_migration;
  }

  if(sample=="signed"){
    utils->splitObject(h_signed_xvar_migration, h_X_xvar_migration);
    delete h_signed_xvar_migration;
  }

  if(sample=="signeddalphat"){
    utils->splitObject(h_signeddalphat_xvar_migration, h_X_xvar_migration);
    delete h_signeddalphat_xvar_migration;
  }

  if(sample=="signeddphit"){
    utils->splitObject(h_signeddphit_xvar_migration, h_X_xvar_migration);
    delete h_signeddphit_xvar_migration;
  }
  //Tejin
  if(sample=="dthetaR"){
    utils->splitObject(h_dthetaR_xvar_migration, h_X_xvar_migration);
    delete h_dthetaR_xvar_migration;
  }
  if(sample=="dthetaP"){
    utils->splitObject(h_dthetaP_xvar_migration, h_X_xvar_migration);
    delete h_dthetaP_xvar_migration;
  }
  if(sample=="dPR"){
    utils->splitObject(h_dPR_xvar_migration, h_X_xvar_migration);
    delete h_dPR_xvar_migration;
  }
  if(sample=="dPP"){
    utils->splitObject(h_dPP_xvar_migration, h_X_xvar_migration);
    delete h_dPP_xvar_migration;
  }
  if(sample=="dPRi"){
    utils->splitObject(h_dPRi_xvar_migration, h_X_xvar_migration);
    delete h_dPRi_xvar_migration;
  }
  if(sample=="dPPi"){
    utils->splitObject(h_dPPi_xvar_migration, h_X_xvar_migration);
    delete h_dPPi_xvar_migration;
  }



  if(sample=="tkitki"){
    utils->splitObject(h_tkitki_xvar_migration, h_X_xvar_migration);
    delete h_tkitki_xvar_migration;
  }

  
  //Get Data POT for this playlist
  //Get data POT for particular playlist
  TChain* tree_data  = utils->getDataTree( "CCQENu", -1 );
  delete tree_data;


  //==================================================================
  // Create ROOT file to store histograms
  // Write to file all the created histograms 
  //==================================================================
  if(xtype!=-1 && ytype!=-1) filename = filename+"_"+axis_name[xtype]+"_"+axis_name[ytype]+".root";
  TFile *f = new TFile( filename.c_str(), "RECREATE" );  
  cout << "write POT" << endl;
  //Write MC and DATA POT to file 
  utils->writePOT( f );
  cout << "write Multiplicity" << endl;
  //Write multiplicity and sample definitions
  TMap *sampleMap = new TMap(2,0);
  string multiplicity_label = (multiplicity==0)? "alltracks": Form("%i_track",multiplicity);
  sampleMap->Add(new TObjString("multiplicity"), new TObjString(multiplicity_label.c_str()) );
  f->WriteTObject( sampleMap, "sample");
  cout << "write n_events" << endl;
  //Write DATA and MC number of events
  TVector2 *evt = new TVector2( data_events, mc_events );
  f->WriteTObject( evt, "n_events" );
    cout << "write fluxParam" << endl;
  //Write out if output ROOT file has the constrained flux or not 
  TParameter<bool> *fluxParam = new TParameter<bool>( "appliedFluxConstraint", applyFluxConstraint ); 
  f->WriteTObject( fluxParam ); 

  f->cd();
  if(sample=="dalphat"){
    cout << "write hists dalphatpt" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_dalphat_xvar->Write();
    h_dalphat_xvar_generated->Write();
  }

  if(sample=="dphit"){
    cout << "write hists dphitpt" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_dphit_xvar->Write();
    h_dphit_xvar_generated->Write();
  }

  if(sample=="pn"){
    cout << "write hists pnpt" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_pn_xvar->Write();
    h_pn_xvar_generated->Write();
  }

  if(sample=="dpt"){
    cout << "write hists dptpt" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_dpt_xvar->Write();
    h_dpt_xvar_generated->Write();
  }

  if(sample=="dptx"){
    cout << "write hists dptxpt" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_dptx_xvar->Write();
    h_dptx_xvar_generated->Write();
  }

  if(sample=="dpty"){
    cout << "write hists dptypt" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_dpty_xvar->Write();
    h_dpty_xvar_generated->Write();
  }

  if(sample=="tp"){
    cout << "write hists tppt" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_tp_xvar->Write();
    h_tp_xvar_generated->Write();
  }

  if(sample=="ptheta"){
    cout << "write hists pthetapt" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_ptheta_xvar->Write();
    h_ptheta_xvar_generated->Write();
  }

  if(sample=="signed"){
    cout << "write hists signedpt" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_signed_xvar->Write();
    h_signed_xvar_generated->Write();
  }

  if(sample=="signeddalphat"){
    cout << "write hists signeddalphatpt" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_signeddalphat_xvar->Write();
    h_signeddalphat_xvar_generated->Write();
  }

  if(sample=="signeddphit"){
    cout << "write hists signeddphitpt" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_signeddphit_xvar->Write();
    h_signeddphit_xvar_generated->Write();
  }


  if(sample=="dthetaR"){
    cout << "write hists dthetaRpt" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_dthetaR_xvar->Write();
    h_dthetaR_xvar_generated->Write();
  }
  if(sample=="dthetaP"){
    cout << "write hists dthetaPpt" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_dthetaP_xvar->Write();
    h_dthetaP_xvar_generated->Write();
  }
  if(sample=="dPR"){
    cout << "write hists dPRpt" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_dPR_xvar->Write();
    h_dPR_xvar_generated->Write();
  }
  if(sample=="dPP"){
    cout << "write hists dPPpt" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_dPP_xvar->Write();
    h_dPP_xvar_generated->Write();
  }
  if(sample=="dPRi"){
    cout << "write hists dPRipt" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_dPRi_xvar->Write();
    h_dPRi_xvar_generated->Write();
  }
  if(sample=="dPPi"){
    cout << "write hists dPPipt" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_dPPi_xvar->Write();
    h_dPPi_xvar_generated->Write();
  }




  if(sample=="tkitki"){
    cout << "write hists tkitki" << endl;
    for(int i=0;i<200;i++){
      cout << "Migration " << i << endl;
      if(h_X_xvar_migration[i]==NULL) break;
      h_X_xvar_migration[i]->Write();
    }
    h_tkitki_xvar->Write();
    h_tkitki_xvar_generated->Write();
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
      "\t-./MigrationMatrixHists  Name_and_path_to_Output_file  Make_Flux_Constraint_Histo  Apply_Flux_Constraint  Multiplicity  Number_MC_files \n\n"<<
      "\t-Name_and_path_to_Output_file\t =\t Name of and path to the Output ROOT file that will be created \n"<<
      "\t-Playlist\t =\t Name of the playlist you want to run over e.g. minerva1 \n"<<
      "\t-Make_Flux_Constraint_Histo\t =\t If TRUE Enter 1; If FALSE Enter 0 \n"<< 
      "\t-Apply_Flux_Constraint\t =\t If TRUE Enter 1; If FALSE Enter 0 \n"<< 
      "\t-Multiplicity\t =\t Enter 0, 1 or 2. Here 0: 1 or 2 tracks; 1: 1-track events only; 2: 2-track events only \n"<<
      "\t-Sample to run\t =\t Enter pzmu,enu,q2 \n"<<
      "\t-Number_MC_files\t =\t Number of MonteCarlo files. To use all files, set this to -1" << std::endl; 
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    return 0; 
  }
  
  //! Default parameters
  std::vector<std::string> par;
  par.push_back("MigrationMatrixHists");
  par.push_back( Form("%s/ana/rootfiles/MigrationMatrixHists.root",getenv("CCQENUROOT") ) );
  par.push_back("minerva1");
  par.push_back("0");
  par.push_back("0");
  par.push_back("0");
  par.push_back("pzmu");
  par.push_back("-1");
  par.push_back("-1");//xtype
  par.push_back("-1");//ytype
  par.push_back("1");//proton, neutron,proton/neutron mode
  par.push_back("PT");//dependent var

  
  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }
  
  bool fluxConstraintHistoNeeded = ( par.at(3) == "1" ) ? true : false; 

  bool applyFluxConstraint = ( par.at(4) == "1" ) ? true : false; 

  if(par.at(11)!="PT" &&
     par.at(11)!="Q2"  )
  {
    cout << "Invalid xvar choice. Please choose PT or Q2. You used: "<<par.at(11) << endl;
    return 1;
  }
  if(par.at(6)!="dalphat" &&
     par.at(6)!="dphit"   &&
     par.at(6)!="pn"      &&
     par.at(6)!="dpt"     &&
     par.at(6)!="dptx"    &&
     par.at(6)!="dpty"    &&
     par.at(6)!="tp"      &&
     par.at(6)!="ptheta"  &&
     par.at(6)!="signed"  &&
     par.at(6)!="signeddalphat" &&
     par.at(6)!="signeddphit" &&
     par.at(6)!="dthetaR" &&
     par.at(6)!="dthetaP" &&
     par.at(6)!="dPR" &&
     par.at(6)!="dPP" &&
     par.at(6)!="dPRi" &&
     par.at(6)!="dPPi" &&
     par.at(6)!="tkitki"){
    cout << "Invalid sample choice. Please choose pzmu, enu, or q2. You used: "<<par.at(6) << endl;
    return 1;
  }



  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;
  
  if(par.at(6)=="tkitki"&& atoi(par[8].c_str())==-1){
    cout << "You need to specify the xtype with tkitki" << endl;
    return 1;
  }
  if(par.at(6)=="tkitki"&& atoi(par[9].c_str())==-1){
    cout << "You need to specify the ytype with tkitki" << endl;
    return 1;
  }


  return MigrationMatrixHists( par[1], par[2], fluxConstraintHistoNeeded, applyFluxConstraint, atoi(par[5].c_str()), par[6], atoi(par[7].c_str()), atoi(par[8].c_str()), atoi(par[9].c_str()), atoi(par[10].c_str()),par[11] );
}
