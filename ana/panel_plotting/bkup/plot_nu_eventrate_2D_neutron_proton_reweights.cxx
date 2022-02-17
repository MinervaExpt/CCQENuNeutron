//#include "myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvH1D.h"
// #include "../util/plot/plot.h"

#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "TEnv.h"
#include "TChain.h"
#include "TF2.h"
#include "Math/DistFunc.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TFile.h"
#include "TROOT.h"
#include "localColor.h"
#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"


#include "include/CCQENuPlotUtils.h"

#include "include/CCQENuUtils.h"
#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/HyperDimLinearizer.h"

#include "include/CommonBins.h"
#include "include/GeneralFunc.h"
#include "include/HDLFunc.h"




using namespace PlotUtils;
//forward declaration
//gROOT->SetBatch();

int nPlaylists = 12;
int nu = 1;
bool fluxHistoExists = true;
void draw2DPlots(string savename, TH2* data, TH2* dataStat, TH2* mc, TH2* mc_qelike_qe_h, TH2* mc_qelike_qe_oth, TH2* mc_qelike_res, TH2* mc_qelike_dis, TH2* mc_qelike_2p2h, string xtitle, string ytitle, string celltitle, string axis="x", double scale =1, double ymin=0, double ymax=2.4, double * multipliers=NULL, bool ratio = false, bool setlogx = false , bool setlogy = false);

void draw1DPlots(CCQENuPlotUtils* utils, string savename, MnvH2D** hists, double pot_data,  double pot_mc, int bmin=0, int bmax=-1, string axis="x", bool bkSub = false, bool logx = false, bool logy=false);

void draw1DPlotsAngle(CCQENuPlotUtils* utils, string savename, MnvH1D** hists, vector<MnvH1D*> &bg_weights, vector<MnvH1D*> &signal_weights, double pot_data,  double pot_mc, bool weigh_bkg, bool weigh_signal, bool bkSub, bool logx, bool logy);


void draw3DAnglePlots(string fname, string fname_bck, bool doBckFitting=true );
void draw3DAnglePlot( string savename, vector<TH2*> &hists, bool ratio, bool doMult, bool doBckSub, string xtitle, string ytitle, string ztitle, string celltitle, string axis, double scale, double ymin=0, double ymax=2, double *multipliers=NULL,  bool setlogx=true, bool setlogy=false, int imin=0, int imax=0);

int getRegion(double dthetaP, double dthetaR);
vector< vector<MnvH1D*>> combineAngularHistograms(MnvH3D** h3, CCQENuUtils* &utils, axis_binning &Q2bins);
void FormatHistos(vector< MnvH1D*> h );
void drawRegions( vector<vector<MnvH1D*>> &hists, bool doRatio, bool doBckSub, int nRegions );


void makePlots(bool doMultipliers, bool doBkgSub, string fname_signal, string fname_bkg_weight, string var, string sideband="Signal", bool drawRCut=false )
{
  string yvar="q2qe";
  bool doCenter = (var == "dthetaP" && drawRCut == true );
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  //three files 2-track, 2-track bkg_fitted
  //CV
  TFile* f1 = new TFile( fname_signal.c_str(), "READ");// original signal histogram
  TFile* f2 = new TFile( fname_bkg_weight.c_str(), "READ");// weighted histogram


  TVector2 *pot = (TVector2*)f1->Get("pot;1");
  //for( int i = 2; i<nPlaylists+1; i++ ) (*pot)+= *((TVector2*)f1->Get(Form( "pot;%d", i) ));


  //TVector2 *pot = (TVector2*)f1->Get("pot");
  double pot_data = pot->X();
  double pot_mc = pot->Y();
  double pot_norm = pot_data/pot_mc;

  CCQENuPlotUtils *putils = new CCQENuPlotUtils( fluxHistoExists );
  CCQENuUtils *utils = new CCQENuUtils( false, fluxHistoExists );
  axis_binning ybins = Q2bins;


  // Get CV histos
  MnvH2D *h_signal[nHistos];
  if (doCenter )
  {
    utils->bookHistos( f1, h_signal, Form( "h_%s_center_%s", var.c_str(), yvar.c_str() ) );
  }
  else 
  {
    utils->bookHistos( f1, h_signal, Form( "h_%s_%s", var.c_str(), yvar.c_str() ) );
  }
  putils->scaleMCHistos( h_signal, pot_norm );

  // Get Background weights
  MnvH2D *h_signal_bkg_weights = (MnvH2D*) f2->Get( Form( "h_weights_%s_yvarbins_qelikenot_%s", var.c_str(), sideband.c_str() ));

  // Get Tuned Backgrounds:
  MnvH2D* h_signal_bkg = (MnvH2D*) h_signal[kQELikeNot]->Clone( Form("h_%s_bkg", var.c_str() ) );
  h_signal_bkg->Multiply( h_signal_bkg, h_signal_bkg_weights );

  MnvH2D* h_signal_bkg_SNP = (MnvH2D*) h_signal[kQELikeNot_SingleNeutralPion]->Clone(Form("h_%s_bkg_snp", var.c_str() ) );

  //h_var_ptmu[kQELikeNot_SingleChargedPion]->Multiply( h_var_ptmu[kQELikeNot_SingleChargedPion], (MnvH2D*) f2->Get(Form("hs_weights_%sptbins_bgType_SingleChargedPion", var.c_str() ) ) );
  //h_var_ptmu[kQELikeNot_SingleNeutralPion]->Multiply( h_var_ptmu[kQELikeNot_SingleNeutralPion], (MnvH2D*) f2->Get(Form("hs_weights_%sptbins_bgType_SingleNeutralPion", var.c_str() ) ) );
  //h_var_ptmu[kQELikeNot_MultiPion]->Multiply( h_var_ptmu[kQELikeNot_MultiPion], (MnvH2D*) f2->Get(Form("hs_weights_%sptbins_bgType_MultiPion", var.c_str() ) ) );



  // Get Tuned Data Signal
  MnvH2D* h_signal_mc_nobkg = (MnvH2D*) h_signal[kMC]->Clone( Form("h_%s_mc_nobkg", var.c_str() ));
  if (doBkgSub ) h_signal_mc_nobkg->Add( h_signal[kQELikeNot], -1 );

  MnvH2D* h_signal_data_nobkg = (MnvH2D*) h_signal[kData]->Clone( Form("h_%s_data_nobkg", var.c_str() ));
  h_signal_data_nobkg->AddMissingErrorBandsAndFillWithCV( *h_signal_bkg );
  if (doBkgSub ) h_signal_data_nobkg->Add( h_signal_bkg, -1 );


  double multipliers[]={ 80, 32, 8, 4, 
                         1 , 1 , 1, 1,
                         1, 1, 2, 4,
                         40};



  //Saving unweighted data, MC into TH2D
  TH2* dataStat = new TH2D( h_signal_data_nobkg->GetCVHistoWithStatError() );
  TH2* data =     new TH2D( h_signal_data_nobkg->GetCVHistoWithError() );
  TH2* mc  =      new TH2D( h_signal_mc_nobkg->GetCVHistoWithStatError() );
  TH2* mc_qelike_qe_h= new TH2D( h_signal[kQELike_QE_H]->GetCVHistoWithStatError() );
  TH2* mc_qelike_qe_oth= new TH2D( h_signal[kQELike_QE_OTH]->GetCVHistoWithStatError() );
  TH2* mc_qelike_res= new TH2D( h_signal[kQELike_RES]->GetCVHistoWithStatError() );
  TH2* mc_qelike_dis= new TH2D(h_signal[kQELike_DIS]->GetCVHistoWithStatError() );
  TH2* mc_qelike_2p2h= new TH2D(h_signal[kQELike_2p2h]->GetCVHistoWithStatError() );

  // ----------------------------------------------------------------------------------
  // First make dthetaP in bins of q2qe
  // ----------------------------------------------------------------------------------

  // scaling histos by width
  string savename = "plots/nu-2d-"+sideband+"-evtrate-q2qe-"+var;
  string ext = ".pdf";
  string doBkg=(doBkgSub)? "_bkgsub": "";
  string xtitle = "#delta#theta_{P}", ytitle = "Event Rate (#times 10^{3}) per Degree";
  if( var == "dthetaR" ) xtitle = "#delta#theta_{R}";
  string celltitle = "Q^{2}_{QE}";
  double scale = 1e-3;
  double ymin = 0, ymax = 2.4;

  draw2DPlots(savename+doBkg+ext, data,  dataStat,  mc,  mc_qelike_qe_h, mc_qelike_qe_oth,  mc_qelike_res,  mc_qelike_dis,  mc_qelike_2p2h, xtitle, ytitle, celltitle, "x", scale, ymin,ymax, multipliers, false );

  savename =  "plots/nu-2d-"+sideband+"-evtrate-ratio-q2qe-"+var;
  ytitle = "Ratio to MnvTune";
  draw2DPlots(savename+doBkg+ext, data,  dataStat,  mc,  mc_qelike_qe_h, mc_qelike_qe_oth,  mc_qelike_res,  mc_qelike_dis,  mc_qelike_2p2h, xtitle, ytitle, celltitle, "x", scale, -.4,2.4, NULL, true );
  //gc->Print("plots/nu-2d-evtrate-model-q2qe-.png");
  //


  //Now let's 

  double multipliers2[]={ 18, 8, 8, 4, 2,2, 2, 1,
                          1, 2, 2, 2, 4, 8,8,18};



  // scaling histos by width
  savename = "plots/nu-2d-"+sideband+"-evtrate-"+var;
  if(doCenter) savename+="-center";
  xtitle = "Q^{2}_{QE} (GeV^{2})"; ytitle = "Event Rate (#times 10^{3}) per GeV^{2}";
  celltitle = "#delta#theta_{P}";
  if( var == "dthetaR" ) celltitle = "#delta#theta_{R}";
  scale = 1e-3;
  ymin = 0; ymax = 2.8;

  draw2DPlots(savename+doBkg+ext, data,  dataStat,  mc,  mc_qelike_qe_h, mc_qelike_qe_oth,  mc_qelike_res,  mc_qelike_dis,  mc_qelike_2p2h, xtitle, ytitle, celltitle, "y", scale, ymin,ymax, multipliers2, false,true );

  savename =  "plots/nu-2d-"+sideband+"-evtrate-ratio-model-"+var; 
  if(doCenter) savename+="-center";
  ytitle = "Ratio to MnvTune";
  draw2DPlots(savename+doBkg+ext, data,  dataStat,  mc,  mc_qelike_qe_h, mc_qelike_qe_oth,  mc_qelike_res,  mc_qelike_dis,  mc_qelike_2p2h, xtitle, ytitle, celltitle, "y", scale, -.4,2.4, NULL, true,true );

}

void draw2DPlots(string savename, TH2* data, TH2* dataStat, TH2* mc, TH2* mc_qelike_qe_h, TH2* mc_qelike_qe_oth, TH2* mc_qelike_res, TH2* mc_qelike_dis, TH2* mc_qelike_2p2h, string xtitle, string ytitle, string celltitle, string axis, double scale , double ymin, double ymax, double * multipliers, bool ratio , bool setlogx  , bool setlogy )
{
  //Plotting unweighed data/MC ratio
  vector<int> mycolors = getColors(2);
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);

  //need to add signal and bkg colors
  mc_qelike_qe_h->SetLineColor( mycolors[2] );
  mc_qelike_qe_oth->SetLineColor( mycolors[3] );
  mc_qelike_res->SetLineColor( mycolors[4] );
  mc_qelike_dis->SetLineColor( mycolors[5] );
  mc_qelike_2p2h->SetLineColor( mycolors[6] );

  data->SetMarkerStyle( 1 );
  data->SetMarkerSize(0.2);
  data->SetLineColor(kBlack);
  data->SetLineWidth(1);

  dataStat->SetMarkerStyle(1);
  dataStat->SetMarkerSize(0.2);
  dataStat->SetLineColor(kBlack);
  dataStat->SetLineWidth(1);

  //set options
  std::vector<std::pair<TH2*, const char*> > histAndOpts;


  histAndOpts.push_back(std::make_pair( (TH2*) dataStat->Clone("dataStat"), "e"));
  histAndOpts.push_back(std::make_pair( (TH2*) data->Clone("data"),     "ep"));
  histAndOpts.push_back(std::make_pair( (TH2*) mc->Clone("mc"),       "hist"));
  histAndOpts.push_back(std::make_pair( (TH2*) mc_qelike_qe_h->Clone("mc_qqeh"),       "histl"));
  histAndOpts.push_back(std::make_pair( (TH2*) mc_qelike_qe_oth->Clone("mc_qqeoth"),   "histl"));
  histAndOpts.push_back(std::make_pair( (TH2*) mc_qelike_res->Clone("mc_qres"),       "histl"));
  histAndOpts.push_back(std::make_pair( (TH2*) mc_qelike_dis->Clone("mc_qdis"),       "histl"));
  histAndOpts.push_back(std::make_pair( (TH2*) mc_qelike_2p2h->Clone("mc_q2p2h"),       "histl"));


  //for( int i = 0; i< histAndOpts.size();i++) histAndOpts[i].first->Scale(1,"width");

  TH2* mcRatio;
  if(ratio) mcRatio = (TH2*) mc->Clone("mcRatio");
 
  for( unsigned int i = 0 ; i< histAndOpts.size(); i++ ) 
  {
    if(ratio) histAndOpts[i].first->Divide(histAndOpts[i].first ,mcRatio);
    else histAndOpts[i].first->Scale(scale, "width");
  }
  GridCanvas* gc = (axis == "x")? plotXAxis1D( histAndOpts, xtitle, celltitle, multipliers ): plotYAxis1D( histAndOpts, xtitle, celltitle, multipliers ) ;
  gc->SetLogx(setlogx);
  gc->SetLogy(setlogy);
  gc->SetGridx(true);
  gc->SetYLimits(ymin, ymax);
  gc->SetYTitle(ytitle.c_str());
  gc->Modified();


  TLegend* leg=new TLegend(0.7, 0.1, 0.9, 0.3);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(data, "MINERvA data", "lpe");
  leg->AddEntry(mc, "MINERvA Tune", "l");
  leg->AddEntry(mc_qelike_qe_h,"QE-H","l");
  leg->AddEntry(mc_qelike_qe_oth,"QE-Oth","l");
  leg->AddEntry(mc_qelike_res,"Resonant","l");
  leg->AddEntry(mc_qelike_dis,"DIS","l");
  leg->AddEntry(mc_qelike_2p2h,"2p2h","l");
  //leg->AddEntry(mc_qelike_2p2h_no_lowrec,"2p2h without fit","l");
  leg->Draw("SAME");

  gc->Print(savename.c_str());
}

void FormatHistos( vector<vector<MnvH2D*>> h )
{
  vector<int> mycolors = getColors(2);
  int nChunks = h[0].size();

  for( int i = 0; i< nChunks; i++ )
  {
    h[kMC][i]->SetLineColor( kRed );
    h[kMC][i]->SetLineWidth(2);

    h[kQELike_QE_H][i]->SetLineColor( mycolors[2] );
    h[kQELike_QE_OTH][i]->SetLineColor( mycolors[3] );
    h[kQELike_RES][i]->SetLineColor( mycolors[4] );
    h[kQELike_DIS][i]->SetLineColor( mycolors[5] );
    h[kQELike_2p2h][i]->SetLineColor( mycolors[6] );
    h[kQELikeNot][i]->SetLineColor( mycolors[7] );

    h[kData][i]->SetMarkerStyle( 1 );
    h[kData][i]->SetMarkerSize(0.2);
    h[kData][i]->SetLineColor(kBlack);
    h[kData][i]->SetLineWidth(1);
  }
}



template<class T>
void DoFit( T** h, vector<T*> &fits )
{
  h[kMC]->Reset();
  h[kQELikeNot]->Reset();
  h[kQELike]->Reset();
  h[kQELike_QE_OTH]->Multiply( h[kQELike_QE_OTH], fits[0] );
  h[kQELike_2p2h]->Multiply( h[kQELike_2p2h], fits[1] );
  h[kQELike_RES]->Multiply( h[kQELike_RES], fits[2] );

  h[kQELikeNot_SingleChargedPion]->Multiply(h[kQELikeNot_SingleChargedPion],fits[3]);
  h[kQELikeNot_SingleNeutralPion]->Multiply(h[kQELikeNot_SingleNeutralPion],fits[4]);
  h[kQELikeNot_MultiPion]->Multiply(h[kQELikeNot_MultiPion],fits[5]);
  
  vector<int> qelike({ kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH });
  for( auto cat: qelike ) h[kQELike]->Add( h[cat] );
  vector<int> qelikenot({ kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } );
  for( auto cat: qelikenot ) h[kQELikeNot]->Add( h[cat] );

  h[kMC]->Add( h[kQELike] );
  h[kMC]->Add( h[kQELikeNot] );
}
void DoFit( MnvH1D **h, vector<MnvH1D*> &fits ){ return DoFit( h, fits ); };
void DoFit( MnvH2D **h, vector<MnvH2D*> &fits ){ return DoFit( h, fits ); };

template<class T>
void DoFit( vector<T*>&h, vector<T*> &fits )
{
  h[kMC]->Reset();
  h[kQELikeNot]->Reset();
  h[kQELike]->Reset();
  h[kQELike_QE_OTH]->Multiply( h[kQELike_QE_OTH], fits[0] );
  h[kQELike_2p2h]->Multiply( h[kQELike_2p2h], fits[1] );
  h[kQELike_RES]->Multiply( h[kQELike_RES], fits[2] );

  h[kQELikeNot_SingleChargedPion]->Multiply(h[kQELikeNot_SingleChargedPion],fits[3]);
  h[kQELikeNot_SingleNeutralPion]->Multiply(h[kQELikeNot_SingleNeutralPion],fits[4]);
  h[kQELikeNot_MultiPion]->Multiply(h[kQELikeNot_MultiPion],fits[5]);
  
  vector<int> qelike({ kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH });
  for( auto cat: qelike ) h[kQELike]->Add( h[cat] );
  vector<int> qelikenot({ kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } );
  for( auto cat: qelikenot ) h[kQELikeNot]->Add( h[cat] );

  h[kMC]->Add( h[kQELike] );
  h[kMC]->Add( h[kQELikeNot] );
}
void DoFit( vector<MnvH1D*>&h, vector<MnvH1D*> &fits ){ return DoFit( h, fits ); };
void DoFit( vector<MnvH2D*>&h, vector<MnvH2D*> &fits ){ return DoFit( h, fits ); };


template<class T> 
vector<T*> GetMCWeights( T**h, vector<T*> &fits )//kqelikenot, kqelike, kmc
{
  vector<T*> ret;
  T* qelikenot = (T*) h[kQELikeNot]->Clone("h_weight_qelikenot");
  T* qelike = (T*) h[kQELike]->Clone("h_weight_qelike");
  T* mc = (T*) h[kMC]->Clone("h_weight_mc");


  vector<int> vqelike({ kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH });
  vector<int> vqelikenot({ kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } );

  vector<T*> tmp;
  for( unsigned int i = 0; i< nHistos; i++ ) tmp.push_back(NULL);
  for( auto c:vqelike ) tmp[c] = (T*) h[c]->Clone( Form("tmp_%s", names[c].c_str() ) );
  for( auto c:vqelikenot ) tmp[c] = (T*) h[c]->Clone( Form("tmp_%s", names[c].c_str() ) );
  
  DoFit( tmp, fits );
  qelikenot->Divide( tmp[kQELikeNot], qelikenot );
  qelike->Divide( tmp[kQELike], qelike);
  mc->Divide( tmp[kMC], mc );
  return vector<T*>( {qelikenot, qelike, mc} );
}
vector<MnvH1D*>GetMCWeights( MnvH1D**h, vector<MnvH1D*> &fits ) { return GetMCWeights(h, fits ); };
vector<MnvH2D*>GetMCWeights( MnvH2D**h, vector<MnvH2D*> &fits ) { return GetMCWeights(h, fits ); };

void makePlots1D(string fname_signal, string fname_bkg_weight, string sideband="Signal",  bool doFit = true)//"#it{p}_{T}^{#mu} (GeV/#it{c})"
{
  cout<<"makePlot1D"<<endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  //three files 2-track, 2-track bkg_fitted
  //CV
  TFile* f1 = new TFile( fname_signal.c_str(), "READ");// original signal histogram
  TFile* f2 = new TFile( fname_bkg_weight.c_str(), "READ");// weighted histogram
  
  MnvPlotter *plotter = new MnvPlotter;

  MnvH2D *h_weight_2d_qe =   (MnvH2D*) f2->Get( Form("hs_weights_dthetaPdthetaR_yvarbins_bgType_qe_oth"));
  MnvH2D *h_weight_2d_2p2h = (MnvH2D*) f2->Get( Form("hs_weights_dthetaPdthetaR_yvarbins_bgType_2p2h"));
  MnvH2D *h_weight_2d_res =  (MnvH2D*) f2->Get( Form("hs_weights_dthetaPdthetaR_yvarbins_bgType_res"));


  MnvH2D *h_weight_2d_qelikenot_scp = (MnvH2D*) f2->Get( Form("hs_weights_dthetaPdthetaR_yvarbins_bgType_qelikenot_scp"));
  MnvH2D *h_weight_2d_qelikenot_snp = (MnvH2D*) f2->Get( Form("hs_weights_dthetaPdthetaR_yvarbins_bgType_qelikenot_snp"));
  MnvH2D *h_weight_2d_qelikenot_mp  = (MnvH2D*) f2->Get( Form("hs_weights_dthetaPdthetaR_yvarbins_bgType_qelikenot_mp"));

  MnvH1D *h_weights_q2qe_qelike_qe = (MnvH1D*) f2->Get( "hs_weights_yvar_bgType_qe_oth");
  MnvH1D *h_weights_q2qe_qelike_2p2h = (MnvH1D*) f2->Get( "hs_weights_yvar_bgType_2p2h");
  MnvH1D *h_weights_q2qe_qelike_res = (MnvH1D*) f2->Get( "hs_weights_yvar_bgType_res");

  MnvH1D *h_weight_1d_qelikenot_scp = (MnvH1D*) f2->Get( "hs_weights_yvar_bgType_qelikenot_scp" );
  MnvH1D *h_weight_1d_qelikenot_snp = (MnvH1D*) f2->Get( "hs_weights_yvar_bgType_qelikenot_snp" );
  MnvH1D *h_weight_1d_qelikenot_mp = (MnvH1D*) f2->Get( "hs_weights_yvar_bgType_qelikenot_mp" );


  vector<MnvH1D*> fits1D({h_weights_q2qe_qelike_qe,h_weights_q2qe_qelike_2p2h,h_weights_q2qe_qelike_res,h_weight_1d_qelikenot_scp,h_weight_1d_qelikenot_snp,h_weight_1d_qelikenot_mp});
  vector<MnvH2D*> fits2D({ h_weight_2d_qe,h_weight_2d_2p2h,h_weight_2d_res,h_weight_2d_qelikenot_scp,h_weight_2d_qelikenot_snp,h_weight_2d_qelikenot_mp});


  TVector2 *pot = (TVector2*)f1->Get("pot;1");
  //for( int i = 2; i<nPlaylists+1; i++ ) (*pot)+= *((TVector2*)f1->Get(Form( "pot;%d", i) ));
  double pot_data = pot->X();
  double pot_mc = pot->Y();
  double pot_norm = pot_data/pot_mc;

  cout<<"===================== POT ==================="<<endl;
  cout<<"Data: "<<pot_data<<endl;
  cout<<"MC: "<<pot_mc<<endl;
  cout<<"Ratio: "<<pot_norm<<endl;

  CCQENuPlotUtils *putils = new CCQENuPlotUtils( fluxHistoExists );
  CCQENuUtils *utils = new CCQENuUtils( false, fluxHistoExists );
  axis_binning xbins = dthetaPerpbins;
  axis_binning ybins = Q2bins;
  axis_binning zbins = dthetaReactbins;
  HyperDimLinearizer* hdl = GetHDL(xbins , ybins, zbins);


  MnvH1D *h_q2qe[nHistos];
  utils->bookHistos( f1, h_q2qe, "h_q2qe" );

  vector<MnvH1D*> mc_weights = GetMCWeights( h_q2qe, fits1D);
  MnvH1D* h_weight_q2qe_qelikenot = mc_weights[0];
  MnvH1D* h_weight_q2qe_qelike= mc_weights[1];
  MnvH1D* h_weight_q2qe_mc = mc_weights[2];

  

  vector<string> histnames({"h_dthetaPdthetaR_q2qe", "h_dthetaPdthetaR_q2qe_nc"});
  int iNC = 0;
  for( auto histname: histnames )
  {
    MnvH2D *h_var_yvar[nHistos];
    utils->bookHistos( f1, h_var_yvar, histname.c_str() );


    if(doFit) DoFit( h_var_yvar, fits2D );

    HDLHistos hdlHistos( hdl, xbins, ybins, zbins, h_var_yvar , true );
    FormatHistos( hdlHistos.GetHistos() );

    int tmin = hdlHistos.GetXaxis()->FindBin(-9.5);
    int tmax = hdlHistos.GetXaxis()->FindBin(9.5);
    vector<vector<MnvH2D*>> slices = hdlHistos.GetHistos();

    vector<MnvH1D*> h_dthetaP = hdlHistos.ProjectionsX( slices, "h_dthetaP" );
    vector<MnvH1D*> h_dthetaR = hdlHistos.ProjectionsZ( slices, "h_dthetaR" );
    vector<MnvH1D*> h_dthetaP_center=hdlHistos.ProjectionsX( slices, "h_dthetaP",tmin,tmax );
    vector<MnvH1D*> h_dthetaR_center=hdlHistos.ProjectionsZ( slices, "h_dthetaR",0,-1,tmin,tmax );

    for( int i = 0; i< nHistos; i++ )
    {
      h_dthetaP[i]->GetXaxis()->SetTitle("#delta#theta_P");
      h_dthetaP_center[i]->GetXaxis()->SetTitle("#delta#theta_P");
      h_dthetaR[i]->GetXaxis()->SetTitle("#delta#theta_R");
      h_dthetaR_center[i]->GetXaxis()->SetTitle("#delta#theta_R");
    }

    TCanvas* c = new TCanvas("c","c", 1200,800 );
    c->cd();
    string savename="";
    string xlabel = "#delta#dtheta_{P}(degree)";
    string ylabel = "evt";

    savename=Form("plots/nu-1d-dthetaP-%s-center_0-fit_%d-nc_%d.pdf",sideband.c_str(), doFit,iNC);
    putils->drawStacked( &(h_dthetaP[0]), "QELike_split_PionInFS", false, pot_data, pot_mc, true, -1, xlabel,ylabel );
    c->Print(savename.c_str());

    savename=Form("plots/nu-1d-dthetaP-%s-center_1-fit_%d-nc_%d.pdf",sideband.c_str(), doFit,iNC);
    putils->drawStacked( &(h_dthetaP_center[0]), "QELike_split_PionInFS", false, pot_data, pot_mc, true, -1, xlabel,ylabel );
    c->Print(savename.c_str());

    xlabel = "#delta#dtheta_{R}(degree)";
    savename=Form("plots/nu-1d-dthetaR-%s-center_0-fit_%d-nc_%d.pdf",sideband.c_str(), doFit,iNC);
    putils->drawStacked( &(h_dthetaR[0]), "QELike_split_PionInFS", false, pot_data, pot_mc, true, -1, xlabel,ylabel );
    c->Print(savename.c_str());

    savename=Form("plots/nu-1d-dthetaR-%s-center_1-fit_%d-nc_%d.pdf",sideband.c_str(), doFit,iNC);
    putils->drawStacked( &(h_dthetaR_center[0]), "QELike_split_PionInFS", false, pot_data, pot_mc, true, -1, xlabel,ylabel );
    c->Print(savename.c_str());

    iNC++;


  }
  //"hs_weights_dpt_yvarbins_bgType_MultiPions"
  //"hs_weights_dpt_yvarbins_bgType_MultiPions"
  cout<<"Done"<<endl;


  //TCanvas* c = new TCanvas("c","c", 1200,800 );
  //c->SetLogx();
  //c->SetGrid(); 
  //savename = Form( "plots/nu-1d-weights-qelikenot-%s.pdf", sideband.c_str() );
  //h_weight_1d_qelikenot->GetXaxis()->SetTitle("Q^{2}_{QE} (GeV^{2})" );
  //h_weight_1d_qelikenot->GetYaxis()->SetTitle("QelikeNot Weights" );
  //h_weight_1d_qelikenot->SetMaximum(2.5);
  //plotter->DrawMCWithErrorBand( h_weight_1d_qelikenot  );
  //c->Print( savename.c_str() );

  //savename = "plots/nu-1d-weights-qelikenot-scp.pdf";
  //h_weight_1d_qelikenot_scp->GetXaxis()->SetTitle("Q^{2}_{QE} (GeV^{2})" );
  //h_weight_1d_qelikenot_scp->GetYaxis()->SetTitle("#pi^{#pm} Weights" );
  //h_weight_1d_qelikenot_scp->SetMaximum(2.5);
  //plotter->DrawMCWithErrorBand( h_weight_1d_qelikenot_scp  );
  //c->Print( savename.c_str() );

  //savename = "plots/nu-1d-weights-qelikenot-snp.pdf";
  //h_weight_1d_qelikenot_snp->GetXaxis()->SetTitle("Q^{2}_{QE} (GeV^{2})" );
  //h_weight_1d_qelikenot_snp->GetYaxis()->SetTitle("#pi^{0} Weights" );
  //h_weight_1d_qelikenot_snp->SetMaximum(2.5);
  //plotter->DrawMCWithErrorBand( h_weight_1d_qelikenot_snp  );
  //c->Print( savename.c_str() );

  //savename = "plots/nu-1d-weights-qelikenot-mp.pdf";
  //h_weight_1d_qelikenot_mp->GetXaxis()->SetTitle("Q^{2}_{QE} (GeV^{2})" );
  //h_weight_1d_qelikenot_mp->GetYaxis()->SetTitle("N#pi Weights" );
  //h_weight_1d_qelikenot_mp->SetMaximum(2.5);
  //plotter->DrawMCWithErrorBand( h_weight_1d_qelikenot_mp );
  //c->Print( savename.c_str() );


  //delete c;

  //savename = "plots/nu-1d-"+sideband+"-evtrate-"+yvar+"-bckweighted";
  //if(doCenter) savename+="-center";
  //cout<<"Drawing "<<savename<<endl;
  //draw1DPlots(putils, savename+ext, h_var_yvar, pot_data,  pot_mc,0,-1, "y", false, logx);
  //savename = "plots/nu-1d-"+sideband+"-evtrate-"+var+"-bckweighted";
  //if(doCenter) savename+="-center";
  //cout<<"Drawing "<<savename<<endl;
  //draw1DPlots(putils, savename+ext, h_var_yvar, pot_data,  pot_mc,0,-1, "x", false);

  //savename = "plots/nu-1d-"+sideband+"-evtrate-"+yvar+"-bckweighted-nobck-";
  //if(doCenter) savename+="-center";
  //cout<<"Drawing "<<savename<<endl;
  //draw1DPlots(putils, savename+ext, h_var_yvar, pot_data,  pot_mc,0,-1, "y", true, logx);
  //savename = "plots/nu-1d-"+sideband+"-evtrate-"+var+"-bckweighted-nobck-";
  //if(doCenter) savename+="-center";
  //cout<<"Drawing "<<savename<<endl;
  //draw1DPlots(putils, savename+ext, h_var_yvar, pot_data,  pot_mc,0,-1, "x", true);




  //MnvPlotter *plotter = new MnvPlotter;
  //savename = "plots/nu-1d-evtrate-ptmu-qelikenot-weight.pdf";
  //TCanvas* c = new TCanvas("c","c", 1200,800 );
  //MnvH1D* h_qelikenot_weights1D = (MnvH1D*) h_qelikenot_weights->ProjectionY("h_qelikenot_weights1D",1,1);
  //h_qelikenot_weights1D->GetXaxis()->SetTitle("p_{T}^{#mu} (GeV/#it{c})" );
  //h_qelikenot_weights1D->GetYaxis()->SetTitle( "MC Weights" );

  //plotter->DrawMCWithErrorBand( h_qelikenot_weights1D );
  //c->Print(savename.c_str() );
}


void makePlotsAngle1D(string fname_signal, string fname_bkg_weight, string fname_signal_weight, string sideband="Signal")//"#it{p}_{T}^{#mu} (GeV/#it{c})"
{
  cout<<"makePlot1D"<<endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  //three files 2-track, 2-track bkg_fitted
  //CV
  TFile* f1 = new TFile( fname_signal.c_str(), "READ");// original signal histogram
  TFile* f2 = new TFile( fname_bkg_weight.c_str(), "READ");// weighted histogram
  TFile* f3 = new TFile( fname_signal_weight.c_str(), "READ");// weighted signal histogram
  

  TVector2 *pot = (TVector2*)f1->Get("pot;1");
  //for( int i = 2; i<nPlaylists+1; i++ ) (*pot)+= *((TVector2*)f1->Get(Form( "pot;%d", i) ));
  double pot_data = pot->X();
  double pot_mc = pot->Y();
  double pot_norm = pot_data/pot_mc;

  CCQENuPlotUtils *putils = new CCQENuPlotUtils( fluxHistoExists );
  CCQENuUtils *utils = new CCQENuUtils( false, fluxHistoExists );

  MnvH1D *h_q2qe_angle_00[nHistos],*h_q2qe_angle_01[nHistos],*h_q2qe_angle_02[nHistos],*h_q2qe_angle_03[nHistos];
  utils->bookHistos( f1, h_q2qe_angle_00, "h_q2qe_angle_00" );
  utils->bookHistos( f1, h_q2qe_angle_01, "h_q2qe_angle_01" );
  utils->bookHistos( f1, h_q2qe_angle_02, "h_q2qe_angle_02" );
  utils->bookHistos( f1, h_q2qe_angle_03, "h_q2qe_angle_03" );

  putils->scaleMCHistos( h_q2qe_angle_00, pot_norm );
  putils->scaleMCHistos( h_q2qe_angle_01, pot_norm );
  putils->scaleMCHistos( h_q2qe_angle_02, pot_norm );
  putils->scaleMCHistos( h_q2qe_angle_03, pot_norm );

  MnvH1D *h_qelike_weights_qe_oth = (MnvH1D*) f3->Get( Form("hs_weights_yvar_bgType_QE_OTH") );
  MnvH1D *h_qelike_weights_res = (MnvH1D*) f3->Get( Form("hs_weights_yvar_bgType_RES") );
  MnvH1D *h_qelike_weights_dis = (MnvH1D*) f3->Get( Form("hs_weights_yvar_bgType_DIS") );
  MnvH1D *h_qelike_weights_2p2h = (MnvH1D*) f3->Get( Form("hs_weights_yvar_bgType_2p2h") );

  MnvH1D *h_qelike_weights_angle_00 = (MnvH1D*) f3->Get("h_weights_q2qe_qelikenot_00");
  MnvH1D *h_qelike_weights_angle_01 = (MnvH1D*) f3->Get("h_weights_q2qe_qelikenot_01");
  MnvH1D *h_qelike_weights_angle_02 = (MnvH1D*) f3->Get("h_weights_q2qe_qelikenot_02");
  MnvH1D *h_qelike_weights_angle_03 = (MnvH1D*) f3->Get("h_weights_q2qe_qelikenot_03");

  MnvH1D *h_qelikenot_weights_scp = (MnvH1D*) f2->Get( "hs_weights_yvar_bgType_SingleChargedPion" );
  MnvH1D *h_qelikenot_weights_snp = (MnvH1D*) f2->Get( "hs_weights_yvar_bgType_SingleNeutralPion" );
  MnvH1D *h_qelikenot_weights_mp = (MnvH1D*) f2->Get( "hs_weights_yvar_bgType_MultiPions" );

  vector<MnvH1D*> bg_weights({ h_qelikenot_weights_scp, h_qelikenot_weights_snp, h_qelikenot_weights_mp} );
  vector<MnvH1D*> signal_weights({ /*h_qelike_weights_qe_oth,*/ h_qelike_weights_res, /*h_qelike_weights_dis,*/ h_qelike_weights_2p2h });


  cout<<"bookHistos done"<<endl;
  bool bkSub, weigh_signal, weigh_bkg;

  bool logx = true;

  string ext = ".pdf";
  string yvar = "q2qe";
  string savename = "plots/nu-angle-00-1d-"+sideband+"-evtrate-"+yvar+"-000";
  bkSub = false; weigh_bkg=false;weigh_signal=false;
  draw1DPlotsAngle(putils, savename+ext, h_q2qe_angle_00, bg_weights, signal_weights, pot_data,  pot_mc,weigh_bkg, weigh_signal, bkSub, logx, false);
  savename = "plots/nu-angle-00-1d-"+sideband+"-evtrate-"+yvar+"-100";
  weigh_bkg=true;weigh_signal=false;bkSub = false; 
  draw1DPlotsAngle(putils, savename+ext, h_q2qe_angle_00,bg_weights, signal_weights,  pot_data,  pot_mc,weigh_bkg, weigh_signal, bkSub, logx, false);
  savename = "plots/nu-angle-00-1d-"+sideband+"-evtrate-"+yvar+"-110";
  weigh_bkg=true;weigh_signal=true;bkSub = false; 
  draw1DPlotsAngle(putils, savename+ext, h_q2qe_angle_00,bg_weights, signal_weights,  pot_data,  pot_mc,weigh_bkg, weigh_signal, bkSub, logx, false);

  savename = "plots/nu-angle-00-1d-"+sideband+"-evtrate-"+yvar+"-111";
  weigh_bkg=true;weigh_signal=true;bkSub = true; 
  draw1DPlotsAngle(putils, savename+ext, h_q2qe_angle_00,bg_weights, signal_weights,  pot_data,  pot_mc,weigh_bkg, weigh_signal, bkSub, logx, false);

  savename = "plots/nu-angle-01-1d-"+sideband+"-evtrate-"+yvar+"-000";
  bkSub = false; weigh_bkg=false;weigh_signal=false;
  draw1DPlotsAngle(putils, savename+ext, h_q2qe_angle_01, bg_weights, signal_weights, pot_data,  pot_mc,weigh_bkg, weigh_signal, bkSub, logx, false);
  savename = "plots/nu-angle-01-1d-"+sideband+"-evtrate-"+yvar+"-100";
  weigh_bkg=true;weigh_signal=false;bkSub = false; 
  draw1DPlotsAngle(putils, savename+ext, h_q2qe_angle_01,bg_weights, signal_weights,  pot_data,  pot_mc,weigh_bkg, weigh_signal, bkSub, logx, false);
  savename = "plots/nu-angle-01-1d-"+sideband+"-evtrate-"+yvar+"-110";
  weigh_bkg=true;weigh_signal=true;bkSub = false; 
  draw1DPlotsAngle(putils, savename+ext, h_q2qe_angle_01,bg_weights, signal_weights,  pot_data,  pot_mc,weigh_bkg, weigh_signal, bkSub, logx, false);

  savename = "plots/nu-angle-01-1d-"+sideband+"-evtrate-"+yvar+"-111";
  weigh_bkg=true;weigh_signal=true;bkSub = true; 
  draw1DPlotsAngle(putils, savename+ext, h_q2qe_angle_01,bg_weights, signal_weights,  pot_data,  pot_mc,weigh_bkg, weigh_signal, bkSub, logx, false);

  savename = "plots/nu-angle-02-1d-"+sideband+"-evtrate-"+yvar+"-000";
  bkSub = false; weigh_bkg=false;weigh_signal=false;
  draw1DPlotsAngle(putils, savename+ext, h_q2qe_angle_02, bg_weights, signal_weights, pot_data,  pot_mc,weigh_bkg, weigh_signal, bkSub, logx, false);
  savename = "plots/nu-angle-02-1d-"+sideband+"-evtrate-"+yvar+"-100";
  weigh_bkg=true;weigh_signal=false;bkSub = false; 
  draw1DPlotsAngle(putils, savename+ext, h_q2qe_angle_02,bg_weights, signal_weights,  pot_data,  pot_mc,weigh_bkg, weigh_signal, bkSub, logx, false);
  savename = "plots/nu-angle-02-1d-"+sideband+"-evtrate-"+yvar+"-110";
  weigh_bkg=true;weigh_signal=true;bkSub = false; 
  draw1DPlotsAngle(putils, savename+ext, h_q2qe_angle_02,bg_weights, signal_weights,  pot_data,  pot_mc,weigh_bkg, weigh_signal, bkSub, logx, false);

  savename = "plots/nu-angle-02-1d-"+sideband+"-evtrate-"+yvar+"-111";
  weigh_bkg=true;weigh_signal=true;bkSub = true; 
  draw1DPlotsAngle(putils, savename+ext, h_q2qe_angle_02,bg_weights, signal_weights,  pot_data,  pot_mc,weigh_bkg, weigh_signal, bkSub, logx, false);

  savename = "plots/nu-angle-03-1d-"+sideband+"-evtrate-"+yvar+"-000";
  bkSub = false; weigh_bkg=false;weigh_signal=false;
  draw1DPlotsAngle(putils, savename+ext, h_q2qe_angle_03, bg_weights, signal_weights, pot_data,  pot_mc,weigh_bkg, weigh_signal, bkSub, logx, false);
  savename = "plots/nu-angle-03-1d-"+sideband+"-evtrate-"+yvar+"-100";
  weigh_bkg=true;weigh_signal=false;bkSub = false; 
  draw1DPlotsAngle(putils, savename+ext, h_q2qe_angle_03,bg_weights, signal_weights,  pot_data,  pot_mc,weigh_bkg, weigh_signal, bkSub, logx, false);
  savename = "plots/nu-angle-03-1d-"+sideband+"-evtrate-"+yvar+"-110";
  weigh_bkg=true;weigh_signal=true;bkSub = false; 
  draw1DPlotsAngle(putils, savename+ext, h_q2qe_angle_03,bg_weights, signal_weights,  pot_data,  pot_mc,weigh_bkg, weigh_signal, bkSub, logx, false);

  savename = "plots/nu-angle-03-1d-"+sideband+"-evtrate-"+yvar+"-111";
  weigh_bkg=true;weigh_signal=true;bkSub = true; 
  draw1DPlotsAngle(putils, savename+ext, h_q2qe_angle_03,bg_weights, signal_weights,  pot_data,  pot_mc,weigh_bkg, weigh_signal, bkSub, logx, false);




  //MnvPlotter *plotter = new MnvPlotter;
  //savename = "plots/nu-angle-1d-evtrate-qelikenot-qelikenot-weight.pdf";
  //TCanvas* c = new TCanvas("c","c", 1200,800 );
  //MnvH1D* h_qelikenot_weights1D = (MnvH1D*) h_qelikenot_weights->ProjectionY("h_qelikenot_weights1D",1,1);
  //h_qelikenot_weights1D->GetXaxis()->SetTitle("p_{T}^{#mu} (GeV/#it{c})" );
  //h_qelikenot_weights1D->GetYaxis()->SetTitle( "MC Weights" );

  //plotter->DrawMCWithErrorBand( h_qelikenot_weights1D );
  //c->Print(savename.c_str() );
}




void draw1DPlots(CCQENuPlotUtils* utils, string savename, MnvH2D** hists, double pot_data,  double pot_mc, int bmin, int bmax, string axis, bool bkSub, bool logx, bool logy)
{
  MnvH1D* h[nHistos];
  for( int i = 0; i< nHistos; i++ )
  {
    if( axis == "x") h[i]= (MnvH1D*) hists[i]->ProjectionX( Form( "%s_projX", hists[i]->GetName() ), bmin, bmax );
    else h[i]= (MnvH1D*) hists[i]->ProjectionY( Form( "%s_projY", hists[i]->GetName() ), bmin, bmax );
  }
  //mainly stacked plots
  //

  if(bkSub)
  {
    h[kData]->Add( h[kQELikeNot], -1 );
    h[kQELikeNot]->Reset();
    h[kQELikeNot_NoPions]->Reset();
    h[kQELikeNot_SingleNeutralPion]->Reset();
    h[kQELikeNot_SingleChargedPion]->Reset();
    h[kQELikeNot_MultiPion]->Reset();
  }
  string xlabel = ( axis == "x")? "#delta#theta_{P} (Degree)": "#it{p}_{T}^{#mu} (GeV/#it{c})";
  string ylabel = ( axis == "x")? "Event Rate (Degree)": "Event Rate (GeV/#it{c})";
  TCanvas* c = new TCanvas("c","c", 1200,800 );
  c->SetLogx( logx );
  c->SetLogy( logy );
  c->cd();
  cout<<"Start DrawStacked"<<endl;
  utils->drawStacked( h, "QELike_split_PionInFS", false, pot_data, pot_mc, true, -1, xlabel,ylabel );
  cout<<"Start Print"<<endl;
  c->Print(savename.c_str());

  delete c;
  for( int i = 0; i< nHistos; i++ )
  {
    delete h[i];
  }

  //utils->Print( c, savename );
}

void draw1DPlotsAngle(CCQENuPlotUtils* utils, string savename, MnvH1D** hists, vector<MnvH1D*> &bg_weights, vector<MnvH1D*> &signal_weights, double pot_data,  double pot_mc, bool weigh_bkg, bool weigh_signal, bool bkSub, bool logx, bool logy)
{
  MnvH1D* h[nHistos];
  for( int i = 0; i< nHistos; i++ )
  {
    h[i] = (MnvH1D*) hists[i]->Clone( Form("%s_copy", hists[i]->GetName() ) );
  }

  h[kData]->AddMissingErrorBandsAndFillWithCV( *h[kMC] );
  if( weigh_bkg )
  {
    h[kQELikeNot_SingleChargedPion]->Multiply(h[kQELikeNot_SingleChargedPion], bg_weights[0] );
    h[kQELikeNot_SingleNeutralPion]->Multiply(h[kQELikeNot_SingleNeutralPion], bg_weights[1] );
    h[kQELikeNot_MultiPion]->Multiply(h[kQELikeNot_MultiPion], bg_weights[2] );
    h[kQELikeNot]->Reset();
    h[kQELikeNot]->Add( h[kQELikeNot_SingleChargedPion] );
    h[kQELikeNot]->Add( h[kQELikeNot_SingleNeutralPion] );
    h[kQELikeNot]->Add( h[kQELikeNot_MultiPion] );
    h[kQELikeNot]->Add( h[kQELikeNot_NoPions] );
  }

  if( weigh_signal )
  {
    //h[kQELike_QE_OTH]->Multiply( h[kQELike_QE_OTH], signal_weights[0] );
    //h[kQELike_QE]->Multiply( h[kQELike_QE], signal_weights[0] );
    h[kQELike_RES]->Multiply( h[kQELike_RES], signal_weights[0] );
    //h[kQELike_DIS]->Multiply( h[kQELike_DIS], signal_weights[1] );
    h[kQELike_2p2h]->Multiply( h[kQELike_2p2h], signal_weights[1] );
    h[kQELike]->Reset();
    h[kQELike]->Add( h[kQELike_QE_OTH] );
    h[kQELike]->Add( h[kQELike_RES] );
    h[kQELike]->Add( h[kQELike_DIS] );
    h[kQELike]->Add( h[kQELike_2p2h] );
    h[kQELike]->Add( h[kQELike_OTH] );
  }



  //mainly stacked plots
  //

  if(bkSub)
  {
    h[kData]->Add( h[kQELikeNot], -1 );
    h[kQELikeNot]->Reset();
    h[kQELikeNot_NoPions]->Reset();
    h[kQELikeNot_SingleNeutralPion]->Reset();
    h[kQELikeNot_SingleChargedPion]->Reset();
    h[kQELikeNot_MultiPion]->Reset();
  }
  string xlabel = "Q^{2}_{QE} (GeV^{2})";
  string ylabel = "Event Rate (GeV^{2})";
  TCanvas* c = new TCanvas("c","c", 1200,800 );
  c->SetLogx( logx );
  c->SetLogy( logy );
  c->cd();
  cout<<"Start DrawStacked"<<endl;
  utils->drawStacked( h, "QELike_split_PionInFS", false, pot_data, pot_mc, true, -1, xlabel,ylabel );
  cout<<"Start Print"<<endl;
  c->Print(savename.c_str());

  delete c;
  for( int i = 0; i< nHistos; i++ )
  {
    delete h[i];
  }

  //utils->Print( c, savename );
}


void draw3DAnglePlots(string fname, string fname_bkg, bool doBckFitting )
{
  cout<<"draw3DAnglePlots"<<endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  //three files 2-track, 2-track bkg_fitted
  //CV
  cout<<"read "<<fname<<endl;
  TFile* f1 = new TFile( fname.c_str(), "READ");// original signal histogram
  TFile* f2 = new TFile( fname_bkg.c_str(), "READ");// original signal histogram
  if(!f1) cout<<"f1 i not read"<<endl;

  TVector2 *pot = (TVector2*)f1->Get("pot;1");
  //for( int i = 2; i<nPlaylists+1; i++ ) (*pot)+= *((TVector2*)f1->Get(Form( "pot;%d", i) ));


  //TVector2 *pot = (TVector2*)f1->Get("pot");
  double pot_data = pot->X();
  double pot_mc = pot->Y();
  double pot_norm = pot_data/pot_mc;
  cout<<"===================== POT ==================="<<endl;
  cout<<"Data: "<<pot_data<<endl;
  cout<<"MC: "<<pot_mc<<endl;
  cout<<"Ratio: "<<pot_norm<<endl;

  CCQENuPlotUtils *putils = new CCQENuPlotUtils( fluxHistoExists );
  CCQENuUtils *utils = new CCQENuUtils( false, fluxHistoExists );
  axis_binning xbins = dthetaPerpbins;
  axis_binning ybins = Q2bins;
  axis_binning zbins = dthetaReactbins;
  HyperDimLinearizer* hdl = GetHDL(xbins , ybins, zbins);



  // Get CV histos
  MnvH2D *h_hdl_signal[nHistos];
  cout<<"booking h2d"<<endl;
  utils->bookHistos( f1, h_hdl_signal, Form( "h_dthetaPdthetaR_q2qe") );
  putils->scaleMCHistos( h_hdl_signal, pot_norm );


  HDLHistos hdlHists( hdl, xbins, ybins, zbins, h_hdl_signal, true );

  MnvH2D *h_bck_weights_qelikenot = (MnvH2D*) f2->Get("h_weights_dthetaPdthetaR_yvarbins_qelikenot_Signal");
  MnvH2D *h_bck_weights_scp = (MnvH2D*) f2->Get("hs_weights_dthetaPdthetaR_yvarbins_bgType_SingleChargedPion");
  MnvH2D *h_bck_weights_snp = (MnvH2D*) f2->Get("hs_weights_dthetaPdthetaR_yvarbins_bgType_SingleNeutralPion");
  MnvH2D *h_bck_weights_mp  = (MnvH2D*) f2->Get("hs_weights_dthetaPdthetaR_yvarbins_bgType_MultiPions");

  if(doBckFitting)
  {
    h_hdl_signal[kMC]->Add(h_hdl_signal[kQELikeNot],-1 );

    h_hdl_signal[kQELikeNot_SingleChargedPion]->Multiply( h_hdl_signal[kQELikeNot_SingleChargedPion], h_bck_weights_scp );
    h_hdl_signal[kQELikeNot_SingleNeutralPion]->Multiply( h_hdl_signal[kQELikeNot_SingleNeutralPion], h_bck_weights_snp );
    h_hdl_signal[kQELikeNot_MultiPion]->Multiply( h_hdl_signal[kQELikeNot_MultiPion], h_bck_weights_mp );

    h_hdl_signal[kQELikeNot]->Reset();
    h_hdl_signal[kQELikeNot]->Add( h_hdl_signal[kQELikeNot_SingleChargedPion] );
    h_hdl_signal[kQELikeNot]->Add( h_hdl_signal[kQELikeNot_SingleNeutralPion] );
    h_hdl_signal[kQELikeNot]->Add( h_hdl_signal[kQELikeNot_MultiPion] );
    h_hdl_signal[kQELikeNot]->Add( h_hdl_signal[kQELikeNot_NoPions] );

    //h_bck_weights_qelikenot->Reset();
    //h_bck_weights_qelikenot->Add( h_hdl_signal[kQELikeNot_SingleChargedPion] );
    //h_bck_weights_qelikenot->Add( h_hdl_signal[kQELikeNot_SingleNeutralPion] );
    //h_bck_weights_qelikenot->Add( h_hdl_signal[kQELikeNot_MultiPion] );
    //h_bck_weights_qelikenot->Add( h_hdl_signal[kQELikeNot_NoPions] );
    //h_bck_weights_qelikenot->Divide( h_bck_weights_qelikenot, h_hdl_signal[kQELikeNot] );

    //h_hdl_signal[kQELikeNot]->Multiply( h_hdl_signal[kQELikeNot], h_bck_weights_qelikenot );


    h_hdl_signal[kMC]->Add(h_hdl_signal[kQELikeNot] );
  }




  MnvH3D *h_dthetaP_q2qe_dthetaR[nHistos];
  cout<<"booking h3d"<<endl;
  utils->bookHistos( h_dthetaP_q2qe_dthetaR, "h_dthetaP_q2qe_dthetaR","#delta#theta_{P}:Q^{2}_{QE}:#delta#theta_{R}", xbins, ybins,zbins );

  cout<<"converting 2D to 3D"<<endl;
  ConvertHDL2DTo3D( hdl, h_hdl_signal, h_dthetaP_q2qe_dthetaR );
  cout<<"converted"<<endl;
  cout<<Form("h_hdl_signal[0] has %.2f events", h_hdl_signal[0]->Integral() )<<endl;
  cout<<Form("h_dthetaP_q2qe_dthetaR[0] has %.2f events", h_dthetaP_q2qe_dthetaR[0]->Integral() )<<endl;

  // 1. Plot dthetaP vs Q2qe in bins of dthetaR
  TH3D *h3[nHistos];
  cout<<"casting h_dthetaP_q2qe_dthetaR"<<endl;
  for( int i = 0; i<nHistos; i++ ) 
  {
    cout<<"working on "<<i<<endl;
    h3[i] = new TH3D(h_dthetaP_q2qe_dthetaR[i]->GetCVHistoWithStatError());
  }
  cout<<"done"<<endl;
  //TH3D *h3_data_sys = new TH3D( h_dthetaP_q2qe_dthetaR[kData]->GetCVHistoWithStatError()); 
  //TH3D *h3_data_sys = new TH3D( h_dthetaP_q2qe_dthetaR[kData]->GetCVHistoWithError()); 
  //TH3D *h3_data_sys = h3[0];
  cout<<"created h3_data_sys"<<endl;

  TH2D *tmp[nHistos];
  TH2D *tmp_data_sys;
  
  double multipliers[]={ 80, 32, 8, 4, 
                         1 , 1 , 1, 1,
                         1, 1, 2, 4,
                         40};
  //double multipliers2[]={ 100, 80, 40, 20, 10,5, 1, 1,
  //                        1, 1, 4, 10, 20, 40,80,100};
  double multipliers2[]={ 100, 80, 40, 1, 1,1, 1, 1,
                          1, 1, 1, 1, 1, 40,80,100};


  cout<<"start plotting"<<endl;
  for( int i = 1; i< 13; i++ ) //dthetaR bins
  {
    int zbin=i+1;
    double dthetaR = h3[0]->GetZaxis()->GetBinCenter(zbin);

    //h3_data_sys->GetZaxis()->SetRange(zbin,zbin);
    //tmp_data_sys =(TH2D*) h3_data_sys->Project3D("xye");

    cout<<"working on "<<i<<"th dthetaR bin at "<<dthetaR<<endl;
    for(int j = 0; j< nHistos; j++ ) 
    {

      cout<<"h3["<<j<<"] has "<<h3[j]->Integral()<<"events"<<endl;
      h3[j]->GetZaxis()->SetRange(zbin,zbin);
      tmp[j]=(TH2D*) h3[j]->Project3D("xye");
      cout<<"tmp["<<j<<"] has "<<tmp[j]->Integral()<<"events"<<endl;
    }
    //dthetaP vs Q2qe

      //Saving unweighted data, MC into TH2D
      TH2* dataStat =        tmp[kData];
      TH2* data =            tmp[kData];
      TH2* mc  =             tmp[kMC];
      TH2* mc_qelike=        tmp[kQELike];
      TH2* mc_qelike_qe_h=     tmp[kQELike_QE_H];
      TH2* mc_qelike_qe_oth=     tmp[kQELike_QE_OTH];
      TH2* mc_qelike_res=    tmp[kQELike_RES];
      TH2* mc_qelike_dis=    tmp[kQELike_DIS];
      TH2* mc_qelike_2p2h=   tmp[kQELike_2p2h];
      TH2* mc_qelikenot =    tmp[kQELikeNot];
      TH2* mc_qelikenot_scp= tmp[kQELikeNot_SingleChargedPion];
      TH2* mc_qelikenot_snp= tmp[kQELikeNot_SingleNeutralPion];
      TH2* mc_qelikenot_mp=  tmp[kQELikeNot_MultiPion];
      TH2* mc_qelikenot_np=  tmp[kQELikeNot_NoPions];

      vector<TH2*> hists({dataStat, data, mc, mc_qelike, mc_qelike_qe_oth, mc_qelike_res,
          mc_qelike_dis, mc_qelike_2p2h,
          mc_qelikenot, mc_qelikenot_scp, mc_qelikenot_snp, mc_qelikenot_mp, mc_qelikenot_np,
          mc_qelike_qe_h});

      string folder="plots/";
      string xtitle = "Q^{2}_{QE} (GeV^{2})";
      string ytitle = "Evt Rate(10^{3})";
      string celltitle="#delta#theta_{P}";
      string ztitle = Form("#delta#theta_{R} = %.2f (degree)",dthetaR );
      double scale = 1e-3;
      double ymax =mc->GetMaximum()*scale*1.2;
      double ymin = -1*ymax/5; 

      //int imin = 0, imax = 0;
      int imin = 3, imax = 12;

      bool ratio = false, doMult = true, doBckSub = false;
      string savename = Form("%s/nu-3d-dthetaP-q2qe-dthetaR_%06.1f-ratio_%d-mult_%d-bckSub_%d",folder.c_str(),dthetaR,ratio,doMult, doBckSub );
      draw3DAnglePlot( savename, hists,ratio,doMult,doBckSub, xtitle, ytitle,ztitle,celltitle,"x",scale,ymin,ymax, multipliers2, true, false,imin,imax  );


      ratio = false; doMult = true; doBckSub = true;
      savename = Form("%s/nu-3d-dthetaP-q2qe-dthetaR_%06.1f-ratio_%d-mult_%d-bckSub_%d",folder.c_str(),dthetaR,ratio,doMult, doBckSub );
      draw3DAnglePlot( savename, hists,ratio,doMult,doBckSub, xtitle, ytitle,ztitle,celltitle,"x",scale,ymin,ymax, multipliers2, true, false,imin,imax  );


      folder="plots/";
      xtitle = "Q^{2}_{QE} (GeV^{2})";
      ytitle = "Ratio to MC";
      celltitle="#delta#theta_{P}";
      ztitle = Form("#delta#theta_{R} = %.2f (degree)",dthetaR );
      scale = 1e-3;

      ymin = -.4; ymax = 2.4;

      ratio = true; doMult = false; doBckSub = false;
      savename = Form("%s/nu-3d-dthetaP-q2qe-dthetaR_%06.1f-ratio_%d-mult_%d-bckSub_%d",folder.c_str(),dthetaR,ratio,doMult, doBckSub );
      draw3DAnglePlot( savename, hists,ratio,doMult,doBckSub, xtitle, ytitle,ztitle,celltitle,"x",scale,ymin,ymax, NULL, true, false,imin,imax  );

      ratio = true; doMult = false; doBckSub = true;
      savename = Form("%s/nu-3d-dthetaP-q2qe-dthetaR_%06.1f-ratio_%d-mult_%d-bckSub_%d",folder.c_str(),dthetaR,ratio,doMult, doBckSub );
      draw3DAnglePlot( savename, hists,ratio,doMult,doBckSub, xtitle, ytitle,ztitle,celltitle,"x",scale,ymin,ymax, NULL, true, false,imin,imax  );



  }
  //now just draw the angular plots dthetaP vs dthetaR

  MnvPlotter *plotter = new MnvPlotter();
  TCanvas* c = new TCanvas("c","c");

  //plotter->SetRedHeatPalette();
  plotter->SetWhiteRainbowPalette();
  vector<int> histoTypes({ kData, kMC,kQELike,kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_DIS, kQELike_2p2h } );
  vector<string> histoNames({ "data","MC","QELike", "QE-H", "QE-C", "RES", "DIS", "2p2h" });
  h3[kQELike]->GetZaxis()->SetRange(-1,0);
  TH2D* hmc = (TH2D*) h3[kQELike]->Project3D("zx");


  int linewidth=10;
  TBox *b0 = new TBox(-10,-10,10,10 );
  b0->SetLineColor(kBlack);
  b0->SetLineWidth(linewidth);
  b0->SetFillStyle(0);

  TBox *b1 = new TBox(-20,-20,20,10 );
  b1->SetLineColor(kBlue);
  b1->SetLineWidth(linewidth);
  b1->SetFillStyle(0);

  TBox *b2 = new TBox(-30,-30,30,10);
  b2->SetLineColor(kAzure);
  b2->SetLineWidth(linewidth);
  b2->SetFillStyle(0);

  TBox *b3 = new TBox(-20,-40,20,-30);
  b3->SetLineColor(kCyan+3);
  b3->SetLineWidth(linewidth);
  b3->SetFillStyle(0);
  TBox *b4 = new TBox(-55,-40,55,-30);
  b4->SetLineColor(kMagenta);
  b4->SetLineWidth(linewidth);
  b4->SetFillStyle(0);
  TBox *b5 = new TBox(-55,-55,55,-40);
  b5->SetLineColor(kRed);
  b5->SetLineWidth(linewidth);
  b5->SetFillStyle(0);






  
  for( int i = 0; i< histoTypes.size();i++ )
  {
    int type = histoTypes[i];
    string name = histoNames[i];
    string savename=Form("plots/nu-2d-dthetaP-dthetaR_%s", name.c_str() );
    h3[type]->GetZaxis()->SetRange(-1,0);
    TH2D* h = (TH2D*) h3[type]->Project3D("zx");
    h->GetXaxis()->SetRangeUser(-55,55);
    h->GetYaxis()->SetRangeUser(-55,55);



    h->GetXaxis()->SetTitle("#delta#theta_{P}");
    h->GetYaxis()->SetTitle("#delta#theta_{R}");

    h->GetZaxis()->SetRangeUser(0,20000);
    h->Draw("colz");
    plotter->AddHistoTitle(name.c_str());

    c->SetLogz();
    c->Print((savename+".pdf").c_str());
    c->Print((savename+".png").c_str());
    b0->Draw();
    b1->Draw();
    b2->Draw();
    b3->Draw();
    b5->Draw();
    b4->Draw();
    //b6->Draw();

    c->Print((savename+"-region.pdf").c_str());



    h->Divide( hmc );
    h->GetZaxis()->SetRangeUser(0,1);
    h->Draw("colz");




    c->SetLogz(false);
    string savename_ratio=Form("plots/nu-2d-ratio-dthetaP-dthetaR_%s", name.c_str() );
    plotter->AddHistoTitle(Form("%s ratio to QELike", name.c_str()));
    c->Print((savename_ratio+".pdf").c_str() );
    //c->Print((savename_ratio+".png").c_str() );
    b0->Draw();
    b1->Draw();
    b2->Draw();
    b3->Draw();
    b5->Draw();
    b4->Draw();
    //b6->Draw();

    c->Print((savename_ratio+"-region.pdf").c_str() );

  }
  c->Close();
  delete c;



  //cout<<"draw regions"<<endl;
  //vector<vector<MnvH1D*>> hvec_qsqqe = combineAngularHistograms( h_dthetaP_q2qe_dthetaR, utils, Q2bins);
  ////cout<<"parsed region"<<hvec_qsqqe[0][kData]->GetEntries()<<endl;
  //int nRegions = 7;
  //drawRegions( hvec_qsqqe, false, false, nRegions );
  //drawRegions( hvec_qsqqe, false, true, nRegions );
  //drawRegions( hvec_qsqqe, true, false, nRegions );
  //drawRegions( hvec_qsqqe, true, true, nRegions );

  f1->Close();
}

void FormatHistos( vector<MnvH1D*> h )
{
  vector<int> mycolors = getColors(2);

  h[kMC]->SetLineColor( kRed );
  h[kMC]->SetLineWidth(2);

  h[kQELike_QE_H]->SetLineColor( mycolors[2] );
  h[kQELike_QE_OTH]->SetLineColor( mycolors[3] );
  h[kQELike_RES]->SetLineColor( mycolors[4] );
  h[kQELike_DIS]->SetLineColor( mycolors[5] );
  h[kQELike_2p2h]->SetLineColor( mycolors[6] );
  h[kQELikeNot]->SetLineColor( mycolors[7] );

  h[kData]->SetMarkerStyle( 1 );
  h[kData]->SetMarkerSize(0.2);
  h[kData]->SetLineColor(kBlack);
  h[kData]->SetLineWidth(1);
}

void drawRegions( vector<vector<MnvH1D*>> &hists, bool doRatio, bool doBckSub, int nRegions )
{
  cout<<"enter drawRegions"<<endl;

  TCanvas *c = new TCanvas("c","c",1200,800);
  c->Divide(4,2);
  c->cd();
  vector<int> histosToDraw({ kMC, kQELikeNot, kQELike_DIS, kQELike_RES, kQELike_2p2h, kQELike_QE_OTH, kQELike_QE_H, kData} );
  vector<string> opts({"hist", "histsame", "histlsame","histlsame","histlsame","histlsame","histlsame","pesame"});
  //vector<TH1D*> th1_histos;
  cout<<"looping over regions"<<endl;
  TLatex *text = new TLatex();
  text->SetTextFont(92);
  text->SetTextSize(0.1);
  text->SetTextColor( kBlue );
  for( int i = 0; i< hists.size(); i++ )
  {
    c->cd(i+1);
    gPad->SetLogx();
    gPad->SetGridx();
    gPad->SetMargin(.15,.00,.15,.1);
    FormatHistos( hists[i] );
    cout<<"Formatted Histos"<<endl;

    TH1D* hmc = new TH1D( hists[i][kMC]->GetCVHistoWithStatError() );
    hmc->SetName("hmcbase");
    TH1D* hqelikenot = new TH1D( hists[i][kQELikeNot]->GetCVHistoWithStatError());
    hqelikenot->SetName("hqelikenotbase");

    if(doBckSub) 
    {
      hmc->Add( hqelikenot, -1 );
    }


    //hists[i][kData]->Draw("pe");
    for( int j = 0; j< histosToDraw.size(); j++ )
    {
      int J = histosToDraw[j];
      cout<<"histo to draw: "<<J<<endl;
      TH1D *h= new TH1D( hists[i][J]->GetCVHistoWithStatError() );
      h->SetName( Form("%d%d", i,j) );
      //th1_histos.push_back(h);
      h->GetXaxis()->SetTitle("Q^{2}_{QE} (GeV^{2})");
      h->SetMinimum(0);
      if( doBckSub && ( J == kData || J == kMC || J == kQELikeNot ) )
      {
        h->Add( hqelikenot, -1 );
      }

      if( doRatio ){
        h->SetMaximum(2.5);
        h->Divide( hmc );
      }
      else
      {
        h->Scale(1,"width");
        h->SetMaximum( h->GetMaximum()*1.5 );
      }

      h->GetXaxis()->SetNdivisions(4);
      h->Draw(opts[j].c_str());
    }
    text->DrawLatexNDC(0.2,0.8, Form("Region %d", i));
    //delete hqelikenot;
    //delete hmc;
  }
  c->cd(8);
  TLegend* leg=new TLegend(0.1, 0.1, 0.9, 0.9);
  leg->SetTextSize(1);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(52);
  leg->SetTextSize(0.05);
  leg->AddEntry(hists[0][kData], "MINERvA data", "lpe");
  leg->AddEntry(hists[0][kMC], "MINERvA Tune", "l");
  leg->AddEntry(hists[0][kQELike_QE_H],"QE-H","l");
  leg->AddEntry(hists[0][kQELike_QE_OTH],"QE-Oth","l");
  leg->AddEntry(hists[0][kQELike_RES],"Resonant","l");
  leg->AddEntry(hists[0][kQELike_DIS],"DIS","l");
  leg->AddEntry(hists[0][kQELike_2p2h],"2p2h","l");
  leg->AddEntry(hists[0][kQELikeNot],"Not QELike","l");
  //leg->AddEntry(mc_qelike_2p2h_no_lowrec,"2p2h without fit","l");
  leg->Draw("SAME");

  
  c->Print(Form("plots/nu-1d-regions-ratio_%d-bcksub-%d.pdf", doRatio, doBckSub));
  delete c;
}

void draw3DAnglePlot( string savename, vector<TH2*> &hists, bool ratio, bool doMult, bool doBckSub, string xtitle, string ytitle, string ztitle, string celltitle, string axis, double scale, double ymin, double ymax, double *multipliers,  bool setlogx, bool setlogy, int imin, int imax)
{
  cout<<"draw3DAnglePlot"<<endl;
  double *mult = (doMult)? multipliers: NULL;
  TH2* dataStat =          (TH2*) hists[0]->Clone("h_0");
  TH2* data =              (TH2*) hists[1]->Clone("h_1");
  TH2* mc  =               (TH2*) hists[2]->Clone("h_2");
  TH2* mc_qelike   =       (TH2*) hists[3]->Clone("h_3");
  TH2* mc_qelike_qe_oth=       (TH2*) hists[4]->Clone("h_4");
  TH2* mc_qelike_res=      (TH2*) hists[5]->Clone("h_5");
  TH2* mc_qelike_dis=      (TH2*) hists[6]->Clone("h_6");
  TH2* mc_qelike_2p2h=     (TH2*) hists[7]->Clone("h_7");
  TH2* mc_qelikenot =      (TH2*) hists[8]->Clone("h_8");
  TH2* mc_qelikenot_scp=   (TH2*) hists[9]->Clone("h_9");
  TH2* mc_qelikenot_snp=   (TH2*) hists[10]->Clone("h_10");
  TH2* mc_qelikenot_mp=    (TH2*) hists[11]->Clone("h_11");
  TH2* mc_qelikenot_np=    (TH2*) hists[12]->Clone("h_12");
  TH2* mc_qelike_qe_h = (TH2*) hists[13]->Clone("h_13");

  if(doBckSub)
  {
    dataStat->Add(mc_qelikenot, -1 );
    data->Add(mc_qelikenot, -1 );
    mc->Add(mc_qelikenot, -1 );
    mc_qelikenot->Reset();
  }

  vector<int> mycolors = getColors(2);
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);

  mc_qelike_qe_oth->SetLineColor( mycolors[3] );
  mc_qelike_qe_h->SetLineColor( mycolors[2] );
  mc_qelike_res->SetLineColor( mycolors[4] );
  mc_qelike_dis->SetLineColor( mycolors[5] );
  mc_qelike_2p2h->SetLineColor( mycolors[6] );
  mc_qelikenot->SetLineColor( mycolors[7] );

  data->SetMarkerStyle( 20 );
  data->SetMarkerSize(.5);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2);
                                  
  dataStat->SetMarkerStyle(20);
  dataStat->SetMarkerSize(.5);
  dataStat->SetLineColor(kBlack);    
  dataStat->SetLineWidth(1);       

  std::vector<std::pair<TH2*, const char*> > histAndOpts;



  histAndOpts.push_back( std::make_pair( mc_qelikenot, "hist"));
  histAndOpts.push_back( std::make_pair( mc_qelike_dis, "histl"));
  histAndOpts.push_back( std::make_pair( mc_qelike_res, "histl"));
  histAndOpts.push_back( std::make_pair( mc_qelike_2p2h, "histl"));
  histAndOpts.push_back( std::make_pair( mc_qelike_qe_oth, "histl"));
  if(!nu) histAndOpts.push_back( std::make_pair( mc_qelike_qe_h, "histl"));
  histAndOpts.push_back( std::make_pair( mc, "hist"));
  histAndOpts.push_back( std::make_pair( data, "ep"));
  histAndOpts.push_back( std::make_pair( dataStat, "e"));

  TH2* mcRatio;
  if(ratio) mcRatio = (TH2*) mc->Clone("mcRatio");
 
  for( unsigned int i = 0 ; i< histAndOpts.size(); i++ ) 
  {
    if(ratio) histAndOpts[i].first->Divide(histAndOpts[i].first ,mcRatio);
    else histAndOpts[i].first->Scale(scale, "width");
  }
  string Xtitle = xtitle+"---"+ztitle;
  GridCanvas* gc = (axis == "x")? plotXAxis1D( histAndOpts, Xtitle, celltitle, mult ,imin,imax): plotYAxis1D( histAndOpts, Xtitle, celltitle, mult,imin,imax ) ;
  gc->SetLogx(setlogx);
  gc->SetLogy(setlogy);
  gc->SetGridx(true);
  gc->SetYLimits(ymin, ymax);
  gc->SetYTitle(ytitle.c_str());
  gc->Modified();


  TLegend* leg=new TLegend(0.7, 0.1, 0.9, 0.3);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(data, "MINERvA data", "lpe");
  leg->AddEntry(mc, "MINERvA Tune", "l");
  if(!nu) leg->AddEntry(mc_qelike_qe_h,"QE-H","l");
  leg->AddEntry(mc_qelike_qe_oth,"QE-C","l");
  leg->AddEntry(mc_qelike_2p2h,"2p2h","l");
  leg->AddEntry(mc_qelike_res,"RES","l");
  leg->AddEntry(mc_qelike_dis,"DIS","l");
  if(!doBckSub) leg->AddEntry(mc_qelikenot,"NotQELike","l");
  leg->Draw("SAME");


  gc->Print((savename+".pdf").c_str());
  gc->Print((savename+".png").c_str());
}



int getRegion(double dthetaP, double dthetaR)
{
  if( abs(dthetaP) > 55 || dthetaR> 10 || -55 > dthetaR )          return -1;
  if( abs(dthetaP) < 10 && 10 > dthetaR && dthetaR > -10 )           return 0; //hydrogen
  else if ( abs(dthetaP) < 20 && abs(dthetaR)<20 )               return 1; //QE-C 1
  else if ( abs(dthetaP) < 30 && abs(dthetaR)<30 )              return 2; //QE-C 3
  else if ( abs(dthetaP) < 20 && dthetaR<-30 && dthetaR >=-40  ) return 3; //2p2h/RES 4
  else if ( abs(dthetaP) > 20 && dthetaR<-30 && dthetaR >=-40  ) return 4; //2p2h/RES 1
  else if ( abs(dthetaP) < 55 && dthetaR>-55 && dthetaR < -40  ) return 5; //2p2h/RES 2
  else return -1;
}
vector<vector<MnvH1D*>> combineAngularHistograms(MnvH3D** h3, CCQENuUtils* &utils,axis_binning &Q2bins)
{
  //this function will combine bins

  const int nRegions = 6;
  vector<vector<MnvH1D*>> ret(nRegions);
  MnvH1D* h1temp = (MnvH1D*) h3[kMC]->ProjectionX();
  for( int i = 0; i< ret.size(); i++ )
  {
    MnvH1D *h1[nHistos];
    utils->bookHistos( h1, Form("h_dthetaP_q2qe_sideband_%02d",i),"#delta#theta_{P}:Q^{2}_{QE}", Q2bins );
    for( int j = 0; j< nHistos; j++ )
    {
      h1[j]->AddMissingErrorBandsAndFillWithCV( *h1temp );
      ret[i].push_back(h1[j]);
    }
  }

  for (int x = 1; x<h3[0]->GetNbinsX()+1; ++x )
  {
    double xmin = h3[0]->GetXaxis()->GetBinLowEdge(x); //dthetaP
    double xmax = h3[0]->GetXaxis()->GetBinUpEdge(x);
    double xc = h3[0]->GetXaxis()->GetBinCenter(x);
    if( xmin < -55 || xmax > 55 ) continue;
    for( int z = 1; z<h3[0]->GetNbinsZ()+1; ++z )
    {
      double zmin = h3[0]->GetZaxis()->GetBinLowEdge(z); //dthetaR
      double zmax = h3[0]->GetZaxis()->GetBinUpEdge(z);
      if( zmin < -55 || zmax > 10 ) continue;
      double zc = h3[0]->GetZaxis()->GetBinCenter(z);
      int region = getRegion(xc,zc);
      cout<<"Region: "<<region<<" at "<<xc<<", "<<zc<<endl;
      if (region == -1) continue;

      for (int h = 0; h<nHistos;h++)
      {
        ret[region][h]->Add( ( h3[h]->ProjectionY("htemp",x,x,z,z)) );
      }
    }
  }
  

  return ret;
}



int main(int argc, char* argv[])
{

  string f_orig = argv[1];
  string f_bkg_weighted = argv[2];
  string f_signal_weighted = argv[3];
  string sideband = argv[4];
  nu = atoi( argv[5] );
  cout<<f_orig<<endl;
  if (argc == 7 ) nPlaylists = atoi(argv[6]);
  //multipliers
  string var = "dthetaP";
  //makePlots(true,false,f_orig, f_bkg_weighted, var, sideband);
  //makePlots(true,true,f_orig, f_bkg_weighted, var, sideband);


  ////draw R cut
  //makePlots(true,false,f_orig, f_bkg_weighted, var, sideband, true);
  //makePlots(true,true,f_orig, f_bkg_weighted, var, sideband, true);

  //var = "dthetaR";
  //makePlots(true,false,f_orig, f_bkg_weighted, var, sideband);
  //makePlots(true,true,f_orig, f_bkg_weighted, var, sideband);

  ////dthetaP and ptmu
  string yvarname = "Q^{2}_{QE} (GeV^{2})";
  string yvar = "q2qe";
  string varname = "#delta#theta_{P} (degree)";
  bool doCenter = false;
  var = "dthetaP";
  makePlots1D( f_orig, f_signal_weighted, sideband, true);
  makePlots1D( f_orig, f_signal_weighted, sideband, false);

  draw3DAnglePlots( f_orig, f_bkg_weighted, true ); //do bck fitting
  //makePlotsAngle1D( f_orig, f_bkg_weighted, f_signal_weighted );
  return 0;
}
