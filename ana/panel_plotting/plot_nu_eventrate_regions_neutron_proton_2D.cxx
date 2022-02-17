//#include "myPlotStyle.h"
#include "TParameter.h"

#include "include/GeneralIncludes.h"


using namespace PlotUtils;
//forward declaration
//gROOT->SetBatch();

int nPlaylists = 12;
int nu = 1;
bool fluxHistoExists = true;
bool drawNuwroOnly = false;
bool hasRESFit = false;
//vector<int> categories_to_fit({ kQELike_QE_OTH, kQELike_2p2h, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion} );
//vector<string> categories_to_fit_names({"qe_oth", "2p2h","qelikenot_scp","qelikenot_snp","qelikenot_mp"});
vector<int> categories_to_fit({ kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion} );
vector<string> categories_to_fit_names({"qelike_qe_oth", "qelike_res","qelike_2p2h","qelikenot_scp", "qelikenot_snp"});
vector<string> categories_to_fit_title_names({"QELike && QE OTH", "QELike && RES", "QELike && 2p2h", "QELikeNot && Single Charged Pion", "QELikeNot && Single Neutral Pion"});


//vector<int> categories_to_fit({ kQELike_QE_OTH, kQELike_2p2h, kQELikeNot_SinglePion} );
//vector<string> categories_to_fit_names({"qelike_qe_oth", "qelike_2p2h","qelikenot_sp"});
//vector<string> categories_to_fit_title_names({"QELike && QE OTH", "QELike && 2p2h", "QELikeNot && Single Pion"});


vector<int> histosUsed({ kData,  kMC, kQELike, kQELike_QE_H, kQELike_QE_OTH, kQELike_2p2h, kQELike_RES, kQELike_DIS, kQELike_OTH, kQELikeNot, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } ); 


//General Func
template<class T>
void FormatHistos(vector< T*> &h );
void FormatHistos( vector<vector<MnvH2D*>> h );
void drawSystematics( CCQENuPlotUtils*utils, MnvH1D* h,  double pot_data, double pot_mc, string name);

TH2D* convert3DToHDL( MnvH2D* h, TH3D* h3, HyperDimLinearizer* hdl, bool resetWeight = false )
{
  cout<<"entering convert3DToHDL"<<endl;
  //h3-> dthetaP vs dthetaR vs Q2qe
  //h -> dthetaP vs Q2qe vs dthetaR
  TH2D* ret = new TH2D( h->GetCVHistoWithStatError() );
  ret->Reset();
  for( int bx = 1; bx < h3->GetNbinsX()+1; bx++)
  {
    double x = h3->GetXaxis()->GetBinCenter( bx );
    for( int by = 1; by < h3->GetNbinsY()+1; by++ )
    {
      double y = h3->GetYaxis()->GetBinCenter( by );
      for( int bz = 1; bz < h3->GetNbinsZ()+1; bz++ )
      {
        double z = h3->GetZaxis()->GetBinCenter( bz );

        double globalx = hdl->GetBin( vector<double>( {x,z,y} ) ).first+0.0001;
        double w = h3->GetBinContent(bx,by,bz);
        int binx = ret->GetXaxis()->FindBin( globalx);
        int biny = ret->GetYaxis()->FindBin( z );
        if(w==0. && resetWeight) 
        {
          cout<<"weights: "<<x<<", "<<y<<", "<<z<<": "<<w<<endl;
          w = 1.;
          cout<<"to fill: "<<x<<", "<<y<<", "<<z<<": "<<w<<endl;
        }
        ret->SetBinContent(binx, biny, w );


        //pair<int,int> xy = hdl->GetBin( vector<double>( {x,z,y} ) );
        //double w = h3->GetBinContent(bx,by,bz);
        //ret->SetBinContent( xy.first, xy.second, w );
      }
    }
  }
  return ret;
}

TH2D* reweighQELike_QE_OTH( MnvH2D** hist, TH3D* h3w, HyperDimLinearizer* hdl )
{
  TH2D* hw = convert3DToHDL( hist[kMC], h3w, hdl, true );
  
  MnvH2D* h = hist[kQELike_QE_OTH];
  hist[kMC]->Add( h, -1 );
  hist[kQELike_QE]->Add( h, -1 );

  cout<<"Multiplying Single"<<endl;
  h->MultiplySingle( h, hw );

  hist[kMC]->Add( h );
  hist[kQELike_QE]->Add( h );
  cout<<"done reweighting"<<endl;

  return hw;
}

template<class MnvHND>
void DoFit( MnvHND**h, vector<MnvHND*>&fits )
{
  for( int i = 0; i< categories_to_fit.size(); i++ )
  {
    int cat = categories_to_fit[i];
    cout<<"apply fit: "<<names[cat]<<endl;
    cout<<"weight: "<<fits[i]->GetName()<<endl;
    h[cat]->Multiply( h[cat],fits[i] );
  }

  //h[kQELikeNot_SingleNeutralPion]->Multiply(h[kQELikeNot_SingleNeutralPion], fits[2]);
  //h[kQELikeNot_SingleChargedPion]->Multiply(h[kQELikeNot_SingleChargedPion], fits[2]);
 

  h[kMC]->Reset();
  h[kQELike]->Reset();
  h[kQELikeNot]->Reset();


  h[kQELike_QE]->Reset();
  h[kQELike_QE]->Add(h[kQELike_QE_OTH]);
  h[kQELike_QE]->Add(h[kQELike_QE_H]);
 
  vector<int> qelike({ kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH });
  for( auto cat: qelike ) h[kQELike]->Add( h[cat] );
  vector<int> qelikenot({ kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } );
  for( auto cat: qelikenot ) h[kQELikeNot]->Add( h[cat] );

  h[kMC]->Add( h[kQELike] );
  h[kMC]->Add( h[kQELikeNot] );
}
template void DoFit( MnvH1D** h, vector<MnvH1D*>&fits );
template void DoFit( MnvH2D** h, vector<MnvH2D*>&fits );
//void DoFit( MnvH1D**h, vector<MnvH1D*>&fits )
//{
//  for( int i = 0; i< categories_to_fit.size(); i++ )
//  {
//    int cat = categories_to_fit[i];
//    cout<<"apply fit: "<<names[cat]<<endl;
//    cout<<"weight: "<<fits[i]->GetName()<<endl;
//    h[cat]->Multiply( h[cat],fits[i] );
//  }
// 
//  h[kQELikeNot_SingleNeutralPion]->Multiply(h[kQELikeNot_SingleNeutralPion], fits[2]);
//  h[kQELikeNot_SingleChargedPion]->Multiply(h[kQELikeNot_SingleChargedPion], fits[2]);
//
//  h[kMC]->Reset();
//  h[kQELike]->Reset();
//  h[kQELikeNot]->Reset();
//
//
//  h[kQELike_QE]->Reset();
//  h[kQELike_QE]->Add(h[kQELike_QE_OTH]);
//  h[kQELike_QE]->Add(h[kQELike_QE_H]);
// 
//  vector<int> qelike({ kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH });
//  for( auto cat: qelike ) h[kQELike]->Add( h[cat] );
//  vector<int> qelikenot({ kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } );
//  for( auto cat: qelikenot ) h[kQELikeNot]->Add( h[cat] );
//
//  h[kMC]->Add( h[kQELike] );
//  h[kMC]->Add( h[kQELikeNot] );
//}

template<class T>
void SubtractBackground( T** h )
{
  h[kData]->AddMissingErrorBandsAndFillWithCV( *h[kMC] );
  h[kData]->Add( h[kQELikeNot], -1 );
  h[kMC]->Add( h[kQELikeNot], -1 );
  h[kQELikeNot]->Reset();
  h[kQELikeNot_SingleNeutralPion]->Reset();
  h[kQELikeNot_SingleChargedPion]->Reset();
  h[kQELikeNot_MultiPion]->Reset();
  h[kQELikeNot_NoPions]->Reset();
}

template<class T>
void RatioToMC( T** h )
{
  T* hmc = (T*) h[kMC]->Clone("hmc");
  vector<int> categories{ kData, kMC, kQELike, kQELike_QE_H, kQELike_QE_OTH, kQELike_QE, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH, kQELikeNot, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions };

  for( auto cat : categories ) h[cat]->Divide(h[cat], hmc);
}

vector<MnvH2D*> GetFitsHisto( TFile* f, string xvar )
{
  vector<MnvH2D*> ret;
  for( unsigned int i = 0; i<categories_to_fit_names.size();i++)
  {
    ret.push_back( (MnvH2D*) f->Get( Form( "hs_weights_%s_yvarbins_bgType_%s", xvar.c_str(), categories_to_fit_names[i].c_str() ) ) );
  }
  return ret;
}
vector<MnvH1D*> GetFitsHisto( TFile* f  )
{
  vector<MnvH1D*> ret;
  for( unsigned int i = 0; i<categories_to_fit_names.size();i++)
  {
    ret.push_back( (MnvH1D*) f->Get( Form( "hs_weights_yvar_bgType_%s", categories_to_fit_names[i].c_str() ) ) );
  }
  return ret;
}

//======================= Draw 1D Distributions  ===============================//

void draw1DDistros( TFile *f_orig, TFile* f_signal_weighted, string sample );
void draw1DDistro( MnvH2D** h, vector<MnvH2D*> &fits, CCQENuPlotUtils* putils, string name, double pot_data, double pot_mc, string axis = "x", string xtitle="x", string ytitle="evt rate", bool doRatio=false, bool doFitting=false, bool doBckSub = false);

void draw1DDistros(  TFile *f_orig, TFile* f_signal_weighted, string sample )
{
  cout<<"draw1DDistros"<<endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  //three files 2-track, 2-track bkg_fitted
  //CV
  //cout<<"read "<<fname<<endl;
  TFile* f1 = f_orig;
  TFile* f2 = f_signal_weighted;
  if(!f1) cout<<"f1 i not read"<<endl;

  TVector2 *pot = (TVector2*)f1->Get("pot");

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
 
  cout<<"==================== Book Histograms ======="<<endl;
  MnvH2D *h_dpty_q2qe[nHistos], *h_dptx_q2qe[nHistos], *h_dpt_q2qe[nHistos], 
         *h_pn_q2qe[nHistos], *h_dalphat_q2qe[nHistos], *h_dphit_q2qe[nHistos],
         *h_dthetaP_q2qe[nHistos], *h_dthetaR_q2qe[nHistos], *h_ptheta_q2qe[nHistos],
         *h_dtheta2D_q2qe[nHistos],
         *h_muontheta_q2qe[nHistos],*h_muonmomentum_q2qe[nHistos];

  utils->bookHistos( f1, h_dtheta2D_q2qe, "h_dtheta2D_q2qe" );
  cout<<"================= Draw Histos ==============="<<endl;

  vector<MnvH2D*> fits =  GetFitsHisto(f2, "dtheta2D");
  draw1DDistro( h_dtheta2D_q2qe, fits, putils, "plot-1d-dtheta2D", pot_data, pot_mc, "x", "#delta#theta_{2D} (degree)", "Event Rate", 0,0,0);
  draw1DDistro( h_dtheta2D_q2qe, fits, putils, "plot-1d-dtheta2D", pot_data, pot_mc, "x", "#delta#theta_{2D} (degree)", "Event Rate", 0,1,0);
  draw1DDistro( h_dtheta2D_q2qe, fits, putils, "plot-1d-dtheta2D", pot_data, pot_mc, "x", "#delta#theta_{2D} (degree)", "Event Rate", 0,1,1);
  draw1DDistro( h_dtheta2D_q2qe, fits, putils, "plot-1d-dtheta2D", pot_data, pot_mc, "x", "#delta#theta_{2D} (degree)", "Ratio To MC", 1,1,1);
}

void draw1DDistro( MnvH2D** h, vector<MnvH2D*> &fits, CCQENuPlotUtils* putils, string name, double pot_data, double pot_mc, string axis, string xtitle, string ytitle, bool doRatio, bool doFitting, bool doBckSub)
{
  bool area_norm = false;
  bool includeData = true;

  MnvH2D* hists[nHistos];
  for(unsigned int i = 0; i<nHistos; i++) hists[i] = new MnvH2D(*( h[i]->Clone(Form("h2d_tmp_%s_%s", h[i]->GetName(), names[i].c_str() ))));

  if(doFitting) DoFit( hists, fits );
  if(doBckSub ) SubtractBackground( hists );

  MnvH1D* hists1d[nHistos];

  for(unsigned int i = 0; i<nHistos; i++ ) 
  {
    MnvH1D* htmp;
    if(axis == "x") hists1d[i] = new MnvH1D(*(hists[i]->ProjectionX(Form("h1d_tmp_%s", hists[i]->GetName()),1,hists[i]->GetNbinsY())) );
    if(axis == "y") hists1d[i] = new MnvH1D(*(hists[i]->ProjectionY(Form("h1d_tmp_%s", hists[i]->GetName()),1,hists[i]->GetNbinsX())) );
    hists1d[i]->GetXaxis()->SetTitle( xtitle.c_str() );
    hists1d[i]->GetYaxis()->SetTitle( ytitle.c_str() );
    if(!doRatio) hists1d[i]->GetYaxis()->SetRangeUser(0,0.02e6);
  }
  //MnvH1D** hists1d = &hists1d_cp[0];
  if(doRatio ) RatioToMC( hists1d );

  TCanvas* c = new TCanvas("cNeutA","Neutron Angulars"); 
  MnvPlotter *plotter = new MnvPlotter;



  double mcscale = -1.; string xaxislabel = xtitle; string yaxislabel = ""; double min_x = -1.; double max_x = -1.; double min_y = -1.; double max_y = -1.;
  bool normalizeHisto = (doRatio)? false: true;
  if (doRatio)
  {
    min_y = 0;
    max_y = 2;
  }
  else
  {
    double value = hists1d[kData]->GetBinNormalizedCopy().GetMaximum();
    max_y = value*1.4;
    min_y = 0;
  }

  putils->drawStacked(hists1d, "QELike_split_PionInFS", area_norm, pot_data, pot_mc, includeData,mcscale,xaxislabel,yaxislabel,min_x, max_x, min_y ,max_y,normalizeHisto );
  TH1D* herr = (TH1D*) hists1d[kMC]->Clone("herr");
  if(!doRatio) herr->Scale(herr->GetBinWidth(1),"width");
  herr->SetFillColorAlpha(kRed, 0.3);
  herr->Draw("e2same");


  //TH1D* mcErr;
  //if(normalizeHisto)  mcErr = (TH1D*) hists1d[kMC]->GetBinNormalizedCopy().GetCVHistoWithError().Clone("err");
  //else mcErr = (TH1D*) hists1d[kMC]->GetCVHistoWithError().Clone("err");
  //mcErr->SetFillStyle(1001);
  //mcErr->SetFillColorAlpha( kRed, 0.35 );
  //mcErr->SetLineWidth(0);
  //mcErr->SetMarkerSize(0);
  //mcErr->Draw("e2same");


  string fname=Form("plots/%s-ratio_%d-bcksub_%d-fit_%d.pdf",name.c_str(), doRatio, doBckSub, doFitting );
  plotter->AddPlotLabel(Form("# Hydrogen: %.2f", hists1d[kQELike_QE_H]->Integral()), 0.3,0.8,0.033,kBlue,32);
  c->Print( fname.c_str() );

  return;
}


//======================= Draw 3D Distributions  ===============================//
// in chunks of dthetaR
void draw3DChunk( string fullinput, string input_bkg, string sample="Signal", bool doFitting=false, string model="" );
void drawChunk(string savename, string chunkname, string axis, HDLHistos &hdlHistos, double *multipliers, bool doRatio, bool doBckSub, bool doFitting );
void drawChunkRegions( HDLHistos &hdlHistos, bool doFitting );

void drawChunkForRegion( vector<vector<MnvH1D*>> &hists, vector<string>&titles, int iRegion, bool doRatio, bool doBckSub,  bool subQELike,bool hasFitting );

void drawChunkForOne(string savename, string chunkname, string axis, vector<TH2D*>&hists, HDLHistos &hdlHistos, double *multipliers, bool doScale );

vector<TH2D*> split3D( TH3D* hist );

void draw3DChunk( string fullinput, string input_bkg, string sample, bool doFitting, string model )
{
  cout<<"draw3DChunk"<<endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  //three files 2-track, 2-track bkg_fitted
  //CV
  cout<<"read "<<fullinput<<endl;
  TFile* f1 = new TFile( fullinput.c_str(), "READ");// original signal histogram
  TFile* f2 = new TFile( input_bkg.c_str(), "READ");// original signal histogram
  if(!f1) cout<<"f1 i not read"<<endl;

  TVector2 *pot = (TVector2*)f1->Get("pot");

  double pot_data = pot->X();
  double pot_mc = pot->Y();
  double pot_norm = pot_data/pot_mc;

  cout<<"===================== POT ==================="<<endl;
  cout<<"Data: "<<pot_data<<endl;
  cout<<"MC: "<<pot_mc<<endl;
  cout<<"Ratio: "<<pot_norm<<endl;

  //CCQENuPlotUtils *putils = new CCQENuPlotUtils( fluxHistoExists );
  CCQENuUtils *utils = new CCQENuUtils( false, fluxHistoExists );
  axis_binning xbins = dthetaPerpbins;
  axis_binning ybins = Q2bins;
  axis_binning zbins = dthetaReactbins;
  HyperDimLinearizer* hdl = GetHDL(xbins , ybins, zbins);

  //============== Get distribution =============
  cout<<"booking histo"<<endl;
  MnvH2D *h_dthetaPdthetaR_q2qe[nHistos];
  utils->bookHistos( f1, h_dthetaPdthetaR_q2qe, "h_dthetaPdthetaR_q2qe");

  TFile *f3 = new TFile("/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/ratios/reweights_2.root", "READ");
  TH3D* ratio_nuwroLFG = (TH3D*) f3->Get("nuwroLFG_mnvGENIEv2_dthetaP_dthetaR_q2qe_qelike_qe_oth_ratio_ratio");
  TH3D* ratio_nuwroSF = (TH3D*) f3->Get("nuwroSF_mnvGENIEv2_dthetaP_dthetaR_q2qe_qelike_qe_oth_ratio_ratio");

  TFile *f4 = new TFile("/minerva/app/users/tejinc/Generators/nuisance/WorkArea/AnalysisCode/MINERVA/reweights/HydrogenAnalysis_Nu_GenieNuWro.root");
  TH3D* nuwroLFG = (TH3D*) f4->Get("nuwroLFG_kPRQsqQE_ccqe");
  TH3D* nuwroSF = (TH3D*) f4->Get("nuwroSF_kPRQsqQE_ccqe");


  MnvH2D* hw;
  vector<TH2D*> hnuwro3DDirect;
  MnvH2D* hnuwroHDL;

  if(model == "nuwroSF") 
  {
    cout<<"reweighting nuwroSF"<<endl;
    hw = (MnvH2D*) reweighQELike_QE_OTH( h_dthetaPdthetaR_q2qe, ratio_nuwroSF, hdl );
    hnuwro3DDirect = split3D( ratio_nuwroSF );
    TH2D* h2dp = convert3DToHDL( h_dthetaPdthetaR_q2qe[kMC], ratio_nuwroSF, hdl );
    hnuwroHDL = new MnvH2D( *h2dp );
  }
  else if(model == "nuwroLFG") 
  {
    cout<<"reweighting nuwroLFG"<<endl;
    hw =(MnvH2D*)  reweighQELike_QE_OTH( h_dthetaPdthetaR_q2qe, ratio_nuwroLFG, hdl );
    hnuwro3DDirect = split3D( ratio_nuwroLFG );
    TH2D* h2dp = convert3DToHDL( h_dthetaPdthetaR_q2qe[kMC], ratio_nuwroLFG, hdl );
    hnuwroHDL = new MnvH2D( *h2dp );
  }
  else
  {
    cout<<"don't reweight"<<endl;
  }


  //============== Get Background Fittings =============
  cout<<"get category fits"<<endl;
  //MnvH2D *h_fit_qelike_qe_oth = (MnvH2D*) f2->Get("hs_weights_dthetaPdthetaR_yvarbins_bgType_qe_oth");
  //MnvH2D *h_fit_qelike_res = (MnvH2D*) f2->Get("hs_weights_dthetaPdthetaR_yvarbins_bgType_res");
  //MnvH2D *h_fit_qelike_2p2h = (MnvH2D*) f2->Get("hs_weights_dthetaPdthetaR_yvarbins_bgType_2p2h");
  //MnvH2D *h_fit_qelikenot_scp = (MnvH2D*) f2->Get("hs_weights_dthetaPdthetaR_yvarbins_bgType_qelikenot_scp");
  //MnvH2D *h_fit_qelikenot_snp = (MnvH2D*) f2->Get("hs_weights_dthetaPdthetaR_yvarbins_bgType_qelikenot_snp");
  //MnvH2D *h_fit_qelikenot_mp = (MnvH2D*) f2->Get("hs_weights_dthetaPdthetaR_yvarbins_bgType_qelikenot_mp");
  //vector<MnvH2D*> fits( {h_fit_qelike_qe_oth,h_fit_qelike_2p2h,h_fit_qelike_res,h_fit_qelikenot_scp,h_fit_qelikenot_snp, h_fit_qelikenot_mp});
  vector<MnvH2D*> fits = GetFitsHisto( f2, "dthetaPdthetaR" );
  if(doFitting) DoFit( h_dthetaPdthetaR_q2qe, fits );
  cout<<"Define hdlHistos"<<endl;
  HDLHistos hdlHistos( hdl, xbins, ybins, zbins, h_dthetaPdthetaR_q2qe , true );
  // vector for category -- vector of MnvH2D per dthetaR
  //vector< vector<MnvH2D*> > fullChunkHists = hdlHistos.GetHistos(); 
  cout<<"Format Histos"<<endl;
  FormatHistos( hdlHistos.GetHistos() );

  double multipliers_Q2qe[]={ 80, 32, 8, 4, 
                         1 , 1 , 1, 1,
                         1, 1, 2, 4,
                         40};
  double multipliers_dthetaP[]={ 100, 80, 40, 1, 1,1, 1, 1,
                          1, 1, 1, 1, 1, 40,80,100};

  //============================================================
  if ( hnuwroHDL )
  {
    bool scale = false;
    cout<<"Get2DHistos"<<endl;
    //vector<TH2D*> hnuwro3DHDL = hdl->Get2DHistos( hnuwroHDL, true );
    cout<<"StartDrawing"<<endl;
    //drawChunkForOne( "nuwro3DHDL", model, "y", hnuwro3DHDL, hdlHistos, NULL, scale );
    drawChunkForOne( "nuwro3DDirect", model, "y", hnuwro3DDirect, hdlHistos, NULL, scale );
  }

  if(drawNuwroOnly) return;
  //============================================================



 
  cout<<"Draw Chunks"<<endl;
  //drawChunk("nu-3d","dthetaR","y", hdlHistos, multipliers_dthetaP, false,false,doFitting);
  //drawChunk("nu-3d","dthetaR","y", hdlHistos, multipliers_dthetaP, false,true,doFitting);
  //drawChunk("nu-3d","dthetaR","y", hdlHistos, NULL, true,false,doFitting);
  //drawChunk("nu-3d","dthetaR","y", hdlHistos, NULL, true,true,doFitting);

  drawChunkRegions(hdlHistos, doFitting);

  //f1->Close();
  //f2->Close();
  cout<<"closed files"<<endl;
  //for(int i =0; i<nHistos;i++) delete h_dthetaPdthetaR_q2qe[nHistos];
  //cout<<"deleted h_dthetaPdthetaR_q2qe"<<endl;
  ////delete[] h_dthetaPdthetaR_q2qe;
  ////cout<<"deleted h_dthetaPdthetaR_q2qe"<<endl;
  //delete hdl;
  //cout<<"deleted hdl"<<endl;

  return;
}


void drawChunk(string savename, string chunkname, string axis, HDLHistos &hdlHistos, double *multipliers, bool doRatio, bool doBckSub, bool doFitting )
{
  string histsname=Form("plots/%s-%s-ratio_%d-bcksub_%d-fit_%d.pdf",savename.c_str(),chunkname.c_str(), doRatio, doBckSub, doFitting) ;
  cout<<"draw "<<histsname<<endl;
  vector<vector<MnvH2D*>> hists = hdlHistos.GetHistos();
  int nChunk = hists[0].size()-2;
  vector<int> histsToUse({kData,kMC, kQELike, kQELikeNot, kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions});

  for( int i = 0; i<nChunk; i++ )
  {
    cout<<"chunk: "<<i<<endl;
    int zbin = i+1;
    double zvalue = hdlHistos.GetZaxis()->GetBinCenter(zbin);
    if(zvalue<-55 || zvalue >55) continue;
    cout<<"zvalue = "<<zvalue<<endl;
    TH2* dataStat = (TH2*) hists[kData][zbin]->GetCVHistoWithStatError().Clone("dataStat");
    //for( int z = 0; z<dataStat->GetNbinsY();z++) cout<<dataStat->GetYaxis()->GetBinCenter(z+1)<<endl;
    TH2* mcErr =(TH2*)  hists[kMC][zbin]->GetCVHistoWithError().Clone("mcErr");
    TH2* mc =(TH2*)  hists[kMC][zbin]->GetCVHistoWithError().Clone("mc");
    mcErr->SetFillStyle(1001);
    mcErr->SetFillColorAlpha(kBlue,0.35);
    mcErr->SetMarkerSize(0 );

    vector<TH2*>histErr;
    cout<<"Creating histErr"<<endl;
    for(unsigned int j = 0; j< nHistos;j++) histErr.push_back( NULL );
    for(auto c: histsToUse) histErr[c] = (TH2*)hists[c][zbin]->GetCVHistoWithError().Clone(Form("h_%s",names[c].c_str() ) ) ;

    //histAndOpts.push_back( std::make_pair( &dataStat, "e" ) );
    cout<<"Creating histAndOpts"<<endl;
    vector<pair<TH2*, const char*> > histAndOpts;
    histAndOpts.push_back( std::make_pair( mcErr, "e2" ) );
    histAndOpts.push_back( std::make_pair( histErr[kMC], "hist" ) );
    histAndOpts.push_back( std::make_pair( histErr[kQELike_QE_H], "histl" ) );
    histAndOpts.push_back( std::make_pair( histErr[kQELike_QE_OTH], "histl" ) );
    histAndOpts.push_back( std::make_pair( histErr[kQELike_RES], "histl" ) );
    histAndOpts.push_back( std::make_pair( histErr[kQELike_2p2h], "histl" ) );
    histAndOpts.push_back( std::make_pair( histErr[kQELike_DIS], "histl" ) );
    histAndOpts.push_back( std::make_pair( histErr[kQELikeNot], "hist" ) );
    histAndOpts.push_back( std::make_pair( histErr[kData], "ep" ) );

    if( doBckSub )
    {
      dataStat->Add( histErr[kQELikeNot], -1 );
      histErr[kData]->Add( histErr[kQELikeNot], -1 );

      mcErr->Add( histErr[kQELikeNot], -1 );
      mc->Add( histErr[kQELikeNot], -1 );
      histErr[kMC]->Add( histErr[kQELikeNot], -1 );
      histErr[kQELikeNot]->Reset();
    }

    for( unsigned int j = 0; j< histAndOpts.size(); j++ )
    {
      if( doRatio ) 
      {
        histAndOpts[j].first->Divide( (histAndOpts[j].first), mc );
      }
      else 
      {
        histAndOpts[j].first->RebinY(2);
        histAndOpts[j].first->Scale(0.001,"width");
      }
    }
    string xtitle="#delta#theta_{P}";
    string ytitle="Q^{2}_{QE} (GeV^{2})";
    GridCanvas* gc = (axis=="x")? plotXAxis1D( histAndOpts, xtitle, ytitle, multipliers ): plotYAxis1D( histAndOpts, ytitle, xtitle, multipliers );
    if(axis=="y") gc->SetLogx(true);
    else gc->SetLogx(false);
    gc->SetGridx(true);
    gc->SetYLimits(-.2,1.2);
    if(doRatio)gc->SetYLimits(-.2,2.2);
    gc->SetYTitle("Evt Rate #times 10^{3}");
    if(doRatio) gc->SetYTitle("Ratio to MnvGENIE");
    gc->Modified();

    if( axis=="y" ) gc->SetXLimits(0,12);

    TLegend* leg=new TLegend(0.7, 0.1, 0.9, 0.3);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry( dataStat, "MINERvA data", "lpe");
    leg->AddEntry( mc, "MINERvA Tune", "l");
    leg->AddEntry( histErr[kQELike_QE_H],"QE-H","l");
    leg->AddEntry( histErr[kQELike_QE_OTH],"QE-Oth","l");
    leg->AddEntry( histErr[kQELike_RES],"Resonant","l");
    leg->AddEntry( histErr[kQELike_2p2h],"2p2h","l");
    leg->AddEntry( histErr[kQELike_DIS],"DIS","l");
    leg->Draw("SAME");

    string fname=Form("plots/%s-%s-%s_%05d-ratio_%d-bcksub_%d-fit_%d.pdf",savename.c_str(),axis.c_str(),chunkname.c_str(), (int)zvalue, doRatio, doBckSub, doFitting) ;
    gc->Print(fname.c_str());
    free(dataStat);
    free(mcErr);
    free(mc);
    for(unsigned int h = 0; h< nHistos;h++) free(histErr[h]);
  }
}



void drawChunkRegions( HDLHistos &hdlHistos, bool doFitting ){
  cout<<"entering draw chunk regions"<<endl;
  
  vector<vector<MnvH2D*>> hists = hdlHistos.GetHistos(); //3D histo dthetaP vs Qsq vs dthetaR, hists[category][chunk]

  int nChunk = hists[0].size()-2;
  vector<int> histsToUse({kData,kMC, kQELike, kQELikeNot, kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions});

  int nRegions = 6;

  //region-- angle bin -- category
  vector< vector< vector<MnvH1D*>>> regioned_histos1D( nRegions);
  vector< vector<string>> regioned_histoAngles( nRegions );
  cout<<"start parsing regions"<<endl;
  for( int iz=0;iz<hists[0].size()+1; iz++ )
  {
    int zbin = iz+1;
    double dthetaR = hdlHistos.GetZaxis()->GetBinCenter(zbin);
    if(dthetaR<-55 || dthetaR >55) continue;
    cout<<"dthetaR = "<<dthetaR<<endl;
    double minR = hdlHistos.GetZaxis()->GetBinLowEdge(zbin);
    double maxR = hdlHistos.GetZaxis()->GetBinUpEdge(zbin);
    for( int ix = 0; ix< hists[0][0]->GetNbinsX()+1; ix++ )
    {
      int xbin = ix+1;
      double dthetaP = hdlHistos.GetXaxis()->GetBinCenter(xbin);
      if(dthetaP<-55 || dthetaP >55) continue;
      double minP = hdlHistos.GetXaxis()->GetBinLowEdge(xbin);
      double maxP = hdlHistos.GetXaxis()->GetBinUpEdge(xbin);
      int iRegion = getRegions(dthetaP, dthetaR);
      if(iRegion < 0 ) continue;

      cout<<"dthetaP = "<<dthetaP<<endl;
      regioned_histos1D[iRegion].push_back( vector<MnvH1D*>(nHistos,NULL) );

      //First add angular text
      string angle_text = Form("%.0f<#delta#theta_{R}<%.0f, %.0f<#delta#theta_{P}<%.0f", minR, maxR, minP,maxP );
      regioned_histoAngles[iRegion].push_back(angle_text);

      for( int icat = 0; icat< histsToUse.size(); icat++ )
      {
        int cat = histsToUse[icat];
        MnvH2D* h2 = hists[cat][zbin];
        MnvH1D* h = (MnvH1D*) h2->ProjectionY( Form("qsq-xbin_%02d-zbin_%02d_%s",xbin,zbin, names[cat].c_str()), xbin,xbin );
        //Second add mnvh1d
        regioned_histos1D[iRegion].back()[cat]=h;
      }
    }
  }
  cout<<"Parsed region, draw region: "<<endl;

  for( int iRegion = 0; iRegion < nRegions; iRegion++ )
  {
    vector<vector<MnvH1D*>> hists = regioned_histos1D[iRegion];
    vector<string> titles = regioned_histoAngles[iRegion];
    //drawChunkForRegion( hists, titles, iRegion,  doRatio,  doBckSub,  subQELike,  doFitting );
    drawChunkForRegion( hists, titles, iRegion,  0,  0,  0,  doFitting );
    drawChunkForRegion( hists, titles, iRegion,  0,  1,  0,  doFitting );
    drawChunkForRegion( hists, titles, iRegion,  0,  1,  1,  doFitting );
    drawChunkForRegion( hists, titles, iRegion,  1,  0,  0,  doFitting );
    drawChunkForRegion( hists, titles, iRegion,  1,  1,  0,  doFitting );
    drawChunkForRegion( hists, titles, iRegion,  1,  1,  1,  doFitting );
  }

}

//region hists
void drawChunkForRegion( vector<vector<MnvH1D*>> &hists, vector<string>&titles, int iRegion, bool doRatio, bool doBckSub,  bool subQELike,bool doFitting )
{
  cout<<"enter drawChunkForRegion"<<endl;
  int nx = 2, ny = 2;
  if( iRegion == 1) { nx=3; ny=3;}
  if( iRegion == 2) { nx=4; ny=3;}
  if( iRegion == 4) { nx=4; ny=3;}

  TCanvas *c = new TCanvas("c","c",1200,800);
  c->Divide(nx,ny);

  vector<int> histosToDraw({ kData, kMC,    kMC,       kQELikeNot, kQELike_DIS, kQELike_RES, kQELike_2p2h, kQELike_QE_OTH, kQELike_QE_H, kData} );
  vector<string>      opts({"pe", "e2same","histsame", "histsame", "histlsame","histlsame","histlsame","histlsame","histlsame","pesame"});

  vector<int> qelikenot_hists({kQELikeNot_SingleNeutralPion, kQELikeNot_SingleChargedPion, kQELikeNot_MultiPion, kQELikeNot_NoPions});

  //vector<TH1D*> th1_histos;
  cout<<"looping over regions"<<endl;
  TLatex *text = new TLatex();
  text->SetTextFont(82);
  text->SetTextSize(0.05);
  text->SetTextColor( kRed );
  for( unsigned int i = 0; i< hists.size(); i++ )
  {
    c->cd(i+1);
    gPad->SetLogx();
    gPad->SetGridx();
    gPad->SetMargin(.15,.00,.15,.1);
    FormatHistos( hists[i] );
    cout<<"Formatted Histos"<<endl;

    MnvH1D* hmc = (MnvH1D*) hists[i][kMC]->Clone("hmcbase");

    MnvH1D* hqelikenot = (MnvH1D*) hists[i][kQELikeNot]->Clone("hqelikenotbase");
    MnvH1D* hqelike = (MnvH1D*) hists[i][kQELike]->Clone("hqelikebase");
    hqelike->Add( hists[i][kQELike_QE_H], -1 );


    THStack* hs = new THStack( Form("%d",i),"");
    for( auto cat : qelikenot_hists ) 
    {
      TH1D* h = (TH1D*) (hists[i][cat]->GetCVHistoWithError()).Clone( Form( "hs_%d_%s",i, names[cat].c_str() ) );
      if(doRatio) h->Divide( hmc );
      else h->Scale(1,"width");
      hs->Add(h);
    }
    if(doBckSub) 
    {
      hmc->Add( hqelikenot, -1 );
      if( subQELike ) hmc->Add( hqelike, -1 );
    }


    //hists[i][kData]->Draw("pe");
    for( unsigned int j = 0; j< histosToDraw.size(); j++ )
    {//Draw histo categories -- begin
      int J = histosToDraw[j];
      cout<<"histo to draw: "<<J<<endl;
      //TH1D *h= new TH1D( hists[i][J]->GetCVHistoWithError() );
      MnvH1D* h = (MnvH1D*) hists[i][J]->Clone( Form("%d%d", i,j) );
      //h->SetName( Form("%d%d", i,j) );
      //th1_histos.push_back(h);
      h->GetXaxis()->SetTitle("Q^{2}_{QE} (GeV^{2})");
      h->GetYaxis()->SetTitle("Events / GeV^{2}");
      h->SetMinimum(0);
      h->GetXaxis()->SetNdivisions(4);


      if(j==1)
      {
        h->SetFillStyle(1001);
        //h->SetFillColor(kRed);
        h->SetFillColorAlpha(kRed,0.35);
      }
      if( doBckSub && ( J == kData || J == kMC || J == kQELikeNot) )
      {
        if( J!=kQELike ) h->Add( hqelikenot, -1 );
      }

      if( doBckSub && subQELike )
      {
        if(J==kData || J==kMC) h->Add( hqelike,-1 );
        if ( J == kQELike_QE_OTH || J == kQELike_RES || J==kQELike_2p2h || J==kQELike_DIS) h->Reset();
      }


      if( doRatio ) h->Divide(h,hmc);

      TH1D* hd = new TH1D( h->GetCVHistoWithError() );
      hd->SetName( Form("th1d_%d%d", i,j) );
      if( doRatio ){
        hd->SetMaximum(2.5);
        hd->SetMinimum(0);
      }
      else
      {
        hd->Scale(1,"width");
        hd->SetMaximum( hd->GetMaximum()*1.5 );
      }

      hd->Draw(opts[j].c_str());
      if(!doBckSub && j == 0) hs->Draw("histsame");
    }//Draw histo categories -- end

    text->DrawLatexNDC(0.2,0.8, titles[i].c_str() );
    //delete hqelikenot;
    //delete hmc;
  }

  
  c->Print(Form("plots/nu-1d-chunkregion_%d-angles-ratio_%d-bcksub_%d-subqelike_%d-fit_%d.pdf",iRegion, doRatio, doBckSub, subQELike, doFitting));

  TLegend* leg=new TLegend(0.5, 0.6, 1, 1);
  leg->SetTextSize(1);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(52);
  leg->SetTextSize(0.03);
  leg->AddEntry(hists[0][kData], "MINERvA data", "lpe");
  leg->AddEntry(hists[0][kMC], "MINERvA Tune", "l");
  leg->AddEntry(hists[0][kQELike_QE_H],"QE-H","l");
  leg->AddEntry(hists[0][kQELike_QE_OTH],"QE-Oth","l");
  leg->AddEntry(hists[0][kQELike_RES],"Resonant","l");
  leg->AddEntry(hists[0][kQELike_DIS],"DIS","l");
  leg->AddEntry(hists[0][kQELike_2p2h],"2p2h","l");
  leg->AddEntry(hists[0][kQELikeNot],"Not QELike","l");
  if( !doBckSub )
  {
    leg->AddEntry(hists[0][kQELikeNot_SingleChargedPion],"SCP","l");
    leg->AddEntry(hists[0][kQELikeNot_SingleNeutralPion],"SNP","l");
    leg->AddEntry(hists[0][kQELikeNot_MultiPion],"MP","l");
    leg->AddEntry(hists[0][kQELikeNot_NoPions],"NP","l");
  }
  //leg->AddEntry(mc_qelike_2p2h_no_lowrec,"2p2h without fit","l");
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->cd();
  leg->Draw();
  c2->Print( Form( "Legend-bck_%d.pdf", doBckSub ) );
}



void drawChunkForOne(string savename, string chunkname, string axis,vector<TH2D*>&hists,  HDLHistos &hdlHistos, double *multipliers, bool doScale )
{
  string histsname=Form("plots/%s-%s.pdf",savename.c_str(),chunkname.c_str()) ;
  cout<<"draw "<<histsname<<endl;
  int nChunk = hists.size()-2;

  for( int i = 0; i<nChunk; i++ )
  {
    cout<<"chunk: "<<i<<endl;
    int zbin = i+1;
    double zvalue = hdlHistos.GetZaxis()->GetBinCenter(zbin);
    if(zvalue<-55 || zvalue >55) continue;
    cout<<"zvalue = "<<zvalue<<endl;
    TH2* mcErr = (TH2*) hists[zbin]->Clone("mcErr");
    TH2* mc = (TH2*) hists[zbin]->Clone("mc");

    mcErr->SetFillStyle(1001);
    mcErr->SetFillColorAlpha(kRed,0.35);
    mcErr->SetMarkerSize(0 );


    //histAndOpts.push_back( std::make_pair( &dataStat, "e" ) );
    cout<<"Creating histAndOpts"<<endl;
    vector<pair<TH2*, const char*> > histAndOpts;
    histAndOpts.push_back( std::make_pair( mcErr, "e2" ) );
    histAndOpts.push_back( std::make_pair( mc, "hist" ) );

    for( unsigned int j = 0; j< histAndOpts.size(); j++ )
    {
      if(doScale) histAndOpts[j].first->Scale(0.001,"width");
    }

    string xtitle="#delta#theta_{P}";
    string ytitle="Q^{2}_{QE} (GeV^{2})";
    GridCanvas* gc = (axis=="x")? plotXAxis1D( histAndOpts, xtitle, ytitle, multipliers ): plotYAxis1D( histAndOpts, ytitle, xtitle, multipliers );
    if(axis=="y") gc->SetLogx(true);
    else gc->SetLogx(false);
    gc->SetGridx(true);
    gc->SetYLimits(0, mc->GetMaximum()*1.5);
    if(!doScale)gc->SetYLimits(-.2,2.2);
    gc->SetYTitle("Evt Rate #times 10^{3}");
    if(!doScale) gc->SetYTitle("Ratio to MnvGENIE");
    gc->Modified();

    if( axis=="y" ) gc->SetXLimits(0,12);

    TLegend* leg=new TLegend(0.7, 0.1, 0.9, 0.3);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry( mc, "nuwro", "l");
    leg->Draw("SAME");

    string fname=Form("plots/%s-%s-%s_%05d-scale_%d.pdf",savename.c_str(),axis.c_str(),chunkname.c_str(), (int)zvalue, doScale );
    gc->Print(fname.c_str());
  }
}


vector<TH2D*> split3D( TH3D* hist )
{
  vector<TH2D*> ret; 
  //xyz ---> y vs xz
  for( int iy = 0; iy < hist->GetNbinsY()+2; iy++ )
  {
    hist->GetYaxis()->SetRange(iy,iy);
    TH2D* h = (TH2D*) hist->Project3D("zx");
    h->SetName( Form("%s_%02d", hist->GetName(), iy ) );
    ret.push_back( h );
  }
  return ret;
}


void FormatHistos( vector<vector<MnvH2D*>> h )
{
  vector<int> mycolors = getColors(2);
  int nChunks = h[0].size();

  for( int i = 0; i< nChunks; i++ )
  {
    cout<<"Formatting "<<i<<endl;
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

    vector<int> qelikenots{kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion,kQELikeNot_NoPions};
    for( unsigned int j = 0; j<qelikenots.size();j++ )
    {
      int cat = qelikenots[j];
      h[cat][i]->SetFillStyle( 1001 );
      h[cat][i]->SetFillColorAlpha( mycolors[7+j], 0.3 );
      h[cat][i]->SetLineColor(mycolors[7+j]);
      h[cat][i]->SetLineStyle( 2 );
      h[cat][i]->SetMarkerSize(0);
    }

  }
}
//===============================================================================


//======================= Region Plots ================================//

void drawRegionAnglePlots(TFile* f_region, TFile* f_orig, TFile* f_signal_weight, bool doBckFitting=true, string sample="Signal", string tag="" );
void drawRegions( vector<vector<MnvH1D*>> &hists, bool doRatio, bool doBckSub, bool subQELike, bool hasFitting, int nRegions, bool angle=false, string tag="" );
void drawRegionsSeparate( vector<vector<MnvH1D*>> &hists, bool doRatio, bool doBckSub, bool subQELike, bool hasFitting, int nRegions, bool angle=false, string tag="" );
void drawOneHist( vector<MnvH1D*> &hist, bool doRatio, bool doBckSub,  bool subQELike,bool hasFitting,string tag, string print="" );

void drawRegionAnglePlots( TFile* f_region, TFile* f_orig,  TFile* f_signal_weight , bool doBckFitting, string sample , string tag)
{
  cout<<"drawRegionAnglePlots"<<endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  //three files 2-track, 2-track bkg_fitted
  //CV
  //cout<<"read "<<fname<<endl;
  TFile* f1 = f_region;
  TFile* f2 = f_signal_weight;
  if(!f1) cout<<"f1 i not read"<<endl;

  TVector2 *pot = (TVector2*)f1->Get("pot");

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
  //HyperDimLinearizer* hdl = GetHDL(xbins , ybins, zbins);

  //vector<MnvH1D*> h_mc_weights = GetFitsHisto( f2 );
  vector<MnvH1D*> fits = GetFitsHisto( f2 );
  //for( int i = 0; i< 24;i++ ) h_mc_weights.push_back( (MnvH1D*) f2->Get( Form( "h_weights_q2qe_mc_%02d", i ) ) );



  // Get CV histos
  //
  //
  vector<MnvH1D*> h_qsq_signal;
  for( unsigned int j = 0; j< nHistos;j ++ )
  {
    string hname = Form("h_q2qe_%s", names[j].c_str()  )  ;
    MnvH1D* h = (MnvH1D*) f_orig->Get(hname.c_str() );
    double norm = (j>0)? pot_norm : 1 ;
    if( j > 0 ) h->Scale(norm);
    h->SetMinimum(0);
    h_qsq_signal.push_back( h );



  }
  h_qsq_signal[0]->AddMissingErrorBandsAndFillWithCV( *h_qsq_signal[1] );



  cout<<"draw regions"<<endl;
  vector<vector<MnvH1D*>> hvec_qsqqe;
  int nRegions = 6;
  for( int i = 0; i< nRegions+1; i ++ )
  {
    int iR = i;
    if (i == nRegions ) iR = 99;
    cout<<"region "<<iR<<endl;
    vector<MnvH1D*> hists;
    for( unsigned int j = 0; j< nHistos;j ++ )
    {
      string hname = Form("h_q2qe_region_%02d_%s",iR, names[j].c_str()  )  ;
      if( tag!="" ) hname = Form("h_q2qe_%s_region_%02d_%s",tag.c_str(), iR, names[j].c_str()  )  ;
      cout<<hname<<endl;
      MnvH1D* h = (MnvH1D*) f1->Get(hname.c_str() );
      if( j == 0 ) h->AddMissingErrorBandsAndFillWithCV( *h_qsq_signal[1] );
      //MnvH1D* h = (MnvH1D*) f1->Get( "h_q2qe_region_00_data" );
      cout<<hname<<endl;
      double norm = (j>0)? pot_norm : 1 ;
      if( j > 0 ) h->Scale(norm);
      h->SetMinimum(0);
      hists.push_back( h );
      cout<<i<<j<<endl;
    }
    hvec_qsqqe.push_back( hists );
  }




  if( doBckFitting )
  {
    cout<<"doBckFitting"<<endl;
    for( unsigned int i = 0; i< hvec_qsqqe.size(); i++ )
    {
      int sideband = 0;
      if(sample=="BlobSideBand") sideband=1;
      else if(sample=="MichelSideBand") sideband=2;
      else if(sample=="MicBlobSideBand") sideband=3;

      int I = i+sideband*(nRegions+1);

      DoFit( &(hvec_qsqqe[i][0]), fits );
    }
    DoFit( &(h_qsq_signal[0]), fits );

    MnvH1D* hdata = (MnvH1D*) hvec_qsqqe[0][kData]->Clone("hdata");
    MnvH1D* hbck = (MnvH1D*)hvec_qsqqe[0][kMC]->Clone("hbck");
    hbck->Add( hvec_qsqqe[0][kQELike_QE_H], -1 );
    hdata->AddMissingErrorBandsAndFillWithCV( *hbck );
    hdata->Add( hbck, -1 );
    drawSystematics( putils, hdata,  pot_data, pot_mc, "region_00");

  }

  cout<<"parsed region "<<hvec_qsqqe[0][kData]->GetEntries()<<endl;
  drawRegions( hvec_qsqqe, false, false, false, doBckFitting,   nRegions, false, tag );
  drawRegions( hvec_qsqqe, false, true,  false, doBckFitting,   nRegions, false, tag );
  drawRegions( hvec_qsqqe, true, false,  false, doBckFitting,   nRegions, false, tag );
  drawRegions( hvec_qsqqe, true, true,   false, doBckFitting,   nRegions, false, tag );

  drawRegions( hvec_qsqqe, false, true,  true, doBckFitting,   nRegions, false, tag );
  drawRegions( hvec_qsqqe, true,  true,  true, doBckFitting,   nRegions, false, tag );

  drawRegionsSeparate( hvec_qsqqe, false, false, false, doBckFitting,   nRegions, false, tag );
  drawRegionsSeparate( hvec_qsqqe, false, true,  false, doBckFitting,   nRegions, false, tag );
  drawRegionsSeparate( hvec_qsqqe, true, false,  false, doBckFitting,   nRegions, false, tag );
  drawRegionsSeparate( hvec_qsqqe, true, true,   false, doBckFitting,   nRegions, false, tag );

  drawRegionsSeparate( hvec_qsqqe, false, true,  true, doBckFitting,   nRegions, false, tag );
  drawRegionsSeparate( hvec_qsqqe, true,  true,  true, doBckFitting,   nRegions, false, tag );

  drawOneHist( h_qsq_signal, false, false, false, doBckFitting, "all", "Full Region" );
  drawOneHist( h_qsq_signal, false, true,  false, doBckFitting, "all", "Full Region" );
  drawOneHist( h_qsq_signal, true, false,  false, doBckFitting, "all", "Full Region" );
  drawOneHist( h_qsq_signal, true, true,   false, doBckFitting, "all", "Full Region" );

  drawOneHist( h_qsq_signal, false, true,  true, doBckFitting, "all", "Full Region"  );
  drawOneHist( h_qsq_signal, true,  true,  true, doBckFitting, "all", "Full Region"  );

  cout<<"End drawRegion"<<endl;


  //f1->Close();
}

template<class T>
void FormatHistos( vector<T*> &h )
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

  vector<int> qelikenots{kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion,kQELikeNot_NoPions};
  for( unsigned int i = 0; i<qelikenots.size();i++ )
  {
    int cat = qelikenots[i];
    h[cat]->SetFillStyle( 1001 );
    h[cat]->SetFillColorAlpha( mycolors[7+i], 0.3 );
    h[cat]->SetLineColor(mycolors[7+i]);
    h[cat]->SetLineStyle( 2 );
    h[cat]->SetLineWidth( 1 );
    h[cat]->SetMarkerSize(0);
  }
}

template void FormatHistos( vector<MnvH1D*> &h );
template void FormatHistos( vector<MnvH2D*> &h );

void drawRegions( vector<vector<MnvH1D*>> &hists, bool doRatio, bool doBckSub,  bool subQELike,bool hasFitting, int nRegions, bool angles, string tag )
{
  cout<<"enter drawRegions"<<endl;

  TCanvas *c = new TCanvas("c","c",1200,800);
  int nCanvas = 8;
  if(!angles)c->Divide(4,2);
  else {c->Divide(3,2); nCanvas=6;}
  c->cd();
  vector<int> histosToDraw({ kData, kMC,    kMC,       kQELikeNot, kQELike_DIS, kQELike_RES, kQELike_2p2h, kQELike_QE_OTH, kQELike_QE_H, kData} );
  vector<string>      opts({"pe", "e2same","histsame", "histsame", "histlsame","histlsame","histlsame","histlsame","histlsame","pesame"});

  vector<int> qelikenot_hists({kQELikeNot_SingleNeutralPion, kQELikeNot_SingleChargedPion, kQELikeNot_MultiPion, kQELikeNot_NoPions});

  //vector<TH1D*> th1_histos;
  cout<<"looping over regions"<<endl;
  TLatex *text = new TLatex();
  text->SetTextFont(92);
  text->SetTextSize(0.1);
  text->SetTextColor( kBlue );
  for( unsigned int i = 0; i< hists.size(); i++ )
  {
    c->cd(i+1);
    gPad->SetLogx();
    gPad->SetGridx();
    gPad->SetMargin(.15,.00,.15,.1);
    FormatHistos( hists[i] );
    cout<<"Formatted Histos"<<endl;

    MnvH1D* hmc = (MnvH1D*) hists[i][kMC]->Clone("hmcbase");

    //TH1D* hmc = new TH1D( hists[i][kMC]->GetCVHistoWithError() );
    //hmc->SetName("hmcbase");
    MnvH1D* hqelikenot = (MnvH1D*) hists[i][kQELikeNot]->Clone("hqelikenotbase");
    MnvH1D* hqelike = (MnvH1D*) hists[i][kQELike]->Clone("hqelikebase");
    hqelike->Add( hists[i][kQELike_QE_H], -1 );
    hqelike->SetMinimum(0);
    //TH1D* hqelikenot = new TH1D( hists[i][kQELikeNot]->GetCVHistoWithError());
    //hqelikenot->SetName("hqelikenotbase");


    THStack* hs = new THStack( Form("%d",i),"");
    for( auto cat : qelikenot_hists ) 
    {
      TH1D* h = (TH1D*) (hists[i][cat]->GetCVHistoWithError()).Clone( Form( "hs_%d_%s",i, names[cat].c_str() ) );
      if(doRatio) h->Divide( hmc );
      else h->Scale(1,"width");
      hs->Add(h);
    }
    if(doBckSub) 
    {
      hmc->Add( hqelikenot, -1 );
      if( subQELike ) hmc->Add( hqelike, -1 );
      hmc->SetMinimum(0);
    }


    //hists[i][kData]->Draw("pe");
    for( unsigned int j = 0; j< histosToDraw.size(); j++ )
    {
      int J = histosToDraw[j];
      cout<<"histo to draw: "<<J<<endl;
      //TH1D *h= new TH1D( hists[i][J]->GetCVHistoWithError() );
      MnvH1D* h = (MnvH1D*) hists[i][J]->Clone( Form("%d%d", i,j) );
      //h->SetName( Form("%d%d", i,j) );
      //th1_histos.push_back(h);
      h->GetXaxis()->SetTitle("Q^{2}_{QE} (GeV^{2})");
      h->GetYaxis()->SetTitle("Events / GeV^{2}");
      h->SetMinimum(0);
      h->GetXaxis()->SetNdivisions(4);


      if(j==1)
      {
        h->SetFillStyle(1001);
        //h->SetFillColor(kRed);
        h->SetFillColorAlpha(kRed,0.35);
      }
      if( doBckSub && ( J == kData || J == kMC || J == kQELikeNot) )
      {
        if( J!=kQELike ) h->Add( hqelikenot, -1 );
        h->SetMinimum(0);
      }

      if( doBckSub && subQELike )
      {
        if(J==kData || J==kMC) h->Add( hqelike,-1 );
        if ( J == kQELike_QE_OTH || J == kQELike_RES || J==kQELike_2p2h || J==kQELike_DIS) h->Reset();
        h->SetMinimum(0);
      }


      if( doRatio ) h->Divide(h,hmc);

      TH1D* hd = new TH1D( h->GetCVHistoWithError() );
      hd->SetName( Form("th1d_%d%d", i,j) );
      if( doRatio ){
        hd->SetMaximum(2.5);
        hd->SetMinimum(0);
      }
      else
      {
        hd->Scale(1,"width");
        hd->SetMaximum( hd->GetMaximum()*1.5 );
      }

      hd->Draw(opts[j].c_str());
      if(!doBckSub && j == 0) hs->Draw("histsame");

    }
    int iR = i;
    if(i==nRegions) iR = 99;
    if(!angles)text->DrawLatexNDC(0.2,0.8, Form("Region %d", iR));
    else text->DrawLatexNDC(0.2,0.8, Form("#delta#theta_{R/P} < %d^{o}", i+1));
    //delete hqelikenot;
    //delete hmc;
  }
  c->cd(nCanvas);
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
  if( !doBckSub )
  {
    leg->AddEntry(hists[0][kQELikeNot_SingleChargedPion],"SCP","l");
    leg->AddEntry(hists[0][kQELikeNot_SingleNeutralPion],"SNP","l");
    leg->AddEntry(hists[0][kQELikeNot_MultiPion],"MP","l");
    leg->AddEntry(hists[0][kQELikeNot_NoPions],"NP","l");


  }
  //leg->AddEntry(mc_qelike_2p2h_no_lowrec,"2p2h without fit","l");
  leg->Draw("SAME");

  
  if(!angles) //c->Print(Form("plots/nu-1d-regions-ratio_%d-bcksub_%d-subqelike_%d-fit_%d.pdf", doRatio, doBckSub, subQELike, hasFitting));
  {
    TString savename = Form("plots/nu-1d-regions-ratio_%d-bcksub_%d-subqelike_%d-fit_%d.pdf", doRatio, doBckSub, subQELike, hasFitting);
    if( tag!="") savename = Form("plots/nu-1d-%s-regions-ratio_%d-bcksub_%d-subqelike_%d-fit_%d.pdf", tag.c_str(), doRatio, doBckSub, subQELike, hasFitting);
    c->Print( savename );
  }
  else c->Print(Form("plots/nu-1d-angles-ratio_%d-bcksub_%d-subqelike_%d-fit_%d.pdf", doRatio, doBckSub, subQELike, hasFitting));
  delete c;
}

void drawRegionsSeparate( vector<vector<MnvH1D*>> &hists, bool doRatio, bool doBckSub,  bool subQELike,bool hasFitting, int nRegions, bool angles, string tag )
{
  cout<<"enter drawRegionsSeparate"<<endl;
  gStyle->SetEndErrorSize(4);
  vector<int> histosToDraw({ kData, kMC,    kMC,       kQELikeNot, kQELike_DIS, kQELike_RES, kQELike_2p2h, kQELike_QE_OTH, kQELike_QE_H, kData} );
  vector<string>      opts({"pe1", "e2same","histsame", "histsame", "histlsame","histlsame","histlsame","histlsame","histlsame","pe1same"});

  vector<int> qelikenot_hists({kQELikeNot_SingleNeutralPion, kQELikeNot_SingleChargedPion, kQELikeNot_MultiPion, kQELikeNot_NoPions});

  //vector<TH1D*> th1_histos;
  cout<<"looping over regions"<<endl;
  TLatex *text = new TLatex();
  text->SetTextFont(92);
  text->SetTextSize(0.1);
  text->SetTextColor( kBlue );

  TLegend* leg=new TLegend(0.7, 0.6, 0.9, 0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(52);
  leg->SetTextSize(0.02);
  leg->AddEntry(hists[0][kData], "MINERvA data", "lpe");
  leg->AddEntry(hists[0][kMC], "MINERvA Tune", "l");
  leg->AddEntry(hists[0][kQELike_QE_H],"QE-H","l");
  leg->AddEntry(hists[0][kQELike_QE_OTH],"QE-Oth","l");
  leg->AddEntry(hists[0][kQELike_RES],"Resonant","l");
  leg->AddEntry(hists[0][kQELike_DIS],"DIS","l");
  leg->AddEntry(hists[0][kQELike_2p2h],"2p2h","l");
  leg->AddEntry(hists[0][kQELikeNot],"Not QELike","l");

  if( !doBckSub )
  {
    leg->AddEntry(hists[0][kQELikeNot_SingleChargedPion],"SCP","l");
    leg->AddEntry(hists[0][kQELikeNot_SingleNeutralPion],"SNP","l");
    leg->AddEntry(hists[0][kQELikeNot_MultiPion],"MP","l");
    leg->AddEntry(hists[0][kQELikeNot_NoPions],"NP","l");
  }




  for( unsigned int i = 0; i< hists.size(); i++ )
  {
    TCanvas *c = new TCanvas("c","c",1200,800);
    gPad->SetLogx();
    gPad->SetGridx();
    gPad->SetMargin(.15,.15,.15,.1);
    FormatHistos( hists[i] );
    cout<<"Formatted Histos"<<endl;

    MnvH1D* hmc = (MnvH1D*) hists[i][kMC]->Clone("hmcbase");

    //TH1D* hmc = new TH1D( hists[i][kMC]->GetCVHistoWithError() );
    //hmc->SetName("hmcbase");
    MnvH1D* hqelikenot = (MnvH1D*) hists[i][kQELikeNot]->Clone("hqelikenotbase");
    MnvH1D* hqelike = (MnvH1D*) hists[i][kQELike]->Clone("hqelikebase");
    hqelike->Add( hists[i][kQELike_QE_H], -1 );
    //TH1D* hqelikenot = new TH1D( hists[i][kQELikeNot]->GetCVHistoWithError());
    //hqelikenot->SetName("hqelikenotbase");


    THStack* hs = new THStack( Form("%d",i),"");
    for( auto cat : qelikenot_hists ) 
    {
      TH1D* h = (TH1D*) (hists[i][cat]->GetCVHistoWithError()).Clone( Form( "hs_%d_%s",i, names[cat].c_str() ) );
      if(doRatio) h->Divide( hmc );
      else h->Scale(1,"width");
      hs->Add(h);
    }
    if(doBckSub) 
    {
      hmc->Add( hqelikenot, -1 );
      if( subQELike ) hmc->Add( hqelike, -1 );
      hmc->SetMinimum(0);
    }


    //hists[i][kData]->Draw("pe");
    for( unsigned int j = 0; j< histosToDraw.size(); j++ )
    {
      int J = histosToDraw[j];
      cout<<"histo to draw: "<<J<<endl;
      //TH1D *h= new TH1D( hists[i][J]->GetCVHistoWithError() );
      MnvH1D* h = (MnvH1D*) hists[i][J]->Clone( Form("%d%d", i,j) );
      //h->SetName( Form("%d%d", i,j) );
      //th1_histos.push_back(h);
      h->GetXaxis()->SetTitle("Q^{2}_{QE} (GeV^{2})");
      h->GetYaxis()->SetTitle("Events / GeV^{2}");
      h->SetMinimum(0);
      h->GetXaxis()->SetNdivisions(4);


      if(j==1)
      {
        h->SetFillStyle(1001);
        //h->SetFillColor(kRed);
        h->SetFillColorAlpha(kRed,0.35);
      }
      if( doBckSub && ( J == kData || J == kMC || J == kQELikeNot) )
      {
        if( J!=kQELike ) h->Add( hqelikenot, -1 );
        h->SetMinimum(0);
      }

      if( doBckSub && subQELike )
      {
        if(J==kData || J==kMC) h->Add( hqelike,-1 );
        if ( J == kQELike_QE_OTH || J == kQELike_RES || J==kQELike_2p2h || J==kQELike_DIS) h->Reset();
        h->SetMinimum(0);
      }


      if( doRatio ) h->Divide(h,hmc);

      TH1D* hd = new TH1D( h->GetCVHistoWithError() );
      TH1D* hdstat = new TH1D( h->GetCVHistoWithStatError() );
      hd->SetName( Form("th1d_%d%d", i,j) );
      hdstat->SetName( Form("th1dstat_%d%d", i,j) );
      if( doRatio ){
        hd->SetMaximum(2.5);
        hd->SetMinimum(0);
      }
      else
      {
        hd->Scale(1,"width");
        hd->SetMaximum( hd->GetMaximum()*1.5 );
        hdstat->Scale(1,"width");
        hdstat->SetMaximum( hd->GetMaximum()*1.5 );
      }

      hd->SetMinimum(0);
      hd->Draw(opts[j].c_str());
      if(J == kData ) 
      {
        //hdstat->SetLineColor(kRed);
        hdstat->Draw("pe1same");
      }
      if(!doBckSub && j == 0) hs->Draw("histsame");

    }
    int iR = i;
    if(i==nRegions) iR = 99;
    if(!angles) text->DrawLatexNDC(0.2,0.8, Form("Region %d", iR));
    else text->DrawLatexNDC(0.2,0.8, Form("#delta#theta_{R/P} < %d^{o}", i+1));
    leg->Draw("same");
    //delete hqelikenot;
    //delete hmc
    if(!angles) 
    {
      TString savename = Form("plots/nu-1d-regions-ratio_%d-bcksub_%d-subqelike_%d-fit_%d-region_%02d.pdf", doRatio, doBckSub, subQELike, hasFitting,iR);
      if( tag!="") savename = Form("plots/nu-1d-%s-regions-ratio_%d-bcksub_%d-subqelike_%d-fit_%d-region_%02d.pdf", tag.c_str(), doRatio, doBckSub, subQELike, hasFitting,iR);
      c->Print( savename );
    }
    else c->Print(Form("plots/nu-1d-angles-ratio_%d-bcksub_%d-subqelike_%d-fit_%d-angle_%02d.pdf", doRatio, doBckSub, subQELike, hasFitting,iR));
    delete c;
  }
}

void drawOneHist( vector<MnvH1D*> &hist, bool doRatio, bool doBckSub,  bool subQELike,bool hasFitting, string tag,string print)
{

  vector<int> histosToDraw({ kData, kMC,    kMC,       kQELikeNot, kQELike_DIS, kQELike_RES, kQELike_2p2h, kQELike_QE_OTH, kQELike_QE_H, kData} );
  vector<string>      opts({"pe", "e2same","histsame", "histsame", "histlsame","histlsame","histlsame","histlsame","histlsame","pesame"});

  vector<int> qelikenot_hists({kQELikeNot_SingleNeutralPion, kQELikeNot_SingleChargedPion, kQELikeNot_MultiPion, kQELikeNot_NoPions});

  //vector<TH1D*> th1_histos;
  cout<<"looping over regions"<<endl;
  TLatex *text = new TLatex();
  text->SetTextFont(92);
  text->SetTextSize(0.1);
  text->SetTextColor( kBlue );

  TLegend* leg=new TLegend(0.7, 0.6, 0.9, 0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(52);
  leg->SetTextSize(0.02);
  leg->AddEntry(hist[kData], "MINERvA data", "lpe");
  leg->AddEntry(hist[kMC], "MINERvA Tune", "l");
  leg->AddEntry(hist[kQELike_QE_H],"QE-H","l");
  leg->AddEntry(hist[kQELike_QE_OTH],"QE-Oth","l");
  leg->AddEntry(hist[kQELike_RES],"Resonant","l");
  leg->AddEntry(hist[kQELike_DIS],"DIS","l");
  leg->AddEntry(hist[kQELike_2p2h],"2p2h","l");
  leg->AddEntry(hist[kQELikeNot],"Not QELike","l");

  if( !doBckSub )
  {
    leg->AddEntry(hist[kQELikeNot_SingleChargedPion],"SCP","l");
    leg->AddEntry(hist[kQELikeNot_SingleNeutralPion],"SNP","l");
    leg->AddEntry(hist[kQELikeNot_MultiPion],"MP","l");
    leg->AddEntry(hist[kQELikeNot_NoPions],"NP","l");
  }




  TCanvas *c = new TCanvas("c","c",1200,800);
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetMargin(.15,.15,.15,.1);
  FormatHistos( hist );
  cout<<"Formatted Histos"<<endl;

  MnvH1D* hmc = (MnvH1D*) hist[kMC]->Clone("hmcbase");

  MnvH1D* hqelikenot = (MnvH1D*) hist[kQELikeNot]->Clone("hqelikenotbase");
  MnvH1D* hqelike = (MnvH1D*) hist[kQELike]->Clone("hqelikebase");
  hqelike->Add( hist[kQELike_QE_H], -1 );


  THStack* hs = new THStack( Form("%d",0),"");
  for( auto cat : qelikenot_hists ) 
  {
    TH1D* h = (TH1D*) (hist[cat]->GetCVHistoWithError()).Clone( Form( "hs_%d_%s",0, names[cat].c_str() ) );
    if(doRatio) h->Divide( hmc );
    else h->Scale(1,"width");
    hs->Add(h);
  }
  if(doBckSub) 
  {
    hmc->Add( hqelikenot, -1 );
    if( subQELike ) hmc->Add( hqelike, -1 );
    hmc->SetMinimum(0);
  }


  //hists[i][kData]->Draw("pe");
  for( unsigned int j = 0; j< histosToDraw.size(); j++ )
  {
    int J = histosToDraw[j];
    cout<<"histo to draw: "<<J<<endl;
    MnvH1D* h = (MnvH1D*) hist[J]->Clone( Form("%d%d", 0,j) );
    h->GetXaxis()->SetTitle("Q^{2}_{QE} (GeV^{2})");
    h->GetYaxis()->SetTitle("Events / GeV^{2}");
    h->SetMinimum(0);
    h->GetXaxis()->SetNdivisions(4);


    if(j==1)
    {
      h->SetFillStyle(1001);
      //h->SetFillColor(kRed);
      h->SetFillColorAlpha(kRed,0.35);
    }
    if( doBckSub && ( J == kData || J == kMC || J == kQELikeNot) )
    {
      if( J!=kQELike ) h->Add( hqelikenot, -1 );
      h->SetMinimum(0);
    }

    if( doBckSub && subQELike )
    {
      if(J==kData || J==kMC) h->Add( hqelike,-1 );
      if ( J == kQELike_QE_OTH || J == kQELike_RES || J==kQELike_2p2h || J==kQELike_DIS) h->Reset();
      h->SetMinimum(0);
    }


    if( doRatio ) h->Divide(h,hmc);

    TH1D* hd = new TH1D( h->GetCVHistoWithError() );
    hd->SetName( Form("th1d_%d%d", 0,j) );
    if( doRatio ){
      hd->SetMaximum(2.5);
      hd->SetMinimum(0);
    }
    else
    {
      hd->Scale(1,"width");
      hd->SetMaximum( hd->GetMaximum()*1.5 );
    }

    hd->SetMinimum(0);
    hd->Draw(opts[j].c_str());
    if(!doBckSub && j == 0) hs->Draw("histsame");

  }
  leg->Draw("same");
  text->DrawLatexNDC(0.2,0.8, print.c_str());
  TString savename = Form("plots/nu-1d-regions-ratio_%d-bcksub_%d-subqelike_%d-fit_%d-region_%s.pdf", doRatio, doBckSub, subQELike, hasFitting, tag.c_str());
  c->Print( savename );
  delete c;
}

void drawSystematics( CCQENuPlotUtils*utils, MnvH1D* h,  double pot_data, double pot_mc, string name)
{

h->GetXaxis()->SetRangeUser(0.05,5);
MnvPlotter *plotter = new MnvPlotter;
//plotter->ApplyStyle(PlotUtils::kCompactStyle);
plotter->width_xspace_per_letter=0.2;
//plotter->legend_offset_y=-0.01;
//plotter->axis_maximum=1;
TCanvas* c = new TCanvas("c","c", 1200,800);
c->SetLogx();
plotter->SetRedHeatPalette();
plotter->error_summary_group_map = utils->getSystematicGroupMap();
plotter->error_color_map         = utils->getSystematicGroupMapColors();
TGaxis::SetMaxDigits(5);
  plotter->axis_maximum_group = 0.2; 
  plotter->headroom = 1.0;
  plotter->legend_n_columns = 2;
  plotter->legend_text_size = 0.02;
  plotter->AddPlotLabel( Form("%s #bullet%s", "#Q^{2}" ,  "Errors"), 0.33, 0.94, 0.025 );

  vector<string> sysGroupNames({"Muon Reconstruction", "Recoil Reconstruction", "Low Recoil Fit", "XSection Models", "FSI Models", "Others", "Flux"});
  vector<string> sysGroupPrintNames({"Muon_Reconstruction", "Recoil_Reconstruction", "Low_Recoil_Fit", "XSection_Models", "FSI_Models", "Others", "Flux"});
  for( int i = 0; i < sysGroupPrintNames.size(); i++ )
  {
    cout<<sysGroupPrintNames[i]<<endl;
    utils->writeNorm( false, pot_data, pot_mc, true ); 
    plotter->DrawErrorSummary(h, "TL", false, true, 0.00001, false, sysGroupNames[i]);
    c->Print(Form("plots/sys_%s_%s.pdf",name.c_str(), sysGroupPrintNames[i].c_str() ) );
  }

    plotter->AddPlotLabel( Form("%s #bullet%s", "#Q^{2}" ,  "Errors"), 0.33, 0.94, 0.025 );
    utils->writeNorm( false, pot_data, pot_mc, true ); 
    plotter->DrawErrorSummary(h, "TL", true, true, 0.00001, false, "");
    c->Print( Form("plots/sys_%s_all.pdf",name.c_str()) );
    delete plotter;
    delete c;
}

//======================= End Region Plots ================================//
//======================= Draw Fits =======================================//

void DrawMCWithErrorBand( TCanvas *c1, MnvH1D* h1d, string title,string xtitle, string ytitle, string savename, bool isRatio=true)
{
  c1->cd();

  TH1D h = h1d->GetCVHistoWithError();
  TH1D herr = h1d->GetCVHistoWithError();
  if (isRatio)
  {
    h.GetYaxis()->SetRangeUser(0,4);
    herr.GetYaxis()->SetRangeUser(0,4);
  }
  h.SetLineColor(kRed);
  h.SetLineWidth(4);
  h.SetMarkerSize(0);

  herr.SetLineWidth(1);
  herr.SetLineColor(kRed);
  herr.SetMarkerSize(0);
  herr.SetFillStyle(1001);
  herr.SetFillColorAlpha(kRed, 0.3);

  h.GetXaxis()->SetTitle(xtitle.c_str() );
  h.GetXaxis()->SetTitleSize(0.05);
  h.GetYaxis()->SetTitle(ytitle.c_str() );
  h.GetYaxis()->SetTitleSize(0.05);

  h.Draw("hist");
  herr.Draw("e2sames");

  TLatex latex;
  latex.SetTextAlign(22);
  latex.SetTextFont(52);
  latex.SetTextSize(0.06);
  latex.DrawLatexNDC(0.5,0.9,title.c_str() );
  c1->SetLogx();
  c1->SetGrid();
  c1->Print( savename.c_str() );

}
void DrawFits( TFile *f )
{
  cout<<"Entering Draw Fits"<<endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  //gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  vector<MnvH1D*> fits = GetFitsHisto( f );
  cout<<"Fit Histo Name: "<<endl;
  for( UInt_t i = 0; i<fits.size(); i++ ) cout<<fits[i]->GetName()<<endl;
  //MnvH1D* h_qe = (MnvH1D*) f->Get("hs_weights_yvar_bgType_qe_oth");
  //MnvH1D* h_2p2h = (MnvH1D*) f->Get("hs_weights_yvar_bgType_2p2h");
  //MnvH1D* h_qelikenot_scp  = (MnvH1D*) f->Get("hs_weights_yvar_bgType_qelikenot_scp");
  //MnvH1D* h_qelikenot_snp  = (MnvH1D*) f->Get("hs_weights_yvar_bgType_qelikenot_snp");
  //MnvH1D* h_qelikenot_mp = (MnvH1D*) f->Get("hs_weights_yvar_bgType_qelikenot_mp");

  TCanvas *c1 = new TCanvas("c1","c1");
  for( UInt_t i = 0; i<categories_to_fit.size(); i++ ) DrawMCWithErrorBand ( c1, fits[i], categories_to_fit_title_names[i], "Q^{2}_{QE} (GeV^{2})", "weight", Form("plots/rw_%s.pdf", categories_to_fit_names[i].c_str()) );
  //DrawMCWithErrorBand ( c1, h_2p2h, "QELike && 2p2h",  "Q^{2}_{QE} (GeV^{2})", "weight", "plots/rw_qelike_2p2h.pdf" );
  //DrawMCWithErrorBand ( c1, h_qelikenot_scp, "QELikeNot && Single Charged Pion",  "Q^{2}_{QE} (GeV^{2})", "weight", "plots/rw_qelikenot_scp.pdf" );
  //DrawMCWithErrorBand ( c1, h_qelikenot_snp, "QELikeNot && Single Neutral Pion",  "Q^{2}_{QE} (GeV^{2})", "weight", "plots/rw_qelikenot_snp.pdf" );
  //DrawMCWithErrorBand ( c1, h_qelikenot_mp,  "QELikeNot && Multi-Pions",  "Q^{2}_{QE} (GeV^{2})", "weight", "plots/rw_qelikenot_mp.pdf" );

 // f->Close();
}

//___________________________________________________________
//___________________________________________________________
//Recoil Plots

vector<TH1D*> convertMnvH1D( vector<MnvH1D*> &hists )
{
  vector<TH1D*> ret(nHistos,NULL );
  for( UInt_t i = 0; i< nHistos; i++ )
  {
    if (!hists[i]) continue;
    ret[i] = new TH1D( hists[i]->GetCVHistoWithError() );
    ret[i]->SetName( Form("h_th1d_%s", names[i].c_str() ) );
  }
  return ret;

}


void drawRecoil( MnvH2D**h , bool doFitting )
{

  cout<<"entering drawRecoil"<<endl;
  //define qsq regions
  vector<double> bin_edges({0,0.025,0.05,0.1,0.3,1,1.5} );
  vector<int>bin_sum_low({0});
  vector<int>bin_sum_high;
  for( int i = 1; i< bin_edges.size(); i++)
  {
    int bin = h[0]->GetYaxis()->FindBin( bin_edges[i] );
    bin_sum_high.push_back(bin);
    bin_sum_low.push_back(bin+1);
  }
  bin_sum_low.pop_back();
  int nRegions = bin_sum_low.size();
  //int bin_low  = h[0]->GetYaxis()->FindBin(0);
  //int bin_mid1 = h[0]->GetYaxis()->FindBin(0.3);
  //int bin_mid2 = h[0]->GetYaxis()->FindBin(1);
  //int bin_high = h[0]->GetYaxis()->FindBin(2);
  //cout<<"got bins: "<< bin_low << "\t" << bin_mid1 <<"\t"<< bin_mid2 <<"\t"<< bin_high<<endl;

  //define histonames
  vector<int> histosToDraw({ kData, kMC,    kMC, kQELikeNot, kQELike_DIS, kQELike_RES, kQELike_2p2h, kQELike_QE_OTH, kQELike_QE_H, kData} );
  vector<string>      opts({"pe", "e2same","histsame", "histsame", "histlsame","histlsame","histlsame","histlsame","histlsame","pesame"});

  //Define Canvas
  TCanvas *c = new TCanvas("c","c",1200,800);
  c->Divide(3,2);
  
  //Define and Formatting Histo1D
  vector< vector<MnvH1D*> > h1_all( nRegions, vector<MnvH1D*>(nHistos,NULL) );
  for( auto cat : histosUsed ) 
  {
    cout<<"Projecting: "<<names[cat]<<endl;
    cout<<"Succeeded?"<<endl;
    for( int i  = 0; i< nRegions; i++ )
    {
      cout<<"Bin: "<<bin_sum_low[i]<<"\t"<<bin_sum_high[i]<<endl;
      h1_all[i][cat] = (MnvH1D*) h[cat]->ProjectionX( Form( "%s_region_%02d",  h[cat]->GetName(),i ), bin_sum_low[i], bin_sum_high[i] );
      cout<<h1_all[i][cat]->GetName()<<"\t";
    }
    cout<<endl;
  }
  for( int i = 0; i< nRegions; i++ ) FormatHistos( h1_all[i] );

  cout<<"Formated Histos"<<endl;

  //vector<TH1D*> th1d_1 = convertMnvH1D( h1_1 );
  //vector<TH1D*> th1d_2 = convertMnvH1D( h1_2 );
  //vector<TH1D*> th1d_3 = convertMnvH1D( h1_3 );

  //Apply draw options

  //TLatex
  TLatex *latex = new TLatex();
  latex->SetTextColor(kBlue);
  latex->SetTextFont(82);
  latex->SetTextSize(0.08);
  double latex_x = .3;
  double latex_y = .8;

  for( int i = 0; i< nRegions; i++ )
  {
    c->cd(i+1);
    for( int j = 0; j<histosToDraw.size(); j++ )
    {
      int cat = histosToDraw[j];
      TH1D* htemp = new TH1D( h1_all[i][cat]->GetCVHistoWithError() );
      if( j == 0 ) htemp->GetYaxis()->SetRangeUser(0, htemp->GetMaximum()*1.3) ;
      if( j == 1 ) 
      {
        htemp->SetMarkerSize(0);
        htemp->SetFillColorAlpha(kRed,0.3);
        htemp->SetLineWidth(0);

      }
      htemp->Draw( opts[j].c_str() );
    }
    latex->DrawLatexNDC(latex_x, latex_y, Form("%.2f<Q^{2}_{QE}<%.2f",bin_edges[i], bin_edges[i+1]));
    gPad->SetMargin(.1, 0.1, 0.1, 0.1 );

  }

  c->cd(6);
  gPad->SetMargin(0,.0,.0,0);
  TLegend *leg = new TLegend(.6,.6,.99,.99);

  vector<int> drawnHisto({kData, kQELike_QE_H, kQELike_QE_OTH, kQELike_2p2h, kQELike_RES, kQELike_DIS, kQELikeNot});
  for( auto cat : drawnHisto ) leg->AddEntry( h1_all[0][cat], names[cat].c_str() );
  leg->Draw();
  c->Update();
  c->Print(Form("plots/recoil-doFitting_%d.pdf", doFitting ));

}

void PlotRecoil( TFile *f_orig, TFile *f_signal_weighted )
{
  cout<<"drawRecoil"<<endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  TFile *f1 = f_orig;
  TFile *f2 = f_signal_weighted;

  //------- POT ----------------
  TVector2 *pot = (TVector2*)f1->Get("pot");
  double pot_data = pot->X();
  double pot_mc = pot->Y();
  double pot_norm = pot_data/pot_mc;

  CCQENuPlotUtils *putils = new CCQENuPlotUtils( fluxHistoExists );
  CCQENuUtils *utils = new CCQENuUtils( false, fluxHistoExists );

  MnvH2D *h_recoil_inc_q2qe[nHistos];
  utils->bookHistos( f1, h_recoil_inc_q2qe, "h_recoil_inc_q2qe");
  putils->scaleMCHistos( h_recoil_inc_q2qe, pot_norm );
  cout<<"Scaled MC Histos"<<endl;

  vector<MnvH2D*> fits = GetFitsHisto( f2, "recoil" );
  cout<<"Got Fits? "<<endl;
  cout<<fits[0]->GetName()<<endl;

  bool doFitting = false;
  drawRecoil( h_recoil_inc_q2qe, doFitting );
  DoFit( h_recoil_inc_q2qe, fits );
  doFitting = true;
  drawRecoil( h_recoil_inc_q2qe, doFitting );
}
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
//Diagnostic Plots

void drawDiagnostic( MnvH2D** hists, string name, string xtitle, string axis="x", bool doRatio=false, bool areaNorm=false )
{
  //Distribution in each Q2
  vector<int> histsToUse({kData,kMC, kQELike, kQELikeNot, kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions});

  double area_norm = (areaNorm)? hists[kData]->Integral()/hists[kMC]->Integral() : 1;

  vector<TH2*>histErr;
  cout<<"Creating histErr"<<endl;
  for(unsigned int j = 0; j< nHistos;j++) histErr.push_back( NULL );
  for(auto c: histsToUse) histErr[c] = (TH2*)hists[c]->GetCVHistoWithError().Clone(Form("h_%s",names[c].c_str() ) ) ;
  for(auto c: histsToUse) 
  {
    if (c==kData) continue;
    histErr[c]->Scale(area_norm);
  }
  TH2* dataStat = (TH2*) histErr[kData]->Clone("dataStat");
  TH2* mcErr =(TH2*)  histErr[kMC]->Clone("mcErr");
  TH2* mc =(TH2*)  histErr[kMC]->Clone("mc");
  mcErr->SetFillStyle(1001);
  mcErr->SetFillColorAlpha(kRed,0.35);
  mcErr->SetMarkerSize(0 );


  //vector<MnvH2D*> vec_hists(hists, hists+nHistos );
  FormatHistos( histErr );

  cout<<"Creating histAndOpts"<<endl;
  vector<pair<TH2*, const char*> > histAndOpts;
  histAndOpts.push_back( std::make_pair( mcErr, "e2" ) );
  histAndOpts.push_back( std::make_pair( histErr[kMC], "hist" ) );
  histAndOpts.push_back( std::make_pair( histErr[kQELike_QE_H], "histl" ) );
  histAndOpts.push_back( std::make_pair( histErr[kQELike_QE_OTH], "histl" ) );
  histAndOpts.push_back( std::make_pair( histErr[kQELike_RES], "histl" ) );
  histAndOpts.push_back( std::make_pair( histErr[kQELike_2p2h], "histl" ) );
  histAndOpts.push_back( std::make_pair( histErr[kQELike_DIS], "histl" ) );
  histAndOpts.push_back( std::make_pair( histErr[kQELikeNot], "hist" ) );
  histAndOpts.push_back( std::make_pair( histErr[kData], "ep" ) );


  for( unsigned int j = 0; j< histAndOpts.size(); j++ )
  {
    if( doRatio ) 
    {
      histAndOpts[j].first->Divide( (histAndOpts[j].first), mc );
    }
    else histAndOpts[j].first->Scale(1);
  }

  //Define Legend First:
  TLegend* leg=new TLegend(0.9, 0.1, 1, 0.3);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry( dataStat, "MINERvA data", "lpe");
  leg->AddEntry( mc, "MINERvA Tune", "l");
  leg->AddEntry( histErr[kQELike_QE_H],"QE-H","l");
  leg->AddEntry( histErr[kQELike_QE_OTH],"QE-Oth","l");
  leg->AddEntry( histErr[kQELike_RES],"Resonant","l");
  leg->AddEntry( histErr[kQELike_2p2h],"2p2h","l");
  leg->AddEntry( histErr[kQELike_DIS],"DIS","l");

  TLatex* text = new TLatex();
  text->SetTextFont(92);
  text->SetTextSize(0.05);
  text->SetTextColor( kBlue );


  double *multipliers = NULL;

  //First project into X axis
  string ytitle="Q^{2}_{QE} (GeV^{2})";

  if(axis == "x")
  {
    GridCanvas* gc = plotXAxis1D( histAndOpts, xtitle, ytitle, multipliers ); 
    gc->SetGridx(true);
    gc->SetYLimits(-.2,mcErr->GetMaximum()*1.5);
    if(doRatio)gc->SetYLimits(-.2,2.2);
    gc->SetYTitle("Evt Rate #times 10^{3}");
    if(doRatio) gc->SetYTitle("Ratio to MnvGENIE");
    gc->Modified();
    leg->Draw("SAME");
    if(areaNorm) text->DrawLatexNDC( 0.7,0.8, "Area Normalized" );

    string fname=Form("plots/diagnostic-%s-axis_x-1d-ratio_%d.pdf",name.c_str(), doRatio );
    gc->Print(fname.c_str());
  }
  else
  {
    GridCanvas* gc = plotYAxis1D( histAndOpts, ytitle, xtitle, multipliers ); 
    gc->SetGridx(true);
    gc->SetYLimits(-.2,1.2);
    if(doRatio)gc->SetYLimits(-.2,2.2);
    gc->SetYTitle("Evt Rate #times 10^{3}");
    if(doRatio) gc->SetYTitle("Ratio to MnvGENIE");
    gc->Modified();
    leg->Draw("SAME");
    if(areaNorm) text->DrawLatexNDC( 0.8,0.8, "Area Normalized" );
    string fname=Form("plots/diagnostic-%s-axis_y-1d-ratio_%d.pdf",name.c_str(), doRatio );
    gc->Print(fname.c_str());
  }

  free(dataStat);
  free(mcErr);
  free(mc);
  for(unsigned int h = 0; h< nHistos;h++) free(histErr[h]);

}


void NeutronDiagnostics( TFile *f_orig, TFile *f_signal_weighted )
{
  cout<<"NeutronDiagnostics"<<endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  TFile *f1 = f_orig;
  //TFile *f2 = f_signal_weighted;

  //------- POT ----------------
  TVector2 *pot = (TVector2*)f1->Get("pot");
  double pot_data = pot->X();
  double pot_mc = pot->Y();
  double pot_norm = pot_data/pot_mc;

  CCQENuPlotUtils *putils = new CCQENuPlotUtils( fluxHistoExists );
  CCQENuUtils *utils = new CCQENuUtils( false, fluxHistoExists );

  MnvH2D *h_blobDist_q2qe[nHistos], *h_blobEnergy_q2qe[nHistos], *h_nBlobs_q2qe[nHistos],*h_n3DBlobs_q2qe[nHistos], *h_n2DBlobs_q2qe[nHistos];
  MnvH2D *h_muonTheta_q2qe[nHistos];
  MnvH2D *h_muonPhi_q2qe[nHistos];
  MnvH2D *h_muonE_q2qe[nHistos];
  MnvH2D *h_protonDEDX_q2qe[nHistos],*h_pionDEDX_q2qe[nHistos];
  MnvH2D *h_nonBlobEnergy_q2qe[nHistos];
  MnvH2D *h_hasTrack_q2qe[nHistos];
  MnvH2D *h_blobMaxE_q2qe[nHistos];
  MnvH2D *h_nClus_q2qe[nHistos];

  utils->bookHistos( f1, h_blobEnergy_q2qe, "h_blobEnergy_q2qe" );
  utils->bookHistos( f1, h_nonBlobEnergy_q2qe, "h_nonBlobEnergy_q2qe" );
  utils->bookHistos( f1, h_blobDist_q2qe, "h_blobDist_q2qe" );
  utils->bookHistos( f1, h_nBlobs_q2qe, "h_nBlobs_q2qe" );
  utils->bookHistos( f1, h_n3DBlobs_q2qe, "h_n3DBlobs_q2qe" );
  utils->bookHistos( f1, h_n2DBlobs_q2qe, "h_n2DBlobs_q2qe" );
  utils->bookHistos( f1, h_muonTheta_q2qe, "h_muonTheta_q2qe" );
  utils->bookHistos( f1, h_muonPhi_q2qe, "h_muonPhi_q2qe" );
  utils->bookHistos( f1, h_muonE_q2qe, "h_muonE_q2qe" );
  utils->bookHistos( f1, h_blobMaxE_q2qe, "h_blobMaxE_q2qe" );
  utils->bookHistos( f1, h_nClus_q2qe, "h_nClus_q2qe" );


  bool areaNorm = true;
  drawDiagnostic( h_blobEnergy_q2qe, "blobE", "E (MeV)", "x", false, areaNorm);
  drawDiagnostic( h_nonBlobEnergy_q2qe, "nonBlobE", "E (MeV)", "x", false,areaNorm );
  drawDiagnostic( h_blobDist_q2qe, "blobDist", "R (m)", "x" , false,areaNorm);
  drawDiagnostic( h_n3DBlobs_q2qe, "n3DBlobs", "n 3DBlobs", "x" , false,areaNorm);
  drawDiagnostic( h_n2DBlobs_q2qe, "n2DBlobs", "n 2DBlobs", "x" , false,areaNorm);
  drawDiagnostic( h_nBlobs_q2qe, "nBlobs", "nBlobs", "x" , false,areaNorm);
  drawDiagnostic( h_nBlobs_q2qe, "nBlobs", "nBlobs", "y" , false,areaNorm);
  drawDiagnostic( h_muonTheta_q2qe, "muonTheta", "#theta_{#mu} (degree)", "x" , false,areaNorm);
  drawDiagnostic( h_muonPhi_q2qe, "muonPhi", "#phi_{#mu} (degree)", "x" , false,areaNorm);
  drawDiagnostic( h_muonE_q2qe, "muonE", "E_{#mu} (GeV)", "x" , false,areaNorm);

  drawDiagnostic( h_blobMaxE_q2qe, "blobMaxE", "E_{blob} (GeV)", "x", false,areaNorm);
  drawDiagnostic( h_nClus_q2qe, "nClus", "#clusters", "x", false ,areaNorm);

  drawDiagnostic( h_blobEnergy_q2qe, "blobE", "E (MeV)", "x", true ,areaNorm);
  drawDiagnostic( h_nonBlobEnergy_q2qe, "nonBlobE", "E (MeV)", "x", true ,areaNorm);
  drawDiagnostic( h_blobDist_q2qe, "blobDist", "R (m)", "x" , true,areaNorm);
  drawDiagnostic( h_n3DBlobs_q2qe, "n3DBlobs", "n 3DBlobs", "x" , true,areaNorm);
  drawDiagnostic( h_n2DBlobs_q2qe, "n2DBlobs", "n 2DBlobs", "x" , true,areaNorm);
  drawDiagnostic( h_nBlobs_q2qe, "nBlobs", "nBlobs", "x" , true,areaNorm);
  drawDiagnostic( h_nBlobs_q2qe, "nBlobs", "nBlobs", "y" , true,areaNorm);
  drawDiagnostic( h_muonTheta_q2qe, "muonTheta", "#theta_{#mu} (degree)", "x" , true,areaNorm);
  drawDiagnostic( h_muonPhi_q2qe, "muonPhi", "#phi_{#mu} (degree)", "x" , true,areaNorm);
  drawDiagnostic( h_muonE_q2qe, "muonE", "E_{#mu} (GeV)", "x" , true,areaNorm);

  drawDiagnostic( h_blobMaxE_q2qe, "blobMaxE", "E_{blob} (GeV)", "x", true,areaNorm);
  drawDiagnostic( h_nClus_q2qe, "nClus", "#clusters", "x", true ,areaNorm);

  vector<MnvH2D*> fit_blobE = GetFitsHisto( f_signal_weighted, "blobEnergy" );
  vector<MnvH2D*> fit_nblob = GetFitsHisto( f_signal_weighted, "nBlobs" );
  vector<MnvH2D*> fit_blobDist = GetFitsHisto( f_signal_weighted, "blobDist" );

  DoFit( h_blobEnergy_q2qe, fit_blobE );
  DoFit( h_nBlobs_q2qe, fit_nblob );
  DoFit( h_n2DBlobs_q2qe, fit_nblob );
  DoFit( h_n3DBlobs_q2qe, fit_nblob );
  DoFit( h_blobDist_q2qe, fit_blobDist );

  drawDiagnostic( h_blobEnergy_q2qe, "blobE-fit_1", "E (MeV)", "x", false );
  drawDiagnostic( h_blobDist_q2qe, "blobDist-fit_1", "R (mm)", "x" , false);
  drawDiagnostic( h_n3DBlobs_q2qe, "n3DBlobs-fit_1", "n 3DBlobs", "x" , false);
  drawDiagnostic( h_n2DBlobs_q2qe, "n2DBlobs-fit_1", "n 2DBlobs", "x" , false);
  drawDiagnostic( h_nBlobs_q2qe, "nBlobs-fit_1", "n Blobs", "x" , false);

  drawDiagnostic( h_blobEnergy_q2qe, "blobE-fit_1", "E (MeV)", "x", true );
  drawDiagnostic( h_blobDist_q2qe, "blobDist-fit_1", "R (mm)", "x" , true);
  drawDiagnostic( h_n3DBlobs_q2qe, "n3DBlobs-fit_1", "n 3DBlobs", "x" , true);
  drawDiagnostic( h_n2DBlobs_q2qe, "n2DBlobs-fit_1", "n 2DBlobs", "x" , true);
  drawDiagnostic( h_nBlobs_q2qe, "nBlobs-fit_1", "n Blobs", "x" , true);





  return;

}

//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
int main(int argc, char* argv[])
{
  string f_orig = argv[1];
  string f_region = argv[2];
  string f_signal_weighted = argv[3];
  string sideband = argv[4];
  nu = atoi( argv[5] );
  cout<<"0: "<<argv[0]<<endl;
  cout<<"f_orig:"<<f_orig<<endl;
  cout<<"f_region:"<<f_region<<endl;
  if (argc == 7 ) nPlaylists = atoi(argv[6]);
  //draw scales with errors
  //
  //Get Files
  TFile* file_orig = new TFile( f_orig.c_str(),"read");
  TFile* file_region = new TFile( f_region.c_str(),"read");
  TFile* file_signal_weighted = new TFile( f_signal_weighted.c_str(),"read");


  //Draw 1D
  //DrawFits( file_signal_weighted );
  //
  string model = "";
  //draw3DChunk( f_orig, f_signal_weighted, sideband, true, model );//do bckfitting
  //draw3DChunk( f_orig, f_signal_weighted, sideband, false , model);//not do bckfitting

  //Draw regions
  drawRegionAnglePlots( file_region, file_orig, file_signal_weighted, true , sideband); //do bck fitting
  drawRegionAnglePlots( file_region, file_orig, file_signal_weighted, false, sideband ); //no bck fitting

  draw1DDistros( file_orig, file_signal_weighted, "" );
  NeutronDiagnostics( file_orig, file_signal_weighted );

  //drawRegionAnglePlots( file_region, file_orig, file_signal_weighted, true , sideband, "vtx"); //do bck fitting
  //drawRegionAnglePlots( file_region, file_orig, file_signal_weighted, false, sideband, "vtx" ); //no bck fitting

  //drawRegionAnglePlots( file_region, file_orig, file_signal_weighted, true , sideband, "nonvtx"); //do bck fitting
  //drawRegionAnglePlots( file_region, file_orig, file_signal_weighted, false, sideband, "nonvtx" ); //no bck fitting

  //PlotRecoil( file_orig, file_signal_weighted );

  //Diagnostics
  return 0;
}
