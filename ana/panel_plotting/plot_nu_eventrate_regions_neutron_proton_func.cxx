//#include "myPlotStyle.h"
#include "TParameter.h"

#include "include/GeneralIncludes.h"

//Modifications to do
//Set 2p2h - qe universe 2 weight to 0
//
//
//
//

using namespace PlotUtils;
//forward declaration
//gROOT->SetBatch();

int nPlaylists = 12;
int nu = 1;
bool fluxHistoExists = true;
bool drawNuwroOnly = false;
bool hasRESFit = false;


int nHistos2= nHistos;
//int nHistos2 = 55;
double Xline=-9999;

//============== Parameters ==============
double sh_lx0 = 0.6, sh_ly0 = 0.45, sh_lx1 = .85, sh_ly1 = .9;
double sh_ltextsize=0.032;

bool hasSystematics = false;

char* scaleOpt = "";

map<int,string> histLegMap={
  {kData, "MINERvA Data"},
  {kMC, "MINERvA Model"},
  {kQELike_QE_H,"CCE Hydrogen"},
  {kQELike_QE_OTH,"QELike CCQE"},
  {kQELike_RES,"QELike Resonant"},
  {kQELike_2p2h,"QELike 2p2h"},
  {kQELike_DIS,"QELike DIS"},
  {kQELikeNot,"Non-QELike"},
  {kQELikeNot_SingleChargedPion,"1 #pi^{+/-}"},
  {kQELikeNot_SingleNeutralPion,"1 #pi^{0}"},
  {kQELikeNot_MultiPion,"N #pi"},
  {kQELikeNot_NoPions,"Others"} };

map<int,string> region_map={
  {0, "CCE Signal"},
  {1, "QE Fit Region"},
  {2, "QE Validation"},
  {3, "Non-QE Validation 1"},
  {4, "Non-QE Validation 2"},
  {5, "Non-QE Fit Region"},
  {99, "Non-QE && Mesons"}
};

map<int,string> region_map_nu={
  {0, "#nu CCE Region"},
  {1, "#nu QE Region"},
  {2, "#nu QE Validation"},
  {3, "#nu Non-QE Validation 1"},
  {4, "#nu Non-QE Validation 2"},
  {5, "#nu Non-QE Fit Region"},
  {99, "#nu Non-QE && Mesons"}
};

map<int,string> region_map_gen = region_map;


//========================================



//vector<int> categories_to_fit({ kQELike_QE_OTH, kQELike_2p2h, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion} );
//vector<string> categories_to_fit_names({"qe_oth", "2p2h","qelikenot_scp","qelikenot_snp","qelikenot_mp"});
//vector<int> categories_to_fit({ kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion} );
//vector<string> categories_to_fit_names({"qelike_qe_oth", "qelike_res","qelike_2p2h","qelikenot_scp", "qelikenot_snp"});
//vector<string> categories_to_fit_title_names({"QELike && QE OTH", "QELike && RES", "QELike && 2p2h", "QELikeNot && Single Charged Pion", "QELikeNot && Single Neutral Pion"});

vector<int> categories_to_fit({ kQELike_QE_OTH,  kQELike_RES, kQELike_2p2h, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion} );
vector<string> categories_to_fit_names({"qelike_qe_oth", "qelike_res", "qelike_2p2h","qelikenot_scp", "qelikenot_snp"});
vector<string> categories_to_fit_title_names({"QELike && QE OTH",  "QELike && RES", "QELike && 2p2h", "QELikeNot && Single Charged Pion", "QELikeNot && Single Neutral Pion"});

//vector<int> categories_to_fit({ kQELike_QE_OTH, kQELike_2p2h, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion} );
//vector<string> categories_to_fit_names({"qelike_qe_oth", "qelike_2p2h","qelikenot_scp", "qelikenot_snp"});
//vector<string> categories_to_fit_title_names({"QELike && QE OTH", "QELike && 2p2h", "QELikeNot && Single Charged Pion", "QElikeNot && Single Neutral Pion"});


vector<int> histosUsed({ kData,  kMC, kQELike, kQELike_QE_H, kQELike_QE_OTH, kQELike_2p2h, kQELike_RES, kQELike_DIS, kQELike_OTH, kQELikeNot, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } ); 


//General Func
template<class T>
void FormatHistos(vector< T*> &h );
void FormatHistos( vector<vector<MnvH2D*>> h );
void drawSystematics( CCQENuPlotUtils*utils, MnvH1D* h,  double pot_data, double pot_mc, string name, bool asfrac=true);
void drawLatSystematicHist( MnvH1D* h,  string sysname, string savename="" ) ;
void drawVertSystematicHist( MnvH1D* h,  string sysname, string savename="" ) ;

int SetLegendEntries1D(TLegend* leg, vector<MnvH1D*> &hist, bool doBackSub, bool subQELike, bool drawNonQELike=true)
{
  leg->AddEntry(hist[kData],        histLegMap[kData].c_str(), "lpe");

  int nEntries = 2;

  if(!subQELike)
  {
    leg->AddEntry(hist[kMC],          histLegMap[kMC].c_str(), "l");
    leg->AddEntry(hist[kQELike_QE_H], histLegMap[kQELike_QE_H].c_str(),"l");
    leg->AddEntry(hist[kQELike_QE_OTH],histLegMap[kQELike_QE_OTH].c_str(),"l");
    leg->AddEntry(hist[kQELike_2p2h],histLegMap[kQELike_2p2h].c_str(),"l");
    leg->AddEntry(hist[kQELike_RES], histLegMap[kQELike_RES].c_str(),"l");
    leg->AddEntry(hist[kQELike_DIS],histLegMap[kQELike_DIS].c_str(),"l");
    nEntries+=5;
  } else leg->AddEntry(hist[kQELike_QE_H], histLegMap[kQELike_QE_H].c_str(),"l");


  if( !doBackSub )
  {
    leg->AddEntry(hist[kQELikeNot], histLegMap[kQELikeNot].c_str(),"l");
    nEntries+=1;
    if (drawNonQELike)
    {
      leg->AddEntry(hist[kQELikeNot_SingleChargedPion],histLegMap[kQELikeNot_SingleChargedPion].c_str(),"l");
      leg->AddEntry(hist[kQELikeNot_SingleNeutralPion],histLegMap[kQELikeNot_SingleNeutralPion].c_str(),"l");
      leg->AddEntry(hist[kQELikeNot_MultiPion],        histLegMap[kQELikeNot_MultiPion].c_str(),"l");
      leg->AddEntry(hist[kQELikeNot_NoPions],          histLegMap[kQELikeNot_NoPions].c_str(),"l");
      nEntries+=4;
    }
  }
  return nEntries;
}


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
void ReconfigureCategories( MnvHND**h)
{

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
template void ReconfigureCategories( MnvH1D** h);
template void ReconfigureCategories( MnvH2D** h);

template<class MnvHND>
vector<MnvHND*> ExtractVertUniverse( MnvHND** h, string errbandName, int iUniv )
{
  vector<MnvHND*> ret(nHistos2, NULL);
  if( h[kMC]->GetVertErrorBand( errbandName ) == NULL )
  {
    throw std::runtime_error( "Errorband does not exist" );
    return ret;
  }
  else if(  h[kMC]->GetVertErrorBand( errbandName )->GetNHists() <= iUniv ) 
  {
    throw std::runtime_error( "Universe does not exist" );
    return ret;
  }

  ret[0] = (MnvHND*) h[kData]->Clone( Form( "%s_%s_%d", h[kData]->GetName(), errbandName.c_str(), iUniv ) );
  ret[0]->ClearAllErrorBands();
  for( auto c: histosUsed )
  {
    if (c == kData) continue;
    ret[c] = new MnvHND( *(h[c]->GetVertErrorBand( errbandName )->GetHist(iUniv)) );
    ret[c]->SetName( Form( "%s_%s_%d", h[c]->GetName(), errbandName.c_str(), iUniv ) );
  }

  return ret;
}
template vector<MnvH1D*> ExtractVertUniverse( MnvH1D**h, string errbandName, int iUniv );
template vector<MnvH2D*> ExtractVertUniverse( MnvH2D**h, string errbandName, int iUniv );

template<class MnvHND>
vector<MnvHND*> ExtractLatUniverse( MnvHND** h, string errbandName, int iUniv )
{
  vector<MnvHND*> ret(nHistos2, NULL);
  if( h[kMC]->GetLatErrorBand( errbandName ) == NULL ) 
  {
    throw std::runtime_error( "Errorband does not exist" );
    return ret;
  }
  else if(  h[kMC]->GetLatErrorBand( errbandName )->GetNHists() >= iUniv ) 
  {
    throw std::runtime_error( "Universe does not exist" );
    return ret;
  }

  ret[0] = (MnvHND*) h[kData]->Clone( Form( "%s_%s_%d", h[kData]->GetName(), errbandName.c_str(), iUniv ) );
  for( auto c: histosUsed )
  {
    if (c == kData) continue;
    ret[c] = new MnvHND( *(h[c]->GetLatErrorBand( errbandName )->GetHist(iUniv)) );
    ret[c]->SetName( Form( "%s_%s_%d", h[c]->GetName(), errbandName.c_str(), iUniv ) );
  }
  return ret;
}
template vector<MnvH1D*> ExtractLatUniverse( MnvH1D**h, string errbandName, int iUniv );
template vector<MnvH2D*> ExtractLatUniverse( MnvH2D**h, string errbandName, int iUniv );


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



void DoFit( MnvH2D**h, vector<MnvH1D*>&fits )
{
  vector<MnvH2D*> fits2D;
  for( int i = 0; i<fits.size();i++ )
  {
    fits2D.push_back( new MnvH2D( h[kData]->GetCVHistoWithStatError() ) );
  }
  ExpandHistos( fits, fits2D, 1 );
  DoFit(h, fits2D);
  for( int i = 0; i<fits2D.size(); i++ ) delete fits2D[i];
}


template<class MnvHND>
vector<MnvHND*> CloneVector( vector<MnvHND*> &hists, bool copyUsedOnly = false )
{
  vector<MnvHND*> ret(hists.size(), NULL);

  if( !copyUsedOnly )
  {
    for( unsigned int i = 0; i< ret.size(); i++ )
    {
      if( hists[i] == NULL ) continue;
      ret[i] = (MnvHND* ) hists[i]->Clone( Form( "%s_clone", hists[i]->GetName() ) );
    }
  }
  else
  {
    for( auto i : histosUsed )
    {
      if( hists[i] == NULL ) continue;
      ret[i] = (MnvHND* ) hists[i]->Clone( Form( "%s_clone", hists[i]->GetName() ) );
    }
  }
  return ret;
}
template vector<MnvH1D*> CloneVector( vector<MnvH1D*> &hists, bool copyUsedOnly);
template vector<MnvH2D*> CloneVector( vector<MnvH2D*> &hists, bool copyUsedOnly);

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
void SubtractBackground( T** h, bool subqelike=false )
{
  h[kData]->AddMissingErrorBandsAndFillWithCV( *h[kMC] );
  h[kData]->Add( h[kQELikeNot], -1 );
  h[kMC]->Add( h[kQELikeNot], -1 );

  h[kQELikeNot]->Reset();
  h[kQELikeNot_SingleNeutralPion]->Reset();
  h[kQELikeNot_SingleChargedPion]->Reset();
  h[kQELikeNot_MultiPion]->Reset();
  h[kQELikeNot_NoPions]->Reset();

  if( subqelike )
  {
    h[kData]->Add( h[kQELike_QE_H] );
    h[kMC]->Add( h[kQELike_QE_H] );
    h[kData]->Add( h[kQELike], -1 );
    h[kMC]->Add( h[kQELike], -1 );

    h[kQELike_QE_OTH]->Reset();
    h[kQELike_2p2h]->Reset();
    h[kQELike_DIS]->Reset();
    h[kQELike_RES]->Reset();

    h[kQELike] = h[kQELike_QE_H];

  }
}



template<class T>
void RatioToMC( T** h )
{
  T* hmc = (T*) h[kMC]->Clone("hmc");
  vector<int> categories{ kData, kMC, kQELike, kQELike_QE_H, kQELike_QE_OTH, kQELike_QE, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH, kQELikeNot, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions };

  for( auto cat : categories ) {
    h[cat]->Divide(h[cat], hmc);
    h[cat]->GetYaxis()->SetTitle("Ratio to Model");
  }
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
void draw1DDistro( MnvH2D** h, vector<MnvH1D*> &fits, CCQENuPlotUtils* putils, string name, double pot_data, double pot_mc, string axis = "x", string xtitle="x", string ytitle="evt rate", bool doRatio=false, bool doFitting=false, bool doBckSub = false, bool subqelike=false, double min=-1., double max=-1.);

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
  MnvH2D *h_dpty_q2qe[nHistos2], *h_dptx_q2qe[nHistos2], *h_dpt_q2qe[nHistos2], 
         *h_pn_q2qe[nHistos2], *h_dalphat_q2qe[nHistos2], *h_dphit_q2qe[nHistos2],
         *h_dthetaP_q2qe[nHistos2], *h_dthetaR_q2qe[nHistos2], *h_ptheta_q2qe[nHistos2],
         *h_dtheta2D_q2qe[nHistos2],
         *h_muontheta_q2qe[nHistos2],*h_muonmomentum_q2qe[nHistos2],
         *h_enu_q2qe[nHistos2];
      

  utils->bookHistos( f1, h_dpt_q2qe, "h_dpt_q2qe" );
  utils->bookHistos( f1, h_dptx_q2qe, "h_dptx_q2qe" );
  utils->bookHistos( f1, h_dpty_q2qe, "h_dpty_q2qe" );

  utils->bookHistos( f1, h_pn_q2qe, "h_pn_q2qe" );
  utils->bookHistos( f1, h_dalphat_q2qe, "h_dalphat_q2qe" );
  utils->bookHistos( f1, h_dphit_q2qe, "h_dphit_q2qe" );

  utils->bookHistos( f1, h_dthetaP_q2qe, "h_dthetaP_q2qe" );
  utils->bookHistos( f1, h_dthetaR_q2qe, "h_dthetaR_q2qe" );
  utils->bookHistos( f1, h_ptheta_q2qe, "h_ptheta_q2qe" );

  utils->bookHistos( f1, h_dtheta2D_q2qe, "h_dtheta2D_q2qe" );
  cout<<"booking muon"<<endl;
  utils->bookHistos( f1, h_muontheta_q2qe, "h_muontheta_q2qe" );
  utils->bookHistos( f1, h_muonmomentum_q2qe, "h_muonmomentum_q2qe" );
  utils->bookHistos( f1, h_enu_q2qe, "h_enu_q2qe" );

  putils->scaleMCHistos(  h_dpt_q2qe,    pot_norm );
  putils->scaleMCHistos(  h_dptx_q2qe,   pot_norm );
  putils->scaleMCHistos(  h_dpty_q2qe,   pot_norm );
  putils->scaleMCHistos(  h_pn_q2qe,     pot_norm );
  putils->scaleMCHistos(  h_dalphat_q2qe,pot_norm );
  putils->scaleMCHistos(  h_dphit_q2qe,  pot_norm );
  putils->scaleMCHistos(  h_dthetaP_q2qe,pot_norm );
  putils->scaleMCHistos(  h_dthetaR_q2qe,pot_norm );
  putils->scaleMCHistos(  h_ptheta_q2qe, pot_norm );
  putils->scaleMCHistos(  h_muontheta_q2qe, pot_norm );
  putils->scaleMCHistos(  h_muonmomentum_q2qe, pot_norm );
  putils->scaleMCHistos(  h_enu_q2qe, pot_norm );
  string axis = "x";
  cout<<"Make Plots"<<endl;
  cout<<"================= Draw Histos ==============="<<endl;

  vector<MnvH1D*> fits = GetFitsHisto(f2);
  //vector<MnvH2D*> fits = GetFitsHisto(f2, "dpt");
  draw1DDistro( h_dpt_q2qe, fits, putils, "plot-1d-dpt", pot_data, pot_mc, "x", "#delta p_{T} (GeV/c)", "Event Rate", 0,0,0);
  draw1DDistro( h_dpt_q2qe, fits, putils, "plot-1d-dpt", pot_data, pot_mc, "x", "#delta p_{T} (GeV/c)", "Event Rate", 0,1,0);
  draw1DDistro( h_dpt_q2qe, fits, putils, "plot-1d-dpt", pot_data, pot_mc, "x", "#delta p_{T} (GeV/c)", "Event Rate", 0,1,1);
  draw1DDistro( h_dpt_q2qe, fits, putils, "plot-1d-dpt", pot_data, pot_mc, "x", "#delta p_{T} (GeV/c)", "Ratio To Model", 1,1,1);
  draw1DDistro( h_dpt_q2qe, fits, putils, "plot-1d-dpt", pot_data, pot_mc, "x", "#delta p_{T} (GeV/c)", "Ratio To Model", 1,1,1,1);
  draw1DDistro( h_dpt_q2qe, fits, putils, "plot-1d-dpt", pot_data, pot_mc, "x", "#delta p_{T} (GeV/c)", "Ratio To Model", 1,1,1,1);

  //fits =  GetFitsHisto(f2, "dptx");
  draw1DDistro( h_dptx_q2qe, fits, putils, "plot-1d-dptx", pot_data, pot_mc, "x", "#delta p_{Tx} (GeV/c)", "Event Rate", 0,0,0);
  draw1DDistro( h_dptx_q2qe, fits, putils, "plot-1d-dptx", pot_data, pot_mc, "x", "#delta p_{Tx} (GeV/c)", "Event Rate", 0,1,0);
  draw1DDistro( h_dptx_q2qe, fits, putils, "plot-1d-dptx", pot_data, pot_mc, "x", "#delta p_{Tx} (GeV/c)", "Event Rate", 0,1,1);
  draw1DDistro( h_dptx_q2qe, fits, putils, "plot-1d-dptx", pot_data, pot_mc, "x", "#delta p_{Tx} (GeV/c)", "Ratio To Model", 1,1,1);
  draw1DDistro( h_dptx_q2qe, fits, putils, "plot-1d-dptx", pot_data, pot_mc, "x", "#delta p_{Tx} (GeV/c)", "Event Rate", 0,1,1,1);
  draw1DDistro( h_dptx_q2qe, fits, putils, "plot-1d-dptx", pot_data, pot_mc, "x", "#delta p_{Tx} (GeV/c)", "Ratio To Model", 1,1,1,1);

  //fits =  GetFitsHisto(f2, "dpty");
  draw1DDistro( h_dpty_q2qe, fits, putils, "plot-1d-dpty", pot_data, pot_mc, "x", "#delta p_{Ty} (GeV/c)", "Event Rate", 0,0,0);
  draw1DDistro( h_dpty_q2qe, fits, putils, "plot-1d-dpty", pot_data, pot_mc, "x", "#delta p_{Ty} (GeV/c)", "Event Rate", 0,1,0);
  draw1DDistro( h_dpty_q2qe, fits, putils, "plot-1d-dpty", pot_data, pot_mc, "x", "#delta p_{Ty} (GeV/c)", "Event Rate", 0,1,1);
  draw1DDistro( h_dpty_q2qe, fits, putils, "plot-1d-dpty", pot_data, pot_mc, "x", "#delta p_{Ty} (GeV/c)", "Ratio To Model", 1,1,1);
  draw1DDistro( h_dpty_q2qe, fits, putils, "plot-1d-dpty", pot_data, pot_mc, "x", "#delta p_{Ty} (GeV/c)", "Event Rate", 0,1,1,1);
  draw1DDistro( h_dpty_q2qe, fits, putils, "plot-1d-dpty", pot_data, pot_mc, "x", "#delta p_{Ty} (GeV/c)", "Ratio To Model", 1,1,1,1);

  //fits =  GetFitsHisto(f2, "pn");
  draw1DDistro( h_pn_q2qe, fits, putils, "plot-1d-pn", pot_data, pot_mc, "x", "p_{n} (GeV/c)", "Event Rate", 0,0,0);
  draw1DDistro( h_pn_q2qe, fits, putils, "plot-1d-pn", pot_data, pot_mc, "x", "p_{n} (GeV/c)", "Event Rate", 0,1,0);
  draw1DDistro( h_pn_q2qe, fits, putils, "plot-1d-pn", pot_data, pot_mc, "x", "p_{n} (GeV/c)", "Event Rate", 0,1,1);
  draw1DDistro( h_pn_q2qe, fits, putils, "plot-1d-pn", pot_data, pot_mc, "x", "p_{n} (GeV/c)", "Ratio To Model", 1,1,1);
  draw1DDistro( h_pn_q2qe, fits, putils, "plot-1d-pn", pot_data, pot_mc, "x", "p_{n} (GeV/c)", "Event Rate", 0,1,1,1);
  draw1DDistro( h_pn_q2qe, fits, putils, "plot-1d-pn", pot_data, pot_mc, "x", "p_{n} (GeV/c)", "Ratio To Model", 1,1,1,1);

  //fits =  GetFitsHisto(f2, "dalphat");
  draw1DDistro( h_dalphat_q2qe, fits, putils, "plot-1d-dalphat", pot_data, pot_mc, "x", "#delta#alpha_{T} (degree)", "Event Rate", 0,0,0);
  draw1DDistro( h_dalphat_q2qe, fits, putils, "plot-1d-dalphat", pot_data, pot_mc, "x", "#delta#alpha_{T} (degree)", "Event Rate", 0,1,0);
  draw1DDistro( h_dalphat_q2qe, fits, putils, "plot-1d-dalphat", pot_data, pot_mc, "x", "#delta#alpha_{T} (degree)", "Event Rate", 0,1,1);
  draw1DDistro( h_dalphat_q2qe, fits, putils, "plot-1d-dalphat", pot_data, pot_mc, "x", "#delta#alpha_{T} (degree)", "Ratio To Model", 1,1,1);
  draw1DDistro( h_dalphat_q2qe, fits, putils, "plot-1d-dalphat", pot_data, pot_mc, "x", "#delta#alpha_{T} (degree)", "Event Rate", 0,1,1,1);
  draw1DDistro( h_dalphat_q2qe, fits, putils, "plot-1d-dalphat", pot_data, pot_mc, "x", "#delta#alpha_{T} (degree)", "Ratio To Model", 1,1,1,1);

  //fits =  GetFitsHisto(f2, "dtheta2D");
  draw1DDistro( h_dtheta2D_q2qe, fits, putils, "plot-1d-dtheta2D", pot_data, pot_mc, "x", "#delta#theta_{2D} (degree)", "Event Rate", 0,0,0);
  draw1DDistro( h_dtheta2D_q2qe, fits, putils, "plot-1d-dtheta2D", pot_data, pot_mc, "x", "#delta#theta_{2D} (degree)", "Event Rate", 0,1,0);
  draw1DDistro( h_dtheta2D_q2qe, fits, putils, "plot-1d-dtheta2D", pot_data, pot_mc, "x", "#delta#theta_{2D} (degree)", "Event Rate", 0,1,1);
  draw1DDistro( h_dtheta2D_q2qe, fits, putils, "plot-1d-dtheta2D", pot_data, pot_mc, "x", "#delta#theta_{2D} (degree)", "Ratio To Model", 1,1,1);
  draw1DDistro( h_dtheta2D_q2qe, fits, putils, "plot-1d-dtheta2D", pot_data, pot_mc, "x", "#delta#theta_{2D} (degree)", "Event Rate", 0,1,1,1);
  draw1DDistro( h_dtheta2D_q2qe, fits, putils, "plot-1d-dtheta2D", pot_data, pot_mc, "x", "#delta#theta_{2D} (degree)", "Ratio To Model", 1,1,1,1);

  //fits =  GetFitsHisto(f2, "muontheta");
  draw1DDistro( h_muontheta_q2qe, fits, putils, "plot-1d-muontheta", pot_data, pot_mc, "x", "#delta#theta (degree)", "Event Rate", 0,0,0);
  draw1DDistro( h_muontheta_q2qe, fits, putils, "plot-1d-muontheta", pot_data, pot_mc, "x", "#delta#theta (degree)", "Event Rate", 0,1,0);
  draw1DDistro( h_muontheta_q2qe, fits, putils, "plot-1d-muontheta", pot_data, pot_mc, "x", "#delta#theta (degree)", "Event Rate", 0,1,1);
  draw1DDistro( h_muontheta_q2qe, fits, putils, "plot-1d-muontheta", pot_data, pot_mc, "x", "#delta#theta (degree)", "Ratio To Model", 1,1,1);
  draw1DDistro( h_muontheta_q2qe, fits, putils, "plot-1d-muontheta", pot_data, pot_mc, "x", "#delta#theta (degree)", "Event Rate", 0,1,1,1);
  draw1DDistro( h_muontheta_q2qe, fits, putils, "plot-1d-muontheta", pot_data, pot_mc, "x", "#delta#theta (degree)", "Ratio To Model", 1,1,1,1);

  //fits =  GetFitsHisto(f2, "muonmomentum");
  draw1DDistro( h_muonmomentum_q2qe, fits, putils, "plot-1d-muonmomentum", pot_data, pot_mc, "x", "p_{#mu} (GeV/c)", "Event Rate", 0,0,0,0, 0,20);
  draw1DDistro( h_muonmomentum_q2qe, fits, putils, "plot-1d-muonmomentum", pot_data, pot_mc, "x", "p_{#mu} (GeV/c)", "Event Rate", 0,1,0,0, 0,20);
  draw1DDistro( h_muonmomentum_q2qe, fits, putils, "plot-1d-muonmomentum", pot_data, pot_mc, "x", "p_{#mu} (GeV/c)", "Event Rate", 0,1,1,0, 0,20);
  draw1DDistro( h_muonmomentum_q2qe, fits, putils, "plot-1d-muonmomentum", pot_data, pot_mc, "x", "p_{#mu} (GeV/c)", "Ratio To Model",1,1,1,0, 0,20);
  draw1DDistro( h_muonmomentum_q2qe, fits, putils, "plot-1d-muonmomentum", pot_data, pot_mc, "x", "p_{#mu} (GeV/c)", "Event Rate", 0,1,1,1, 0,20);
  draw1DDistro( h_muonmomentum_q2qe, fits, putils, "plot-1d-muonmomentum", pot_data, pot_mc, "x", "p_{#mu} (GeV/c)", "Ratio To Model",1,1,1,1, 0,20);

  draw1DDistro( h_enu_q2qe, fits, putils, "plot-1d-enu", pot_data, pot_mc, "x", "E_{#nu} (GeV)", "Event Rate", 0,0,0,0, 0,20);
  draw1DDistro( h_enu_q2qe, fits, putils, "plot-1d-enu", pot_data, pot_mc, "x", "E_{#nu} (GeV)", "Event Rate", 0,1,0,0, 0,20);
  draw1DDistro( h_enu_q2qe, fits, putils, "plot-1d-enu", pot_data, pot_mc, "x", "E_{#nu} (GeV)", "Event Rate", 0,1,1,0, 0,20);
  draw1DDistro( h_enu_q2qe, fits, putils, "plot-1d-enu", pot_data, pot_mc, "x", "E_{#nu} (GeV)", "Ratio To Model",1,1,1,0, 0,20);
  draw1DDistro( h_enu_q2qe, fits, putils, "plot-1d-enu", pot_data, pot_mc, "x", "E_{#nu} (GeV)", "Event Rate", 0,1,1,1, 0,20);
  draw1DDistro( h_enu_q2qe, fits, putils, "plot-1d-enu", pot_data, pot_mc, "x", "E_{#nu} (GeV)", "Ratio To Model",1,1,1,1, 0,20);



  //fits =  GetFitsHisto(f2, "dphit");
  //TCanvas *c = new TCanvas("c","c");
  //for( auto fit: fits )
  //{
  //  fit->Draw( "colz" );
  //  c->Print( Form("plots/fit_%s.pdf", fit->GetName()));
  //}
  draw1DDistro( h_dphit_q2qe, fits, putils, "plot-1d-dphit", pot_data, pot_mc, "x", "#delta#phi_{T} (degree)", "Event Rate", 0,0,0);
  draw1DDistro( h_dphit_q2qe, fits, putils, "plot-1d-dphit", pot_data, pot_mc, "x", "#delta#phi_{T} (degree)", "Event Rate", 0,1,0);
  draw1DDistro( h_dphit_q2qe, fits, putils, "plot-1d-dphit", pot_data, pot_mc, "x", "#delta#phi_{T} (degree)", "Event Rate", 0,1,1);
  draw1DDistro( h_dphit_q2qe, fits, putils, "plot-1d-dphit", pot_data, pot_mc, "x", "#delta#phi_{T} (degree)", "Ratio To Model", 1,1,1);
  draw1DDistro( h_dphit_q2qe, fits, putils, "plot-1d-dphit", pot_data, pot_mc, "x", "#delta#phi_{T} (degree)", "Event Rate", 0,1,1,1);
  draw1DDistro( h_dphit_q2qe, fits, putils, "plot-1d-dphit", pot_data, pot_mc, "x", "#delta#phi_{T} (degree)", "Ratio To Model", 1,1,1,1);

  //fits =  GetFitsHisto(f2, "dthetaP");
  draw1DDistro( h_dthetaP_q2qe, fits, putils, "plot-1d-dthetaP", pot_data, pot_mc, "x", "#delta#theta_{P} (degree)", "Event Rate", 0,0,0);
  draw1DDistro( h_dthetaP_q2qe, fits, putils, "plot-1d-dthetaP", pot_data, pot_mc, "x", "#delta#theta_{P} (degree)", "Event Rate", 0,1,0);
  draw1DDistro( h_dthetaP_q2qe, fits, putils, "plot-1d-dthetaP", pot_data, pot_mc, "x", "#delta#theta_{P} (degree)", "Event Rate", 0,1,1);
  draw1DDistro( h_dthetaP_q2qe, fits, putils, "plot-1d-dthetaP", pot_data, pot_mc, "x", "#delta#theta_{P} (degree)", "Event Rate", 0,1,1,1);
  draw1DDistro( h_dthetaP_q2qe, fits, putils, "plot-1d-dthetaP", pot_data, pot_mc, "x", "#delta#theta_{P} (degree)", "Ratio To Model", 1,1,1);
  draw1DDistro( h_dthetaP_q2qe, fits, putils, "plot-1d-dthetaP", pot_data, pot_mc, "x", "#delta#theta_{P} (degree)", "Ratio To Model", 1,1,1,1);
  draw1DDistro( h_dthetaP_q2qe, fits, putils, "plot-1d-dthetaP", pot_data, pot_mc, "x", "#delta#theta_{P} (degree)", "Ratio To Model", 1,1,0,0);
  draw1DDistro( h_dthetaP_q2qe, fits, putils, "plot-1d-dthetaP", pot_data, pot_mc, "x", "#delta#theta_{P} (degree)", "Ratio To Model", 1,0,0,0);

  //fits =  GetFitsHisto(f2, "dthetaR");
  draw1DDistro( h_dthetaR_q2qe, fits, putils, "plot-1d-dthetaR", pot_data, pot_mc, "x", "#delta#theta_{R} (degree)", "Event Rate", 0,0,0);
  draw1DDistro( h_dthetaR_q2qe, fits, putils, "plot-1d-dthetaR", pot_data, pot_mc, "x", "#delta#theta_{R} (degree)", "Event Rate", 0,1,0);
  draw1DDistro( h_dthetaR_q2qe, fits, putils, "plot-1d-dthetaR", pot_data, pot_mc, "x", "#delta#theta_{R} (degree)", "Event Rate", 0,1,1);
  draw1DDistro( h_dthetaR_q2qe, fits, putils, "plot-1d-dthetaR", pot_data, pot_mc, "x", "#delta#theta_{R} (degree)", "Ratio To Model", 1,1,1);
  draw1DDistro( h_dthetaR_q2qe, fits, putils, "plot-1d-dthetaR", pot_data, pot_mc, "x", "#delta#theta_{R} (degree)", "Event Rate", 0,1,1,1);
  draw1DDistro( h_dthetaR_q2qe, fits, putils, "plot-1d-dthetaR", pot_data, pot_mc, "x", "#delta#theta_{R} (degree)", "Ratio To Model", 1,1,1,1);
  draw1DDistro( h_dthetaR_q2qe, fits, putils, "plot-1d-dthetaR", pot_data, pot_mc, "x", "#delta#theta_{R} (degree)", "Ratio To Model", 1,1,0,0);
  draw1DDistro( h_dthetaR_q2qe, fits, putils, "plot-1d-dthetaR", pot_data, pot_mc, "x", "#delta#theta_{R} (degree)", "Ratio To Model", 1,0,0,0);

  //fits =  GetFitsHisto(f2, "ptheta");
  draw1DDistro( h_ptheta_q2qe, fits, putils, "plot-1d-ptheta", pot_data, pot_mc, "x", "#theta_{proton} (degree)", "Event Rate", 0,0,0);
  draw1DDistro( h_ptheta_q2qe, fits, putils, "plot-1d-ptheta", pot_data, pot_mc, "x", "#theta_{proton} (degree)", "Event Rate", 0,1,0);
  draw1DDistro( h_ptheta_q2qe, fits, putils, "plot-1d-ptheta", pot_data, pot_mc, "x", "#theta_{proton} (degree)", "Event Rate", 0,1,1);
  draw1DDistro( h_ptheta_q2qe, fits, putils, "plot-1d-ptheta", pot_data, pot_mc, "x", "#theta_{proton} (degree)", "Ratio To Model", 1,1,1);
  draw1DDistro( h_ptheta_q2qe, fits, putils, "plot-1d-ptheta", pot_data, pot_mc, "x", "#theta_{proton} (degree)", "Event Rate", 0,1,1,1);
  draw1DDistro( h_ptheta_q2qe, fits, putils, "plot-1d-ptheta", pot_data, pot_mc, "x", "#theta_{proton} (degree)", "Ratio To Model", 1,1,1,1);
}

void draw1DDistro( MnvH2D** h, vector<MnvH1D*> &fits, CCQENuPlotUtils* putils, string name, double pot_data, double pot_mc, string axis, string xtitle, string ytitle, bool doRatio, bool doFitting, bool doBckSub, bool subqelike, double min_x, double max_x)
{
  bool area_norm = false;
  bool includeData = true;

  MnvH2D* hists[nHistos2];
  for(unsigned int i = 0; i<nHistos2; i++) hists[i] = new MnvH2D(*( h[i]->Clone(Form("h2d_tmp_%s_%s", h[i]->GetName(), names[i].c_str() ))));

  if(doFitting) DoFit( hists, fits );
  if(doBckSub ) SubtractBackground( hists, subqelike );

  MnvH1D* hists1d[nHistos2];

  for(unsigned int i = 0; i<nHistos2; i++ ) 
  {
    MnvH1D* htmp;
    int minbin = 1;
    if(axis == "x") hists1d[i] = new MnvH1D(*(hists[i]->ProjectionX(Form("h1d_tmp_%s", hists[i]->GetName()),1,hists[i]->GetNbinsY())) );
    if(axis == "y") hists1d[i] = new MnvH1D(*(hists[i]->ProjectionY(Form("h1d_tmp_%s", hists[i]->GetName()),1,hists[i]->GetNbinsX())) );
    hists1d[i]->GetXaxis()->SetTitle( xtitle.c_str() );
    hists1d[i]->GetYaxis()->SetTitle( ytitle.c_str() );
  }
  //MnvH1D** hists1d = &hists1d_cp[0];
  if(doRatio ) RatioToMC( hists1d );

  TCanvas* c = new TCanvas("cNeutA","Neutron Angulars"); 
  MnvPlotter *plotter = new MnvPlotter;



  double mcscale = -1.; string xaxislabel = xtitle; string yaxislabel = ""; double min_y = -1.; double max_y = -1.;
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

  string plotType = "QELike_split_PionInFS";
  if( doBckSub )
  {
    plotType = "QELikeOnly";
    if (subqelike) plotType = "HOnly";
  }
  hists1d[kQELike_QE_H]->Scale(1.08);
  putils->drawStacked(hists1d, "QELike_split_PionInFS", area_norm, pot_data, pot_mc, includeData,mcscale,xaxislabel,yaxislabel,min_x, max_x, min_y ,max_y,normalizeHisto );
  hists1d[kQELike_QE_H]->Scale(1/1.08);
  if( min_x != -1 ) cout<<min_x<<"------------"<<max_x<<endl;
  MnvH1D* herr = (MnvH1D*) hists1d[kMC]->Clone("herr");
  if(!doRatio) herr->Scale(herr->GetBinWidth(1), scaleOpt);
  herr->SetFillColorAlpha(kBlue, 0.3);
  if( normalizeHisto ) herr->GetBinNormalizedCopy().Draw("e2same");
  else herr->Draw("e2same");

  string fname=Form("plots/%s-ratio_%d-bcksub_%d-subqelike_%d-fit_%d.pdf",name.c_str(), doRatio, doBckSub,subqelike, doFitting );
  plotter->AddPlotLabel(Form("# Hydrogen: %.2f", hists1d[kQELike_QE_H]->Integral()), 0.3,0.8,0.033,kBlue,32);
  c->Print( fname.c_str() );

  return;
}
void draw2DDistro( MnvH2D** h, vector<MnvH1D*> &fits, CCQENuPlotUtils* putils, string name, double pot_data, double pot_mc, string axis, string xtitle, string ytitle, bool doRatio, bool doFitting, bool doBckSub, bool subqelike, double min_x, double max_x)
{
  bool area_norm = false;
  bool includeData = true;

  MnvH2D* hists[nHistos2];
  for(unsigned int i = 0; i<nHistos2; i++) hists[i] = new MnvH2D(*( h[i]->Clone(Form("h2d_tmp_%s_%s", h[i]->GetName(), names[i].c_str() ))));

  if(doFitting) DoFit( hists, fits );
  if(doBckSub ) SubtractBackground( hists, subqelike );

  MnvH1D* hists1d[nHistos2];

  for(unsigned int i = 0; i<nHistos2; i++ ) 
  {
    MnvH1D* htmp;
    int minbin = 1;
    if(axis == "x") hists1d[i] = new MnvH1D(*(hists[i]->ProjectionX(Form("h1d_tmp_%s", hists[i]->GetName()),1,hists[i]->GetNbinsY())) );
    if(axis == "y") hists1d[i] = new MnvH1D(*(hists[i]->ProjectionY(Form("h1d_tmp_%s", hists[i]->GetName()),1,hists[i]->GetNbinsX())) );
    hists1d[i]->GetXaxis()->SetTitle( xtitle.c_str() );
    hists1d[i]->GetYaxis()->SetTitle( ytitle.c_str() );
  }
  //MnvH1D** hists1d = &hists1d_cp[0];
  if(doRatio ) RatioToMC( hists1d );

  TCanvas* c = new TCanvas("cNeutA","Neutron Angulars"); 
  MnvPlotter *plotter = new MnvPlotter;



  double mcscale = -1.; string xaxislabel = xtitle; string yaxislabel = ""; double min_y = -1.; double max_y = -1.;
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

  string plotType = "QELike_split_PionInFS";
  if( doBckSub )
  {
    plotType = "QELikeOnly";
    if (subqelike) plotType = "HOnly";
  }
  putils->drawStacked(hists1d, "QELike_split_PionInFS", area_norm, pot_data, pot_mc, includeData,mcscale,xaxislabel,yaxislabel,min_x, max_x, min_y ,max_y,normalizeHisto );
  if( min_x != -1 ) cout<<min_x<<"------------"<<max_x<<endl;
  TH1D* herr = (TH1D*) hists1d[kMC]->Clone("herr");
  if(!doRatio) herr->Scale(herr->GetBinWidth(1), scaleOpt);
  herr->SetFillColorAlpha(kBlue, 0.3);
  herr->Draw("e2same");

  string fname=Form("plots/%s-ratio_%d-bcksub_%d-subqelike_%d-fit_%d.pdf",name.c_str(), doRatio, doBckSub,subqelike, doFitting );
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
  MnvH2D *h_dthetaPdthetaR_q2qe[nHistos2];
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
  //for(int i =0; i<nHistos2;i++) delete h_dthetaPdthetaR_q2qe[nHistos2];
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
    for(unsigned int j = 0; j< nHistos2;j++) histErr.push_back( NULL );
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
        histAndOpts[j].first->Scale(0.001, scaleOpt);
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
    leg->AddEntry( mc, "MINERvA Model", "l");
    leg->AddEntry( histErr[kQELike_QE_H],"CCE","l");
    leg->AddEntry( histErr[kQELike_QE_OTH],"QELike CCQE","l");
    leg->AddEntry( histErr[kQELike_2p2h],"QELike 2p2h","l");
    leg->AddEntry( histErr[kQELike_RES],"QELike Resonant","l");
    leg->AddEntry( histErr[kQELike_DIS],"QELike DIS","l");
    leg->Draw("SAME");

    string fname=Form("plots/%s-%s-%s_%05d-ratio_%d-bcksub_%d-fit_%d.pdf",savename.c_str(),axis.c_str(),chunkname.c_str(), (int)zvalue, doRatio, doBckSub, doFitting) ;
    gc->Print(fname.c_str());
    free(dataStat);
    free(mcErr);
    free(mc);
    for(unsigned int h = 0; h< nHistos2;h++) free(histErr[h]);
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
      regioned_histos1D[iRegion].push_back( vector<MnvH1D*>(nHistos2,NULL) );

      //First add angular text
      string angle_text = Form("%.0f<#delta#theta_{R}.0f, %.0f<#delta#theta_{P}<%.0f", minR, maxR, minP,maxP );
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
  cout<<"drawChunkForRegion: looping over regions"<<endl;
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
      else h->Scale(1, scaleOpt);
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
        hd->Scale(1, scaleOpt);
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
  leg->AddEntry(hists[0][kMC], "MINERvA Model", "l");
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
      if(doScale) histAndOpts[j].first->Scale(0.001, scaleOpt);
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
    h[kMC][i]->SetLineWidth(3);

    h[kQELike_QE_H][i]->SetLineColor( mycolors[2] );
    //h[kQELike_QE_OTH][i]->SetLineColor( mycolors[3] );
    h[kQELike_QE_OTH][i]->SetLineColor( kGreen+1 );
    h[kQELike_RES][i]->SetLineColor( mycolors[4] );
    h[kQELike_DIS][i]->SetLineColor( mycolors[5] );
    h[kQELike_2p2h][i]->SetLineColor( mycolors[6] );
    int kQELikeNotColorIndex = 8;
    h[kQELikeNot][i]->SetLineColor( mycolors[kQELikeNotColorIndex] );

    h[kData][i]->SetMarkerStyle( 1 );
    h[kData][i]->SetMarkerSize(0.2);
    h[kData][i]->SetLineColor(kBlack);
    h[kData][i]->SetLineWidth(3);

    vector<int> qelike{kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_DIS, kQELike_2p2h, kQELikeNot};
    for( unsigned int j = 0; j<qelike.size();j++ )
    {
      int cat = qelike[j];
      h[cat][i]->SetLineWidth(2);
    }

    vector<int> qelikenots{kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion,kQELikeNot_NoPions};
    for( unsigned int j = 0; j<qelikenots.size();j++ )
    {
      int cat = qelikenots[j];
      h[cat][i]->SetFillStyle( 1001 );
      h[cat][i]->SetFillColorAlpha( mycolors[kQELikeNotColorIndex+j], 0.3 );
      h[cat][i]->SetLineColor(mycolors[kQELikeNotColorIndex+j]);
      h[cat][i]->SetLineStyle( 2 );
      h[cat][i]->SetMarkerSize(0);
    }

  }
}
//===============================================================================


//======================= Region Plots ================================//

void drawRegionAnglePlots(TFile* f_region, TFile* f_signal_weight, bool doBckFitting=true, string sample="Signal", string tag="", double xmin=1,double xmax=-1 );
void drawRegions( vector<vector<MnvH1D*>> &hists, bool doRatio, bool doBckSub, bool subQELike, bool hasFitting, int nRegions, bool angle=false, string tag="", double minx=1, double maxx=-1 );
void drawRegionsSeparate( vector<vector<MnvH1D*>> &hists, bool doRatio, bool doBckSub, bool subQELike, bool hasFitting, int nRegions, bool angle=false, string tag="", double xmin=1, double xmax=-1, double xline=-99999999, bool drawNonQELike = true );
void drawOneHist( vector<MnvH1D*> &hist, bool doRatio, bool doBckSub,  bool subQELike,bool hasFitting,string tag, string print="", double xmin=1,double xmax=-1, double xline=-9999999., bool drawNonQELike=true );

void drawRegionAnglePlots( TFile* f_region, TFile* f_signal_weight , bool doBckFitting, string sample , string tag,double xmin, double xmax )
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



  cout<<"--- Get CV histos ---"<<endl;
  // Get CV histos
  //
  //
  vector<MnvH1D*> h_qsq_signal;
  vector<MnvH1D*> h_qsq_non99;
  cout<<"------ nHistos2: "<<nHistos2<<endl;
  for( unsigned int j = 0; j< nHistos2;j ++ )
  {
    cout<<j<<endl;


    string hname = Form("h_q2qe_%s", names[j].c_str()  )  ;
    MnvH1D* h = (MnvH1D*) f1->Get(hname.c_str() );
    double norm = (j>0)? pot_norm : 1 ;
    if( j > 0 ) h->Scale(norm);
    h->SetMinimum(0);
    h_qsq_signal.push_back( h );

    hname = Form("h_q2qe_non99_%s", names[j].c_str()  )  ;
    MnvH1D* h2 = (MnvH1D*) f1->Get(hname.c_str() );
    if( j > 0 ) h2->Scale(norm);
    h2->SetMinimum(0);
    h_qsq_non99.push_back( h2 );


  }
  h_qsq_signal[0]->AddMissingErrorBandsAndFillWithCV( *h_qsq_signal[1] );
  h_qsq_non99[0]->AddMissingErrorBandsAndFillWithCV( *h_qsq_signal[1] );



  cout<<"draw regions"<<endl;
  vector<vector<MnvH1D*>> hvec_qsqqe;
  vector<vector<MnvH1D*>> hvec_qsqqe_sys_2;
  vector<vector<MnvH1D*>> hvec_qsqqe_sys_1;
  vector<vector<MnvH1D*>> hvec_qsqqe_sys_0;

  vector<vector<MnvH1D*>> hvec_qsqqe_sys_rw_neutron_0;
  vector<vector<MnvH1D*>> hvec_qsqqe_sys_rw_neutron_1;


  vector<vector<MnvH1D*>> hvec_qsqqe_sys_neutInt_0;
  vector<vector<MnvH1D*>> hvec_qsqqe_sys_neutInt_1;

  vector<vector<MnvH1D*>> hvec_qsqqe_sys_mares_0,hvec_qsqqe_sys_mares_1;
  int nRegions = 6;
  for( int i = 0; i< nRegions+1; i ++ )
  {
    int iR = i;
    if (i == nRegions ) iR = 99;
    vector<MnvH1D*> hists;
    for( unsigned int j = 0; j< nHistos2;j ++ )
    {
      string hname = Form("h_q2qe_region_%02d_%s",iR, names[j].c_str()  )  ;
      if( tag!="" ) hname = Form("h_q2qe_%s_region_%02d_%s",tag.c_str(), iR, names[j].c_str()  )  ;
      MnvH1D* h = (MnvH1D*) f1->Get(hname.c_str() );
      if( j == 0 ) h->AddMissingErrorBandsAndFillWithCV( *h_qsq_signal[1] );
      //MnvH1D* h = (MnvH1D*) f1->Get( "h_q2qe_region_00_data" );
      cout<<hname<<endl;
      //double norm = (j>0)? pot_norm : 1 ;
      //if( j > 0 ) h->Scale(norm);
      h->SetMinimum(0);
      hists.push_back( h );
      cout<<i<<j<<endl;
    }
    hvec_qsqqe.push_back( hists );
  }




  if( doBckFitting )
  {
    for( unsigned int i = 0; i< hvec_qsqqe.size(); i++ )
    {
      DoFit( &(hvec_qsqqe[i][0]), fits );
    }
    DoFit( &(h_qsq_signal[0]), fits );
    DoFit( &(h_qsq_non99[0]), fits );

    if(hasSystematics)
    {
      string latName="Muon_Energy_Resolution";
      string latName2="BeamAngleX";
      string vertName = "Neutron_Interaction";
      string vertName2 = "Low_Recoil_2p2h_Tune";
      string vertName3 = "Reweight_Neutron";

      MnvH1D* hdata = (MnvH1D*) hvec_qsqqe[0][kData]->Clone("hdata");
      drawLatSystematicHist( hdata, latName  , "dataRaw" ) ;
      drawLatSystematicHist( hdata, latName2  , "dataRaw" ) ;
      drawVertSystematicHist( hdata, vertName  , "dataRaw" ) ;
      drawVertSystematicHist( hdata, vertName2  , "dataRaw" ) ;
      drawVertSystematicHist( hdata, vertName3  , "dataRaw" ) ;

      MnvH1D* hqelikeh= (MnvH1D*) hvec_qsqqe[0][kQELike_QE_H]->Clone("hqelike_qe_h");
      drawLatSystematicHist( hqelikeh, latName  , "mch" ) ;
      drawLatSystematicHist( hqelikeh, latName2  , "mch" ) ;
      drawVertSystematicHist( hqelikeh, vertName  , "mch" ) ;
      drawVertSystematicHist( hqelikeh, vertName2 , "mch" ) ;
      drawVertSystematicHist( hqelikeh, vertName3 , "mch" ) ;

      MnvH1D* hbck = (MnvH1D*)hvec_qsqqe[0][kMC]->Clone("hbck");
      drawLatSystematicHist( hbck, latName  , "mc" ) ;
      drawLatSystematicHist( hbck, latName2  , "mc" ) ;
      drawVertSystematicHist( hbck, vertName  , "mc" ) ;
      drawVertSystematicHist( hbck, vertName2  , "mc" ) ;
      drawVertSystematicHist( hbck, vertName3  , "mc" ) ;

      hbck->Add( hvec_qsqqe[0][kQELike_QE_H], -1 );
      drawLatSystematicHist( hbck, latName  , "mcbck" ) ;
      drawLatSystematicHist( hbck, latName2  , "mcbck" ) ;
      drawVertSystematicHist( hbck, vertName  , "mcbck" ) ;
      drawVertSystematicHist( hbck, vertName2 , "mcbck" ) ;
      drawVertSystematicHist( hbck, vertName3 , "mcbck" ) ;
      
      hdata->AddMissingErrorBandsAndFillWithCV( *hbck );
      hdata->Add( hbck, -1 );
      drawLatSystematicHist( hdata, latName  , "data" ) ;
      drawLatSystematicHist( hdata, latName2  , "data" ) ;
      drawVertSystematicHist( hdata, vertName  , "data" ) ;
      drawVertSystematicHist( hdata, vertName2  , "data" ) ;
      drawVertSystematicHist( hdata, vertName3  , "data" ) ;
      bool asfrac=true;
      drawSystematics( putils, hdata,  pot_data, pot_mc, "region_data_00", asfrac);
      drawSystematics( putils, hqelikeh,  pot_data, pot_mc, "region_mc_00", asfrac);
      drawSystematics( putils, hbck,  pot_data, pot_mc, "region_mcbck_00", asfrac);


      hdata->Divide(hdata,hqelikeh);
      drawSystematics( putils, hdata,  pot_data, pot_mc, "region_ratio_00", asfrac);
    }

  }

  if( hasSystematics )
  {
    for( unsigned int i = 0; i < hvec_qsqqe.size(); i++ )
    {
      vector<MnvH1D*> hists0 = ExtractVertUniverse( &hvec_qsqqe[i][0], "Low_Recoil_2p2h_Tune", 0 );
      FormatHistos( hists0 );
      hvec_qsqqe_sys_0.push_back( hists0 );

      vector<MnvH1D*> hists1 = ExtractVertUniverse( &hvec_qsqqe[i][0], "Low_Recoil_2p2h_Tune", 1 );
      FormatHistos( hists1 );
      hvec_qsqqe_sys_1.push_back( hists1 );

      vector<MnvH1D*> hists2 = ExtractVertUniverse( &hvec_qsqqe[i][0], "Low_Recoil_2p2h_Tune", 2 );
      FormatHistos( hists2 );
      hvec_qsqqe_sys_2.push_back( hists2 );

      vector<MnvH1D*> histsr0 = ExtractVertUniverse( &hvec_qsqqe[i][0], "Reweight_Neutron", 0 );
      FormatHistos( histsr0 );
      hvec_qsqqe_sys_rw_neutron_0.push_back(histsr0);
      vector<MnvH1D*> histsr1 = ExtractVertUniverse( &hvec_qsqqe[i][0], "Reweight_Neutron", 1 );
      FormatHistos( histsr1 );
      hvec_qsqqe_sys_rw_neutron_1.push_back(histsr1);

      vector<MnvH1D*> histsni0 = ExtractVertUniverse( &hvec_qsqqe[i][0], "Neutron_Interaction", 0 );
      FormatHistos( histsr0 );
      hvec_qsqqe_sys_neutInt_0.push_back(histsni0);
      vector<MnvH1D*> histsni1 = ExtractVertUniverse( &hvec_qsqqe[i][0], "Neutron_Interaction", 1 );
      FormatHistos( histsr1 );
      hvec_qsqqe_sys_neutInt_1.push_back(histsni1);

      vector<MnvH1D*> histsmares0 = ExtractVertUniverse( &hvec_qsqqe[i][0], "GENIE_MaRES", 0 );
      FormatHistos( histsmares0 );
      hvec_qsqqe_sys_mares_0.push_back(histsmares0);
      vector<MnvH1D*> histsmares1 = ExtractVertUniverse( &hvec_qsqqe[i][0], "GENIE_MaRES", 1 );
      FormatHistos( histsmares1 );
      hvec_qsqqe_sys_mares_1.push_back(histsmares1);
    }


    cout<<"Drawing hvec_qsqqe_sys_2"<<endl;
    drawRegions( hvec_qsqqe_sys_0, false, false, false, doBckFitting,   nRegions, false, "sys_2p2h_0" );
    drawRegions( hvec_qsqqe_sys_0, false, true,  false, doBckFitting,   nRegions, false, "sys_2p2h_0" );
    drawRegions( hvec_qsqqe_sys_0, true, false,  false, doBckFitting,   nRegions, false, "sys_2p2h_0" );
    drawRegions( hvec_qsqqe_sys_0, true, true,   false, doBckFitting,   nRegions, false, "sys_2p2h_0" );
    drawRegions( hvec_qsqqe_sys_0, false, true,  true, doBckFitting,   nRegions, false,  "sys_2p2h_0" );
    drawRegions( hvec_qsqqe_sys_0, true,  true,  true, doBckFitting,   nRegions, false,  "sys_2p2h_0" );

    drawRegions( hvec_qsqqe_sys_1, false, false, false, doBckFitting,   nRegions, false, "sys_2p2h_1" );
    drawRegions( hvec_qsqqe_sys_1, false, true,  false, doBckFitting,   nRegions, false, "sys_2p2h_1" );
    drawRegions( hvec_qsqqe_sys_1, true, false,  false, doBckFitting,   nRegions, false, "sys_2p2h_1" );
    drawRegions( hvec_qsqqe_sys_1, true, true,   false, doBckFitting,   nRegions, false, "sys_2p2h_1" );
    drawRegions( hvec_qsqqe_sys_1, false, true,  true, doBckFitting,   nRegions, false,  "sys_2p2h_1" );
    drawRegions( hvec_qsqqe_sys_1, true,  true,  true, doBckFitting,   nRegions, false,  "sys_2p2h_1" );

    drawRegions( hvec_qsqqe_sys_2, false, false, false, doBckFitting,   nRegions, false, "sys_2p2h_2" );
    drawRegions( hvec_qsqqe_sys_2, false, true,  false, doBckFitting,   nRegions, false, "sys_2p2h_2" );
    drawRegions( hvec_qsqqe_sys_2, true, false,  false, doBckFitting,   nRegions, false, "sys_2p2h_2" );
    drawRegions( hvec_qsqqe_sys_2, true, true,   false, doBckFitting,   nRegions, false, "sys_2p2h_2" );
    drawRegions( hvec_qsqqe_sys_2, false, true,  true, doBckFitting,   nRegions, false,  "sys_2p2h_2" );
    drawRegions( hvec_qsqqe_sys_2, true,  true,  true, doBckFitting,   nRegions, false,  "sys_2p2h_2" );

    drawRegions( hvec_qsqqe_sys_rw_neutron_0, false, false, false, doBckFitting,   nRegions, false, "sys_rwn_0" );
    drawRegions( hvec_qsqqe_sys_rw_neutron_0, false, true,  false, doBckFitting,   nRegions, false, "sys_rwn_0" );
    drawRegions( hvec_qsqqe_sys_rw_neutron_0, true, false,  false, doBckFitting,   nRegions, false, "sys_rwn_0" );
    drawRegions( hvec_qsqqe_sys_rw_neutron_0, true, true,   false, doBckFitting,   nRegions, false, "sys_rwn_0" );
    drawRegions( hvec_qsqqe_sys_rw_neutron_0, false, true,  true, doBckFitting,   nRegions, false,  "sys_rwn_0" );
    drawRegions( hvec_qsqqe_sys_rw_neutron_0, true,  true,  true, doBckFitting,   nRegions, false,  "sys_rwn_0" );

    drawRegions( hvec_qsqqe_sys_rw_neutron_1, false, false, false, doBckFitting,   nRegions, false, "sys_rwn_1" );
    drawRegions( hvec_qsqqe_sys_rw_neutron_1, false, true,  false, doBckFitting,   nRegions, false, "sys_rwn_1" );
    drawRegions( hvec_qsqqe_sys_rw_neutron_1, true, false,  false, doBckFitting,   nRegions, false, "sys_rwn_1" );
    drawRegions( hvec_qsqqe_sys_rw_neutron_1, true, true,   false, doBckFitting,   nRegions, false, "sys_rwn_1" );
    drawRegions( hvec_qsqqe_sys_rw_neutron_1, false, true,  true, doBckFitting,   nRegions, false,  "sys_rwn_1" );
    drawRegions( hvec_qsqqe_sys_rw_neutron_1, true,  true,  true, doBckFitting,   nRegions, false,  "sys_rwn_1" );

    drawRegions( hvec_qsqqe_sys_neutInt_0, false, false, false, doBckFitting,   nRegions, false, "sys_ni_0" );
    drawRegions( hvec_qsqqe_sys_neutInt_0, false, true,  false, doBckFitting,   nRegions, false, "sys_ni_0" );
    drawRegions( hvec_qsqqe_sys_neutInt_0, true, false,  false, doBckFitting,   nRegions, false, "sys_ni_0" );
    drawRegions( hvec_qsqqe_sys_neutInt_0, true, true,   false, doBckFitting,   nRegions, false, "sys_ni_0" );
    drawRegions( hvec_qsqqe_sys_neutInt_0, false, true,  true, doBckFitting,   nRegions, false,  "sys_ni_0" );
    drawRegions( hvec_qsqqe_sys_neutInt_0, true,  true,  true, doBckFitting,   nRegions, false,  "sys_ni_0" );

    drawRegions( hvec_qsqqe_sys_neutInt_1, false, false, false, doBckFitting,   nRegions, false, "sys_ni_1" );
    drawRegions( hvec_qsqqe_sys_neutInt_1, false, true,  false, doBckFitting,   nRegions, false, "sys_ni_1" );
    drawRegions( hvec_qsqqe_sys_neutInt_1, true, false,  false, doBckFitting,   nRegions, false, "sys_ni_1" );
    drawRegions( hvec_qsqqe_sys_neutInt_1, true, true,   false, doBckFitting,   nRegions, false, "sys_ni_1" );
    drawRegions( hvec_qsqqe_sys_neutInt_1, false, true,  true, doBckFitting,   nRegions, false,  "sys_ni_1" );
    drawRegions( hvec_qsqqe_sys_neutInt_1, true,  true,  true, doBckFitting,   nRegions, false,  "sys_ni_1" );


    drawRegions( hvec_qsqqe_sys_mares_0, false, false, false, doBckFitting,   nRegions, false, "sys_mares_0" );
    drawRegions( hvec_qsqqe_sys_mares_0, false, true,  false, doBckFitting,   nRegions, false, "sys_mares_0" );
    drawRegions( hvec_qsqqe_sys_mares_0, true, false,  false, doBckFitting,   nRegions, false, "sys_mares_0" );
    drawRegions( hvec_qsqqe_sys_mares_0, true, true,   false, doBckFitting,   nRegions, false, "sys_mares_0" );
    drawRegions( hvec_qsqqe_sys_mares_0, false, true,  true, doBckFitting,   nRegions, false,  "sys_mares_0" );
    drawRegions( hvec_qsqqe_sys_mares_0, true,  true,  true, doBckFitting,   nRegions, false,  "sys_mares_0" );

    drawRegions( hvec_qsqqe_sys_mares_1, false, false, false, doBckFitting,   nRegions, false, "sys_mares_1" );
    drawRegions( hvec_qsqqe_sys_mares_1, false, true,  false, doBckFitting,   nRegions, false, "sys_mares_1" );
    drawRegions( hvec_qsqqe_sys_mares_1, true, false,  false, doBckFitting,   nRegions, false, "sys_mares_1" );
    drawRegions( hvec_qsqqe_sys_mares_1, true, true,   false, doBckFitting,   nRegions, false, "sys_mares_1" );
    drawRegions( hvec_qsqqe_sys_mares_1, false, true,  true, doBckFitting,   nRegions, false,  "sys_mares_1" );
    drawRegions( hvec_qsqqe_sys_mares_1, true,  true,  true, doBckFitting,   nRegions, false,  "sys_mares_1" );
  }


  cout<<"parsed region "<<hvec_qsqqe[0][kData]->GetEntries()<<endl;
  drawRegions( hvec_qsqqe, false, false, false, doBckFitting,   nRegions, false, tag, xmin, xmax );
  drawRegions( hvec_qsqqe, false, true,  false, doBckFitting,   nRegions, false, tag, xmin, xmax );
  drawRegions( hvec_qsqqe, true, false,  false, doBckFitting,   nRegions, false, tag, xmin, xmax );
  drawRegions( hvec_qsqqe, true, true,   false, doBckFitting,   nRegions, false, tag, xmin, xmax );

  drawRegions( hvec_qsqqe, false, true,  true, doBckFitting,   nRegions, false, tag, xmin, xmax );
  drawRegions( hvec_qsqqe, true,  true,  true, doBckFitting,   nRegions, false, tag, xmin, xmax );

  bool drawLeg;
  drawRegionsSeparate( hvec_qsqqe, false, false, false, doBckFitting,   nRegions, false, tag, xmin, xmax, Xline, false );
  drawRegionsSeparate( hvec_qsqqe, false, true,  false, doBckFitting,   nRegions, false, tag, xmin, xmax, Xline, false );
  drawRegionsSeparate( hvec_qsqqe, true, false,  false, doBckFitting,   nRegions, false, tag, xmin, xmax, Xline, false );
  drawRegionsSeparate( hvec_qsqqe, true, true,   false, doBckFitting,   nRegions, false, tag, xmin, xmax, Xline, false );

  drawRegionsSeparate( hvec_qsqqe, false, true,  true, doBckFitting,   nRegions, false, tag, xmin, xmax, Xline, false );
  drawRegionsSeparate( hvec_qsqqe, true,  true,  true, doBckFitting,   nRegions, false, tag, xmin, xmax, Xline, false );


  drawOneHist( h_qsq_signal, false, false, false, doBckFitting, "all", "Full Region", xmin, xmax, Xline, false );
  drawOneHist( h_qsq_signal, false, true,  false, doBckFitting, "all", "Full Region", xmin, xmax, Xline, false );
  drawOneHist( h_qsq_signal, true, false,  false, doBckFitting, "all", "Full Region", xmin, xmax, Xline, false );
  drawOneHist( h_qsq_signal, true, true,   false, doBckFitting, "all", "Full Region", xmin, xmax, Xline, false );

  drawOneHist( h_qsq_signal, false, true,  true, doBckFitting, "all", "Full Region" , xmin, xmax, Xline, false );
  drawOneHist( h_qsq_signal, true,  true,  true, doBckFitting, "all", "Full Region" , xmin, xmax, Xline, false );

  drawOneHist( h_qsq_non99, false, false, false, doBckFitting, "not99", "Not Region 99", xmin, xmax, Xline, false );
  drawOneHist( h_qsq_non99, false, true,  false, doBckFitting, "not99", "Not Region 99", xmin, xmax, Xline, false );
  drawOneHist( h_qsq_non99, true, false,  false, doBckFitting, "not99", "Not Region 99", xmin, xmax, Xline, false );
  drawOneHist( h_qsq_non99, true, true,   false, doBckFitting, "not99", "Not Region 99", xmin, xmax, Xline, false );

  drawOneHist( h_qsq_non99, false, true,  true, doBckFitting, "not99", "Not Region 99" , xmin, xmax, Xline, false );
  drawOneHist( h_qsq_non99, true,  true,  true, doBckFitting, "not99", "Not Region 99" , xmin, xmax, Xline, false );





  //_______________________________________________________
  // Draw Inner Angles
  if( tag == "" )
  {
    cout<<"draw angles"<<endl;
    vector<vector<MnvH1D*>> hvecA_qsqqe;
    int nAngles = 4;
    for( int i = 0; i< nAngles; i ++ )
    {
      int iR = i+1;
      vector<MnvH1D*> hists;
      for( unsigned int j = 0; j< nHistos2;j ++ )
      {
        string hname = Form("h_q2qe_angle_%02d_%s",iR, names[j].c_str()  )  ;
        MnvH1D* h = (MnvH1D*) f1->Get(hname.c_str() );
        if( j == 0 ) h->AddMissingErrorBandsAndFillWithCV( *h_qsq_signal[1] );
        //MnvH1D* h = (MnvH1D*) f1->Get( "h_q2qe_region_00_data" );
        cout<<hname<<endl;
        double norm = (j>0)? pot_norm : 1 ;
        if( j > 0 ) h->Scale(norm);
        h->SetMinimum(0);
        hists.push_back( h );
        cout<<i<<j<<endl;
        cout<<h<<endl;
      }
      hvecA_qsqqe.push_back( hists );
    }

    if( doBckFitting )
    {
      for( unsigned int i = 0; i< hvecA_qsqqe.size(); i++ )
      {
        int sideband = 0;
        if(sample=="BlobSideBand") sideband=1;
        else if(sample=="MichelSideBand") sideband=2;
        else if(sample=="MicBlobSideBand") sideband=3;

        DoFit( &(hvecA_qsqqe[i][0]), fits );
      }
      MnvH1D* hdata = (MnvH1D*) hvecA_qsqqe[2][kData]->Clone("hdata");
      MnvH1D* hbck = (MnvH1D*)hvec_qsqqe[2][kMC]->Clone("hbck");
      hbck->Add( hvec_qsqqe[2][kQELike_QE_H], -1 );
      hdata->AddMissingErrorBandsAndFillWithCV( *hbck );
      hdata->Add( hbck, -1 );
      drawSystematics( putils, hdata,  pot_data, pot_mc, "angle_02");

    }
    cout<<"parsed angle "<<hvecA_qsqqe[0][kData]->GetEntries()<<endl;
    drawRegions( hvecA_qsqqe, false, false, false, doBckFitting,   nAngles, true );
    drawRegions( hvecA_qsqqe, false, true,  false, doBckFitting,   nAngles , true);
    drawRegions( hvecA_qsqqe, true, false,  false, doBckFitting,   nAngles, true );
    drawRegions( hvecA_qsqqe, true, true,   false, doBckFitting,   nAngles, true );

    drawRegions( hvecA_qsqqe, false, true,  true, doBckFitting,   nAngles, true );
    drawRegions( hvecA_qsqqe, true,  true,  true, doBckFitting,   nAngles, true );

    drawRegionsSeparate( hvecA_qsqqe, false, false, false, doBckFitting,   nAngles, true );
    drawRegionsSeparate( hvecA_qsqqe, false, true,  false, doBckFitting,   nAngles, true );
    drawRegionsSeparate( hvecA_qsqqe, true, false,  false, doBckFitting,   nAngles, true );
    drawRegionsSeparate( hvecA_qsqqe, true, true,   false, doBckFitting,   nAngles, true );

    drawRegionsSeparate( hvecA_qsqqe, false, true,  true, doBckFitting,   nAngles, true );
    drawRegionsSeparate( hvecA_qsqqe, true,  true,  true, doBckFitting,   nAngles, true );
  }




  //f1->Close();
}

template<class T>
void FormatHistos(vector<T*> &h )
{
  int size=3;
  gStyle->SetEndErrorSize(size);
  vector<int> mycolors = getColors(2);
  for( auto c: histosUsed )
  {
    h[c]->SetLineStyle(1);
    h[c]->SetLineWidth(size);
  }
  h[kMC]->SetLineWidth(3);
  h[kData]->SetLineWidth(2);

  h[kMC]->SetLineColor( kRed );

  h[kQELike_QE_H]->SetLineColor( mycolors[2] );
  //h[kQELike_QE_OTH]->SetLineColor( mycolors[3] );
  h[kQELike_QE_OTH]->SetLineColor( kGreen+1 );
  h[kQELike_RES]->SetLineColor( mycolors[4] );
  h[kQELike_DIS]->SetLineColor( mycolors[5] );
  h[kQELike_2p2h]->SetLineColor( mycolors[6] );
  int co=8;
  //h[kQELikeNot]->SetLineColor( mycolors[co] );
  h[kQELikeNot]->SetLineColor( kGray );

  h[kData]->SetMarkerStyle( 1 );
  h[kData]->SetMarkerSize(0.2);
  h[kData]->SetLineColor(kBlack);


  vector<int> qelikenots{kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion,kQELikeNot_NoPions};
  for( unsigned int i = 0; i<qelikenots.size();i++ )
  {
    int cat = qelikenots[i];
    h[cat]->SetFillStyle( 1001 );
    h[cat]->SetFillColorAlpha( mycolors[co+i], 0.3 );
    h[cat]->SetLineColor(mycolors[co+i]);
    h[cat]->SetLineStyle( 2 );
    h[cat]->SetLineWidth( 1 );
    h[cat]->SetMarkerSize(0);
  }
}

template void FormatHistos( vector<MnvH1D*> &h );
template void FormatHistos( vector<MnvH2D*> &h );

void drawRegions( vector<vector<MnvH1D*>> &hists, bool doRatio, bool doBckSub,  bool subQELike,bool hasFitting, int nRegions, bool angles, string tag, double minx, double maxx )
{
  cout<<"enter drawRegions"<<endl;
  cout<<hists.size()<<" regions with "<<hists[0].size()<<" histos"<<endl;
  cout<<"first histo: "<<hists[0][0]->GetName()<<endl;

  string xtitle = "Q^{2}_{QE} (GeV^{2})";
  string ytitle = "Events / GeV^{2}";
  if( scaleOpt == "" ) ytitle = "Events/Bin";
  if( doRatio ) ytitle = "Ratio to Model";

  TCanvas *c = new TCanvas("c","c",1200,800);
  int nCanvas = 8;
  if(!angles)c->Divide(4,2);
  else {c->Divide(3,2); nCanvas=6;}
  c->cd();
  vector<int> histosToDraw({ kData, kMC,    kMC,       kQELikeNot, kQELike_DIS, kQELike_RES, kQELike_2p2h, kQELike_QE_OTH, kQELike_QE_H, kData} );
  vector<string>      opts({"pe1", "e2same","histsame", "histsame", "histlsame","histlsame","histlsame","histlsame","histlsame","pe1same"});

  vector<int> qelikenot_hists({kQELikeNot_SingleNeutralPion, kQELikeNot_SingleChargedPion, kQELikeNot_MultiPion, kQELikeNot_NoPions});

  //vector<TH1D*> th1_histos;
  cout<<"drawRegions: looping over regions"<<endl;
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

    if(doBckSub) 
    {
      hmc->Add( hqelikenot, -1 );
      if( subQELike ) hmc->Add( hqelike, -1 );
      hmc->SetMinimum(0);
    }



    THStack* hs = new THStack( Form("%d",i),"");
    for( auto cat : qelikenot_hists ) 
    {
      TH1D* h = (TH1D*) (hists[i][cat]->GetCVHistoWithError()).Clone( Form( "hs_%d_%s",i, names[cat].c_str() ) );
      if(doRatio) {
        h->Divide( hmc );
      }
      else h->Scale(1, scaleOpt);
      hs->Add(h);
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
      h->GetXaxis()->SetTitle(xtitle.c_str() );
      h->GetYaxis()->SetTitle(ytitle.c_str() );
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
        if ( J == kQELike_QE_OTH || J == kQELike_RES || J==kQELike_2p2h || J==kQELike_DIS || J==kMC) h->Reset();
        h->SetMinimum(0);
      }


      if( doRatio ) 
      {
        //if(J!=0) h->Divide(h,hmc);
        //else h->DivideSingle(h,  hmc );
        h->Divide(h,hmc);
        if(J==0 && doBckSub && subQELike &&i==0) 
        {
          TH1D*htest = new TH1D(h->GetCVHistoWithError() );
          cout<<"========Get Error"<<endl;
          for( int b=1; b<h->GetNbinsX()+1;b++) cout<<htest->GetBinError(b)<<"\t";
          cout<<endl;
        }
      }

      TH1D* hd = new TH1D( h->GetCVHistoWithError() );
      TH1D* hds = new TH1D( h->GetCVHistoWithStatError() );
      hd->SetName( Form("th1d_%d%d", i,j) );
      if( doRatio ){
        hd->SetMaximum(2.5);
        hd->SetMinimum(0);
      }
      else
      {
        hd->Scale(1, scaleOpt);
        hds->Scale(1, scaleOpt);
        hd->SetMaximum( hd->GetMaximum()*1.5 );
      }

      hd->GetXaxis()->SetRangeUser(minx,maxx);
      hd->Draw(opts[j].c_str());
      if( J==0 )
      {
        hds->Draw("e1same");

      }
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
  int nLegEntries = SetLegendEntries1D(leg, hists[0], doBckSub, subQELike );
  //leg->AddEntry(hists[0][kData], "MINERvA data", "lpe");
  //leg->AddEntry(hists[0][kMC], "MINERvA Tune", "l");
  //leg->AddEntry(hists[0][kQELike_QE_H],"QE-H","l");
  //leg->AddEntry(hists[0][kQELike_QE_OTH],"QE-Oth","l");
  //leg->AddEntry(hists[0][kQELike_RES],"Resonant","l");
  //leg->AddEntry(hists[0][kQELike_DIS],"DIS","l");
  //leg->AddEntry(hists[0][kQELike_2p2h],"2p2h","l");
  //leg->AddEntry(hists[0][kQELikeNot],"Not QELike","l");
  //if( !doBckSub )
  //{
  //  leg->AddEntry(hists[0][kQELikeNot_SingleChargedPion],"SCP","l");
  //  leg->AddEntry(hists[0][kQELikeNot_SingleNeutralPion],"SNP","l");
  //  leg->AddEntry(hists[0][kQELikeNot_MultiPion],"MP","l");
  //  leg->AddEntry(hists[0][kQELikeNot_NoPions],"NP","l");


  //}
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

void drawRegionsSeparate( vector<vector<MnvH1D*>> &hists, bool doRatio, bool doBckSub,  bool subQELike,bool hasFitting, int nRegions, bool angles, string tag, double xmin, double xmax, double xline, bool drawNonQELike )
{
  cout<<"enter drawRegionsSeparate"<<endl;
  myPlotStyle();
  gStyle->SetEndErrorSize(4);
  vector<int> histosToDraw({ kData, kMC,    kMC,       kQELikeNot, kQELike_DIS, kQELike_RES, kQELike_2p2h, kQELike_QE_OTH, kQELike_QE_H, kData} );
  vector<string>      opts({"pe1", "e2same","histsame", "histsame", "histlsame","histlsame","histlsame","histlsame","histlsame","pe1same"});

  string xtitle = "Q^{2}_{QE} (GeV^{2})";
  string ytitle = "Events / GeV^{2}";
  if( scaleOpt == "" ) ytitle = "Events / Bin";
  if( doRatio ) ytitle = "Ratio to Model";

  vector<int> qelikenot_hists({kQELikeNot_SingleNeutralPion, kQELikeNot_SingleChargedPion, kQELikeNot_MultiPion, kQELikeNot_NoPions});

  //vector<TH1D*> th1_histos;
  cout<<"drawRegionsSeparate: looping over regions"<<endl;
  TLatex *text = new TLatex();
  text->SetTextFont(92);
  text->SetTextSize(0.05);
  text->SetTextColor( kBlue );
  double x0 = sh_lx0, x1 = sh_lx1, y0 = sh_ly0, y1 = sh_ly1;
  int nColumns = 1;
  if (doRatio)
  {
    x0 = 0.4; y0 = 0.7;
    y1 = 0.88;
    nColumns=2;
  }

  TLegend* leg=new TLegend(x0, y0, x1,y1);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(52);
  leg->SetTextSize(sh_ltextsize);
  leg->SetNColumns(nColumns);
  SetLegendEntries1D( leg, hists[0], doBckSub, subQELike, drawNonQELike );
  //leg->AddEntry(hists[0][kData], "MINERvA data", "lpe");
  //leg->AddEntry(hists[0][kMC], "MINERvA Tune", "l");
  //leg->AddEntry(hists[0][kQELike_QE_H],"Hydrogen CCQE","l");
  //if( !subQELike )
  //{
  //  leg->AddEntry(hists[0][kQELike_QE_OTH],"QELike QE Carbon","l");
  //  leg->AddEntry(hists[0][kQELike_RES],"QELike RES","l");
  //  leg->AddEntry(hists[0][kQELike_DIS],"QELike DIS","l");
  //  leg->AddEntry(hists[0][kQELike_2p2h],"QELike 2p2h","l");
  //  leg->AddEntry(hists[0][kQELikeNot],"Not QELike","l");
  //}

  //if( !doBckSub )
  //{
  //  leg->AddEntry(hists[0][kQELikeNot_SingleChargedPion],"1 #pi^{+/-}","l");
  //  leg->AddEntry(hists[0][kQELikeNot_SingleNeutralPion],"1 #pi^{0}","l");
  //  leg->AddEntry(hists[0][kQELikeNot_MultiPion],"N #pi","l");
  //  leg->AddEntry(hists[0][kQELikeNot_NoPions],"Others","l");
  //}




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
      if(doRatio) 
      {
        h->Divide( hmc );
      }
      else h->Scale(1, scaleOpt);
      hs->Add(h);
    }
    if(doBckSub) 
    {
      hmc->Add( hqelikenot, -1 );
      if( subQELike ) hmc->Add( hqelike, -1 );
      hmc->SetMinimum(0);
    }

    double ymin = 1;
    double ymax = -1;

    //hists[i][kData]->Draw("pe");
    for( unsigned int j = 0; j< histosToDraw.size(); j++ )
    {
      int J = histosToDraw[j];
      cout<<"histo to draw: "<<J<<endl;
      //TH1D *h= new TH1D( hists[i][J]->GetCVHistoWithError() );
      MnvH1D* h = (MnvH1D*) hists[i][J]->Clone( Form("%d%d", i,j) );
      //h->SetName( Form("%d%d", i,j) );
      //th1_histos.push_back(h);
      h->GetXaxis()->SetTitle(xtitle.c_str() );
      h->GetYaxis()->SetTitle(ytitle.c_str() );
      h->SetMinimum(0);
      h->GetXaxis()->SetNdivisions(4);

      double size = 0.055;
      double size2 = 0.055;
      h->GetXaxis()->SetTitleSize(size);
      h->GetXaxis()->SetTitleOffset(1.3);
      h->GetXaxis()->SetLabelSize(size2);
      h->GetYaxis()->SetTitleSize(size);
      h->GetYaxis()->SetTitleOffset(1.3);
      h->GetYaxis()->SetLabelSize(size2);


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
        if ( J == kQELike_QE_OTH || J == kQELike_RES || J==kQELike_2p2h || J==kQELike_DIS || J==kMC) h->Reset();
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
        hd->Scale(1, scaleOpt);
        hd->SetMaximum( hd->GetMaximum()*1.5 );
        hdstat->Scale(1, scaleOpt);
        hdstat->SetMaximum( hd->GetMaximum()*1.5 );
      }

      hd->SetMinimum(0);
      hd->GetXaxis()->SetRangeUser(xmin, xmax);
      hd->Draw(opts[j].c_str());
      if(J == kData ) 
      {
        //hdstat->SetLineColor(kRed);
        hdstat->Draw("pe1same");
        ymin = hd->GetMinimum();
        ymax = hd->GetMaximum();
      }
      if(!doBckSub && j == 0 && drawNonQELike) hs->Draw("histsame");

    }
    int iR = i;
    if(i==nRegions) iR = 99;
    if(!angles)
    { 
      //if( i == 0 ) text->DrawLatexNDC(0.2,0.8, Form("Signal"));
      //else text->DrawLatexNDC(0.2,0.8, Form("Region %d", iR));
      text->DrawLatexNDC(0.2,0.8, Form("%s",region_map_gen[iR].c_str()));
    }
    else text->DrawLatexNDC(0.2,0.8, Form("#delta#theta_{R/P} < %d^{o}", i+1));
    //leg->Draw("same");


    if(xline>0)
    {
      double height = 0.8;
      if (doRatio) height=0.6;
      MnvPlotter *plotter = new MnvPlotter;
      plotter->arrow_line_width=3;
      plotter->arrow_line_color=kGray+2;
      plotter->AddCutArrow( xline, 0, (ymax-0)*height+ymin,.2,"R");
    }


    //delete hqelikenot;
    //delete hmc
    if(!angles) 
    {
      TString savename = Form("plots/nu-1d-regions-ratio_%d-bcksub_%d-subqelike_%d-fit_%d-region_%02d.pdf", doRatio, doBckSub, subQELike, hasFitting,iR);
      if( tag!="") savename = Form("plots/nu-1d-%s-regions-ratio_%d-bcksub_%d-subqelike_%d-fit_%d-region_%02d.pdf", tag.c_str(), doRatio, doBckSub, subQELike, hasFitting,iR);
      c->Print( savename );
    }
    else c->Print(Form("plots/nu-1d-angles-ratio_%d-bcksub_%d-subqelike_%d-fit_%d-angle_%02d.pdf", doRatio, doBckSub, subQELike, hasFitting,iR));


  TCanvas *ccanvas;
  if( nColumns == 1 ) ccanvas = new TCanvas("ccanvas","ccanvas",100,200);
  else ccanvas = new TCanvas("ccanvas","ccanvas",200,150);
  ccanvas->cd();
  leg->SetX1NDC(.0);
  leg->SetY1NDC(.0);
  leg->SetX2NDC(1);
  leg->SetY2NDC(1);
  leg->Draw();
  ccanvas->Print(Form("plots/one_hist_leg_ratio_%d-bcksub_%d-subqelike_%d.pdf", doRatio, doBckSub, subQELike));
  delete ccanvas;


    delete c;
  }
}

void drawOneHist( vector<MnvH1D*> &hist, bool doRatio, bool doBckSub,  bool subQELike,bool hasFitting, string tag,string print, double xmin, double xmax, double xline, bool drawNonQELike)
{

  MnvPlotter *plotter = new MnvPlotter;
  vector<int> histosToDraw({ kData, kMC,    kMC,       kQELikeNot, kQELike_DIS, kQELike_RES, kQELike_2p2h, kQELike_QE_OTH, kQELike_QE_H, kData} );
  vector<string>      opts({"pe", "e2same","histsame", "histsame", "histlsame","histlsame","histlsame","histlsame","histlsame","pesame"});

  vector<int> qelikenot_hists({kQELikeNot_SingleNeutralPion, kQELikeNot_SingleChargedPion, kQELikeNot_MultiPion, kQELikeNot_NoPions});


  string xtitle = "Q^{2}_{QE} (GeV^{2})";
  string ytitle = "Events / GeV^{2}";
  if( scaleOpt == "" ) ytitle = "Events/Bin";
  if( doRatio ) ytitle = "Ratio to Model";




  //vector<TH1D*> th1_histos;
  cout<<"drawOneHist: looping over regions"<<endl;
  TLatex *text = new TLatex();
  text->SetTextFont(92);
  text->SetTextSize(0.05);
  text->SetTextColor( kBlue );
  double x0 = sh_lx0, x1 = sh_lx1, y0 = sh_ly0, y1 = sh_ly1;
  int nColumns = 1;
  if (doRatio)
  {
    x0 = 0.4; y0 = 0.7;
    y1 = 0.88;
    nColumns=2;
  }

  TLegend* leg=new TLegend(x0,y0,x1,y1);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(52);
  leg->SetTextSize(sh_ltextsize);
  leg->SetNColumns(nColumns);
  SetLegendEntries1D( leg, hist, doBckSub, subQELike , drawNonQELike);



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
    if(doRatio)
    { 
      h->Divide( hmc );
    }
    else h->Scale(1, scaleOpt);
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
    h->GetXaxis()->SetTitle(xtitle.c_str());
    h->GetYaxis()->SetTitle(ytitle.c_str());
    h->SetMinimum(0);
    h->GetXaxis()->SetNdivisions(4);

    h->GetXaxis()->SetTitleOffset(1.3);
    h->GetYaxis()->SetTitleOffset(1.3);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);


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
      hd->Scale(1, scaleOpt);
      hd->SetMaximum( hd->GetMaximum()*1.5 );
    }

    hd->SetMinimum(0);
    hd->GetXaxis()->SetRangeUser(xmin, xmax);
    hd->Draw(opts[j].c_str());
    if(!doBckSub && j == 0 && drawNonQELike) hs->Draw("histsame");

    if(xline> 0)
    {
      double ymin = hd->GetMinimum();
      double ymax = hd->GetMaximum();
      plotter->AddCutArrow( xline, ymin, (ymax-ymin)*0.8+ymin,2,"R");
    }

  }
  //leg->Draw("same");
  text->DrawLatexNDC(0.2,0.8, print.c_str());

  //TLine line(xline, ymin ,xline, ymax ); 
  //line->SetLineColor(kRed);




  TString savename = Form("plots/nu-1d-regions-ratio_%d-bcksub_%d-subqelike_%d-fit_%d-region_%s.pdf", doRatio, doBckSub, subQELike, hasFitting, tag.c_str());
  c->Print( savename );
  //Creating new legend canvas
  TCanvas *ccanvas;
  if( nColumns == 1 ) ccanvas = new TCanvas("ccanvas","ccanvas",100,200);
  else ccanvas = new TCanvas("ccanvas","ccanvas",200,150);
  ccanvas->cd();
  leg->SetX1NDC(.0);
  leg->SetY1NDC(.0);
  leg->SetX2NDC(1);
  leg->SetY2NDC(1);
  leg->Draw();
  ccanvas->Print(Form("plots/one_hist_leg_ratio_%d-bcksub_%d-subqelike_%d.pdf", doRatio, doBckSub, subQELike));
  delete c;
  delete ccanvas;
}

void drawSystematics( CCQENuPlotUtils*utils, MnvH1D* h,  double pot_data, double pot_mc, string name, bool asfrac)
{
  if(!hasSystematics) return;
  h->GetVertErrorBand("Low_Recoil_2p2h_Tune")->SetUnivWgt(2,0.5);
  h->GetXaxis()->SetRangeUser(0.05,5);
  MnvPlotter *plotter = new MnvPlotter;
  plotter->ApplyStyle(PlotUtils::kCCQENuStyle);
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

  //plotter->axis_maximum = 1.;
  if( !asfrac )plotter->axis_maximum = 1000;
  plotter->headroom = 1.0;
  plotter->legend_n_columns = 2;
  plotter->legend_text_size = 0.02;
  plotter->AddPlotLabel( Form("%s #bullet%s", "#Q^{2}" ,  "Errors"), 0.33, 0.94, 0.025 );

  vector<string> sysGroupNames({"Muon Reconstruction", "Recoil Reconstruction", "Low Recoil Fit", "XSection Models", "FSI Models", "Others", "Flux"});
  vector<string> sysGroupPrintNames({"Muon_Reconstruction", "Recoil_Reconstruction", "Low_Recoil_Fit", "XSection_Models", "FSI_Models", "Others", "Flux"});
  for( int i = 0; i < sysGroupPrintNames.size(); i++ )
  {
    //plotter->axis_maximum = 1000;
    cout<<sysGroupPrintNames[i]<<endl;
    utils->writeNorm( false, pot_data, pot_mc, true ); 
    plotter->DrawErrorSummary(h, "TL", false, true, 0.00001, false, sysGroupNames[i],asfrac);
    c->Print(Form("plots/sys_%s_%s.pdf",name.c_str(), sysGroupPrintNames[i].c_str() ) );
  }

    //plotter->axis_maximum = 1000;
    plotter->AddPlotLabel( Form("%s #bullet%s", "#Q^{2}" ,  "Errors"), 0.33, 0.94, 0.025 );
    utils->writeNorm( false, pot_data, pot_mc, true ); 
    plotter->DrawErrorSummary(h, "TL", true, true, 0.00001, false, "", asfrac);
    c->Print( Form("plots/sys_%s_all.pdf",name.c_str()) );
    delete plotter;
    delete c;
}

void drawLatSystematicHist( MnvH1D* h,  string sysname, string savename ) 
{
  h->GetXaxis()->SetRangeUser(0.05,5);
  MnvLatErrorBand *eb = h->GetLatErrorBand( sysname );
  TH1D* h0 = eb->GetHist(0);
  TH1D* h1 = eb->GetHist(1);
  TCanvas c;
  c.cd();
  eb->Draw("pe");
  h0->Draw("histsame");
  h1->Draw("histsame");
  TLegend* leg = new TLegend(.7,.7,.9,.9);
  leg->AddEntry(eb, "eb");
  leg->AddEntry(h0, "e0");
  leg->AddEntry(h1, "e1");
  leg->Draw();
  c.Print(Form("plots/0___lat_%s_%s.pdf", sysname.c_str(), savename.c_str() ));
}

void drawVertSystematicHist( MnvH1D* h,  string sysname, string savename ) 
{
  h->GetXaxis()->SetRangeUser(0.05,5);
  MnvVertErrorBand *eb = h->GetVertErrorBand( sysname );
  int nHists = eb->GetNHists();

  TLegend* leg = new TLegend(.7,.7,.9,.9);
  TCanvas c;
  c.cd();
  leg->AddEntry(eb, "CV hist");
  eb->Draw("pe");

  for( int i = 0; i< nHists; i++ )
  {
    TH1D* h = eb->GetHist(i);
    h->Draw("histsame");
    leg->AddEntry(h, Form("Universe %d",i) );

  }
  c.SetLogx();
  leg->Draw();
  c.Print(Form("plots/0___vert_%s_%s.pdf", sysname.c_str(), savename.c_str() ));
}

//======================= End Region Plots ================================//
//======================= Draw Fits =======================================//

void DrawMCWithErrorBand( TCanvas *c1, MnvH1D* h1d, string title,string xtitle, string ytitle, string savename, bool isRatio=true, bool drawSys = true )
{
  c1->cd();

  TH1D h = (drawSys)? h1d->GetCVHistoWithError(): h1d->GetCVHistoWithStatError();
  TH1D herr = (drawSys)? h1d->GetCVHistoWithError(): h1d->GetCVHistoWithStatError();
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
  (drawSys)?  latex.DrawLatexNDC(0.5,0.8, "systematic error"):latex.DrawLatexNDC(0.5,0.8, "fitting error");
  c1->SetLogx();
  c1->SetGrid();
  c1->Print( savename.c_str() );

}
void DrawFits( TFile *f, bool drawSys = true )
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
  for( UInt_t i = 0; i<categories_to_fit.size(); i++ ) 
  {
    fits[i]->GetXaxis()->SetTitle("Reconstructed Q^{2}_{QE}");
    DrawMCWithErrorBand ( c1, fits[i], categories_to_fit_title_names[i], "Q^{2}_{QE} (GeV^{2})", "weight", Form("plots/rw_%s.pdf", categories_to_fit_names[i].c_str()), true,drawSys );
  }
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
  vector<TH1D*> ret(nHistos2,NULL );
  for( UInt_t i = 0; i< nHistos2; i++ )
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
  vector< vector<MnvH1D*> > h1_all( nRegions, vector<MnvH1D*>(nHistos2,NULL) );
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
    latex->DrawLatexNDC(latex_x, latex_y, Form("%.2f<Q^{2}_{QE}.2f",bin_edges[i], bin_edges[i+1]));
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


void FormatPartType( MnvH2D* h, string xtitle )
{
  h->GetYaxis()->SetLabelSize(0.1);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitle(xtitle.c_str() );
  h->GetXaxis()->SetTitleSize(0.05);

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

  MnvH2D *h_recoil_inc_q2qe[nHistos2];
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

  MnvH2D* h_blobMaxE_PartType[nHistos2], *h_blobE_PartType[nHistos2];
  MnvH2D* h_mainBlobMaxE_PartType[nHistos2], *h_mainBlobE_PartType[nHistos2];

  utils->bookHistos( f1, h_blobMaxE_PartType, "h_blobMaxE_PartType" );
  utils->bookHistos( f1, h_blobE_PartType, "h_blobE_PartType" );
  utils->bookHistos( f1, h_mainBlobMaxE_PartType, "h_mainBlobMaxE_PartType" );
  utils->bookHistos( f1, h_mainBlobE_PartType, "h_mainBlobE_PartType" );

  FormatPartType( h_blobMaxE_PartType[kMC], "Leading Cluster E (MeV)" );
  FormatPartType( h_mainBlobMaxE_PartType[kMC], "Leading Cluster E (MeV)" );
  FormatPartType( h_blobE_PartType[kMC], "Blob Total E (MeV)" );
  FormatPartType( h_blobE_PartType[kQELike_QE_OTH], "Blob Total E (MeV)" );
  FormatPartType( h_blobE_PartType[kMC], "Blob Total E (MeV)" );
  FormatPartType( h_mainBlobE_PartType[kMC], "Blob Total E (MeV)" );
  setBlackbodyPalette();
  TCanvas *c = new TCanvas();
  c->cd();
  h_blobMaxE_PartType[kMC]->Draw("colz");
  c->Print("./plots/particles-all-blobMaxE.pdf");

  h_blobE_PartType[kMC]->Draw("colz");
  c->Print("./plots/particles-all-blobE.pdf");
  h_blobE_PartType[kMC]->GetXaxis()->SetRangeUser(0,10);
  c->Print("./plots/particles-all-blobE_10MeV.pdf");




  h_blobE_PartType[kQELike_QE_OTH]->Draw("colz");
  c->Print("./plots/particles-qelikeoth-blobE.pdf");
  h_blobE_PartType[kQELike_QE_OTH]->GetXaxis()->SetRangeUser(0,10);
  c->Print("./plots/particles-qelikeoth-blobE_10MeV.pdf");



  h_mainBlobMaxE_PartType[kMC]->Draw("colz");
  c->Print("./plots/particles-mainCandidate-blobMaxE.pdf");
  h_mainBlobE_PartType[kMC]->Draw("colz");
  c->Print("./plots/particles-mainCandidate-blobE.pdf");
  h_mainBlobE_PartType[kMC]->GetXaxis()->SetRangeUser(0,10);
  h_mainBlobE_PartType[kMC]->Draw("colz");
  c->Print("./plots/particles-mainCandidate-blobE_10MeV.pdf");



}
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
//Diagnostic Plots

void drawDiagnostic( MnvH2D** hists, string name, string xtitle, string axis="x", bool doRatio=false, double xmin = 0, double xmax=-1, bool logx = false, string ytitle="Q^{2}_{QE} (GeV^{2})", bool drawInt=false )
{
  //Distribution in each Q2
  vector<int> histsToUse({kData,kMC, kQELike, kQELikeNot, kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions});

  TH2* dataStat = (TH2*) hists[kData]->GetCVHistoWithStatError().Clone("dataStat");
  TH2* mcErr =(TH2*)  hists[kMC]->GetCVHistoWithError().Clone("mcErr");
  TH2* mc =(TH2*)  hists[kMC]->GetCVHistoWithError().Clone("mc");
  mcErr->SetFillStyle(1001);
  mcErr->SetFillColorAlpha(kRed,0.35);
  mcErr->SetMarkerSize(0 );

  vector<TH2*>histErr;
  cout<<"Creating histErr"<<endl;
  for(unsigned int j = 0; j< nHistos2;j++) histErr.push_back( NULL );
  for(auto c: histsToUse) histErr[c] = (TH2*)hists[c]->GetCVHistoWithError().Clone(Form("h_%s",names[c].c_str() ) ) ;

  //vector<MnvH2D*> vec_hists(hists, hists+nHistos2 );
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
  
  vector<TH2D*> hInts = hists[kMC]->GetVertErrorBand("Neutron_Interaction")->GetHists();
  if(drawInt)
  {
    histAndOpts.push_back( std::make_pair( hInts[0], "hist" ) );
    histAndOpts.push_back( std::make_pair( hInts[1], "hist" ) );
  }


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
  if( drawInt )
  {
    leg->AddEntry( hInts[0], "N Int 0" );
    leg->AddEntry( hInts[1], "N Int 1" );
  }


  double *multipliers = NULL;

  //First project into X axis
  //string ytitle="Q^{2}_{QE} (GeV^{2})";

  if(axis == "x")
  {
    GridCanvas* gc = plotXAxis1D( histAndOpts, xtitle, ytitle, multipliers ); 
    gc->SetGridx(true);
    gc->SetLogx(logx);
    gc->SetYLimits(-.2,mcErr->GetMaximum()*1.5);




    if( xmin<xmax) gc->SetXLimits(xmin,xmax);
    if(doRatio)gc->SetYLimits(-.2,2.2);
    else gc->SetYLimits(0, histErr[kData]->GetMaximum()*1.3 );
    gc->SetYTitle("Evt Rate #times 10^{3}");
    if(doRatio) gc->SetYTitle("Ratio to MnvGENIE");
    gc->Modified();
    leg->Draw("SAME");

    string fname=Form("plots/diagnostic-%s-axis_x-1d-ratio_%d.pdf",name.c_str(), doRatio );
    gc->Print(fname.c_str());
  }
  else
  {
    GridCanvas* gc = plotYAxis1D( histAndOpts, ytitle, xtitle, multipliers ); 
    gc->SetGridx(true);
    gc->SetLogx(logx);
    gc->SetYLimits(-.2,1.2);
    if( xmin<xmax) gc->SetXLimits(xmin,xmax);
    if(doRatio)gc->SetYLimits(-.2,2.2);
    else gc->SetYLimits(0, histErr[kData]->GetMaximum()*1.3 );
    gc->SetYTitle("Evt Rate #times 10^{3}");
    if(doRatio) gc->SetYTitle("Ratio to MnvGENIE");
    gc->Modified();
    leg->Draw("SAME");
    string fname=Form("plots/diagnostic-%s-axis_y-1d-ratio_%d.pdf",name.c_str(), doRatio );
    gc->Print(fname.c_str());
  }

  free(dataStat);
  free(mcErr);
  free(mc);
  for(unsigned int h = 0; h< nHistos2;h++) free(histErr[h]);

}

//void drawDiagnostic2( MnvH2D** hists, string name, string xtitle, string axis="x", bool doRatio=false, double xmin = 0, double xmax=-1, bool logx = false )
//{
//  //Distribution in each Q2
//  vector<int> histsToUse({kData,kMC, kQELike, kQELikeNot, kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions});
//
//  vector<MnvH2D*> vhists(nHistos2, NULL);
//  for( auto h: histsToUse ) vhists[h] = (MnvH2D*) hists[h]->Clone( Form( "%s_cl", hists[h]->GetName() ) );
//
//  MnvH2D* data = (MnvH2D*) hists[kData]->Clone("data");
//  MnvH2D* mcErr = (MnvH2D*) hists[kMC]->Clone("mcErr");
//  mcErr->SetFillStyle(1001);
//  mcErr->SetFillColorAlpha(kRed,0.35);
//  mcErr->SetMarkerSize(0 );
//
//  vector<TH2*>histErr;
//  cout<<"Creating histErr"<<endl;
//  for(unsigned int j = 0; j< nHistos2;j++) histErr.push_back( NULL );
//  for(auto c: histsToUse) histErr[c] = (TH2*)hists[c]->GetCVHistoWithError().Clone(Form("h_%s",names[c].c_str() ) ) ;
//
//  //vector<MnvH2D*> vec_hists(hists, hists+nHistos2 );
//  FormatHistos( histErr );
//
//  cout<<"Creating histAndOpts"<<endl;
//  vector<pair<TH2*, const char*> > histAndOpts;
//  histAndOpts.push_back( std::make_pair( mcErr, "e2" ) );
//  histAndOpts.push_back( std::make_pair( histErr[kMC], "hist" ) );
//  histAndOpts.push_back( std::make_pair( histErr[kQELike_QE_H], "histl" ) );
//  histAndOpts.push_back( std::make_pair( histErr[kQELike_QE_OTH], "histl" ) );
//  histAndOpts.push_back( std::make_pair( histErr[kQELike_RES], "histl" ) );
//  histAndOpts.push_back( std::make_pair( histErr[kQELike_2p2h], "histl" ) );
//  histAndOpts.push_back( std::make_pair( histErr[kQELike_DIS], "histl" ) );
//  histAndOpts.push_back( std::make_pair( histErr[kQELikeNot], "hist" ) );
//  histAndOpts.push_back( std::make_pair( histErr[kData], "ep" ) );
//
//
//  for( unsigned int j = 0; j< histAndOpts.size(); j++ )
//  {
//    if( doRatio ) 
//    {
//      histAndOpts[j].first->Divide( (histAndOpts[j].first), mc );
//    }
//    else histAndOpts[j].first->Scale(1);
//  }
//
//  //Define Legend First:
//  TLegend* leg=new TLegend(0.9, 0.1, 1, 0.3);
//  leg->SetFillStyle(0);
//  leg->SetBorderSize(0);
//  leg->SetTextSize(0.03);
//  leg->AddEntry( dataStat, "MINERvA data", "lpe");
//  leg->AddEntry( mc, "MINERvA Tune", "l");
//  leg->AddEntry( histErr[kQELike_QE_H],"QE-H","l");
//  leg->AddEntry( histErr[kQELike_QE_OTH],"QE-Oth","l");
//  leg->AddEntry( histErr[kQELike_RES],"Resonant","l");
//  leg->AddEntry( histErr[kQELike_2p2h],"2p2h","l");
//  leg->AddEntry( histErr[kQELike_DIS],"DIS","l");
//
//
//  double *multipliers = NULL;
//
//  //First project into X axis
//  string ytitle="Q^{2}_{QE} (GeV^{2})";
//
//  if(axis == "x")
//  {
//    GridCanvas* gc = plotXAxis1D( histAndOpts, xtitle, ytitle, multipliers ); 
//    gc->SetGridx(true);
//    gc->SetLogx(logx);
//    gc->SetYLimits(-.2,mcErr->GetMaximum()*1.5);
//    if( xmin<xmax) gc->SetXLimits(xmin,xmax);
//    if(doRatio)gc->SetYLimits(-.2,2.2);
//    else gc->SetYLimits(0, histErr[kData]->GetMaximum()*1.3 );
//    gc->SetYTitle("Evt Rate #times 10^{3}");
//    if(doRatio) gc->SetYTitle("Ratio to MnvGENIE");
//    gc->Modified();
//    leg->Draw("SAME");
//
//    string fname=Form("plots/diagnostic-%s-axis_x-1d-ratio_%d.pdf",name.c_str(), doRatio );
//    gc->Print(fname.c_str());
//  }
//  else
//  {
//    GridCanvas* gc = plotYAxis1D( histAndOpts, ytitle, xtitle, multipliers ); 
//    gc->SetGridx(true);
//    gc->SetLogx(logx);
//    gc->SetYLimits(-.2,1.2);
//    if(doRatio)gc->SetYLimits(-.2,2.2);
//    else gc->SetYLimits(0, histErr[kData]->GetMaximum()*1.3 );
//    gc->SetYTitle("Evt Rate #times 10^{3}");
//    if(doRatio) gc->SetYTitle("Ratio to MnvGENIE");
//    gc->Modified();
//    leg->Draw("SAME");
//    string fname=Form("plots/diagnostic-%s-axis_y-1d-ratio_%d.pdf",name.c_str(), doRatio );
//    gc->Print(fname.c_str());
//  }
//
//  free(dataStat);
//  free(mcErr);
//  free(mc);
//  for(unsigned int h = 0; h< nHistos2;h++) free(histErr[h]);
//
//}


void NeutronDiagnostics( TFile *f_orig, TFile *f_signal_weighted, string tag = "" )
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

  MnvH2D *h_blobDist_q2qe[nHistos2], *h_blobEnergy_q2qe[nHistos2], *h_nBlobs_q2qe[nHistos2],*h_n3DBlobs_q2qe[nHistos2], *h_n2DBlobs_q2qe[nHistos2];
  MnvH2D *h_blobEnergyLow_q2qe[nHistos2];
  MnvH2D *h_muonTheta_q2qe[nHistos2];
  MnvH2D *h_muonPhi_q2qe[nHistos2];
  MnvH2D *h_muonE_q2qe[nHistos2];
  MnvH2D *h_protonDEDX_q2qe[nHistos2],*h_pionDEDX_q2qe[nHistos2];
  MnvH2D *h_nonBlobEnergy_q2qe[nHistos2];
  MnvH2D *h_hasTrack_q2qe[nHistos2];
  MnvH2D *h_blobMaxE_q2qe[nHistos2];
  MnvH2D *h_nClus_q2qe[nHistos2];
  MnvH2D *h_enu_q2qe[nHistos2];
  utils->bookHistos( f1, h_blobEnergy_q2qe, "h_blobEnergy_q2qe"+tag );
  utils->bookHistos( f1, h_blobEnergyLow_q2qe, "h_blobEnergyLow_q2qe"+tag );
  utils->bookHistos( f1, h_nonBlobEnergy_q2qe, "h_nonBlobEnergy_q2qe"+tag );
  utils->bookHistos( f1, h_blobDist_q2qe, "h_blobDist_q2qe"+tag );
  utils->bookHistos( f1, h_nBlobs_q2qe, "h_nBlobs_q2qe"+tag );
  utils->bookHistos( f1, h_n3DBlobs_q2qe, "h_n3DBlobs_q2qe"+tag );
  utils->bookHistos( f1, h_n2DBlobs_q2qe, "h_n2DBlobs_q2qe"+tag );
  utils->bookHistos( f1, h_muonTheta_q2qe, "h_muonTheta_q2qe"+tag );
  utils->bookHistos( f1, h_muonPhi_q2qe, "h_muonPhi_q2qe"+tag );
  utils->bookHistos( f1, h_muonE_q2qe, "h_muonE_q2qe"+tag );
  utils->bookHistos( f1, h_blobMaxE_q2qe, "h_blobMaxE_q2qe"+tag );
  utils->bookHistos( f1, h_nClus_q2qe, "h_nClus_q2qe"+tag );
  utils->bookHistos( f1, h_enu_q2qe, "h_enu_q2qe"+tag );

  drawDiagnostic( h_blobEnergy_q2qe, "blobE"+tag, "E (MeV)", "x", false, 0,-1, false, "Q^{2}_{QE} (GeV^{2})", false );
  //drawDiagnostic( h_blobEnergy_q2qe, "blobEInt"+tag, "E (MeV)", "x", false, 0,-1, false, "Q^{2}_{QE} (GeV^{2})", true );
  //drawDiagnostic( h_blobEnergyLow_q2qe, "blobLowE"+tag, "E (MeV)", "x", false );
  //drawDiagnostic( h_nonBlobEnergy_q2qe, "nonBlobE"+tag, "E (MeV)", "x", false );
  //drawDiagnostic( h_blobDist_q2qe, "blobDist"+tag, "R (mm)", "x" , false, 0,2000);
  //drawDiagnostic( h_n3DBlobs_q2qe, "n3DBlobs"+tag, "n 3DBlobs", "x" , false);
  //drawDiagnostic( h_n2DBlobs_q2qe, "n2DBlobs"+tag, "n 2DBlobs", "x" , false);
  //drawDiagnostic( h_nBlobs_q2qe, "nBlobs"+tag, "nBlobs", "x" , false);
  //drawDiagnostic( h_nBlobs_q2qe, "nBlobs"+tag, "nBlobs", "y" , false);
  //drawDiagnostic( h_muonTheta_q2qe, "muonTheta"+tag, "#theta_{#mu} (degree)", "x" , false);
  //drawDiagnostic( h_muonPhi_q2qe, "muonPhi"+tag, "#phi_{#mu} (degree)", "x" , false);
  //drawDiagnostic( h_muonE_q2qe, "muonE"+tag, "E_{#mu} (GeV)", "x" , false);
  //drawDiagnostic( h_enu_q2qe, "enu"+tag, "E_{#nu} (GeV)", "x" , false);


  //drawDiagnostic( h_blobMaxE_q2qe, "blobMaxE"+tag, "E_{blob} (GeV)", "x", false);
  //drawDiagnostic( h_nClus_q2qe, "nClus"+tag, "#clusters", "x", false );

  //drawDiagnostic( h_blobEnergy_q2qe, "blobE"+tag, "E (MeV)", "x", true );

  drawDiagnostic( h_blobEnergy_q2qe, "blobE"+tag, "E (MeV)", "x", true, 0,-1, false, "Q^{2}_{QE} (GeV^{2})", false );
  //drawDiagnostic( h_blobEnergy_q2qe, "blobEInt"+tag, "E (MeV)", "x", true, 0,-1, false, "Q^{2}_{QE} (GeV^{2})", true );
  //drawDiagnostic( h_blobEnergyLow_q2qe, "blobLowE"+tag, "E (MeV)", "x", true );
  //drawDiagnostic( h_nonBlobEnergy_q2qe, "nonBlobE"+tag, "E (MeV)", "x", true );
  //drawDiagnostic( h_blobDist_q2qe, "blobDist"+tag, "R (mm)", "x" , true, 0,2000);
  //drawDiagnostic( h_n3DBlobs_q2qe, "n3DBlobs"+tag, "n 3DBlobs", "x" , true);
  //drawDiagnostic( h_n2DBlobs_q2qe, "n2DBlobs"+tag, "n 2DBlobs", "x" , true);
  //drawDiagnostic( h_nBlobs_q2qe, "nBlobs"+tag, "nBlobs", "x" , true);
  //drawDiagnostic( h_nBlobs_q2qe, "nBlobs"+tag, "nBlobs", "y" , true);
  //drawDiagnostic( h_muonTheta_q2qe, "muonTheta"+tag, "#theta_{#mu} (degree)", "x" , true);
  //drawDiagnostic( h_muonPhi_q2qe, "muonPhi"+tag, "#phi_{#mu} (degree)", "x" , true);
  //drawDiagnostic( h_muonE_q2qe, "muonE"+tag, "E_{#mu} (GeV)", "x" , true);
  //drawDiagnostic( h_enu_q2qe, "enu"+tag, "E_{#nu} (GeV)", "x" , true);

  //drawDiagnostic( h_blobMaxE_q2qe, "blobMaxE"+tag, "E_{leading cluster} (GeV)", "x", true);
  //drawDiagnostic( h_nClus_q2qe, "nClus"+tag, "#clusters", "x", true );

  vector<MnvH2D*> fit_blobE = GetFitsHisto( f_signal_weighted, "blobEnergy" );
  //vector<MnvH2D*> fit_blobDist = GetFitsHisto( f_signal_weighted, "blobDist" );
  //vector<MnvH1D*> fit_q2 = GetFitsHisto( f_signal_weighted );
  //vector<MnvH2D*> fit_nblob = GetFitsHisto( f_signal_weighted, "nBlobs" );

  DoFit( h_blobEnergy_q2qe, fit_blobE );
  //DoFit( h_blobEnergyLow_q2qe, fit_blobE );
  //DoFit( h_blobDist_q2qe, fit_blobDist );
  //DoFit( h_nBlobs_q2qe, fit_q2 );
  //DoFit( h_n2DBlobs_q2qe, fit_q2 );
  //DoFit( h_n3DBlobs_q2qe, fit_q2 );
  //DoFit( h_blobMaxE_q2qe, fit_q2 );
  //DoFit( h_nClus_q2qe, fit_q2 );
  //DoFit( h_enu_q2qe, fit_q2 );

  drawDiagnostic( h_blobEnergy_q2qe, "blobE-fit_1"+tag, "E (MeV)", "x", false );
  //drawDiagnostic( h_blobEnergyLow_q2qe, "blobLowE-fit_1"+tag, "E (MeV)", "x", false );
  //drawDiagnostic( h_blobDist_q2qe, "blobDist-fit_1"+tag, "R (mm)", "x" , false, 0,2000);
  //drawDiagnostic( h_n3DBlobs_q2qe, "n3DBlobs-fit_1"+tag, "n 3DBlobs", "x" , false);
  //drawDiagnostic( h_n2DBlobs_q2qe, "n2DBlobs-fit_1"+tag, "n 2DBlobs", "x" , false);
  //drawDiagnostic( h_nBlobs_q2qe, "nBlobs-fit_1"+tag, "n Blobs", "x" , false);

  //drawDiagnostic( h_blobMaxE_q2qe, "blobMaxE-fit_1"+tag, "E_{leading cluster} (MeV)","x",false );
  //drawDiagnostic( h_nClus_q2qe, "nClus-fit_1"+tag, "#clusters", "x",false );

  //drawDiagnostic( h_enu_q2qe, "enu-fit_1"+tag, "E_{#nu} (GeV)", "x",false );

  drawDiagnostic( h_blobEnergy_q2qe, "blobE-fit_1"+tag, "E (MeV)", "x", true );
  //drawDiagnostic( h_blobEnergyLow_q2qe, "blobLowE-fit_1"+tag, "E (MeV)", "x", true );
  //drawDiagnostic( h_blobDist_q2qe, "blobDist-fit_1"+tag, "R (mm)", "x" , true, 0,2000);
  //drawDiagnostic( h_n3DBlobs_q2qe, "n3DBlobs-fit_1"+tag, "n 3DBlobs", "x" , true);
  //drawDiagnostic( h_n2DBlobs_q2qe, "n2DBlobs-fit_1"+tag, "n 2DBlobs", "x" , true);
  //drawDiagnostic( h_nBlobs_q2qe, "nBlobs-fit_1"+tag, "n Blobs", "x" , true);

  //drawDiagnostic( h_blobMaxE_q2qe, "blobMaxE-fit_1"+tag, "E_{leading cluster} (MeV)","x",true);
  //drawDiagnostic( h_nClus_q2qe, "nClus-fit_1"+tag, "#clusters", "x",true );
  //drawDiagnostic( h_enu_q2qe, "enu-fit_1"+tag, "E_{#nu} (GeV)", "x",true );



  return;
}



MnvH2D* ComponentToDataFit( MnvH2D** h, int iHisto )
{
  MnvH2D* data = (MnvH2D*) h[kData]->Clone("data");
  vector<int> component_histos({ kQELike_QE_H, kQELike_QE_OTH, kQELike_2p2h, kQELike_RES, kQELike_DIS, kQELike_OTH, kQELikeNot} );

  for( auto i: component_histos ) // subtract other component
  {
    if (i == iHisto) continue;
    data->Add( h[i],-1 );
  }
  MnvH2D* fit = (MnvH2D*) h[iHisto]->Clone( Form("%s_fit", h[iHisto]->GetName() ) );
  int xfirst =1, yfirst=1;
  int xlast = fit->GetNbinsX(), ylast = fit->GetNbinsY();
  //MnvH1D* scale_q2_1d = (MnvH1D*) fit->ProjectionY("scale1D", xfirst,xlast);
  //MnvH2D* scale_q2_2d = (MnvH2D*) fit->Clone("scale2D");


  //scale_q2_2d->ClearAllErrorBands();
  //scale_q2_2d->Reset();
  //scale_q2_1d->Divide( scale_q2_1d, (MnvH1D*) data->ProjectionY("dataY", xfirst,xlast) );//scale data to mc


  //ExpandHisto(scale_q2_1d, scale_q2_2d, 1 );//expand y into x

  //data->Multiply(data, scale_q2_2d); //scale data to mc

  fit->Divide( data, fit );//fit mc to scaled data.

  return fit;
}

void NeutronLepton2D( TFile *f_orig, TFile *f_signal_weighted, string tag = "" )
{
  cout<<"NeutronLepton2D"<<endl;
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

  MnvH2D *h_muonTheta_q2qe[nHistos2];
  MnvH2D *h_muonPhi_q2qe[nHistos2];
  MnvH2D *h_muonE_q2qe[nHistos2];

  vector<int> regions({0,1,2,3,4,5,99});
  vector<vector<MnvH2D*>> enuHists(regions.size(), vector<MnvH2D*>(nHistos2,NULL) );
  vector<vector<MnvH2D*>> muonthetaHists(regions.size(), vector<MnvH2D*>(nHistos2,NULL) );

  utils->bookHistos( f1, h_muonTheta_q2qe, "h_muonTheta_q2qe"+tag );
  utils->bookHistos( f1, h_muonPhi_q2qe, "h_muonPhi_q2qe"+tag );
  utils->bookHistos( f1, h_muonE_q2qe, "h_muonE_q2qe"+tag );

  for( unsigned int i = 0; i< regions.size(); i++ )
  {
    int r = regions[i];
    utils->bookHistos( f1, &enuHists[i][0], Form("h_enu_q2qe_region_%02d", r) );
    utils->bookHistos( f1, &muonthetaHists[i][0], Form("h_muontheta_q2qe_region_%02d", r) );
  }

  drawDiagnostic( h_muonTheta_q2qe, "muonTheta"+tag, "#theta_{#mu} (degree)", "x" , false);
  drawDiagnostic( h_muonPhi_q2qe, "muonPhi"+tag, "#phi_{#mu} (degree)", "x" , false);
  drawDiagnostic( h_muonE_q2qe, "muonE"+tag, "E_{#mu} (GeV)", "x" , false);

  drawDiagnostic( h_muonTheta_q2qe, "muonTheta"+tag, "#theta_{#mu} (degree)", "x" , true);
  drawDiagnostic( h_muonPhi_q2qe, "muonPhi"+tag, "#phi_{#mu} (degree)", "x" , true);
  drawDiagnostic( h_muonE_q2qe, "muonE"+tag, "E_{#mu} (GeV)", "x" , true);


  //drawDiagnostic( h_muonTheta_q2qe, "q2qeMuonTheta"+tag, "Q^{2}_{QE} (GeV^{2})", "y" , false, 0, -1, true);
  //drawDiagnostic( h_muonTheta_q2qe, "q2qeMuonTheta"+tag, "Q^{2}_{QE} (GeV^{2})", "y" , true, 0, -1, true);

  //drawDiagnostic( h_enu_q2qe, "q2qeEnu"+tag, "Q^{2}_{QE} (GeV^{2})", "y" , false,0.15,50,true);
  //drawDiagnostic( h_enu_q2qe, "q2qeEnu"+tag, "Q^{2}_{QE} (GeV^{2})", "y" , true,0.15,50,true);

  vector<MnvH1D*> fit_q2 = GetFitsHisto( f_signal_weighted );
  for( unsigned int i = 0; i< regions.size(); i++ )
  {
    //continue;
    int r = regions[i];
    drawDiagnostic( &enuHists[i][0], Form("enu-region_%02d%s",r,tag.c_str()), "E_{#nu} (GeV)", "x" , false,0.15,50,true);
    drawDiagnostic( &enuHists[i][0], Form("enu-region_%02d%s",r,tag.c_str()), "E_{#nu} (GeV)", "x" , true,0.15,50,true);

    drawDiagnostic( &enuHists[i][0], Form("q2qeEnu-region_%02d%s",r,tag.c_str()), "Q^{2}_{QE} (GeV^2)", "y" , false,0.15,50,true, "E_{#nu}^{reco} (GeV)");
    drawDiagnostic( &enuHists[i][0], Form("q2qeEnu-region_%02d%s",r,tag.c_str()), "Q^{2}_{QE} (GeV^2)", "y" , true,0.15,50,true, "E_{#nu}^{reco} (GeV)");

    DoFit( &enuHists[i][0], fit_q2 );


    drawDiagnostic( &muonthetaHists[i][0], Form("muontheta-region_%02d%s",r,tag.c_str()), "#theta_{#mu} (degree)", "x" , false,0.15,50,true);
    drawDiagnostic( &muonthetaHists[i][0], Form("muontheta-region_%02d%s",r,tag.c_str()), "#theta_{#mu} (degree)", "x" , true,0.15,50,true);

    drawDiagnostic( &muonthetaHists[i][0], Form("q2qeMuonTheta-region_%02d%s",r,tag.c_str()), "Q^{2}_{QE} (GeV^2)", "y" , false,0.15,50,true, "#theta_{#mu}");
    drawDiagnostic( &muonthetaHists[i][0], Form("q2qeMuonTheta-region_%02d%s",r,tag.c_str()), "Q^{2}_{QE} (GeV^2)", "y" , true,0.15,50,true, "#theta_{#mu}");

    DoFit( &muonthetaHists[i][0], fit_q2 );

  }


  vector<MnvH2D*> enuHistsC0 = CloneVector( enuHists[0], true );
  vector<MnvH2D*> enuHistsC1 = CloneVector( enuHists[1], true );
  vector< vector<MnvH2D*> > enuHistsC;
  enuHistsC.push_back( enuHistsC0 );
  enuHistsC.push_back( enuHistsC1 );

  //DoFit( h_enu_q2qe, fit_q2 );
  DoFit( h_muonTheta_q2qe, fit_q2 );
  DoFit( h_muonPhi_q2qe, fit_q2 );
  DoFit( h_muonE_q2qe, fit_q2 );

  drawDiagnostic( h_muonTheta_q2qe, "muonTheta-fit_1"+tag, "#theta_{#mu} (degree)", "x",false );
  drawDiagnostic( h_muonPhi_q2qe, "muonPhi-fit_1"+tag, "#phi_{#mu} (degree)", "x",false );
  drawDiagnostic( h_muonE_q2qe, "muonE-fit_1"+tag, "E_{#mu} (GeV)", "x",false );
  //drawDiagnostic( h_enu_q2qe, "enu-fit_1"+tag, "E_{#nu} (GeV)", "x",false,0.15,50,true );

  drawDiagnostic( h_muonTheta_q2qe, "muonTheta-fit_1"+tag, "#theta_{#mu} (degree)", "x",true );
  drawDiagnostic( h_muonPhi_q2qe, "muonPhi-fit_1"+tag, "#phi_{#mu} (degree)", "x",true );
  drawDiagnostic( h_muonE_q2qe, "muonE-fit_1"+tag, "E_{#mu} (GeV)", "x",true );
  //drawDiagnostic( h_enu_q2qe, "enu-fit_1"+tag, "E_{#nu} (GeV)", "x",true,0.15,50,true );


  drawDiagnostic( h_muonTheta_q2qe, "q2qeMuonTheta-fit_1"+tag, "Q^{2}_{QE} (GeV^{2})", "y" , false, 0, -1, true);
  drawDiagnostic( h_muonTheta_q2qe, "q2qeMuonTheta-fit_1"+tag, "Q^{2}_{QE} (GeV^{2})", "y" , true, 0, -1, true);




  for( unsigned int i = 0; i< regions.size(); i++ )
  {
    //continue;
    int r = regions[i];
    drawDiagnostic( &enuHists[i][0], Form("enu-region_%02d-fit_1%s",r,tag.c_str()), "E_{#nu} (GeV)", "x" , false,0.15,50,true);
    drawDiagnostic( &enuHists[i][0], Form("enu-region_%02d-fit_1%s",r,tag.c_str()), "E_{#nu} (GeV)", "x" , true,0.15,50,true);
    drawDiagnostic( &enuHists[i][0], Form("q2qeEnu-region_%02d-fit_1%s",r,tag.c_str()), "E_{#nu}^{reco} (GeV)", "y" , false,0.15,50,true, "Q^{2}_{QE} (GeV^2)");
    drawDiagnostic( &enuHists[i][0], Form("q2qeEnu-region_%02d-fit_1%s",r,tag.c_str()), "E_{#nu}^{reco} (GeV)", "y" , true,0.15,50,true,  "Q^{2}_{QE} (GeV^2)");

    drawDiagnostic( &muonthetaHists[i][0], Form("muontheta-region_%02d-fit_1%s",r,tag.c_str()), "#theta_{#mu} (degree)", "x" , false,0.15,50,true);
    drawDiagnostic( &muonthetaHists[i][0], Form("muontheta-region_%02d-fit_1%s",r,tag.c_str()), "#theta_{#mu} (degree)", "x" , true,0.15,50,true);

    drawDiagnostic( &muonthetaHists[i][0], Form("q2qeMuonTheta-region_%02d-fit_1%s",r,tag.c_str()),"#theta_{#mu}" , "y" , false,0.15,50,true, "Q^{2}_{QE} (GeV^2)");
    drawDiagnostic( &muonthetaHists[i][0], Form("q2qeMuonTheta-region_%02d-fit_1%s",r,tag.c_str()),"#theta_{#mu}" , "y" , true,0.15,50,true,  "Q^{2}_{QE} (GeV^2)");

    SubtractBackground(&enuHists[i][0], true);
    SubtractBackground(&muonthetaHists[i][0], true);

    drawDiagnostic( &enuHists[i][0], Form("enu-region_%02d-fit_1-nobck_1%s",r,tag.c_str()), "E_{#nu} (GeV)", "x" , false,0.15,50,true);
    drawDiagnostic( &enuHists[i][0], Form("enu-region_%02d-fit_1-nobck_1%s",r,tag.c_str()), "E_{#nu} (GeV)", "x" , true,0.15,50,true);
    drawDiagnostic( &enuHists[i][0], Form("q2qeEnu-region_%02d-fit_1-nobck_1%s",r,tag.c_str()), "E_{#nu}^{reco} (GeV)", "y" , false,0.15,50,true, "Q^{2}_{QE} (GeV^2)");
    drawDiagnostic( &enuHists[i][0], Form("q2qeEnu-region_%02d-fit_1-nobck_1%s",r,tag.c_str()), "E_{#nu}^{reco} (GeV)", "y" , true,0.15,50,true,  "Q^{2}_{QE} (GeV^2)");

    drawDiagnostic( &muonthetaHists[i][0], Form("muontheta-region_%02d-fit_1-nobck_1%s",r,tag.c_str()), "#theta_{#mu} (degree)", "x" , false,0.15,50,true);
    drawDiagnostic( &muonthetaHists[i][0], Form("muontheta-region_%02d-fit_1-nobck_1%s",r,tag.c_str()), "#theta_{#mu} (degree)", "x" , true,0.15,50,true);

    drawDiagnostic( &muonthetaHists[i][0], Form("q2qeMuonTheta-region_%02d-fit_1-nobck_1%s",r,tag.c_str()), "#theta_{#mu}", "y" , false,0.15,50,true, "Q^{2}_{QE} (GeV^2)");
    drawDiagnostic( &muonthetaHists[i][0], Form("q2qeMuonTheta-region_%02d-fit_1-nobck_1%s",r,tag.c_str()), "#theta_{#mu}", "y" , true,0.15,50,true,  "Q^{2}_{QE} (GeV^2)");


  }

  MnvH2D* enu01_fit = ComponentToDataFit( &enuHistsC[1][0], kQELike_QE_OTH );

  TCanvas *c = new TCanvas("c1","c1");
  enu01_fit->Draw("colz");
  c->Print("plots/enu01_fit.pdf");
  enuHistsC[1][0]->Draw("colz");
  c->Print("plots/enu01_data.pdf");
  enuHistsC[1][kQELike_QE_OTH]->Draw("colz");
  c->Print("plots/enu01_qelike_qe_oth.pdf");


  enuHistsC[0][kQELike_QE_OTH]->MultiplySingle(enuHistsC[0][kQELike_QE_OTH], (TH2*) enu01_fit );
  enuHistsC[1][kQELike_QE_OTH]->MultiplySingle(enuHistsC[1][kQELike_QE_OTH], (TH2*) enu01_fit );


  for( int i = 0; i<2; i++ )
  {

    int r = regions[i];
    drawDiagnostic( &enuHistsC[i][0], Form("fittest-enu%02d-fit_0-scaleEnu%s",r,tag.c_str()), "E_{#nu} (GeV)", "x" , false,0.15,50,true);
    drawDiagnostic( &enuHistsC[i][0], Form("fittest-enu%02d-fit_0-scaleEnu%s",r,tag.c_str()), "E_{#nu} (GeV)", "x" , true,0.15,50,true);
    drawDiagnostic( &enuHistsC[i][0], Form("fittest-q2qeEnu%02d-fit_0-scaleEnu%s",r,tag.c_str()), "E_{#nu} (GeV)", "y" , false,0.15,50,true,"E_{#nu}^{reco} (GeV)" );
    drawDiagnostic( &enuHistsC[i][0], Form("fittest-q2qeEnu%02d-fit_0-scaleEnu%s",r,tag.c_str()), "E_{#nu} (GeV)", "y" , true,0.15,50,true,"E_{#nu}^{reco} (GeV)" );

    enuHistsC[i][kQELike_QE_OTH]->MultiplySingle(enuHistsC[i][kQELike_QE_OTH], (TH2*) enu01_fit );
    ReconfigureCategories( &enuHistsC[i][0] );
    drawDiagnostic( &enuHistsC[i][0], Form("fittest-enu%02d-fit_1-scaleEnu%s",r,tag.c_str()), "E_{#nu} (GeV)", "x" , false,0.15,50,true);
    drawDiagnostic( &enuHistsC[i][0], Form("fittest-enu%02d-fit_1-scaleEnu%s",r,tag.c_str()), "E_{#nu} (GeV)", "x" , true,0.15,50,true);
    drawDiagnostic( &enuHistsC[i][0], Form("fittest-q2qeEnu%02d-fit_1-scaleEnu%s",r,tag.c_str()), "E_{#nu} (GeV)", "y" , false,0.15,50,true,"E_{#nu}^{reco} (GeV)" );
    drawDiagnostic( &enuHistsC[i][0], Form("fittest-q2qeEnu%02d-fit_1-scaleEnu%s",r,tag.c_str()), "E_{#nu} (GeV)", "y" , true,0.15,50,true,"E_{#nu}^{reco} (GeV)" );

    SubtractBackground( &enuHistsC[i][0], true );
    drawDiagnostic( &enuHistsC[i][0], Form("fittest-enu%02d-fit_1-nobck_1-scaleEnu%s",r,tag.c_str()), "E_{#nu} (GeV)", "x" , false,0.15,50,true);
    drawDiagnostic( &enuHistsC[i][0], Form("fittest-enu%02d-fit_1-nobck_1-scaleEnu%s",r,tag.c_str()), "E_{#nu} (GeV)", "x" , true,0.15,50,true);
    drawDiagnostic( &enuHistsC[i][0], Form("fittest-q2qeEnu%02d-fit_1-nobck_1-scaleEnu%s",r,tag.c_str()), "E_{#nu} (GeV)", "y" , false,0.15,50,true,"E_{#nu}^{reco} (GeV)" );
    drawDiagnostic( &enuHistsC[i][0], Form("fittest-q2qeEnu%02d-fit_1-nobck_1-scaleEnu%s",r,tag.c_str()), "E_{#nu} (GeV)", "y" , true,0.15,50,true,"E_{#nu}^{reco} (GeV)" );

  }



  return;
}

void NeutronDiagnostics2( TFile *f_orig, TFile *f_signal_weighted, string tag = "" )
{
  cout<<"NeutronDiagnostics2"<<endl;
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

  MnvH2D *h_blobDist_q2qe[nHistos2], *h_blobEnergy_q2qe[nHistos2], *h_nBlobs_q2qe[nHistos2],*h_n3DBlobs_q2qe[nHistos2], *h_n2DBlobs_q2qe[nHistos2];
  MnvH2D *h_blobEnergyLow_q2qe[nHistos2];
  MnvH2D *h_nonBlobEnergy_q2qe[nHistos2];
  MnvH2D *h_blobMaxE_q2qe[nHistos2];
  MnvH2D *h_nClus_q2qe[nHistos2];
  utils->bookHistos( f1, h_blobEnergy_q2qe, "h_blobEnergy_q2qe"+tag );
  utils->bookHistos( f1, h_blobEnergyLow_q2qe, "h_blobEnergyLow_q2qe"+tag );
  utils->bookHistos( f1, h_nonBlobEnergy_q2qe, "h_nonBlobEnergy_q2qe"+tag );
  utils->bookHistos( f1, h_blobDist_q2qe, "h_blobDist_q2qe"+tag );
  utils->bookHistos( f1, h_nBlobs_q2qe, "h_nBlobs_q2qe"+tag );
  utils->bookHistos( f1, h_n3DBlobs_q2qe, "h_n3DBlobs_q2qe"+tag );
  utils->bookHistos( f1, h_n2DBlobs_q2qe, "h_n2DBlobs_q2qe"+tag );
  utils->bookHistos( f1, h_blobMaxE_q2qe, "h_blobMaxE_q2qe"+tag );
  utils->bookHistos( f1, h_nClus_q2qe, "h_nClus_q2qe"+tag );

  cout<<"blobMaxE: "<<h_blobMaxE_q2qe[kMC]<<endl;

  drawDiagnostic( h_blobEnergy_q2qe, "blobE"+tag, "E (MeV)", "x", false );
  //drawDiagnostic( h_blobEnergyLow_q2qe, "blobLowE"+tag, "E (MeV)", "x", false );
  //drawDiagnostic( h_nonBlobEnergy_q2qe, "nonBlobE"+tag, "E (MeV)", "x", false );
  //drawDiagnostic( h_blobDist_q2qe, "blobDist"+tag, "R (mm)", "x" , false, 0,2000);
  //drawDiagnostic( h_n3DBlobs_q2qe, "n3DBlobs"+tag, "n 3DBlobs", "x" , false);
  //drawDiagnostic( h_n2DBlobs_q2qe, "n2DBlobs"+tag, "n 2DBlobs", "x" , false);
  //drawDiagnostic( h_nBlobs_q2qe, "nBlobs"+tag, "nBlobs", "x" , false);
  //drawDiagnostic( h_nBlobs_q2qe, "nBlobs"+tag, "nBlobs", "y" , false);
  //drawDiagnostic( h_blobMaxE_q2qe, "blobMaxE"+tag, "E_{blob} (GeV)", "x", false);
  //drawDiagnostic( h_nClus_q2qe, "nClus"+tag, "#clusters", "x", false );

  drawDiagnostic( h_blobEnergy_q2qe, "blobE"+tag, "E (MeV)", "x", true );
  //drawDiagnostic( h_blobEnergyLow_q2qe, "blobLowE"+tag, "E (MeV)", "x", true );
  //drawDiagnostic( h_nonBlobEnergy_q2qe, "nonBlobE"+tag, "E (MeV)", "x", true );
  //drawDiagnostic( h_blobDist_q2qe, "blobDist"+tag, "R (mm)", "x" , true, 0,2000);
  //drawDiagnostic( h_n3DBlobs_q2qe, "n3DBlobs"+tag, "n 3DBlobs", "x" , true);
  //drawDiagnostic( h_n2DBlobs_q2qe, "n2DBlobs"+tag, "n 2DBlobs", "x" , true);
  //drawDiagnostic( h_nBlobs_q2qe, "nBlobs"+tag, "nBlobs", "x" , true);
  //drawDiagnostic( h_nBlobs_q2qe, "nBlobs"+tag, "nBlobs", "y" , true);

  //drawDiagnostic( h_blobMaxE_q2qe, "blobMaxE"+tag, "E_{leading cluster} (GeV)", "x", true);
  //drawDiagnostic( h_nClus_q2qe, "nClus"+tag, "#clusters", "x", true );

  vector<MnvH2D*> fit_blobE = GetFitsHisto( f_signal_weighted, "blobEnergy" );
  vector<MnvH2D*> fit_blobDist = GetFitsHisto( f_signal_weighted, "blobDist" );
  vector<MnvH1D*> fit_q2 = GetFitsHisto( f_signal_weighted );
  //vector<MnvH2D*> fit_nblob = GetFitsHisto( f_signal_weighted, "nBlobs" );

  DoFit( h_blobEnergy_q2qe, fit_blobE );
  //DoFit( h_blobDist_q2qe, fit_blobDist );
  //DoFit( h_nBlobs_q2qe, fit_q2 );
  //DoFit( h_n2DBlobs_q2qe, fit_q2 );
  //DoFit( h_n3DBlobs_q2qe, fit_q2 );
  //DoFit( h_blobMaxE_q2qe, fit_q2 );
  //DoFit( h_nClus_q2qe, fit_q2 );
  //DoFit( h_blobEnergyLow_q2qe, fit_q2 );

  drawDiagnostic( h_blobEnergy_q2qe, "blobE-fit_1"+tag, "E (MeV)", "x", false );
  //drawDiagnostic( h_blobEnergyLow_q2qe, "blobLowE-fit_1"+tag, "E (MeV)", "x", false );
  //drawDiagnostic( h_blobDist_q2qe, "blobDist-fit_1"+tag, "R (mm)", "x" , false, 0,2000);
  //drawDiagnostic( h_n3DBlobs_q2qe, "n3DBlobs-fit_1"+tag, "n 3DBlobs", "x" , false);
  //drawDiagnostic( h_n2DBlobs_q2qe, "n2DBlobs-fit_1"+tag, "n 2DBlobs", "x" , false);
  //drawDiagnostic( h_nBlobs_q2qe, "nBlobs-fit_1"+tag, "n Blobs", "x" , false);

  //drawDiagnostic( h_blobMaxE_q2qe, "blobMaxE-fit_1"+tag, "E_{leading cluster} (MeV)","x",false );
  //drawDiagnostic( h_nClus_q2qe, "nClus-fit_1"+tag, "#clusters", "x",false );

  drawDiagnostic( h_blobEnergy_q2qe, "blobE-fit_1"+tag, "E (MeV)", "x", true );
  //drawDiagnostic( h_blobEnergyLow_q2qe, "blobLowE-fit_1"+tag, "E (MeV)", "x", true );
  //drawDiagnostic( h_blobDist_q2qe, "blobDist-fit_1"+tag, "R (mm)", "x" , true, 0,2000);
  //drawDiagnostic( h_n3DBlobs_q2qe, "n3DBlobs-fit_1"+tag, "n 3DBlobs", "x" , true);
  //drawDiagnostic( h_n2DBlobs_q2qe, "n2DBlobs-fit_1"+tag, "n 2DBlobs", "x" , true);
  //drawDiagnostic( h_nBlobs_q2qe, "nBlobs-fit_1"+tag, "n Blobs", "x" , true);

  //drawDiagnostic( h_blobMaxE_q2qe, "blobMaxE-fit_1"+tag, "E_{leading cluster} (MeV)","x",true);
  //drawDiagnostic( h_nClus_q2qe, "nClus-fit_1"+tag, "#clusters", "x",true );



  return;
}




//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
//int main(int argc, char* argv[])
//{
//  string f_orig = argv[1];
//  string f_region = argv[2];
//  string f_signal_weighted = argv[3];
//  string sideband = argv[4];
//  nu = atoi( argv[5] );
//  cout<<f_orig<<endl;
//  if (argc == 7 ) nPlaylists = atoi(argv[6]);
//  //draw scales with errors
//  //
//  //Get Files
//  TFile* file_orig = new TFile( f_orig.c_str(),"read");
//  TFile* file_region = new TFile( f_region.c_str(),"read");
//  TFile* file_signal_weighted = new TFile( f_signal_weighted.c_str(),"read");
//
//
//  //Draw 1D
//  DrawFits( file_signal_weighted );
//  //
//  string model = "";
//  draw3DChunk( f_orig, f_signal_weighted, sideband, true, model );//do bckfitting
//  draw3DChunk( f_orig, f_signal_weighted, sideband, false , model);//not do bckfitting
//
//  //Draw regions
//  drawRegionAnglePlots( file_region, file_signal_weighted, true , sideband); //do bck fitting
//  drawRegionAnglePlots( file_region, file_signal_weighted, false, sideband ); //no bck fitting
//
//  drawRegionAnglePlots( file_region, file_signal_weighted, true , sideband, "vtx"); //do bck fitting
//  drawRegionAnglePlots( file_region, file_signal_weighted, false, sideband, "vtx" ); //no bck fitting
//
//  drawRegionAnglePlots( file_region, file_signal_weighted, true , sideband, "nonvtx"); //do bck fitting
//  drawRegionAnglePlots( file_region, file_signal_weighted, false, sideband, "nonvtx" ); //no bck fitting
//
//  draw1DDistros( file_orig, file_signal_weighted, "" );
//  PlotRecoil( file_orig, file_signal_weighted );
//
//  //Diagnostics
//  NeutronDiagnostics( file_orig, file_signal_weighted );
//  return 0;
//}
