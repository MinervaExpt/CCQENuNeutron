//#include "myPlotStyle.h"
#include "TParameter.h"

#include "include/GeneralIncludes.h"
#include "PlotUtils/HyperDimLinearizer.h"
#include "include/GeneralFunc.h"

#include "TPaletteAxis.h"


using namespace PlotUtils;
//forward declaration
//gROOT->SetBatch();

bool fluxHistoExists = true;
bool drawNuwroOnly = false;
bool hasRESFit = false;
string scaleOpt = "width";
string xvar = "ptmu";
//vector<int> categories_to_fit({ kQELike_QE_OTH, kQELike_2p2h, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion} );
//vector<string> categories_to_fit_names({"qe_oth", "2p2h","qelikenot_scp","qelikenot_snp","qelikenot_mp"});
vector<int> categories_to_fit({ kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion} );
vector<string> categories_to_fit_names({"qelikenot_scp", "qelikenot_snp", "qelikenot_mp"});
vector<string> categories_to_fit_title_names({"Single Charged Pion", "Single Neutral Pion", "Multi Pion"});


//vector<int> categories_to_fit({ kQELike_QE_OTH, kQELike_2p2h, kQELikeNot_SinglePion} );
//vector<string> categories_to_fit_names({"qelike_qe_oth", "qelike_2p2h","qelikenot_sp"});
//vector<string> categories_to_fit_title_names({"QELike && QE OTH", "QELike && 2p2h", "QELikeNot && Single Pion"});
//



vector<int> histosUsed({ kData,  kMC, kQELike_PN, kQELike_QE_PN, kQELike_2p2h_PN, kQELike_RES_PN, kQELike_DIS_PN, kQELike_OTH_PN, kQELike_NotPN, kQELikeNot, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } ); 



vector<int> histosUsedAll({ kData,  kMC,kQELike, kQELike_QE, kQELike_2p2h, kQELike_RES, kQELike_DIS, kQELike_OTH, kQELike_PN, kQELike_QE_PN, kQELike_2p2h_PN, kQELike_RES_PN, kQELike_DIS_PN, kQELike_OTH_PN, kQELike_NotPN, 
kQELike_QE_NotPN, kQELike_2p2h_NotPN, kQELike_RES_NotPN, kQELike_DIS_NotPN, kQELike_OTH_NotPN, kQELikeNot, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions, kQELike_NoProton_NotPN, kQELike_NoNeutron_NotPN, kQELike_NoNucleon_NotPN } ); 

vector<int> histosQELike({ kQELike, kQELike_QE, kQELike_2p2h, kQELike_RES, kQELike_DIS, kQELike_OTH });
vector<int> histosQELike_PN({ kQELike_PN, kQELike_QE_PN, kQELike_2p2h_PN, kQELike_RES_PN, kQELike_DIS_PN, kQELike_OTH_PN });
vector<int> histosQELike_NotPN({ kQELike_NotPN, kQELike_QE_NotPN, kQELike_2p2h_NotPN, kQELike_RES_NotPN, kQELike_DIS_NotPN, kQELike_OTH_NotPN });

vector<int> QELikeNotHistos({ kQELikeNot, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } ); 


//General Func
template<class T>
void FormatHistos(vector< T*> &h );
void FormatHistos( vector<vector<MnvH2D*>> h );
void drawSystematics( CCQENuPlotUtils*utils, MnvH1D* h,  double pot_data, double pot_mc, string name, bool asfrac=true);
void drawLatSystematicHist( MnvH1D* h,  string sysname, string savename="" ) ;


void draw1DDistro2( MnvH2D** h, vector<MnvH1D*> &fits, CCQENuPlotUtils* putils, string name, double pot_data, double pot_mc, string axis = "x", string xtitle="x", string ytitle="evt rate", int combineBin=1, bool doRatio=false, bool doFitting=false, bool doBckSub = false, bool subqelike=false, double min=-9999999, double max=-9999999, string tag="");
void draw1DDistro3( MnvH1D** h, vector<MnvH1D*> &fits, CCQENuPlotUtils* putils, string name, double pot_data, double pot_mc, string axis = "x", string xtitle="x", string ytitle="evt rate", bool doRatio=false, bool doFitting=false, bool doBckSub = false, bool subqelike=false, double min=-1., double max=-1.);

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

void ApplyFit( MnvH2D**h, vector<MnvH1D*>&fits )
{
  for( int i = 0; i< categories_to_fit.size(); i++ )
  {
    int cat = categories_to_fit[i];
    cout<<"apply fit: "<<names[cat]<<endl;
    cout<<"weight: "<<fits[i]->GetName()<<endl;
    MnvH2D* fit = (MnvH2D*) h[kData]->Clone("fit");
    cout<<"fit nx: "<<fit->GetNbinsX()<<", ny: "<<fit->GetNbinsY()<<endl;
    fit->Reset();
    fit->ClearAllErrorBands();

    ExpandHisto( fits[i], fit, 1 );
    //for( int j = 0; j< fits[i]->GetNbinsX(); j++ )
    //{
    //  cout<<fits[i]->GetBinContent(j+1);
    //  cout<<", "<<fit->GetBinContent(1, j+1)<<endl;
    //}
    h[cat]->Multiply( h[cat],fit );
    delete fit;
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
  //h[kQELike]->Reset();
  h[kQELikeNot]->Reset();


  //h[kQELike_QE]->Reset();
  //h[kQELike_QE]->Add(h[kQELike_QE_OTH]);
  //h[kQELike_QE]->Add(h[kQELike_QE_H]);
 
  //vector<int> qelike({ kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH });
  //for( auto cat: qelike ) h[kQELike]->Add( h[cat] );
  
  vector<int> qelikenot({ kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } );
  for( auto cat: qelikenot ) h[kQELikeNot]->Add( h[cat] );

  h[kMC]->Add( h[kQELike] );
  h[kMC]->Add( h[kQELikeNot] );
}
template void DoFit( MnvH1D** h, vector<MnvH1D*>&fits );
template void DoFit( MnvH2D** h, vector<MnvH2D*>&fits );



void DoFit( MnvH2D**h, vector<MnvH1D*>&fits )
{
  for( int i = 0; i< categories_to_fit.size(); i++ )
  {
    int cat = categories_to_fit[i];
    cout<<"apply fit: "<<names[cat]<<endl;
    cout<<"weight: "<<fits[i]->GetName()<<endl;
    MnvH2D* fit = (MnvH2D*) h[kData]->Clone("fit");
    cout<<"fit nx: "<<fit->GetNbinsX()<<", ny: "<<fit->GetNbinsY()<<endl;
    fit->Reset();
    fit->ClearAllErrorBands();

    ExpandHisto( fits[i], fit, 1 );
    //for( int j = 0; j< fits[i]->GetNbinsX(); j++ )
    //{
    //  cout<<fits[i]->GetBinContent(j+1);
    //  cout<<", "<<fit->GetBinContent(1, j+1)<<endl;
    //}
    h[cat]->Multiply( h[cat],fit );
    delete fit;
  }

  //h[kQELikeNot_SingleNeutralPion]->Multiply(h[kQELikeNot_SingleNeutralPion], fits[2]);
  //h[kQELikeNot_SingleChargedPion]->Multiply(h[kQELikeNot_SingleChargedPion], fits[2]);
 

  h[kMC]->Reset();
  //h[kQELike]->Reset();
  h[kQELikeNot]->Reset();


  //h[kQELike_QE]->Reset();
  //h[kQELike_QE]->Add(h[kQELike_QE_OTH]);
  //h[kQELike_QE]->Add(h[kQELike_QE_H]);
 
  //vector<int> qelike({ kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH });
  //for( auto cat: qelike ) h[kQELike]->Add( h[cat] );
  
  vector<int> qelikenot({ kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } );
  for( auto cat: qelikenot ) h[kQELikeNot]->Add( h[cat] );

  h[kMC]->Add( h[kQELike] );
  h[kMC]->Add( h[kQELikeNot] );
}


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

  for( auto cat : histosUsedAll ) h[cat]->Divide(h[cat], hmc);
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


void Draw2DDistro( TFile *f_orig );

void draw1DDistrosQuadrant( TFile *f_orig, TFile* f_signal_weighted, string histname = "h_pn_q2qe", string savename="plot-1d-pn", string xtitle="p_{n} (GeV)" );
void draw1DDistro( vector<MnvH2D**> h, vector<MnvH1D*> &fits, CCQENuPlotUtils* putils, string name, double pot_data, double pot_mc, string axis = "x", string xtitle="x", string ytitle="evt rate", bool doRatio=false, bool doFitting=false, bool doBckSub = false, bool subqelike=false, double min=-1., double max=-1.);

void plot1DInPad( MnvH1D** h, bool ratio, bool subbkg,bool scale,bool drawLegend, string xtitle, string ytitle, double xmin, double xmax, double ymin, double ymax, bool inPad=true )
{
  //Set style
  vector<TH1D*> histos(nHistos,NULL);
  for( auto id : histosUsedAll ) 
  {
    cout<<"ID: "<<names[id]<<endl;
    histos[id] =(TH1D*) h[id]->Clone( Form("%s_clone", h[id]->GetName() ) );

    if(scale) histos[id]->Scale(1,"width");
    histos[id]->SetLineColor(CCQENuPlotUtils::colors[id] );
    histos[id]->SetLineStyle( 0 );
    histos[id]->SetLineWidth( 0 );
    histos[id]->GetXaxis()->SetRangeUser(xmin,xmax);
    histos[id]->GetYaxis()->SetRangeUser(ymin,ymax);
    histos[id]->GetXaxis()->SetTitle( xtitle.c_str() );
    histos[id]->GetYaxis()->SetTitle( ytitle.c_str() );
    //histos[id]->SetFillColorAlpha( CCQENuPlotUtils::colors[id], 1 );
  }
  cout<<"Set Styles"<<endl;
  histos[kData]->SetLineWidth( 1 );
  histos[kData]->SetLineStyle( 1 );

  cout<<"histos[kQELike_QE_PN]"<<histos[kQELike_QE_PN]<<endl;
  for( int i = 0; i< histosQELike.size(); i++ )
  {
    int qelike = histosQELike[i];
    int qelike_pn = histosQELike_PN[i];
    int qelike_notpn = histosQELike_NotPN[i];
    histos[qelike_pn]->SetFillColor( CCQENuPlotUtils::colors[ qelike ] );
    histos[qelike_notpn]->SetFillColorAlpha( CCQENuPlotUtils::colors[ qelike ], 0.3 );

    histos[qelike_pn]->SetLineColor( CCQENuPlotUtils::colors[ qelike ] );
    histos[qelike_notpn]->SetLineColor( CCQENuPlotUtils::colors[ qelike ] );


    histos[qelike_pn]->SetFillStyle( 1001 );
    histos[qelike_notpn]->SetFillStyle( 1001 );


  }



  cout<<"Set PN Styles"<<endl;
  histos[kMC]->SetFillStyle(0);

  for( auto id : QELikeNotHistos )
  {
    //if( id != QELikeNotHistos[0] )
    //{
      histos[id]->SetLineStyle( 2 );
      histos[id]->SetFillColor( CCQENuPlotUtils::colors[id] );
      histos[id]->SetFillColorAlpha( CCQENuPlotUtils::colors[id], 0.3 );
      histos[id]->SetFillStyle( 3001 );
    //}
  }

  if(subbkg)
  {
    histos[kData]->Add( histos[kQELikeNot], -1 );
    histos[kMC]->Add( histos[kQELikeNot], -1 );

    histos[kQELikeNot]->Reset();
    histos[kQELikeNot_SingleNeutralPion]->Reset();
    histos[kQELikeNot_SingleChargedPion]->Reset();
    histos[kQELikeNot_MultiPion]->Reset();
    histos[kQELikeNot_NoPions]->Reset();

  }

  if( ratio )
  {
    MnvH1D* mc = (MnvH1D*) histos[kMC]->Clone("denom");
    for( auto id: histosUsedAll ) 
    {
      histos[id]->Divide( histos[id], mc );
      histos[id]->GetYaxis()->SetRangeUser(0,2);
    }
  }

  //draw
  gPad->cd();
  histos[kData]->Draw();
  THStack *hs = new THStack("hs","");
  if(!subbkg)
  {
    hs->Add( histos[kQELikeNot_NoPions] );
    hs->Add( histos[kQELikeNot_MultiPion] );
    hs->Add( histos[kQELikeNot_SingleChargedPion] );
    hs->Add( histos[kQELikeNot_SingleNeutralPion] );
  }
  hs->Add(histos[kQELike_OTH_NotPN] );
  hs->Add(histos[kQELike_DIS_NotPN] );
  hs->Add(histos[kQELike_RES_NotPN] );
  hs->Add(histos[kQELike_2p2h_NotPN] );
  hs->Add(histos[kQELike_QE_NotPN] );
  hs->Add(histos[kQELike_OTH_PN] );
  hs->Add(histos[kQELike_DIS_PN] );
  hs->Add(histos[kQELike_RES_PN] );
  hs->Add(histos[kQELike_2p2h_PN] );
  hs->Add(histos[kQELike_QE_PN] );

  hs->Draw("histsame");
  histos[kMC]->Draw("histsame");
  histos[kData]->Draw("same");

  histos[kQELike_NoProton_NotPN]->SetLineColor(kBlack);
  histos[kQELike_NoNeutron_NotPN]->SetLineColor(kRed);
  histos[kQELike_NoNucleon_NotPN]->SetLineColor(kBlue);

  //histos[kQELike_NoProton_NotPN]->Draw("histsame");
  //histos[kQELike_NoNeutron_NotPN]->Draw("histlsame");
  //histos[kQELike_NoNucleon_NotPN]->Draw("same");

  //TLatex* latex = new TLatex();
  //latex->DrawLatexNDC(0,0.8, label.c_str() );

  if(drawLegend)
  {
    TLegend* leg = NULL;
    if(inPad)
    {
      if(!ratio) leg = new TLegend(0.7,0.5,.99,.99);
      else 
      {
        leg = new TLegend(0.3,0.8,.99,.99);
        leg->SetNColumns(3);
        leg->SetFillStyle(0);
      }
    }
    else
    {
      if(!ratio) leg = new TLegend(0.7,0.6,.85,.9);
      else 
      {
        leg = new TLegend(0.15,0.65,.85,.9);
        leg->SetNColumns(3);
        leg->SetFillStyle(0);
      }
    }
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry( histos[kData], "Data", "p");
    leg->AddEntry( histos[kMC], "MC","l");
    leg->AddEntry( histos[kQELike_QE_PN], "QE","lf");
    leg->AddEntry( histos[kQELike_2p2h_PN], "2p2h","lf");
    leg->AddEntry( histos[kQELike_RES_PN], "RES","lf");
    leg->AddEntry( histos[kQELike_DIS_PN], "DIS","lf");
    leg->AddEntry( histos[kQELike_OTH_PN], "OTH","lf");
    leg->AddEntry( histos[kQELike_QE_NotPN], "Not PN && QE","lf");
    leg->AddEntry( histos[kQELike_2p2h_NotPN], "2p2h","lf");
    leg->AddEntry( histos[kQELike_RES_NotPN], "RES","lf");
    leg->AddEntry( histos[kQELike_DIS_NotPN], "DIS","lf");
    leg->AddEntry( histos[kQELike_OTH_NotPN], "OTH","lf");
    if(!subbkg)
    {
      //leg->AddEntry( histos[kQELikeNot], "QELikeNot","lf");
      leg->AddEntry( histos[kQELikeNot_SingleNeutralPion], "1#pi^{0}","lf");
      leg->AddEntry( histos[kQELikeNot_SingleChargedPion], "1#pi^{+/-}","lf");
      leg->AddEntry( histos[kQELikeNot_MultiPion], "N#pi","lf");
      leg->AddEntry( histos[kQELikeNot_NoPions], "0#pi","lf");
    }
    leg->Draw();

  }


}


void draw1DDistro2( MnvH2D** h, vector<MnvH1D*> &fits, CCQENuPlotUtils* putils, string name, double pot_data, double pot_mc, string axis, string xtitle, string ytitle, int combineBin, bool doRatio, bool doFitting, bool doBckSub, bool subqelike, double min_x, double max_x, string tag)
{
  bool area_norm = false;
  bool includeData = true;

  MnvH2D* hists[nHistos];
  for(unsigned int i = 0; i<nHistos; i++) hists[i] = new MnvH2D(*( h[i]->Clone(Form("h2d_tmp_%s_%s", h[i]->GetName(), names[i].c_str() ))));

  if(doFitting) DoFit( hists, fits );
  if(doBckSub ) SubtractBackground( hists, subqelike );

  MnvH1D* hists1d[nHistos];

  for(unsigned int i = 0; i<nHistos; i++ ) 
  {
    int minbin = 1;
    if(axis == "x") hists1d[i] = new MnvH1D(*(hists[i]->ProjectionX(Form("h1d_tmp_%s", hists[i]->GetName()),1,hists[i]->GetNbinsY())) );
    if(axis == "y") hists1d[i] = new MnvH1D(*(hists[i]->ProjectionY(Form("h1d_tmp_%s", hists[i]->GetName()),1,hists[i]->GetNbinsX())) );
    hists1d[i]->GetXaxis()->SetTitle( xtitle.c_str() );
    hists1d[i]->GetYaxis()->SetTitle( ytitle.c_str() );
    if(combineBin > 1) hists1d[i]=(MnvH1D*) hists1d[i]->Rebin(combineBin, hists1d[i]->GetName() );
  }
  //MnvH1D** hists1d = &hists1d_cp[0];
  if(doRatio ) RatioToMC( hists1d );

  TCanvas* c = new TCanvas("cNeutA","Neutron Angulars"); 
  MnvPlotter *plotter = new MnvPlotter;



  double mcscale = -1.; string xaxislabel = xtitle; string yaxislabel = ""; double min_y = -1.; double max_y = -1.;
  bool normalizeHisto = (doRatio)? false: true;
  bool scaleHisto =false;
  if (doRatio)
  {
    min_y = 0;
    max_y = 2;
  }
  else
  {
    //double value = hists1d[kData]->GetBinNormalizedCopy().GetMaximum();
    double value = hists1d[kData]->GetMaximum();
    max_y = value*1.4;
    min_y = 0;
  }

  string plotType = "QELike_split_PN";
  if( doBckSub )
  {
    plotType = "QELikeOnly";
    if (subqelike) plotType = "HOnly";
  }
  //putils->drawStacked(hists1d, "QELike_split_PN", area_norm, pot_data, pot_mc, includeData,mcscale,xaxislabel,yaxislabel,min_x, max_x, min_y ,max_y,normalizeHisto );

  bool drawLegend = true;
  plot1DInPad( hists1d, doRatio, doBckSub, scaleHisto, drawLegend, xaxislabel, yaxislabel, min_x, max_x, min_y, max_y ,false);


  if( min_x != -1 ) cout<<min_x<<"------------"<<max_x<<endl;
  TH1D* herr = (TH1D*) hists1d[kMC]->Clone("herr");
  if(!doRatio && scaleHisto) herr->Scale(herr->GetBinWidth(1),scaleOpt.c_str() );
  herr->SetFillColorAlpha(kBlue, 0.3);
  herr->Draw("e2same");
  if( tag!="" )
  {
    TLatex latex;
    latex.SetTextColor(kBlue);
    latex.SetTextFont(82);
    latex.SetTextSize(0.1);
    latex.DrawLatexNDC( 0.2,0.85, tag.c_str() );
  }


  string fname=Form("plots/%s-ratio_%d-bcksub_%d-subqelike_%d-fit_%d.pdf",name.c_str(), doRatio, doBckSub,subqelike, doFitting );
  //plotter->AddPlotLabel(Form("# Hydrogen: %.2f", hists1d[kQELike_QE_H]->Integral()), 0.3,0.8,0.033,kBlue,32);
  c->Print( fname.c_str() );

  return;
}


void draw1DDistro3( MnvH1D** h, vector<MnvH1D*> &fits, CCQENuPlotUtils* putils, string name, double pot_data, double pot_mc, string axis, string xtitle, string ytitle, bool doRatio, bool doFitting, bool doBckSub, bool subqelike, double min_x, double max_x)
{
  bool area_norm = false;
  bool includeData = true;

  MnvH1D* hists1d[nHistos];
  for(unsigned int i = 0; i<nHistos; i++) 
  {
    hists1d[i] = new MnvH1D(*( h[i]->Clone(Form("h1d_tmp_%s_%s", h[i]->GetName(), names[i].c_str() ))));
    hists1d[i]->GetXaxis()->SetTitle( xtitle.c_str() );
    hists1d[i]->GetYaxis()->SetTitle( ytitle.c_str() );
  }

  //if(doFitting) DoFit( hists1d, fits );
  if(doBckSub ) SubtractBackground( hists1d, subqelike );


  //MnvH1D** hists1d = &hists1d_cp[0];
  if(doRatio ) RatioToMC( hists1d );

  TCanvas* c = new TCanvas("cNeutA","Neutron Angulars"); 
  MnvPlotter *plotter = new MnvPlotter;



  double mcscale = -1.; string xaxislabel = xtitle; string yaxislabel = ""; double min_y = -1.; double max_y = -1.;
  bool normalizeHisto = (doRatio)? false: true;
  bool scaleHisto =false;
  if (doRatio)
  {
    min_y = 0;
    max_y = 2;
  }
  else
  {
    //double value = hists1d[kData]->GetBinNormalizedCopy().GetMaximum();
    double value = hists1d[kData]->GetMaximum();
    max_y = value*1.4;
    min_y = 0;
  }

  string plotType = "QELike_split_PN";
  if( doBckSub )
  {
    plotType = "QELikeOnly";
    if (subqelike) plotType = "HOnly";
  }
  //putils->drawStacked(hists1d, "QELike_split_PN", area_norm, pot_data, pot_mc, includeData,mcscale,xaxislabel,yaxislabel,min_x, max_x, min_y ,max_y,normalizeHisto );

  bool drawLegend = true;
  plot1DInPad( hists1d, doRatio, doBckSub, scaleHisto, drawLegend, xaxislabel, yaxislabel, min_x, max_x, min_y, max_y ,false);


  if( min_x != -1 ) cout<<min_x<<"------------"<<max_x<<endl;
  TH1D* herr = (TH1D*) hists1d[kMC]->Clone("herr");
  if(!doRatio) herr->Scale(herr->GetBinWidth(1),scaleOpt.c_str() );
  herr->SetFillColorAlpha(kBlue, 0.3);
  herr->Draw("e2same");


  string fname=Form("plots/%s-ratio_%d-bcksub_%d-subqelike_%d-fit_%d.pdf",name.c_str(), doRatio, doBckSub,subqelike, doFitting );
  //plotter->AddPlotLabel(Form("# Hydrogen: %.2f", hists1d[kQELike_QE_H]->Integral()), 0.3,0.8,0.033,kBlue,32);
  c->Print( fname.c_str() );

  return;
}


void Draw1DHistograms(  TFile *f_orig, TFile* f_signal_weighted, string yvar="ptmu" )
{
  cout<<"Draw1DHistograms"<<endl;
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

  cout<<"==================== Book Histograms ======="<<endl;
  MnvH2D *h_dphitt_yvar[nHistos];
  MnvH2D *h_dphitt_left_yvar[nHistos];
  MnvH2D *h_dphitt_right_yvar[nHistos];
  MnvH2D *h_dphitt_sign_yvar[nHistos];
  MnvH2D *h_pnAngle_yvar[nHistos];
  MnvH2D *h_pnEndDist_yvar[nHistos];
  MnvH2D *h_blobDist_yvar[nHistos];
  MnvH1D *h_vtx_energy[nHistos];
      
//h_dphitt_sign_ptmu_left_mc
  utils->bookHistos( f1, h_dphitt_yvar,     Form( "h_dphitt_%s", yvar.c_str() ) );
  utils->bookHistos( f1, h_dphitt_left_yvar,     Form( "h_dphitt_sign_%s_left", yvar.c_str() ) );
  utils->bookHistos( f1, h_dphitt_right_yvar,     Form( "h_dphitt_sign_%s_right", yvar.c_str() ) );
  utils->bookHistos( f1, h_dphitt_sign_yvar,     Form( "h_dphitt_sign_%s", yvar.c_str() ) );
  utils->bookHistos( f1, h_pnAngle_yvar,     Form( "h_pnAngle_%s", yvar.c_str() ) );
  utils->bookHistos( f1, h_pnEndDist_yvar,     Form( "h_pnEndDist_%s", yvar.c_str() ) );
  utils->bookHistos( f1, h_blobDist_yvar,     Form( "h_blobDist_%s", yvar.c_str() ) );
  utils->bookHistos( f1, h_vtx_energy,     Form( "h_vtx_energy" ) );

  putils->scaleMCHistos(  h_dphitt_yvar,    pot_norm );
  putils->scaleMCHistos(  h_pnAngle_yvar,    pot_norm );
  string axis = "x";
  cout<<"Make Plots"<<endl;
  cout<<"================= Draw Histos ==============="<<endl;

  vector<MnvH1D*> fits = GetFitsHisto(f2);
  //vector<MnvH2D*> fits = GetFitsHisto(f2, "dpt");
  //bool doRatio, bool doFitting, bool doBckSub
  int rebin=4;
  draw1DDistro2( h_dphitt_yvar, fits, putils, "plot-1d-dphitt", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Event Rate",rebin, 0,0,0);
  draw1DDistro2( h_dphitt_yvar, fits, putils, "plot-1d-dphitt", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Event Rate",rebin, 0,1,0);
  draw1DDistro2( h_dphitt_yvar, fits, putils, "plot-1d-dphitt", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Ratio To MC",rebin, 1,0,0);
  draw1DDistro2( h_dphitt_yvar, fits, putils, "plot-1d-dphitt", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Ratio To MC",rebin, 1,1,0);
  draw1DDistro2( h_dphitt_yvar, fits, putils, "plot-1d-dphitt", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Event Rate",rebin, 0,1,1);
  draw1DDistro2( h_dphitt_yvar, fits, putils, "plot-1d-dphitt", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Ratio To MC",rebin, 1,1,1);

  draw1DDistro2( h_dphitt_left_yvar, fits, putils, "plot-1d-dphitt-left", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Event Rate",rebin, 0,0,0,0,-180,180, "-");
  draw1DDistro2( h_dphitt_left_yvar, fits, putils, "plot-1d-dphitt-left", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Event Rate",rebin, 0,1,0,0,-180,180, "-");
  draw1DDistro2( h_dphitt_left_yvar, fits, putils, "plot-1d-dphitt-left", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Ratio To MC",rebin, 1,0,0,0,-180,180, "-");
  draw1DDistro2( h_dphitt_left_yvar, fits, putils, "plot-1d-dphitt-left", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Ratio To MC",rebin, 1,1,0,0,-180,180, "-");
  draw1DDistro2( h_dphitt_left_yvar, fits, putils, "plot-1d-dphitt-left", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Event Rate", rebin,0,1,1,0,-180,180, "-");
  draw1DDistro2( h_dphitt_left_yvar, fits, putils, "plot-1d-dphitt-left", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Ratio To MC",rebin, 1,1,1,0,-180,180, "-");

  draw1DDistro2( h_dphitt_right_yvar, fits, putils, "plot-1d-dphitt-right", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Event Rate",rebin, 0,0,0,0,-180,180, "+");
  draw1DDistro2( h_dphitt_right_yvar, fits, putils, "plot-1d-dphitt-right", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Event Rate",rebin, 0,1,0,0,-180,180, "+");
  draw1DDistro2( h_dphitt_right_yvar, fits, putils, "plot-1d-dphitt-right", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Ratio To MC",rebin, 1,0,0,0,-180,180, "+");
  draw1DDistro2( h_dphitt_right_yvar, fits, putils, "plot-1d-dphitt-right", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Ratio To MC",rebin, 1,1,0,0,-180,180, "+");
  draw1DDistro2( h_dphitt_right_yvar, fits, putils, "plot-1d-dphitt-right", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Event Rate", rebin,0,1,1,0,-180,180, "+");
  draw1DDistro2( h_dphitt_right_yvar, fits, putils, "plot-1d-dphitt-right", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Ratio To MC",rebin, 1,1,1,0,-180,180, "+");

  draw1DDistro2( h_dphitt_sign_yvar, fits, putils, "plot-1d-dphitt-sign", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Event Rate",rebin, 0,0,0,0,-90,90);
  draw1DDistro2( h_dphitt_sign_yvar, fits, putils, "plot-1d-dphitt-sign", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Event Rate",rebin, 0,1,0,0,-90,90);
  draw1DDistro2( h_dphitt_sign_yvar, fits, putils, "plot-1d-dphitt-sign", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Ratio To MC",rebin, 1,0,0,0,-90,90);
  draw1DDistro2( h_dphitt_sign_yvar, fits, putils, "plot-1d-dphitt-sign", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Ratio To MC",rebin, 1,1,0,0,-90,90);
  draw1DDistro2( h_dphitt_sign_yvar, fits, putils, "plot-1d-dphitt-sign", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Event Rate", rebin,0,1,1,0,-90,90);
  draw1DDistro2( h_dphitt_sign_yvar, fits, putils, "plot-1d-dphitt-sign", pot_data, pot_mc, "x", "#delta #phi_{TT} (GeV/c)", "Ratio To MC",rebin, 1,1,1,0,-90,90);


  rebin=1;

  draw1DDistro2( h_pnAngle_yvar, fits, putils, "plot-1d-pnAngle", pot_data, pot_mc, "x", "#theta_{p,n} (GeV/c)", "Event Rate",rebin, 0,0,0);
  draw1DDistro2( h_pnAngle_yvar, fits, putils, "plot-1d-pnAngle", pot_data, pot_mc, "x", "#theta_{p,n} (GeV/c)", "Event Rate",rebin, 0,1,0);
  draw1DDistro2( h_pnAngle_yvar, fits, putils, "plot-1d-pnAngle", pot_data, pot_mc, "x", "#theta_{p,n} (GeV/c)", "Ratio To MC",rebin, 1,0,0);
  draw1DDistro2( h_pnAngle_yvar, fits, putils, "plot-1d-pnAngle", pot_data, pot_mc, "x", "#theta_{p,n} (GeV/c)", "Ratio To MC",rebin, 1,1,0);
  draw1DDistro2( h_pnAngle_yvar, fits, putils, "plot-1d-pnAngle", pot_data, pot_mc, "x", "#theta_{p,n} (GeV/c)", "Event Rate", rebin,0,1,1);
  draw1DDistro2( h_pnAngle_yvar, fits, putils, "plot-1d-pnAngle", pot_data, pot_mc, "x", "#theta_{p,n} (GeV/c)", "Ratio To MC",rebin, 1,1,1);

  draw1DDistro2( h_blobDist_yvar, fits, putils, "plot-1d-ptmu", pot_data, pot_mc, "y", "Muon p_{T} (GeV/c)", "Event Rate",rebin, 0,0,0);
  draw1DDistro2( h_blobDist_yvar, fits, putils, "plot-1d-ptmu", pot_data, pot_mc, "y", "Muon p_{T} (GeV/c)", "Event Rate",rebin, 0,1,0);
  draw1DDistro2( h_blobDist_yvar, fits, putils, "plot-1d-ptmu", pot_data, pot_mc, "y", "Muon p_{T} (GeV/c)", "Ratio To MC",rebin, 1,0,0);
  draw1DDistro2( h_blobDist_yvar, fits, putils, "plot-1d-ptmu", pot_data, pot_mc, "y", "Muon p_{T} (GeV/c)", "Ratio To MC",rebin, 1,1,0);
  draw1DDistro2( h_blobDist_yvar, fits, putils, "plot-1d-ptmu", pot_data, pot_mc, "y", "Muon p_{T} (GeV/c)", "Event Rate", rebin,0,1,1);
  draw1DDistro2( h_blobDist_yvar, fits, putils, "plot-1d-ptmu", pot_data, pot_mc, "y", "Muon p_{T} (GeV/c)", "Ratio To MC",rebin, 1,1,1);


  draw1DDistro2( h_pnEndDist_yvar, fits, putils, "plot-1d-pnEndDist", pot_data, pot_mc, "x", "Dist EndPts (mm)", "Event Rate", rebin,  0,0,0);
  draw1DDistro2( h_pnEndDist_yvar, fits, putils, "plot-1d-pnEndDist", pot_data, pot_mc, "x", "Dist EndPts (mm)", "Event Rate", rebin,  0,1,0);
  draw1DDistro2( h_pnEndDist_yvar, fits, putils, "plot-1d-pnEndDist", pot_data, pot_mc, "x", "Dist EndPts (mm)", "Ratio To MC",rebin,  1,0,0);
  draw1DDistro2( h_pnEndDist_yvar, fits, putils, "plot-1d-pnEndDist", pot_data, pot_mc, "x", "Dist EndPts (mm)", "Ratio To MC",rebin,  1,1,0);
  draw1DDistro2( h_pnEndDist_yvar, fits, putils, "plot-1d-pnEndDist", pot_data, pot_mc, "x", "Dist EndPts (mm)", "Event Rate", rebin,  0,1,1);
  draw1DDistro2( h_pnEndDist_yvar, fits, putils, "plot-1d-pnEndDist", pot_data, pot_mc, "x", "Dist EndPts (mm)", "Ratio To MC",rebin,  1,1,1);

  draw1DDistro2( h_blobDist_yvar, fits, putils, "plot-1d-blobDist", pot_data, pot_mc, "x", "Blob Vtx Dist (mm)", "Event Rate", rebin,  0,0,0);
  draw1DDistro2( h_blobDist_yvar, fits, putils, "plot-1d-blobDist", pot_data, pot_mc, "x", "Blob Vtx Dist (mm)", "Event Rate", rebin,  0,1,0);
  draw1DDistro2( h_blobDist_yvar, fits, putils, "plot-1d-blobDist", pot_data, pot_mc, "x", "Blob Vtx Dist (mm)", "Ratio To MC",rebin,  1,0,0);
  draw1DDistro2( h_blobDist_yvar, fits, putils, "plot-1d-blobDist", pot_data, pot_mc, "x", "Blob Vtx Dist (mm)", "Ratio To MC",rebin,  1,1,0);
  draw1DDistro2( h_blobDist_yvar, fits, putils, "plot-1d-blobDist", pot_data, pot_mc, "x", "Blob Vtx Dist (mm)", "Event Rate", rebin,  0,1,1);
  draw1DDistro2( h_blobDist_yvar, fits, putils, "plot-1d-blobDist", pot_data, pot_mc, "x", "Blob Vtx Dist (mm)", "Ratio To MC",rebin,  1,1,1);






  draw1DDistro3( h_vtx_energy, fits, putils, "plot-1d-vtxEnergy", pot_data, pot_mc, "x", "E_{vtx} (MeV/c)", "Event Rate", 0,0,0);
  draw1DDistro3( h_vtx_energy, fits, putils, "plot-1d-vtxEnergy", pot_data, pot_mc, "x", "E_{vtx} (MeV/c)", "Ratio To MC", 1,0,0);



}


void draw1DDistrosQuadrant(  TFile *f_orig, TFile* f_signal_weighted, string histname, string savename, string xtitle )
{
  cout<<"draw1DDistrosQuadrant"<<endl;
  cout<<"histo name: "<<histname<<endl;
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
  if(!f1) cout<<"f1 is not read"<<endl;

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
  MnvH2D *h_ll[nHistos], *h_lr[nHistos], *h_rl[nHistos], *h_rr[nHistos];
    

  utils->bookHistos( f1, h_ll, Form( "%s_ll", histname.c_str()) );
  utils->bookHistos( f1, h_lr, Form( "%s_lr", histname.c_str()) );
  utils->bookHistos( f1, h_rr, Form( "%s_rr", histname.c_str()) );
  utils->bookHistos( f1, h_rl, Form( "%s_rl", histname.c_str()) );
  cout<<"Booked quadrant"<<endl;
  cout<<h_ll[kMC]<<endl;
  cout<<h_lr[kMC]<<endl;
  cout<<h_rl[kMC]<<endl;
  cout<<h_rr[kMC]<<endl;

  vector<MnvH2D**> histos({ h_ll, h_lr, h_rl, h_rr } );
  for( auto h : histos )
  {
    putils->scaleMCHistos( h,    pot_norm );
  }

  string axis = "x";
  cout<<"Make Plots"<<endl;
  cout<<"================= Draw Histos ==============="<<endl;

  vector<MnvH1D*> fits = GetFitsHisto(f2);
  cout<<"================= Print Fits ================"<<endl;
  for(unsigned int i = 0; i< fits.size(); i++ )cout<<"Fit histo address: "<<fits[i]<<endl;

  draw1DDistro( histos, fits, putils, savename,  pot_data, pot_mc, "x", xtitle, "Event Rate ( #times 10^{3} )", 0,0,0);
  draw1DDistro( histos, fits, putils, savename,  pot_data, pot_mc, "x", xtitle, "Event Rate ( #times 10^{3} )", 0,1,0);

  draw1DDistro( histos, fits, putils, savename,  pot_data, pot_mc, "x", xtitle, "Ratio To MC", 1,0,0);
  draw1DDistro( histos, fits, putils, savename,  pot_data, pot_mc, "x", xtitle, "Ratio To MC", 1,1,0);
}

void draw1DDistro( vector<MnvH2D**> h, vector<MnvH1D*> &fits, CCQENuPlotUtils* putils, string name, double pot_data, double pot_mc, string axis, string xtitle, string ytitle, bool doRatio, bool doFitting, bool doBckSub, bool subqelike, double min_x, double max_x)
{
  cout<<"------------drawing: "<<name<<endl;
  bool area_norm = false;
  bool includeData = true;
  int nRegions = 4;

  vector< vector<MnvH2D*>> h2ds(nRegions, vector<MnvH2D*>(nHistos) );
  cout<<"transfer histos to new container"<<endl;
  for( int i = 0; i < nRegions; i++ )
  {
    for( unsigned int j = 0 ; j< nHistos; j++ )
    {
      h2ds[i][j] = new MnvH2D(*( h[i][j]->Clone(Form("h2d_tmp_%d_%s_%s", i, h[i][j]->GetName(), names[j].c_str() ))));
      h2ds[i][j]->Scale(0.001);
    }
    cout<<"N Data for 2D: "<<i<<" : "<<h2ds[i][kData]->Integral()<<endl;
  }
  for( int i = 0; i<nRegions; i++ )
  {
    if(doFitting) DoFit( &h2ds[i][0],fits );
    if(doBckSub ) SubtractBackground( &h2ds[i][0], subqelike );
  }

  cout<<"Define 1D histo and containers"<<endl;
  vector<vector<MnvH1D*>> h1ds(nRegions, vector<MnvH1D*>(nHistos) );


  for( int j = 0; j< nRegions; j++ )
  {
    cout<<"Projecting "<<j<<endl;
    for(unsigned int i = 0; i<nHistos; i++ ) 
    {
      //cout<<"do "<<j<<", "<<i<<" ... ";
      if(axis == "x") h1ds[j][i] = new MnvH1D(*(h2ds[j][i]->ProjectionX(Form("h1d_tmp_%d_%d",j,i),1,h2ds[j][i]->GetNbinsY())) );
      if(axis == "y") h1ds[j][i] = new MnvH1D(*(h2ds[j][i]->ProjectionY(Form("h1d_tmp_%d_%d",j,i),1,h2ds[j][i]->GetNbinsX())) );
      //cout<<"done"<<endl;
    }
    cout<<"N Data for "<<j<<" : "<<h1ds[j][kData]->Integral()<<endl;
  }
  //MnvH1D** hists1d = &hists1d_cp[0];
  //if(doRatio ) 
  //{
  //  for( auto h : h1ds ) RatioToMC( &h[0] );
  //}


  MnvPlotter *plotter = new MnvPlotter;
  TCanvas* c = new TCanvas("cNeutA","Neutron Angulars"); 
  c->Divide(2,2, 0, 0);



  double mcscale = -1.; string xaxislabel = ""; string yaxislabel = ""; double min_y = -1.; double max_y = -1.;
  bool normalizeHisto = (doRatio)? false: true;
  bool scaleHisto = true;
  if (doRatio)
  {
    min_y = 0;
    max_y = 2;
  }
  else
  {
    vector<TH1D> HISTS;
    HISTS.push_back(h1ds[0][kData]->GetCVHistoWithStatError() );
    HISTS.push_back(h1ds[1][kData]->GetCVHistoWithStatError() );
    HISTS.push_back(h1ds[2][kData]->GetCVHistoWithStatError() );
    HISTS.push_back(h1ds[3][kData]->GetCVHistoWithStatError() );
    if(scaleHisto) 
    {
      HISTS[0].Scale(1,"width");
      HISTS[1].Scale(1,"width");
      HISTS[2].Scale(1,"width");
      HISTS[3].Scale(1,"width");
    }
    
    double value = -999;
    for( int i = 0; i<4; i++ ) 
    {
      double tmp = HISTS[i].GetMaximum();
      if( tmp> value ) value = tmp;
    }
    max_y = value*1.4;
    cout<<"max_y = "<<max_y<<endl;
    min_y = 0;
  }

  string plotType = "QELike_split_PN";
  if( doBckSub )
  {
    plotType = "QELikeOnly";
    if (subqelike) plotType = "HOnly";
  }

  string label = "--";
  for( int i = 0; i< 4; i++ )
  {
    c->cd(i+1);
    gPad->SetTicky(2);

    if(i==1) label = "-+";
    if(i==2) label = "+-";
    if(i==3) label = "++";
    bool drawLegend = (i==3);
    //putils->drawStacked(&(h1ds[i][0]), "QELike_split_PionInFS", area_norm, pot_data, pot_mc, includeData,mcscale,xaxislabel,yaxislabel,min_x, max_x, min_y ,max_y,normalizeHisto );
    plot1DInPad( &(h1ds[i][0]), doRatio, false, scaleHisto, drawLegend, xaxislabel, yaxislabel, min_x, max_x, min_y, max_y );
    //if( min_x != -1 ) cout<<min_x<<"------------"<<max_x<<endl;
    //TH1D* herr = (TH1D*) h1ds[i][kMC]->Clone("herr");
    //if(!doRatio) herr->Scale(herr->GetBinWidth(1),"width");
    //herr->SetFillColorAlpha(kBlue, 0.3);
    //herr->Draw("e2same");
  }

  c->cd();
  TLatex* latex = new TLatex();
  latex->SetTextFont(62);
  latex->SetTextSize(0.04);
  latex->SetTextAlign(22);
  latex->DrawLatexNDC(0.5,0.04, xtitle.c_str() );
  latex->SetTextAngle(90);
  latex->DrawLatexNDC(0.02,0.5, ytitle.c_str() );

  latex->SetTextAngle(0);
  latex->SetTextFont(82);
  latex->SetTextSize(0.04);
  latex->SetTextColor(kBlue);
  latex->DrawLatexNDC(0.1,0.93,"--");
  latex->DrawLatexNDC(0.54,0.93,"-+");
  latex->DrawLatexNDC(0.1,0.49,"+-");
  latex->DrawLatexNDC(0.54,0.49,"++");

  string fname=Form("plots/%s-ratio_%d-bcksub_%d-subqelike_%d-fit_%d.pdf",name.c_str(), doRatio, doBckSub,subqelike, doFitting );
  //plotter->AddPlotLabel( xtitle.c_str(), 0.4,0.
  //plotter->AddPlotLabel(Form("# Hydrogen: %.2f", hists1d[kQELike_QE_H]->Integral()), 0.3,0.8,0.033,kBlue,32);
  c->Print( fname.c_str() );

  return;
}


void DrawCuts(TCanvas *c1, int i, string axis, double x0, double y0, double x1, double y1, int set=1)
{
  // 0:pn angle, 1: pnEndDist, 2: blobvtxDist
  vector<double> values({30,1500,400});
  vector<double> values2({30,1500,1800});
  if(set==2) values = values2;
  double val = values[i];
  if( axis == "x" )
  {
    x0 = val;
    x1 = val;
  }
  else
  {
    y0 = val;
    y1 = val;
  }
  cout<<"==============="<<endl;
  //cout<<X0<<", "<<X1<<", "<<Y0<<", "<<Y1<<endl;
  cout<<x0<<", "<<x1<<", "<<y0<<", "<<y1<<endl;

  TLine* line = new TLine(x0,y0,x1,y1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->SetLineColor( kRed );
  c1->cd();
  line->Draw("same");
}

void Draw2DDistro( TFile *f_orig )
{
  cout<<"Draw2DDistro"<<endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();

  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  //three files 2-track, 2-track bkg_fitted
  //CV
  CCQENuPlotUtils *putils = new CCQENuPlotUtils( fluxHistoExists );

  vector<int> types_to_draw( {kQELike_QE_PN, kQELike_2p2h_PN,kQELike_RES_PN, kQELike_DIS_PN,  kQELike_QE_NotPN, kQELike_2p2h_NotPN,kQELike_RES_NotPN} );
  vector<string> hnames({ "h_pnEndDist_pnAngle","h_blobDist_pnAngle", "h_blobDist_pnEndDist" });
  vector<string> xtitles({ "EndDist (mm)", "Blob Dist (mm)", "Blob Dist (mm)" });
  vector<string> ytitles({ "#theta_{p,n} (degree)", "#theta_{p,n} (degree)", "EndDist (mm)"});

  // 0:pn angle, 1: pnEndDist, 2: blobvtxDist
  vector<int> xcuts({1,2,2});
  vector<int> ycuts({0,0,1});

  vector<vector<MnvH2D*>> histos;
  for( int i = 0; i<hnames.size(); i++ )
  {
    string name = hnames[i];
    vector<MnvH2D*> hs(nHistos,NULL);
    putils->bookHistos(f_orig,  &hs[0], name.c_str() );
    for( auto h : hs )
    {
      h->GetXaxis()->SetTitle(  xtitles[i].c_str() );
      h->GetYaxis()->SetTitle(  ytitles[i].c_str() );
    }
    //actually, I can just draw them here..
    //histos.push_back( hs );
    
    int xcut = xcuts[i];
    int ycut = ycuts[i];
    
    TCanvas *c1 = new TCanvas("c1","c1");
    for( auto id : types_to_draw )
    {
      c1->cd();

      MnvH2D* h = (MnvH2D*) hs[id]->Clone("tmp");
      double x0 = h->GetXaxis()->GetXmin();
      double x1 = h->GetXaxis()->GetXmax();
      double y0 = h->GetYaxis()->GetXmin();
      double y1 = h->GetYaxis()->GetXmax();

      h->Draw("colz");

      gPad->Update();
      TPaletteAxis *palette = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
      palette->SetY1NDC(0.25);
      gPad->Update();
      
      for(int i = 1; i<=2; i++)
      {
        if(xcut>=0) DrawCuts(c1, xcut, "x", x0, y0, x1, y1, i);
        if(ycut>=0) DrawCuts(c1, ycut, "y", x0, y0, x1, y1, i);
      }



      c1->Print( Form( "plots/2D_EvtRate_%s_%s.pdf",name.c_str(), names[id].c_str() ) );
      h->Divide( h, hs[kQELike] );
      h->GetZaxis()->SetRangeUser(0,1);
      h->Draw("colz");

      for(int i = 1; i<=2; i++)
      {
        if(xcut>=0) DrawCuts(c1, xcut, "x", x0, y0, x1, y1, i);
        if(ycut>=0) DrawCuts(c1, ycut, "y", x0, y0, x1, y1, i);
      }



      c1->Print( Form( "plots/2D_Ratio_%s_%s.pdf",name.c_str(), names[id].c_str() ) );
    }

      c1->cd();

      MnvH2D* h = (MnvH2D*) hs[kQELike]->Clone("qelike_signal");
      h->Add( hs[kQELike_QE_NotPN], -1);
      h->Add( hs[kQELike_2p2h_NotPN], -1);
      h->Add( hs[kQELike_RES_NotPN], -1);
      h->Add( hs[kQELike_DIS_NotPN], -1);
      h->Add( hs[kQELike_OTH_NotPN], -1);
      double x0 = h->GetXaxis()->GetXmin();
      double x1 = h->GetXaxis()->GetXmax();
      double y0 = h->GetYaxis()->GetXmin();
      double y1 = h->GetYaxis()->GetXmax();

      h->Draw("colz");

      gPad->Update();
      TPaletteAxis *palette = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
      palette->SetY1NDC(0.25);
      gPad->Update();
      for(int i = 1; i<=2; i++)
      {
        if(xcut>=0) DrawCuts(c1, xcut, "x", x0, y0, x1, y1, i);
        if(ycut>=0) DrawCuts(c1, ycut, "y", x0, y0, x1, y1, i);
      }
     

      c1->Print( Form( "plots/2D_EvtRate_%s_%s.pdf",name.c_str(), "qelikeSignal" ) );
      h->Divide( h, hs[kQELike] );
      h->GetZaxis()->SetRangeUser(0,1);
      h->Draw("colz");

      for(int i = 1; i<=2; i++)
      {
        if(xcut>=0) DrawCuts(c1, xcut, "x", x0, y0, x1, y1, i);
        if(ycut>=0) DrawCuts(c1, ycut, "y", x0, y0, x1, y1, i);
      }


      c1->Print( Form( "plots/2D_Ratio_%s_%s.pdf",name.c_str(), "qelikeSignal" ) );
  }
}





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
  //c1->SetLogx();
  c1->SetGrid();
  c1->Print( savename.c_str() );

}

vector<MnvH1D*> GetRatio( MnvH1D** h)
{
  vector<MnvH1D*> ret( nHistos);
  for( uint i = 0; i< nHistos; i++ )
  {
    ret[i] = (MnvH1D*) h[i]->Clone( Form( "h_tmp_%s", names[i].c_str() ) );
    ret[i]->Divide( ret[i], h[kMC] );
  }
  return ret;
}

void DrawFits( TFile* f, TFile *fweight )
{
  cout<<"Entering Draw Fits"<<endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  //gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  vector<MnvH1D*> fits = GetFitsHisto( fweight );
  cout<<"Fit Histo Name: "<<endl;
  for( UInt_t i = 0; i<fits.size(); i++ ) cout<<fits[i]->GetName()<<endl;

  TCanvas *c1 = new TCanvas("c1","c1");
  for( UInt_t i = 0; i<categories_to_fit.size(); i++ ) DrawMCWithErrorBand ( c1, fits[i], categories_to_fit_title_names[i], "p_{T,#mu}(GeV)", "weight", Form("plots/rw_%s.pdf", categories_to_fit_names[i].c_str()) );


  TVector2 *pot = (TVector2*)f->Get("pot");
  double pot_data = pot->X();
  double pot_mc = pot->Y();
  double pot_norm = pot_data/pot_mc;

  bool includeData= true;
  bool area_norm = false;
  double mcscale = -1;
  string xaxislabel = "#p_{T,#mu} (GeV)";
  string yaxislabel = "Event Rate";
  double min_x = 0, min_y = 0;
  double max_x = -1, max_y = -1;
  bool normalizeHisto = false;


  CCQENuPlotUtils *putils = new CCQENuPlotUtils( fluxHistoExists );

  MnvH1D* h_xvar[nHistos];
  putils->bookHistos(f,  h_xvar, Form("h_%s",xvar.c_str() ) );
  vector<MnvH1D*> h_xvar_ratio = GetRatio(h_xvar);

  c1->cd();

  max_y = h_xvar[kData]->GetMaximum()*1.3;
  putils->drawStacked(h_xvar, "QELike_split_PionInFS", area_norm, pot_data, pot_mc, includeData,mcscale,xaxislabel,yaxislabel,min_x, max_x, min_y ,max_y,normalizeHisto );
  c1->Print("plots/h-ptmu-fit_0-ratio_0.pdf");
  normalizeHisto = false;
  putils->drawStacked(&h_xvar_ratio[0], "QELike_split_PionInFS", area_norm, pot_data, pot_mc, includeData,mcscale,xaxislabel,yaxislabel,min_x, max_x, min_y ,2,normalizeHisto );
  c1->Print("plots/h-ptmu-fit_0-ratio_1.pdf");


  DoFit(h_xvar, fits );
  normalizeHisto = false;
  putils->drawStacked(h_xvar, "QELike_split_PionInFS", area_norm, pot_data, pot_mc, includeData,mcscale,xaxislabel,yaxislabel,min_x, max_x, min_y ,max_y,normalizeHisto );
  c1->Print("plots/h-ptmu-fit_1-ratio_0.pdf");
  normalizeHisto = false;
  h_xvar_ratio = GetRatio(h_xvar);
  putils->drawStacked(&h_xvar_ratio[0], "QELike_split_PionInFS", area_norm, pot_data, pot_mc, includeData,mcscale,xaxislabel,yaxislabel,min_x, max_x, min_y ,2,normalizeHisto );
  c1->Print("plots/h-ptmu-fit_1-ratio_1.pdf");

  MnvH2D* h_enuFit_xvar_mc = (MnvH2D*) f->Get("h_fitEnu_ptmu_mc");
  MnvH2D* h_enu_xvar_mc = (MnvH2D*) f->Get("h_enu_ptmu_mc");

  h_enuFit_xvar_mc->SetLineColor(kRed);
  h_enuFit_xvar_mc->ProjectionX()->Draw("hist");
  h_enu_xvar_mc->ProjectionX()->Draw("histsame");
  TLatex latex;
  latex.SetTextAlign(22);
  latex.SetTextFont(62);
  latex.SetTextSize(0.06);
  latex.DrawLatexNDC(0.5,0.04, "Reconstructed E_{#nu} (GeV)" );
  c1->Print("plots/h-enu-test.pdf");

  
}



//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
int main(int argc, char* argv[])
{
  //string f_signal = argv[1];
  //string f_blob = argv[2];
  //string f_micblob = argv[3];
  //string f_michel = argv[4];
  //string f_weighted = argv[5];
  //Get Files

  for (int i = 1; i<argc;i++ ) cout<<"Parameter "<<i<<" :"<<argv[i]<<endl;
  TFile* f_input = new TFile( argv[1],"read");
  TFile* f_weighted = new TFile( argv[2],"read");

 
  draw1DDistrosQuadrant( f_input, f_weighted, "h_pn_ptmu", "pn", "p_{n} (GeV/c)" );
  draw1DDistrosQuadrant( f_input, f_weighted, "h_dpt_ptmu", "dpt", "#delta p_{T} (GeV/c)" );
  draw1DDistrosQuadrant( f_input, f_weighted, "h_dptx_ptmu", "dptx", "#delta p_{Tx} (GeV/c)" );
  draw1DDistrosQuadrant( f_input, f_weighted, "h_dpty_ptmu", "dpty", "#delta p_{Ty} (GeV/c)" );
  draw1DDistrosQuadrant( f_input, f_weighted, "h_dalphat_ptmu", "dalphat", "#delta#alpha_{T} (degree)" );
  draw1DDistrosQuadrant( f_input, f_weighted, "h_dphit_ptmu", "dphit", "#delta#phi_{T} (degree)" );
  Draw1DHistograms( f_input, f_weighted, "ptmu" );
  DrawFits( f_input, f_weighted );
  Draw2DDistro( f_input );

  return 0;
}
