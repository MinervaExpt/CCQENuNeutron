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

vector<int> categories_to_fit({ kQELike_QE_OTH,  kQELike_RES, kQELike_2p2h, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion} );
vector<string> categories_to_fit_names({"qelike_qe_oth", "qelike_res", "qelike_2p2h","qelikenot_scp", "qelikenot_snp"});
vector<string> categories_to_fit_title_names({"QELike && QE OTH",  "QELike && RES", "QELike && 2p2h", "QELikeNot && Single Charged Pion", "QELikeNot && Single Neutral Pion"});

vector<int> histosUsed({ kData,  kMC, kQELike, kQELike_QE_H, kQELike_QE_OTH, kQELike_2p2h, kQELike_RES, kQELike_DIS, kQELike_OTH, kQELikeNot, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } ); 


//General Func
template<class T>
void FormatHistos(vector< T*> &h );
void FormatHistos( vector<vector<MnvH2D*>> h );
void drawSystematics( CCQENuPlotUtils*utils, MnvH1D* h,  double pot_data, double pot_mc, string name, bool asfrac=true);
void drawLatSystematicHist( MnvH1D* h,  string sysname, string savename="" ) ;

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
    h[kQELikeNot][i]->SetLineColor( kBlack );

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

  MnvH2D* h_blobMaxE_PartType[nHistos], *h_blobE_PartType[nHistos];
  MnvH2D* h_mainBlobMaxE_PartType[nHistos], *h_mainBlobE_PartType[nHistos];

  utils->bookHistos( f1, h_blobMaxE_PartType, "h_blobMaxE_PartType" );
  utils->bookHistos( f1, h_blobE_PartType, "h_blobE_PartType" );
  utils->bookHistos( f1, h_mainBlobMaxE_PartType, "h_mainBlobMaxE_PartType" );
  utils->bookHistos( f1, h_mainBlobE_PartType, "h_mainBlobE_PartType" );

  FormatPartType( h_blobMaxE_PartType[kMC], "Leading Cluster E (MeV)" );
  FormatPartType( h_mainBlobMaxE_PartType[kMC], "Leading Cluster E (MeV)" );
  FormatPartType( h_blobE_PartType[kMC], "Blob Total E (MeV)" );
  FormatPartType( h_mainBlobE_PartType[kMC], "Blob Total E (MeV)" );
  setBlackbodyPalette();
  TCanvas *c = new TCanvas();
  c->cd();
  h_blobMaxE_PartType[kMC]->Draw("colz");
  c->Print("./plots/particles-all-blobMaxE.pdf");
  h_blobE_PartType[kMC]->Draw("colz");
  c->Print("./plots/particles-all-blobE.pdf");

  h_mainBlobMaxE_PartType[kMC]->Draw("colz");
  c->Print("./plots/particles-mainCandidate-blobMaxE.pdf");
  h_mainBlobE_PartType[kMC]->Draw("colz");
  c->Print("./plots/particles-mainCandidate-blobE.pdf");



}
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
//Diagnostic Plots

void drawDiagnostic( MnvH2D** hists, string name, string xtitle, string axis="x", bool doRatio=false, double xmin = 0, double xmax=-1, bool doLog=false )
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
  for(unsigned int j = 0; j< nHistos;j++) histErr.push_back( NULL );
  for(auto c: histsToUse) histErr[c] = (TH2*)hists[c]->GetCVHistoWithError().Clone(Form("h_%s",names[c].c_str() ) ) ;

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
  histAndOpts.push_back( std::make_pair( histErr[kQELikeNot], "hist" ) );
  histAndOpts.push_back( std::make_pair( histErr[kQELikeNot_SingleNeutralPion], "hist" ) );
  histAndOpts.push_back( std::make_pair( histErr[kQELikeNot_SingleChargedPion], "hist" ) );
  histAndOpts.push_back( std::make_pair( histErr[kQELikeNot_MultiPion], "hist" ) );
  histAndOpts.push_back( std::make_pair( histErr[kQELikeNot_NoPions], "hist" ) );
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
  TLegend* leg=new TLegend(0.75, 0.1, .85, 0.3);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);

  TLegend* leg2=new TLegend(0.9, 0.1, 1, 0.3);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);

  leg->AddEntry( dataStat, "MINERvA data", "lpe");
  leg->AddEntry( mc, "MINERvA Tune", "l");
  leg->AddEntry( histErr[kQELike_QE_H],"QE-H","l");
  leg->AddEntry( histErr[kQELike_QE_OTH],"QE-Oth","l");
  leg->AddEntry( histErr[kQELike_RES],"QEL Resonant","l");
  leg->AddEntry( histErr[kQELike_2p2h],"QEL 2p2h","l");
  leg->AddEntry( histErr[kQELike_DIS],"QEL DIS","l");
  leg2->AddEntry( histErr[kQELikeNot_SingleNeutralPion],"1#pi^{0}","l");
  leg2->AddEntry( histErr[kQELikeNot_SingleChargedPion],"1#pi^{+/-}","l");
  leg2->AddEntry( histErr[kQELikeNot_MultiPion],"N#pi","l");
  leg2->AddEntry( histErr[kQELikeNot_NoPions],"other","l");



  double *multipliers = NULL;

  //First project into X axis
  string ytitle="Q^{2}_{QE} (GeV^{2})";

  if(axis == "x")
  {
    GridCanvas* gc = plotXAxis1D( histAndOpts, xtitle, ytitle, multipliers ); 
    gc->SetGridx(true);
    gc->SetYLimits(-.2,mcErr->GetMaximum()*1.5);
    if( xmin<xmax) gc->SetXLimits(xmin,xmax);
    if(doRatio)gc->SetYLimits(-.2,2.2);
    gc->SetYTitle("Evt Rate #times 10^{3}");
    if(doRatio) gc->SetYTitle("Ratio to MnvGENIE");
    if(!doRatio ) gc->SetLogy(doLog);
    gc->Modified();
    leg->Draw("SAME");
    leg2->Draw("SAME");

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
    if(!doRatio ) gc->SetLogy(doLog);
    gc->Modified();
    leg->Draw("SAME");
    leg2->Draw("SAME");
    string fname=Form("plots/diagnostic-%s-axis_y-1d-ratio_%d.pdf",name.c_str(), doRatio );
    gc->Print(fname.c_str());
  }

  free(dataStat);
  free(mcErr);
  free(mc);
  for(unsigned int h = 0; h< nHistos;h++) free(histErr[h]);

}
//void drawDiagnosticPar( map<pair<double,double>, vector<MnvH2D*>>hmap, string name, string xtitle, string axis="x", bool doRatio=false, double xmin = 0, double xmax=-1 )
//{
//  //Distribution in each Q2
//  vector<int> histsToUse({kData,kMC, kQELike, kQELikeNot, kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions});
//
//  TH2* dataStat = (TH2*) hists[kData]->GetCVHistoWithStatError().Clone("dataStat");
//  TH2* mcErr =(TH2*)  hists[kMC]->GetCVHistoWithError().Clone("mcErr");
//  TH2* mc =(TH2*)  hists[kMC]->GetCVHistoWithError().Clone("mc");
//  mcErr->SetFillStyle(1001);
//  mcErr->SetFillColorAlpha(kRed,0.35);
//  mcErr->SetMarkerSize(0 );
//
//  vector<TH2*>histErr;
//  cout<<"Creating histErr"<<endl;
//  for(unsigned int j = 0; j< nHistos;j++) histErr.push_back( NULL );
//  for(auto c: histsToUse) histErr[c] = (TH2*)hists[c]->GetCVHistoWithError().Clone(Form("h_%s",names[c].c_str() ) ) ;
//
//  //vector<MnvH2D*> vec_hists(hists, hists+nHistos );
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
//    gc->SetYLimits(-.2,mcErr->GetMaximum()*1.5);
//    if( xmin<xmax) gc->SetXLimits(xmin,xmax);
//    if(doRatio)gc->SetYLimits(-.2,2.2);
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
//    gc->SetYLimits(-.2,1.2);
//    if(doRatio)gc->SetYLimits(-.2,2.2);
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
//  for(unsigned int h = 0; h< nHistos;h++) free(histErr[h]);
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
  TVector2* pot = (TVector2*)f1->Get("pot");
  double pot_data = pot->X();
  double pot_mc = pot->Y();
  double pot_norm = pot_data/pot_mc;
  //double pot_norm = 1;


  vector<double> dropFrac({0.1,0.2,0.3});
  vector<double> dropEnergy({5, 10, 15,20});

  CCQENuPlotUtils *putils = new CCQENuPlotUtils( fluxHistoExists );
  CCQENuUtils *utils = new CCQENuUtils( false, fluxHistoExists );


  vector<MnvH2D*> fit_blobE = GetFitsHisto( f_signal_weighted, "blobEnergy" );
  vector<MnvH2D*> fit_blobN = GetFitsHisto( f_signal_weighted, "blobN" );
  vector<MnvH2D*> fit_blobType = GetFitsHisto( f_signal_weighted, "blobType" );
  vector<MnvH2D*> fit_blobDist = GetFitsHisto( f_signal_weighted, "blobDist" );

  map<pair<double,double>, vector<MnvH2D*>> 
    h_blobE_q2qe, 
    h_blobDist_q2qe, 
    h_blobN_q2qe,
    h_blobN2D_q2qe,
    h_blobN3D_q2qe,
    h_blobType_q2qe;
  map<pair<double,double>, vector<MnvH2D*>> 
    h_inner_blobE_q2qe, 
    h_inner_blobDist_q2qe, 
    h_inner_blobN_q2qe,
    h_inner_blobN2D_q2qe,
    h_inner_blobN3D_q2qe,
    h_inner_blobType_q2qe;

  vector<pair<double,double>> pars;
  pars.push_back({0,0});

  for( auto frac: dropFrac )
  {
    for( auto ene : dropEnergy )
    {
      pair<double,double> par = make_pair(frac,ene);
      pars.push_back(par);
      h_blobE_q2qe[par].resize(nHistos);
      h_blobDist_q2qe[par].resize(nHistos);
      h_blobN_q2qe[par].resize(nHistos);
      h_blobN2D_q2qe[par].resize(nHistos);
      h_blobN3D_q2qe[par].resize(nHistos);
      h_blobType_q2qe[par].resize(nHistos);

      h_inner_blobE_q2qe[par].resize(nHistos);
      h_inner_blobDist_q2qe[par].resize(nHistos);
      h_inner_blobN_q2qe[par].resize(nHistos);
      h_inner_blobN2D_q2qe[par].resize(nHistos);
      h_inner_blobN3D_q2qe[par].resize(nHistos);
      h_inner_blobType_q2qe[par].resize(nHistos);

      utils->bookHistos( f1, &(h_blobE_q2qe[par][0]), Form( "h_blobE_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_blobDist_q2qe[par][0]), Form( "h_blobDist_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_blobN_q2qe[par][0]), Form( "h_blobN_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_blobN2D_q2qe[par][0]), Form( "h_blobN2D_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_blobN3D_q2qe[par][0]), Form( "h_blobN3D_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_blobType_q2qe[par][0]), Form( "h_blobType_q2qe_%.2f_%.1f", frac,ene ) );

      utils->bookHistos( f1, &(h_inner_blobE_q2qe[par][0]), Form( "h_inner_blobE_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_inner_blobDist_q2qe[par][0]), Form( "h_inner_blobDist_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_inner_blobN_q2qe[par][0]), Form( "h_inner_blobN_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_inner_blobN2D_q2qe[par][0]), Form( "h_inner_blobN2D_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_inner_blobN3D_q2qe[par][0]), Form( "h_inner_blobN3D_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_inner_blobType_q2qe[par][0]), Form( "h_inner_blobType_q2qe_%.2f_%.1f", frac,ene ) );
    }
  }
      pair<double,double> par = make_pair(0,0);
      double frac=0,ene=0;
      h_blobE_q2qe[par].resize(nHistos);
      h_blobDist_q2qe[par].resize(nHistos);
      h_blobN_q2qe[par].resize(nHistos);
      h_blobN2D_q2qe[par].resize(nHistos);
      h_blobN3D_q2qe[par].resize(nHistos);
      h_blobType_q2qe[par].resize(nHistos);

      h_inner_blobE_q2qe[par].resize(nHistos);
      h_inner_blobDist_q2qe[par].resize(nHistos);
      h_inner_blobN_q2qe[par].resize(nHistos);
      h_inner_blobN2D_q2qe[par].resize(nHistos);
      h_inner_blobN3D_q2qe[par].resize(nHistos);
      h_inner_blobType_q2qe[par].resize(nHistos);

      utils->bookHistos( f1, &(h_blobE_q2qe[par][0]), Form( "h_blobE_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_blobDist_q2qe[par][0]), Form( "h_blobDist_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_blobN_q2qe[par][0]), Form( "h_blobN_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_blobN2D_q2qe[par][0]), Form( "h_blobN2D_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_blobN3D_q2qe[par][0]), Form( "h_blobN3D_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_blobType_q2qe[par][0]), Form( "h_blobType_q2qe_%.2f_%.1f", frac,ene ) );

      utils->bookHistos( f1, &(h_inner_blobE_q2qe[par][0]), Form( "h_inner_blobE_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_inner_blobDist_q2qe[par][0]), Form( "h_inner_blobDist_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_inner_blobN_q2qe[par][0]), Form( "h_inner_blobN_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_inner_blobN2D_q2qe[par][0]), Form( "h_inner_blobN2D_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_inner_blobN3D_q2qe[par][0]), Form( "h_inner_blobN3D_q2qe_%.2f_%.1f", frac,ene ) );
      utils->bookHistos( f1, &(h_inner_blobType_q2qe[par][0]), Form( "h_inner_blobType_q2qe_%.2f_%.1f", frac,ene ) );

      cout<<"Print histogram pointer address: "<<h_blobE_q2qe[par][0]<<endl;
      cout<<"Print histogram address: "<<&h_blobE_q2qe[par][0]<<endl;
  
  for(auto par: pars)
  {
    string parStr = Form( "%.2f_%.1f", par.first, par.second );
    for( int i = 0; i< parStr.size(); i++ )
    {
      if( parStr[i] == '.' ) parStr[i] = 'd';
    }

    drawDiagnostic( &h_blobE_q2qe[par][0], Form("blobE_%s", parStr.c_str() ), "E (MeV)", "x", false );
    drawDiagnostic( &h_blobDist_q2qe[par][0], Form("blobDist_%s", parStr.c_str() ), "R (mm)", "x", false,0,2500 );
    drawDiagnostic( &h_blobType_q2qe[par][0], Form("blobType_%s", parStr.c_str() ), "Type", "x", false,1,2e4,true );
    drawDiagnostic( &h_blobN_q2qe[par][0], Form("blobN_%s", parStr.c_str() ), "N Blobs", "x", false );
    drawDiagnostic( &h_blobN2D_q2qe[par][0], Form("blobN2D_%s", parStr.c_str() ), "N Blobs", "x", false );
    drawDiagnostic( &h_blobN3D_q2qe[par][0], Form("blobN3D_%s", parStr.c_str() ), "N Blobs", "x", false,0,8 );
    drawDiagnostic( &h_inner_blobN3D_q2qe[par][0], Form("blobN3Di_%s", parStr.c_str() ), "N Blobs", "x", false,0,8 );

    drawDiagnostic( &h_blobE_q2qe[par][0], Form("blobE_%s", parStr.c_str() ), "E (MeV)", "x", true );
    drawDiagnostic( &h_blobDist_q2qe[par][0], Form("blobDist_%s", parStr.c_str() ), "R (mm)", "x", true,0,2500 );
    drawDiagnostic( &h_blobType_q2qe[par][0], Form("blobType_%s", parStr.c_str() ), "Type", "x", true );
    drawDiagnostic( &h_blobN_q2qe[par][0], Form("blobN_%s", parStr.c_str() ), "N Blobs", "x", true );
    drawDiagnostic( &h_blobN2D_q2qe[par][0], Form("blobN2D_%s", parStr.c_str() ), "N Blobs", "x", true );
    drawDiagnostic( &h_blobN3D_q2qe[par][0], Form("blobN3D_%s", parStr.c_str() ), "N Blobs", "x", true,0,8 );
    drawDiagnostic( &h_inner_blobN3D_q2qe[par][0], Form("blobN3Di_%s", parStr.c_str() ), "N Blobs", "x", true,0,8);

    DoFit( &h_blobE_q2qe[par][0], fit_blobE );
    DoFit( &h_blobDist_q2qe[par][0], fit_blobDist );
    DoFit( &h_blobType_q2qe[par][0], fit_blobType );
    DoFit( &h_blobN_q2qe[par][0], fit_blobN );
    DoFit( &h_blobN2D_q2qe[par][0], fit_blobN );
    DoFit( &h_blobN3D_q2qe[par][0], fit_blobN );
    DoFit( &h_inner_blobN3D_q2qe[par][0], fit_blobN );

    drawDiagnostic( &h_blobE_q2qe[par][0],    Form("fit_blobE_%s", parStr.c_str() ), "E (MeV)", "x", false );
    drawDiagnostic( &h_blobDist_q2qe[par][0], Form("fit_blobDist_%s", parStr.c_str() ), "R (mm)", "x", false,0,2500 );
    drawDiagnostic( &h_blobType_q2qe[par][0], Form("fit_blobType_%s", parStr.c_str() ), "Type", "x", false,1,2e4,true );
    drawDiagnostic( &h_blobN_q2qe[par][0],    Form("fit_blobN_%s", parStr.c_str() ), "N Blobs", "x", false );
    drawDiagnostic( &h_blobN2D_q2qe[par][0],    Form("fit_blobN2D_%s", parStr.c_str() ), "N Blobs", "x", false );
    drawDiagnostic( &h_blobN3D_q2qe[par][0],    Form("fit_blobN3D_%s", parStr.c_str() ), "N Blobs", "x", false,0,8 );
    drawDiagnostic( &h_inner_blobN3D_q2qe[par][0],    Form("fit_blobN3Di_%s", parStr.c_str() ), "N Blobs", "x", false,0,8 );

    drawDiagnostic( &h_blobE_q2qe[par][0],    Form("fit_blobE_%s", parStr.c_str() ), "E (MeV)", "x", true );
    drawDiagnostic( &h_blobDist_q2qe[par][0], Form("fit_blobDist_%s", parStr.c_str() ), "R (mm)", "x", true,0,2500 );
    drawDiagnostic( &h_blobType_q2qe[par][0], Form("fit_blobType_%s", parStr.c_str() ), "Type", "x", true );
    drawDiagnostic( &h_blobN_q2qe[par][0],    Form("fit_blobN_%s", parStr.c_str() ), "N Blobs", "x", true );
    drawDiagnostic( &h_blobN2D_q2qe[par][0],    Form("fit_blobN2D_%s", parStr.c_str() ), "N Blobs", "x", true );
    drawDiagnostic( &h_blobN3D_q2qe[par][0],    Form("fit_blobN3D_%s", parStr.c_str() ), "N Blobs", "x", true,0,8 );
    drawDiagnostic( &h_inner_blobN3D_q2qe[par][0],    Form("fit_blobN3Di_%s", parStr.c_str() ), "N Blobs", "x", true,0,8 );
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

  MnvH2D *h_blobDist_q2qe[nHistos], *h_blobEnergy_q2qe[nHistos], *h_nBlobs_q2qe[nHistos],*h_n3DBlobs_q2qe[nHistos], *h_n2DBlobs_q2qe[nHistos];
  MnvH2D *h_blobEnergyLow_q2qe[nHistos];
  MnvH2D *h_nonBlobEnergy_q2qe[nHistos];
  MnvH2D *h_blobMaxE_q2qe[nHistos];
  MnvH2D *h_nClus_q2qe[nHistos];
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
  drawDiagnostic( h_blobEnergyLow_q2qe, "blobLowE"+tag, "E (MeV)", "x", false );
  drawDiagnostic( h_nonBlobEnergy_q2qe, "nonBlobE"+tag, "E (MeV)", "x", false );
  drawDiagnostic( h_blobDist_q2qe, "blobDist"+tag, "R (mm)", "x" , false, 0,2000);
  drawDiagnostic( h_n3DBlobs_q2qe, "n3DBlobs"+tag, "n 3DBlobs", "x" , false);
  drawDiagnostic( h_n2DBlobs_q2qe, "n2DBlobs"+tag, "n 2DBlobs", "x" , false);
  drawDiagnostic( h_nBlobs_q2qe, "nBlobs"+tag, "nBlobs", "x" , false);
  drawDiagnostic( h_nBlobs_q2qe, "nBlobs"+tag, "nBlobs", "y" , false);
  drawDiagnostic( h_blobMaxE_q2qe, "blobMaxE"+tag, "E_{blob} (GeV)", "x", false);
  drawDiagnostic( h_nClus_q2qe, "nClus"+tag, "#clusters", "x", false );

  drawDiagnostic( h_blobEnergy_q2qe, "blobE"+tag, "E (MeV)", "x", true );
  drawDiagnostic( h_blobEnergyLow_q2qe, "blobLowE"+tag, "E (MeV)", "x", true );
  drawDiagnostic( h_nonBlobEnergy_q2qe, "nonBlobE"+tag, "E (MeV)", "x", true );
  drawDiagnostic( h_blobDist_q2qe, "blobDist"+tag, "R (mm)", "x" , true, 0,2000);
  drawDiagnostic( h_n3DBlobs_q2qe, "n3DBlobs"+tag, "n 3DBlobs", "x" , true);
  drawDiagnostic( h_n2DBlobs_q2qe, "n2DBlobs"+tag, "n 2DBlobs", "x" , true);
  drawDiagnostic( h_nBlobs_q2qe, "nBlobs"+tag, "nBlobs", "x" , true);
  drawDiagnostic( h_nBlobs_q2qe, "nBlobs"+tag, "nBlobs", "y" , true);

  drawDiagnostic( h_blobMaxE_q2qe, "blobMaxE"+tag, "E_{leading cluster} (GeV)", "x", true);
  drawDiagnostic( h_nClus_q2qe, "nClus"+tag, "#clusters", "x", true );

  vector<MnvH2D*> fit_blobE = GetFitsHisto( f_signal_weighted, "blobEnergy" );
  vector<MnvH2D*> fit_blobDist = GetFitsHisto( f_signal_weighted, "blobDist" );
  vector<MnvH1D*> fit_q2 = GetFitsHisto( f_signal_weighted );
  //vector<MnvH2D*> fit_nblob = GetFitsHisto( f_signal_weighted, "nBlobs" );

  DoFit( h_blobEnergy_q2qe, fit_blobE );
  DoFit( h_blobDist_q2qe, fit_blobDist );
  DoFit( h_nBlobs_q2qe, fit_q2 );
  DoFit( h_n2DBlobs_q2qe, fit_q2 );
  DoFit( h_n3DBlobs_q2qe, fit_q2 );
  DoFit( h_blobMaxE_q2qe, fit_q2 );
  DoFit( h_nClus_q2qe, fit_q2 );
  DoFit( h_blobEnergyLow_q2qe, fit_q2 );

  drawDiagnostic( h_blobEnergy_q2qe, "blobE-fit_1"+tag, "E (MeV)", "x", false );
  drawDiagnostic( h_blobEnergyLow_q2qe, "blobLowE-fit_1"+tag, "E (MeV)", "x", false );
  drawDiagnostic( h_blobDist_q2qe, "blobDist-fit_1"+tag, "R (mm)", "x" , false, 0,2000);
  drawDiagnostic( h_n3DBlobs_q2qe, "n3DBlobs-fit_1"+tag, "n 3DBlobs", "x" , false);
  drawDiagnostic( h_n2DBlobs_q2qe, "n2DBlobs-fit_1"+tag, "n 2DBlobs", "x" , false);
  drawDiagnostic( h_nBlobs_q2qe, "nBlobs-fit_1"+tag, "n Blobs", "x" , false);

  drawDiagnostic( h_blobMaxE_q2qe, "blobMaxE-fit_1"+tag, "E_{leading cluster} (MeV)","x",false );
  drawDiagnostic( h_nClus_q2qe, "nClus-fit_1"+tag, "#clusters", "x",false );

  drawDiagnostic( h_blobEnergy_q2qe, "blobE-fit_1"+tag, "E (MeV)", "x", true );
  drawDiagnostic( h_blobEnergyLow_q2qe, "blobLowE-fit_1"+tag, "E (MeV)", "x", true );
  drawDiagnostic( h_blobDist_q2qe, "blobDist-fit_1"+tag, "R (mm)", "x" , true, 0,2000);
  drawDiagnostic( h_n3DBlobs_q2qe, "n3DBlobs-fit_1"+tag, "n 3DBlobs", "x" , true);
  drawDiagnostic( h_n2DBlobs_q2qe, "n2DBlobs-fit_1"+tag, "n 2DBlobs", "x" , true);
  drawDiagnostic( h_nBlobs_q2qe, "nBlobs-fit_1"+tag, "n Blobs", "x" , true);

  drawDiagnostic( h_blobMaxE_q2qe, "blobMaxE-fit_1"+tag, "E_{leading cluster} (MeV)","x",true);
  drawDiagnostic( h_nClus_q2qe, "nClus-fit_1"+tag, "#clusters", "x",true );



  return;

}

//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
int main(int argc, char* argv[])
{

  cout<<"main"<<endl;
  string f_orig = argv[1];
  string f_signal_weighted = argv[2];
  TFile* file_orig = new TFile( f_orig.c_str(),"read");
  TFile* file_signal_weighted = new TFile( f_signal_weighted.c_str(),"read");
  cout<<"acquired files"<<endl;


  //Diagnostics
  NeutronDiagnostics( file_orig, file_signal_weighted );
  return 0;
}
