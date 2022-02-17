//#include "myPlotStyle.h"
#include "TParameter.h"
#include "include/GeneralIncludes.h"

#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/MnvPlotter.h"

#include "TMinuit.h"
#include "TFitter.h"

#include "FA.h"
#include "FA05.h"
#include "Zexp.h"

#include "TMultiGraph.h"

#include "TSpline.h"


using namespace PlotUtils;

using namespace BBA05;
//forward declaration
//gROOT->SetBatch();



double xs_xmin = 0.0;
double xs_xmax = 10.0;

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


TH1D* weightedQ2;

map<double,double> WeightedQ2({
{ 0.003125 , 0.003125 },
{ 0.009375 , 0.009375 },
{ 0.01875  , 0.01875  },
{ 0.03125  , 0.03125  },
{ 0.04375  , 0.04375  },
{ 0.075    , 0.075    },
{ 0.125    , 0.125    },
{ 0.175    , 0.175    },
{ 0.25     , 0.25     },
{ 0.35     , 0.3499   },
{ 0.5      , 0.4995   },
{ 0.7      , 0.6993   },
{ 0.9      , 0.8991   },
{ 1.1      , 1.099    },
{ 1.6      , 1.582    },
{ 3        , 2.815    },
{ 5        , 4.768    },
{ 8        , 7.596   } } );


double GetWeightedQ2Center( double Q2 )
{
  int ibin = weightedQ2->FindBin(Q2);
  return weightedQ2->GetBinContent(ibin);
}

MnvH1D* xs_data;

//vector<int> categories_to_fit({ kQELike_QE_OTH, kQELike_2p2h, kQELikeNot_SinglePion} );
//vector<string> categories_to_fit_names({"qelike_qe_oth", "qelike_2p2h","qelikenot_sp"});
//vector<string> categories_to_fit_title_names({"QELike && QE OTH", "QELike && 2p2h", "QELikeNot && Single Pion"});


vector<int> histosUsed({ kData,  kMC, kQELike, kQELike_QE_H, kQELike_QE_OTH, kQELike_2p2h, kQELike_RES, kQELike_DIS, kQELike_OTH, kQELikeNot, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } ); 


void DrawHorizontalLine(double y, double x0, double x1, TCanvas *c)
{
  c->cd();
  TLine *line = new TLine();
  line->SetLineColor( kRed );
  line->SetLineWidth( 2 );
  line->DrawLine(x0,y,x1,y);
}

//----Migration Matrix-----
TH2D* GetMigHistos(TFile*f, string name )
{
  TH2D* ret = (TH2D*) f->Get( (name+"_qelike_qe_h").c_str() );
  return ret;
}

vector<TH2D*> GetMigHistos(TFile*f, string xname, string yname, string tname )
{
  TH2D* migration = (TH2D*) f->Get( Form("h_%s_%s_%s_migration", xname.c_str(), yname.c_str(), tname.c_str() ) );
  TH2D* reco = (TH2D*) f->Get( Form("h_%s_%s_%s_reco", xname.c_str(), yname.c_str(), tname.c_str() ) );
  TH2D* truth = (TH2D*) f->Get( Form("h_%s_%s_%s_truth", xname.c_str(), yname.c_str(), tname.c_str() ) );
  return vector<TH2D*>({migration,reco,truth});
}

void drawMig( TCanvas* c, TH2D* h, string name, string var, bool log=true )
{

  c->cd();
  h->GetXaxis()->SetTitle( Form(" Reconstructed %s", var.c_str() ) );
  h->GetYaxis()->SetTitle( Form(" Truth %s", var.c_str()) );
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetTitleOffset(1.3);
  h->GetYaxis()->SetTitleOffset(1.3);

  h->Draw("colz");

  c->SetLogx( log );
  c->SetLogy( log );
  //c->SetGrid();
  c->Print( Form("./plots/migration-%s.pdf",name.c_str() ) );
}

TH2D* combineMigration( TH2D* migration, TH2D *reco, bool xaxis = true )
{
  int nX = reco->GetNbinsX();
  int nY = reco->GetNbinsY();
  int nBin = (xaxis)? nX : nY;
  TH2D* ret = new TH2D("migration","migration", nBin, 0, nBin, nBin,0,nBin );
  int sizeChunk = nX+2;

  for( int i = 0; i< migration->GetNbinsX(); i++ )
  {
    int xbin = i%sizeChunk;
    for( int j = 0; j< migration->GetNbinsY(); j++ )
    {
      int ybin = j%sizeChunk;
      double value = ret->GetBinContent( xbin, ybin );
      value+= migration->GetBinContent( i+1, j+1 );
      ret->SetBinContent( xbin, ybin, value );
    }
  }
  return ret;
}

void draw2DMig( TFile*f, TCanvas* c, string xvar, string yvar, string var="qelike_qe_h", bool log=true )
{
  vector<TH2D*> histos = GetMigHistos(f, xvar, yvar, var);
  drawMig( c, histos[0], Form("%s-%s",xvar.c_str(), yvar.c_str() ),"",false );
  TH2D* mig_row = Normalize( histos[0], false );
  drawMig( c, mig_row, Form( "%s-%s-rowNorm", xvar.c_str(), yvar.c_str() ),"",false );

  TH2D* mig_1D = combineMigration( histos[0], histos[1] );
  drawMig( c, mig_1D, Form("%s",xvar.c_str() ),"",false );
  TH2D* mig_1D_row = Normalize( mig_1D, false );
  drawMig( c, mig_1D_row, Form("%s-rowNorm",xvar.c_str() ),"",false );
}

void DrawMig( TFile *f )
{
  cout<<"Entering Draw Mig"<<endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  //setRedPalette();
  setBlackbodyPalette();
  //setCorrelationPalette(0.5);
  //setRainbowToWhitePalette();
  
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.15);

  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetEndErrorSize(2);

  TCanvas *c = new TCanvas("c","c",400,400);

  TH2D* mig_q2qe = GetMigHistos(f, "h_q2qe_angle_10_migration");
  TH2D* mig_q2qe_row = Normalize( mig_q2qe, false );

  mig_q2qe_row->SetMaximum(1);

  MnvPlotter plotter;
  plotter.SetROOT6Palette(55);
  drawMig( c, mig_q2qe_row, "q2qe-rowNorm", "Q^{2} (GeV^{2})" );

  draw2DMig(f, c, "muonmomentum", "q2qe" );
  draw2DMig(f, c, "muontheta", "q2qe" );
  draw2DMig(f, c, "enu", "q2qe" );
  draw2DMig(f, c, "dthetaR", "q2qe", "qelike_qe_h" );
  draw2DMig(f, c, "dthetaP", "q2qe", "qelike_qe_h" );
  
  
  //vector<TH2D*> histos_muonmomentum_q2qe = GetMigHistos(f, "muonmomentum", "q2qe", "qelike_qe_h");
  //drawMig( c, histos_muonmomentum_q2qe[0], "muonmomentum-q2qe","",false );
  //TH2D* mig_muonmomentum_q2qe_row = Normalize( histos_muonmomentum_q2qe[0], false );
  //drawMig( c, mig_muonmomentum_q2qe_row, "muonmomentum-q2qe-rowNorm","",false );

  //vector<TH2D*> histos_muontheta_q2qe = GetMigHistos(f, "muontheta", "q2qe", "qelike_qe_h");
  //drawMig( c, histos_muontheta_q2qe[0], "muontheta-q2qe","",false );
  //TH2D* mig_muontheta_q2qe_row = Normalize( histos_muontheta_q2qe[0], false );
  //drawMig( c, mig_muontheta_q2qe_row, "muontheta-q2qe-rowNorm","",false );

  //vector<TH2D*> histos_dthetaR_q2qe = GetMigHistos(f, "dthetaR", "q2qe", "qelike_qe_h");
  //drawMig( c, histos_dthetaR_q2qe[0], "dthetaR-q2qe","",false );
  //TH2D* mig_dthetaR_q2qe_row = Normalize( histos_dthetaR_q2qe[0], false );
  //drawMig( c, mig_dthetaR_q2qe_row, "dthetaR-q2qe-rowNorm","",false );

  //vector<TH2D*> histos_dthetaP_q2qe = GetMigHistos(f, "dthetaP", "q2qe", "qelike_qe_h");
  //drawMig( c, histos_dthetaP_q2qe[0], "dthetaP-q2qe","",false );
  //TH2D* mig_dthetaP_q2qe_row = Normalize( histos_dthetaP_q2qe[0], false );
  //drawMig( c, mig_dthetaP_q2qe_row, "dthetaP-q2qe-rowNorm","",false );

  c->Close();
}


//================= Efficiency ==============

MnvH1D* GetEfficiency( TFile* f )
{
  cout<<"Getting Efficiency 1D"<<endl;
  MnvH1D* num = (MnvH1D*) f->Get("h_q2qe_angle_10_qelike_qe_h");
  MnvH1D* den = (MnvH1D*) f->Get("h_truth_q2qe_qelike_qe_h");
  cout<<num->GetName()<<endl;
  cout<<den->GetName()<<endl;
  den->AddMissingErrorBandsAndFillWithCV( *num );
  MnvH1D* eff = (MnvH1D*) num->Clone("eff");
  eff->Divide( eff, den );
  return eff;
}
MnvH2D* GetEfficiency( TFile* f, string xvar, string yvar, string type )
{
  MnvH2D* num = (MnvH2D*) f->Get( Form("h_%s_%s_%s", xvar.c_str(),  yvar.c_str(), type.c_str() ));
  MnvH2D* den = (MnvH2D*) f->Get( Form("h_truth_%s_%s_%s", xvar.c_str(),  yvar.c_str(), type.c_str() ));
  den->AddMissingErrorBandsAndFillWithCV( *num );
  MnvH2D* eff = (MnvH2D*) num->Clone("eff");
  eff->Divide( eff, den );
  return eff;
}

void Draw1D( TCanvas*c, MnvH1D* h, string savename, string xtitle, bool logx=true, bool drawAvg=false )
{
  c->cd();
  TH1D hsys = h->GetCVHistoWithError();
  TH1D hstat = h->GetCVHistoWithStatError();
  hsys.GetXaxis()->SetTitle( xtitle.c_str() );
  hsys.GetYaxis()->SetTitle( "Efficiency" );

  hsys.GetXaxis()->SetTitleSize( 0.05 );
  hsys.GetYaxis()->SetTitleSize( 0.05 );
  hsys.GetXaxis()->SetTitleOffset( 1.3 ); 
  hsys.GetYaxis()->SetTitleOffset( 1.3 );


  hsys.SetMinimum(0);

  hsys.SetLineWidth(2);
  hstat.SetLineWidth(2);


  hsys.Draw("e1");
  hstat.Draw("e1same");


  if(drawAvg)
  {
    hstat.Fit("pol0","","",0.05,2);
    hstat.Draw("e1same");

  }




  c->SetLogx( logx );
  c->SetGrid();
  c->Print( Form("plots/eff-%s.pdf", savename.c_str() ) );
}

void Draw2D( TCanvas*c, MnvH2D* h, string savename, string xtitle, string celltitle, double xmin, double xmax)
{
  c->cd();
  int binmin = h->GetYaxis()->FindBin(xmin);
  int binmax = h->GetYaxis()->FindBin(xmax);
  TH2* hsys =  (TH2*) h->GetCVHistoWithError().Clone("effErr");
  TH2* hstat = (TH2*) h->GetCVHistoWithStatError().Clone("effStat");
  vector<pair<TH2*, const char*> > histAndOpts;
  histAndOpts.push_back( std::make_pair( hsys, "hist" ) );
  histAndOpts.push_back( std::make_pair( hstat, "hist" ) );
  GridCanvas *gc = plotXAxis1D( histAndOpts, xtitle, celltitle,NULL,binmin,binmax );
  gc->SetGridx(true);
  gc->SetYTitle( "Efficiency" );
  gc->Print( Form("plots/eff-%s.pdf", savename.c_str() ) );
}

void DrawEff( TFile *f, string tag="" )
{
  cout<<"Entering Draw Eff"<<endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  //setRedPalette();
  setBlackbodyPalette();
  //setCorrelationPalette(0.5);
  //setRainbowToWhitePalette();
  
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetEndErrorSize(2);
  gStyle->SetTitleOffset(1.2,"Y");

  cout<<"creating new canvas"<<endl;
  TCanvas *c = new TCanvas("c","c",400,400);
  MnvH1D* eff_q2qe = GetEfficiency( f );
  cout<<"got efficiency, now draw 1D"<<endl;
  Draw1D(c, eff_q2qe, "q2qe"+tag, "Q^{2}_{QE} (GeV^{2})", true);
  Draw1D(c, eff_q2qe, "q2qe_avg"+tag, "Q^{2}_{QE} (GeV^{2})", true, true);

  string celltitle = "Q^{2}";
  MnvH2D* eff_dthetaP_q2qe =  GetEfficiency( f, "dthetaP", "q2qe", "qelike_qe_h" );
  Draw2D( c, eff_dthetaP_q2qe, "dthetaP-q2qe"+tag,"#delta #theta_{R}", celltitle,0.05,1.4 );

  MnvH2D* eff_dthetaR_q2qe =  GetEfficiency( f, "dthetaR", "q2qe", "qelike_qe_h" );
  Draw2D( c, eff_dthetaR_q2qe, "dthetaR-q2qe"+tag,"#delta #theta_{P}", celltitle,0.05,1.4 );

  MnvH2D* eff_muonmomentum_q2qe =  GetEfficiency( f, "muonmomentum", "q2qe", "qelike_qe_h" );
  Draw2D( c, eff_muonmomentum_q2qe, "muonmomentum-q2qe"+tag,"p_{#mu} (GeV)", celltitle,0.05,1.4 );

  MnvH2D* eff_muontheta_q2qe =  GetEfficiency( f, "muontheta", "q2qe", "qelike_qe_h" );
  Draw2D( c, eff_muontheta_q2qe, "muontheta-q2qe"+tag,"#theta_{#mu}", celltitle,0.05,1.4 );
}

//=========================================================
// Unfolding and Xs
vector<MnvH1D*> GetAllHistos( TFile* f )
{
  MnvH1D* data = (MnvH1D*) f->Get("h_q2qe_region_00_data_nobck_hydrogen");
  MnvH1D* mc = (MnvH1D*) f->Get("h_q2qe_region_00_qelike_qe_h_nobck_hydrogen");
  MnvH1D* data_unfold = (MnvH1D*) f->Get("h_q2qe_region_00_data_nobck_hydrogen_unfold");
  MnvH1D* mc_unfold = (MnvH1D*) f->Get("h_q2qe_region_00_qelike_qe_h_nobck_hydrogen_unfold");
  MnvH1D* data_unfold_effcor = (MnvH1D*) f->Get("h_q2qe_region_00_data_nobck_hydrogen_unfold_effcor");
  MnvH1D* mc_unfold_effcor = (MnvH1D*) f->Get("h_q2qe_region_00_qelike_qe_h_nobck_hydrogen_unfold_effcor");
  MnvH1D* data_unfold_effcor_xs = (MnvH1D*) f->Get("h_q2qe_region_00_data_nobck_hydrogen_unfold_effcor_cross_section");
  MnvH1D* mc_unfold_effcor_xs = (MnvH1D*) f->Get("h_q2qe_region_00_qelike_qe_h_nobck_hydrogen_unfold_effcor_cross_section");
  MnvH1D* eff_denum_xs = (MnvH1D*) f->Get("h_effhist_demH_cross_section");

  eff_denum_xs->SetName("h_eff_denum_xs");
  vector<MnvH1D*> ret({ data,mc, data_unfold, mc_unfold, data_unfold_effcor, mc_unfold_effcor, data_unfold_effcor_xs, mc_unfold_effcor_xs, eff_denum_xs} );

  for ( auto h : ret )
  {
    h->GetVertErrorBand("Low_Recoil_2p2h_Tune")->SetUnivWgt(2,0.5);
  }

  for ( int i = 2; i<ret.size(); i++ )
  {
    ret[i]->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
  }
  data_unfold_effcor_xs->GetYaxis()->SetTitle("d#sigma/dQ^{2} (cm^{2}/GeV^{2}/Hydrogen)");
  mc_unfold_effcor_xs->GetYaxis()->SetTitle("d#sigma/dQ^{2} (cm^{2}/GeV^{2}/Hydrogen)");

  return ret;
}

vector<MnvH2D*> GetAllHistos2D( TFile* f )
{
  MnvH2D* data = (MnvH2D*) f->Get("h_enu_q2qe_data_nobck_hydrogen");
  MnvH2D* mc = (MnvH2D*) f->Get("h_enu_q2qe_qelike_qe_h_nobck_hydrogen");
  MnvH2D* data_unfold = (MnvH2D*) f->Get("h_enu_q2qe_data_nobck_hydrogen_unfold");
  MnvH2D* mc_unfold = (MnvH2D*) f->Get("h_enu_q2qe_qelike_qe_h_nobck_hydrogen_unfold");
  MnvH2D* data_unfold_effcor = (MnvH2D*) f->Get("h_enu_q2qe_data_nobck_hydrogen_unfold_effcor");
  MnvH2D* mc_unfold_effcor = (MnvH2D*) f->Get("h_enu_q2qe_qelike_qe_h_nobck_hydrogen_unfold_effcor");
  MnvH2D* data_unfold_effcor_xs = (MnvH2D*) f->Get("h_enu_q2qe_data_nobck_hydrogen_unfold_effcor_cross_section");
  MnvH2D* mc_unfold_effcor_xs = (MnvH2D*) f->Get("h_enu_q2qe_qelike_qe_h_nobck_hydrogen_unfold_effcor_cross_section");

  MnvH2D* eff_denum_xs = NULL; //(MnvH2D*) f->Get("h_effhist_demH_cross_section");
  //eff_denum_xs->SetName("h_eff_denum_xs");
  return vector<MnvH2D*>({ data,mc, data_unfold, mc_unfold, data_unfold_effcor, mc_unfold_effcor, data_unfold_effcor_xs, mc_unfold_effcor_xs, eff_denum_xs} );
}



template<class T>
void FormatHistos( vector<T *> &hists )
{
  for( uint i = 0; i< hists.size(); i++ )
  {
    if( i%2 == 0 ) 
    {
      hists[i]->SetLineColor( kBlack );
      hists[i]->SetMarkerStyle( 23 );
      hists[i]->SetMarkerColor( kBlack );
      hists[i]->SetMarkerSize( 0.03 );
    }
    else hists[i]->SetLineColor( kBlue );
  }

  hists[0]->SetMarkerStyle(22);

  hists[2]->SetLineColor(kRed);
  hists[2]->SetMarkerColor(kRed);
  hists[2]->SetLineStyle(1);
  hists[2]->SetMarkerStyle(32);

  hists[3]->SetLineStyle(2);

}
template void FormatHistos(vector<MnvH1D*> &hists );
template void FormatHistos(vector<MnvH2D*> &hists );

void drawUnfold( TCanvas* c, vector<MnvH1D*> &hists, int iteration )
{
  gStyle->SetErrorX(0.001);
  c->cd();
  TLegend *leg = new TLegend(.7,.7,.9,.9);
  leg->AddEntry( hists[0], "Data" );
  leg->AddEntry( hists[1], "MC" );
  leg->AddEntry( hists[2], "Data Unfolded" );
  leg->AddEntry( hists[3], "MC Unfolded" );

  hists[0]->Draw("pe1");
  hists[1]->Draw("histsame");
  hists[2]->Draw("pe1same");
  hists[3]->Draw("histsame");
  leg->Draw();
  c->SetLogx();
  c->SetGrid();
  c->Print(Form("plots/xs-unfold-q2qe-it_%d.pdf", iteration));
}

void drawCorrelation( TCanvas* c, vector<MnvH1D*> &hists, int iteration )
{
  gStyle->SetErrorX(0.001);
  c->cd();
  c->SetLogx(false);
  c->SetLogy(false);
  c->SetLogz(true);
  MnvH1D* data = (MnvH1D*) hists[6]->Clone("errormatrix");
  data->Scale(1,"width");

  TH2D errMatrix_all =   TH2D(data->GetTotalErrorMatrix( true,   false, false));
  //6-14
  for( int i =0 ; i< 20; i ++ )
  {
    bool doClean = (i<6 || i>14);
    if( doClean )
    {
      for( int j=0;j<20;j++)
      {
       errMatrix_all.SetBinContent(i,j,0);
       errMatrix_all.SetBinContent(j,i,0);
      }

    }
  }

  TH2D errMatrix_sys =   TH2D(data->GetTotalErrorMatrix( false,  false, false ));
  TH2D errMatrix_stat =  TH2D(data->GetStatErrorMatrix( false ));


  double max = errMatrix_all.GetMaximum();
  double factor = 1;
  double min = pow(10,-84);
  errMatrix_all.GetZaxis()->SetRangeUser(min, max*factor);
  errMatrix_all.Draw("colz");
  c->Print(Form("plots/ErrorMatrix-all-it_%d.pdf", iteration));
  errMatrix_sys.GetZaxis()->SetRangeUser(min, max*factor);
  errMatrix_sys.Draw("colz");
  c->Print(Form("plots/ErrorMatrix-sys-it_%d.pdf", iteration));
  errMatrix_stat.GetZaxis()->SetRangeUser(min, max*factor);
  errMatrix_stat.Draw("colz");
  c->Print(Form("plots/ErrorMatrix-stat-it_%d.pdf", iteration));

  TMatrixD corrMatrix_all = hists[6]->GetTotalCorrelationMatrix( false,  true );
  TMatrixD corrMatrix_sys = hists[6]->GetTotalCorrelationMatrix( false,  false);
  TMatrixD corrMatrix_stat = corrMatrix_all- corrMatrix_sys;

  c->SetLogz(false);
  corrMatrix_all.Draw("colz");
  c->Print(Form("plots/CorrelationMatrix-all-it_%d.pdf", iteration));
  corrMatrix_sys.Draw("colz");
  c->Print(Form("plots/CorrelationMatrix-sys-it_%d.pdf", iteration));
}

void drawEnuQ2( TCanvas* c  )
{
  gStyle->SetErrorX(0.001);
  c->cd();
  c->SetLogx(false);
  c->SetLogy(false);
  //c->SetLogz(true);

  TH2D enuQ2 = BBA05::h_enu_q2qe_mc->GetCVHistoWithStatError();
  enuQ2.Scale(1,"width");
  enuQ2.GetXaxis()->SetTitle("E_{#nu} (GeV)");
  enuQ2.GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
  enuQ2.GetZaxis()->SetTitle("N/GeV^{3}");
  enuQ2.GetXaxis()->SetRangeUser(0,10);
  enuQ2.GetYaxis()->SetRangeUser(xs_xmin, xs_xmax);


  double max = enuQ2.GetMaximum();
  double factor = 1;
  double min = 0.;
  enuQ2.GetZaxis()->SetRangeUser(min,max*1.1);
  enuQ2.Draw("colz");
  c->Print(Form("plots/EnuQ2.pdf"));
}




TH1D* GetXsHisto( double MA, bool useKevin, int form = 0 )//0: BBA05, 1: BBBA07
{
  TH1D* h = new TH1D(BBA05::h_q2qe_base->GetCVHistoWithError());
  h->Reset();
  h->SetName(Form("h_%d",form) );

  for( int i = 1; i< h->GetNbinsX()+1; i++ )
  {
    double q2 = h->GetBinCenter( i );
    bool isNeutrino = false;

    double fa = BBA05::FA(q2,MA);
    double xs = (form == 0)? BBA05::XsIntegral(q2,fa,isNeutrino, useKevin) : BBBA07::XsIntegral(q2,fa,isNeutrino, useKevin);
    h->SetBinContent(i, xs );
  }
  return h;
}

TH1D* GetZexpXsecHisto( int mode )//0: Meyer, 1: Hydrogen
{
  TH1D* h = new TH1D(BBA05::h_q2qe_base->GetCVHistoWithError());
  h->Reset();
  h->SetName(Form("h_zexp") );

  for( int i = 1; i< h->GetNbinsX()+1; i++ )
  {
    double q2 = h->GetBinCenter( i );
    bool isNeutrino = false;

    double fa = Zexp::MeyerFitCV(q2);
    double faerr = Zexp::MeyerFitError(q2);
    if( mode == 1 )
    {
      fa = Zexp::TejinFitCV(q2);
      faerr = Zexp::TejinFitError(q2);
    }
    //double xs = (form == 0)? BBA05::XsIntegral(q2,fa,isNeutrino, useKevin) : BBBA07::XsIntegral(q2,fa,isNeutrino, useKevin);
    double xs = BBA05::XsIntegral(q2,fa,isNeutrino, false);
    double xserr1=BBA05::XsIntegral(q2,fa+faerr,isNeutrino, false);
    double xserr2=BBA05::XsIntegral(q2,fa-faerr,isNeutrino, false);
    double dxsec = abs( xserr1 -xs) + abs(xserr2-xs);
    dxsec/=2;

    h->SetBinContent(i, xs );
    h->SetBinError(i, dxsec );
  }
  return h;
}

TGraphErrors GetZexpXsecGraph( int mode, bool ratio=false, double Ma = 1.015 )//0: Meyer, 1: Hydrogen
{
  //TGraphErrors Fa_BBA05G, Fa_BBBA07G;
  //for( int i = 1; i< Fa_BBA05->GetNbinsX()+1; i++ )
  //{
  //  double Q2 = Fa_BBA05->GetBinCenter(i);
  //  Fa_BBA05G.SetPoint(i,GetWeightedQ2Center(Q2), Fa_BBA05->GetBinContent(i) );
  //  Fa_BBA05G.SetPointError(i, 0, Fa_BBA05->GetBinError(i) );

  //  Fa_BBBA07G.SetPoint(i, GetWeightedQ2Center(Q2), Fa_BBBA07->GetBinContent(i) );
  //  Fa_BBBA07G.SetPointError(i, 0, Fa_BBBA07->GetBinError(i) );
  //}

  TGraphErrors graph;

  double q2Max=10;
  double q2Min=0.01;
  int nSteps = 1000;
  double step = q2Max/nSteps;
  for( int i = 1; i< nSteps; i++ )
  {
    double q2 = q2Min+step*i;
    bool isNeutrino = false;

    double fa, faerr;
    if( mode == 0 )
    {
      fa = Zexp::MeyerFitCV(q2);
      faerr = Zexp::MeyerFitError(q2);
    }
    else if( mode == 1 )
    {
      fa = Zexp::TejinFitCV(q2);
      faerr = Zexp::TejinFitError(q2);
    }
    //double xs = (form == 0)? BBA05::XsIntegral(q2,fa,isNeutrino, useKevin) : BBBA07::XsIntegral(q2,fa,isNeutrino, useKevin);
    double xs = BBA05::XsIntegral(q2,fa,isNeutrino, false);
    double xserr1=BBA05::XsIntegral(q2,fa+faerr,isNeutrino, false);
    double xserr2=BBA05::XsIntegral(q2,fa-faerr,isNeutrino, false);
    double dxsec = abs( xserr1 -xs) + abs(xserr2-xs);
    dxsec/=2;

    if( ratio )
    {
      double dipole = BBA05::FA(q2,Ma);
      bool isNeutrino=false, useKevin = false;
      double dipolexs = BBA05::XsIntegral(q2,dipole,isNeutrino, useKevin);
      xs/=dipolexs;
      dxsec/=dipolexs;
    }

    graph.SetPoint(i, q2, xs );
    graph.SetPointError(i, 0, dxsec );
  }
  return graph;
}

double GetXsValueBBA05( double Q2, double Fa)//0: BBA05, 1: BBBA07
{
  bool isNeutrino = false, useKevin = false;
  double xs = BBA05::XsIntegral(Q2,Fa,isNeutrino, useKevin);
  return xs;
}

double GetXsValueBBBA07( double Q2, double Fa)//0: BBA05, 1: BBBA07
{
  bool isNeutrino = false, useKevin = false;
  double xs = BBBA07::XsIntegral(Q2,Fa,isNeutrino, useKevin);
  return xs;
}

double GetXsValueBBA05Wrapper( double *x, double *par )//solving for Ma using fixed Q2
{
  double Q2 = par[0];
  double Fa = x[0];
  return GetXsValueBBA05(Q2,Fa)*1e40;
}
double GetXsValueBBBA07Wrapper( double *x, double *par )//solving for Ma using fixed Q2
{
  double Q2 = par[0];
  double Fa = x[0];
  return GetXsValueBBBA07(Q2,Fa)*1e40;
}

TF1* FAFuncBBA05 = new TF1("MaBBA05", GetXsValueBBA05Wrapper, BBA05::gA*1.5,0,1 );
TF1* FAFuncBBBA07 = new TF1("MaBBBA07", GetXsValueBBBA07Wrapper, BBA05::gA*1.5,0,1 );

double GetFaJacobianBBA05Wrapper(double*x, double*par)
{
  double Q2 = x[0];
  double Fa = x[1];
  double ret=BBA05::XsIntegral(Q2, Fa, false, false);
  return 1/ret;
}
double GetFaJacobianBBBA07Wrapper(double*x, double*par)
{
  double Q2 = x[0];
  double Fa = x[1];
  double ret=BBBA07::XsIntegral(Q2, Fa, false, false);
  return 1/ret;
}

TF2* FAJacobianBBA05 = new TF2("FAJacobianBBA05",GetFaJacobianBBA05Wrapper,0,20,-2,2,0);
TF2* FAJacobianBBBA07 = new TF2("FAJacobianBBBA07",GetFaJacobianBBBA07Wrapper,0,20,-2,2,0);


double CalcFA( double Q2, double fa )
{
  return fa/BBA05::gA;
}

MnvH1D* SolveFA( MnvH1D* data_xsec, bool UseBBA05 )
{
  cout<<"Enetering SolveFA"<<endl;
  TH1D hdata = data_xsec->GetCVHistoWithError();
  hdata.Scale(1,"width");
  TH1D hfa = hdata; hfa.SetName("hfa"); hfa.Reset();
  TH1D hfaE1 = hfa; hfaE1.SetName("hfa_e1");
  TH1D hfaE2 = hfa; hfaE2.SetName("hfa_e2");

  TF1* f = (UseBBA05)? FAFuncBBA05 : FAFuncBBBA07;

  MnvH1D* ret = (MnvH1D*) data_xsec->Clone("FA");
  ret->ClearAllErrorBands();
  ret->Reset();
  ret->AddVertErrorBand("Error",2);
  vector<TH1D*> errbandHists = ret->GetVertErrorBand("Error")->GetHists();

  //cout<<"Print MA"<<endl;
  for( int i = 0; i< hdata.GetNbinsX()+1; i ++ )
  {

    int ibin = i+1;
    double Q2 = hdata.GetBinCenter(ibin);
    Q2 = GetWeightedQ2Center( Q2 );

    //if( Q2<0.05  ) continue;
    f->SetParameter(0, Q2);

    double xs = hdata.GetBinContent(ibin)*1e40;
    double xsErr = hdata.GetBinError(ibin)*1e40;
    //cout<<"Check xs: "<<xs<<"+/-"<<xsErr<<endl;

    vector<double> XsValues({ xs-xsErr, xs, xs+xsErr } );
    vector<double> FaValues;
    for( int j = 0; j < 3; j++ )
    {
      double this_xs = XsValues[j];//scale this up
      double fa = f->GetX(this_xs, BBA05::gA*1.5,0,1e-10,1000  );
      //cout<<j<<", Ma: "<<ma<<", "<<f->Eval(ma)<<", "<<this_xs<<endl;
      FaValues.push_back(fa);
    }
    //cout<<"Q2: "<<Q2<<". Ma: "<<MaValues[1]<<endl;
    //cout<<"Xs: "<<XsValues[1]<<", "<<GetXsValueBBA05(Q2,MaValues[1])*1e40<<endl;
    //cout<<"Test: "<<f->Eval( MaValues[1] )<<endl;
    cout<<FaValues[1]<<endl;
    double fa_cv = FaValues[1];
    double fa_1 = FaValues[0];
    double fa_2 = FaValues[2];
    ret->SetBinContent(ibin, fa_cv );
    ret->GetVertErrorBand("Error")->SetBinContent(ibin, fa_cv );
    errbandHists[0]->SetBinContent(ibin, fa_1);
    errbandHists[1]->SetBinContent(ibin, fa_2);
    double avgerr = (abs(fa_1-fa_cv) + abs(fa_2-fa_cv))/2;
    ret->SetBinError(ibin, avgerr);

  }

  //ret->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");

  ret->GetYaxis()->SetTitle("#font[32]{F}_{A}");
  return ret;
}

//TH1D* GetFAHisto( MnvH1D* data_xsec, bool UseBBA05 )
//{
//  vector<TH1D> hists = SolveFa( data_xsec, UseBBA05 );
//  TH1D hfa_cv = hists[1];
//  TH1D hfa_E1 = hists[0];
//  TH1D hfa_E2 = hists[2];
//
//  TH1D* hfa = (TH1D*) hfa_cv.Clone("hfa");
//  hfa->Reset();
//
//  for(int i = 0; i< hfa->GetNbinsX()+1; i++ )
//  {
//    int ibin = i+1;
//    double Q2 = hfa_cv.GetBinCenter(ibin);
//    if( Q2<0.05 || Q2>2 ) continue;
//    double fa_cv = hfa_cv.GetBinContent(ibin);
//    double fa_e1 = hfa_E1.GetBinContent(ibin);
//    double fa_e2 = hfa_E1.GetBinContent(ibin);
//    hfa->SetBinContent(ibin, fa_cv);
//    double err = (abs(fa_e1-fa_cv)+abs(fa_e2-fa_cv))/2.;
//    hfa->SetBinError( ibin, err );
//    //cout<<Q2<<" "<<fa_cv<<", "<<err<<endl;
//  }
//  hfa->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
//  hfa->GetYaxis()->SetTitle("F_{A}/g_{A}");
//  return hfa;
//}


bool useKevin = true;
int form = 0;

double GetScale( TH1D* h )
{
  //int b0 = h->FindBin(0.25);
  //int b1 = h->FindBin(1.4);
  //int b07 = h->FindBin(0.7);
  ////b0 = b1 = b07;
  //b0 = b1;
  //double dataInt = xs_data->Integral(b0,b1);
  double hInt = h->Integral(1,1);
  double scale = 22*pow(10,-39)/hInt;
  return scale;
}

double chi2Func( double *par )
{
  double MA = par[0];
  TH1D* h = GetXsHisto( MA, useKevin, form );
  //double scale = GetScale( h );
  //h->Scale(scale);
  double chi2 = 0;
  int b0 = h->FindBin(0.25);
  int b1 = h->FindBin(1.4);
  for( int i = b0; i<=b1; i++ )
  {
    double d = xs_data->GetBinContent(i);//*pow(10,38);
    double derr = xs_data->GetBinError(i);//*pow(10,38);
    double f =  h->GetBinContent(i);//*pow(10,38);
    chi2+= pow(d-f,2)/derr;
  }
  double residual = abs(chi2/(b1-b0+1)-1);
  //cout<<dataInt<<", "<<hInt<<", "<<MA<<", "<<chi2<<", "<<residual<<endl;
  delete h;
  return chi2;
  //return residual;
}


double result = 0;
void minuitFunction( int &nDim, double* gout, double &result, double par[], int flg)
{
  result= chi2Func( par );
}

TFitter* FitMA( int f )
{
  useKevin = false;
  form = f;
  TFitter* minimizer = new TFitter(1);
  double print_val = -1;
  minimizer->ExecuteCommand("SET PRINTOUT",&print_val,1);
  minimizer->SetFCN(minuitFunction);
  minimizer->SetParameter(0,"ma",1.4,0.00001,0.5,2);
  //minimizer->ExecuteCommand("SIMPLEX",0,0);
  minimizer->ExecuteCommand("MIGRAD",0,0);
  double fitMA = minimizer->GetParameter(0);
  double fitMAErr = minimizer->GetParError(0);
  double fitMAResidual = minimizer->GetParError(0);
  cout<<"FitMA: "<<fitMA<<" "<<fitMAErr<<endl;
  //return vector<double>{fitMA, fitMAErr };
  return minimizer;

}




void drawXs( TCanvas* c, vector<MnvH1D*> &hists,int iteration, string tag="" )
{
  gStyle->SetEndErrorSize(10);
  gStyle->SetErrorX(0);
  //TFile f_genie("/minerva/data/users/tejinc/Documents/Neutrons/figures/cross_sections/GENIEXSECEXTRACT_CCQENuNeutron_Q2QE_minervame6D.root");
  TFile f_genie("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/GENIEXSecExtract/cmt/GENIEXSECEXTRACT_CCQENuNeutron_Q2QE_minervame6D.root");
  MnvH1D *h_ds_dq2_genie = (MnvH1D*) f_genie.Get("ds_dq2_xsec");

  MnvPlotter plotter;
  gStyle->SetHistLineWidth(4);
  gStyle->SetErrorX(0);
  c->cd();
  c->SetLogx(false);
  c->SetLogy(false);
  c->SetLogz(false);

  TLegend* leg = new TLegend(.6,.6,.8,.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  MnvH1D* data = (MnvH1D*) hists[6]->Clone("data");
  MnvH1D* mc = (MnvH1D*) hists[7]->Clone("mc");
  MnvH1D* eff_denum_xs = (MnvH1D*) hists[8]->Clone("eff_denum_xs");

  data->SetLineWidth(2);
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(1);
  data->SetMarkerColor( kBlack );
  data->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");

  mc->SetLineWidth(2);

  data->Scale(1,"width");
  mc->Scale(1,"width");
  //vector<TH1D*> MaErrBands = mc->GetVertErrorBand("GENIE_MaCCQE")->GetHists();
  eff_denum_xs->Scale(1,"width");


  xs_data = data;//not sure why it's here..
  //mc->Scale( GetScale( mc ) );


  TH1D h_sys = data->GetCVHistoWithError();
  TH1D h_stat = data->GetCVHistoWithStatError();
  TH1D h_mc = mc->GetCVHistoWithError();

  data->GetXaxis()->SetRangeUser(xs_xmin, xs_xmax);
  mc->GetXaxis()->SetRangeUser(xs_xmin, xs_xmax);


  TH1D* h_mc_err = (TH1D*) mc->GetCVHistoWithError().Clone("h_mc_er");
  h_mc_err->SetFillStyle(1001);
  h_mc_err->SetFillColorAlpha(kRed,0.35);
  h_mc_err->SetMarkerSize(0);

  //h_sys.GetXaxis()->SetRangeUser(0.05,1.4);
  h_sys.GetYaxis()->SetRangeUser(0,0.3*pow(10,-37));
  h_sys.GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
  h_sys.Draw("e1");
  h_stat.Draw("e1same");
  h_mc.Draw("histlsame");

  //h_mc_err->GetXaxis()->SetRangeUser(0.05,2);
  h_mc_err->Draw("e5same");
  //for( auto h: MaErrBands )
  //{
  //  h->Draw("histlsame");
  //}
  //h_sys.Draw("e1same");
  //h_stat.Draw("pe1same");

  h_ds_dq2_genie->SetLineColor(kOrange);
  h_ds_dq2_genie->SetLineWidth(2);
  //h_ds_dq2_genie->Draw("histlsame");

  double ymax = 25*pow(10,-39);

  leg->AddEntry( data, "Data");
  leg->AddEntry( mc, "MC");
  //leg->AddEntry( MaErrBands[0], "M_{A} - 1#sigma");
  //leg->AddEntry( MaErrBands[1], "M_{A} + 1#sigma");
  //leg->AddEntry( h_ds_dq2_genie, "GENIE XS (6D)");



  c->SetLogx(true);
  plotter.arrow_line_color = kRed;
  plotter.arrow_line_width = 2;
  plotter.AddCutArrow(0.0125, 0, ymax, .02, "R");
  plotter.AddCutArrow(6, 0, ymax, 1, "L");

  leg->Draw();
  
  c->Print(Form("plots/xs-final-q2qe-it_%d-%s.pdf", iteration, tag.c_str()));



  double testMA = 0.99;
  useKevin = false;
  TH1D* h_testMA_BBA05 = GetXsHisto( testMA, useKevin, 0 );
  TH1D* h_testMA_BBBA07 = GetXsHisto( testMA, useKevin, 1 );



  int binBegin = mc->FindBin(0.05+0.00001);

  h_testMA_BBA05->SetLineColor(kRed);
  h_testMA_BBA05->Draw("histlsame");

  h_testMA_BBBA07->SetLineColor(kMagenta);
  h_testMA_BBBA07->Draw("histlsame");


  leg->AddEntry(h_testMA_BBA05, Form("Ma=%.2f, BBBA05", testMA));
  leg->AddEntry(h_testMA_BBBA07, Form("Ma=%.2f, BBBA07", testMA));
  leg->Draw();
  c->Print(Form("plots/xs-final-q2qe-it_%d-fit-%s.pdf",iteration, tag.c_str()));


  TH1D* h_Meyer_BBA05 = GetZexpXsecHisto( 0 );
  TH1D* h_Hydrogen_BBA05 = GetZexpXsecHisto( 1 );
  h_Meyer_BBA05->SetLineColor( kOrange );
  h_Hydrogen_BBA05->SetLineColor( kRed );
  h_Hydrogen_BBA05->SetLineStyle( 2 );

  h_Meyer_BBA05->Draw("histlsame");
  h_Hydrogen_BBA05->Draw("histlsame");

  leg->AddEntry(h_Meyer_BBA05, Form("Meyer Fit" ));
  leg->AddEntry(h_Hydrogen_BBA05, Form("Hydrogen Fit"));

  c->Print(Form("plots/xs-final-q2qe-it_%d-fit-%s-zexp.pdf",iteration, tag.c_str()));

  /*

  TH1D* h_Meyer_BBA05_Ratio = (TH1D*) h_Meyer_BBA05->Clone("h_Meyer_BBA05_Ratio");
  TH1D* h_Hydrogen_BBA05_Ratio = (TH1D*) h_Hydrogen_BBA05->Clone("h_Hydrogen_BBA05_Ratio");
  h_Meyer_BBA05_Ratio->Divide(h_testMA_BBA05 );
  h_Hydrogen_BBA05_Ratio->Divide(h_testMA_BBA05); 

  TH1D* h_sys_ratio = (TH1D*) h_sys.Clone("h_sys_ratio");
  TH1D* h_stat_ratio= (TH1D*) h_stat.Clone("h_stat_ratio");
  h_sys_ratio->Divide( h_testMA_BBA05 );
  h_stat_ratio->Divide( h_testMA_BBA05 );

  h_Meyer_BBA05_Ratio->GetYaxis()->SetRangeUser(0,2);
  h_Meyer_BBA05_Ratio->Draw(

  c->Print(Form("plots/xs-final-q2qe-it_%d-fit-%s-zexp-ratio.pdf",iteration, tag.c_str()));
  */


  TFitter* fitBBA05 = FitMA( 0 );
  TFitter* fitBBBA07 = FitMA(1 );



  TH1D* h_ma_fit_BBA05 = GetXsHisto( fitBBA05->GetParameter(0), false,0 );
  h_ma_fit_BBA05->SetLineColor(kRed);
  h_ma_fit_BBA05->SetLineStyle( 2 );
  h_ma_fit_BBA05->Draw("histlsame");

  TH1D* h_ma_fit_BBBA07 = GetXsHisto( fitBBBA07->GetParameter(0), false,1 );
  h_ma_fit_BBBA07->SetLineColor(kMagenta);
  h_ma_fit_BBBA07->SetLineStyle( 2 );
  h_ma_fit_BBBA07->Draw("histlsame");



  leg->AddEntry( h_ma_fit_BBA05, Form("Ma=%.3f#pm%.3f,BBBA05", fitBBA05->GetParameter(0),fitBBA05->GetParError(0) ) );
  leg->AddEntry( h_ma_fit_BBBA07, Form("Ma=%.3f#pm%.3f,BBBA07", fitBBBA07->GetParameter(0),fitBBBA07->GetParError(0) ) );

  leg->Draw();
  c->Print(Form("plots/xs-final-q2qe-it_%d-fit-%s-new.pdf",iteration,tag.c_str()));
  //Start Fitting: 
  //
  //Draw ratio --

  TLegend* leg2 = new TLegend(.8,.8,.95,.95);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  mc->AddMissingErrorBandsAndFillWithCV(*data);

  data->Divide(data, mc);
  data->GetYaxis()->SetTitle("Data/MC");
  //data->GetXaxis()->SetRangeUser(0.05,6);

  MnvH1D *mc_ratio = (MnvH1D*) mc->Clone("mc_ratio");
  mc_ratio->DivideSingle(mc_ratio, &h_mc);

  cout<<"mc_ratio->DivideSingle"<<endl;

  MnvH1D *mc_ratio_err = (MnvH1D*) mc_ratio->Clone("mc_ratio_err");
  cout<<" MnvH1D *mc_ratio_err = (MnvH1D*) mc_ratio->Clone(mc_ratio_err);"<<endl;

  mc_ratio_err->SetFillStyle(1001);
  mc_ratio_err->SetFillColorAlpha(kRed,0.35);
  mc_ratio_err->SetMarkerSize(0);
  //eff_denum_xs->Divide(eff_denum_xs, h_ds_dq2_genie);
  h_testMA_BBBA07->Divide(h_testMA_BBBA07, mc);
  h_testMA_BBA05->Divide(h_testMA_BBBA07, mc);

  data->GetYaxis()->SetRangeUser(0,2);
  //data->GetXaxis()->SetRangeUser(0.05,4);

  data->Draw("e1");

  mc_ratio->Draw("histlsame");
  mc_ratio_err->Draw("e5same");

  data->Draw("e1same");
  //leg2->AddEntry(data,"data xs");
  //leg2->AddEntry(mc_ratio,"mc xs");
  //leg2->AddEntry(eff_denum_xs,"eff denum xs");
  //leg2->AddEntry( h_testMA_BBA05, "BBA05, Ma=.99");
  //leg2->AddEntry( h_testMA_BBBA07, "BBBA07, Ma=.99");
  //leg2->Draw();
  c->SetLogx();
  c->Print(Form("plots/xs-ratio-final-q2qe-it_%d-%s.pdf",iteration,tag.c_str()));
  cout<<"Printed"<<endl;
}

void drawXsZexp( TCanvas* c, vector<MnvH1D*> &hists,int iteration,  bool ratio = false, string tag="")
{
  gStyle->SetEndErrorSize(10);
  gStyle->SetErrorX(0);
  //TFile f_genie("/minerva/data/users/tejinc/Documents/Neutrons/figures/cross_sections/GENIEXSECEXTRACT_CCQENuNeutron_Q2QE_minervame6D.root");
  TFile f_genie("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/GENIEXSecExtract/cmt/GENIEXSECEXTRACT_CCQENuNeutron_Q2QE_minervame6D.root");
  MnvH1D *h_ds_dq2_genie = (MnvH1D*) f_genie.Get("ds_dq2_xsec");

  MnvPlotter plotter;
  gStyle->SetHistLineWidth(4);
  gStyle->SetErrorX(0);
  c->cd();
  c->SetLogx(false);
  c->SetLogy(false);
  c->SetLogz(false);


  MnvH1D* data = (MnvH1D*) hists[6]->Clone("data");
  MnvH1D* mc = (MnvH1D*) hists[7]->Clone("mc");
  MnvH1D* eff_denum_xs = (MnvH1D*) hists[8]->Clone("eff_denum_xs");

  data->SetLineWidth(2);
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(1);
  data->SetMarkerColor( kBlack );
  data->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");

  mc->SetLineWidth(2);

  data->Scale(1,"width");
  mc->Scale(1,"width");
  //vector<TH1D*> MaErrBands = mc->GetVertErrorBand("GENIE_MaCCQE")->GetHists();
  eff_denum_xs->Scale(1,"width");


  xs_data = data;//not sure why it's here..
  //mc->Scale( GetScale( mc ) );


  TH1D h_sys = data->GetCVHistoWithError();
  TH1D h_stat = data->GetCVHistoWithStatError();
  TH1D h_mc = mc->GetCVHistoWithError();

  data->GetXaxis()->SetRangeUser(xs_xmin, xs_xmax);
  mc->GetXaxis()->SetRangeUser(xs_xmin, xs_xmax);


  TH1D* h_mc_err = (TH1D*) mc->GetCVHistoWithError().Clone("h_mc_er");
  h_mc_err->SetFillStyle(1001);
  h_mc_err->SetFillColorAlpha(kRed,0.35);
  h_mc_err->SetMarkerSize(0);

  //h_sys.GetXaxis()->SetRangeUser(0.05,1.4);

  h_sys.GetYaxis()->SetRangeUser(-0.001*pow(10,-37),0.3*pow(10,-37));
  h_sys.GetXaxis()->SetTitle("Q^{2} (GeV^{2})");

  
  TGraphErrors h_Meyer_BBA05 = GetZexpXsecGraph( 0, ratio );
  TGraphErrors h_Hydrogen_BBA05 = GetZexpXsecGraph( 1, ratio );
  h_Meyer_BBA05.SetLineColor( kOrange+1 );
  h_Meyer_BBA05.SetLineWidth(2);
  h_Meyer_BBA05.SetFillColorAlpha( kOrange+1,0.3);
  h_Meyer_BBA05.SetFillStyle( 1001 );
  h_Meyer_BBA05.SetMarkerStyle( 0 );
  h_Meyer_BBA05.SetMarkerSize( 0 );

  h_Hydrogen_BBA05.SetLineColor( kMagenta-1 );
  h_Hydrogen_BBA05.SetLineWidth(2);
  h_Hydrogen_BBA05.SetFillColorAlpha( kMagenta-1,0.2);
  h_Hydrogen_BBA05.SetFillStyle( 1001 );
  h_Hydrogen_BBA05.SetMarkerStyle( 0 );
  h_Hydrogen_BBA05.SetMarkerSize( 0 );


  double testMA = 1.015;
  useKevin = false;
  TH1D* h_testMA_BBA05 = GetXsHisto( testMA, useKevin, 0 );
  h_testMA_BBA05->SetLineWidth(2);
  h_testMA_BBA05->SetLineColor( kBlack );

  if( ratio )
  {
    h_sys.Divide( h_testMA_BBA05 );
    h_stat.Divide( h_testMA_BBA05 );
    h_testMA_BBA05->Divide( h_testMA_BBA05 );
    h_sys.GetYaxis()->SetRangeUser(0,8);
    h_sys.GetYaxis()->SetTitle( "Ratio to Dipole F_{A}" );
  }

  h_sys.Draw("e1");
  h_stat.Draw("e1same");
  h_testMA_BBA05->Draw("histlsame");

  h_Meyer_BBA05.Draw("l3");
  h_Hydrogen_BBA05.Draw("l3");

  h_sys.Draw("e1same");
  h_stat.Draw("e1same");

  //h_Meyer_BBA05_err->Draw("e5same");
  //h_Hydrogen_BBA05_err->Draw("e5same");
  TLegend* leg; 
  if(ratio) leg= new TLegend(.4,.7,.7,.9);
  else leg= new TLegend(.5,.7,.8,.9);
  //leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->AddEntry( &h_sys, "Data" ,"pl");
  leg->AddEntry( &h_Meyer_BBA05, "Meyer Fit" ,"lf");
  leg->AddEntry( &h_Hydrogen_BBA05, "Hydrogen Fit" ,"lf");
  leg->AddEntry( h_testMA_BBA05, "Dipole M_{A}=1.015 GeV/c^{2}" ,"l");
  leg->Draw();
  c->SetLogx();
  c->Print(Form("plots/xs-final-zexp-it_%d-ratio_%d-fit-%s.pdf",iteration,int(ratio),tag.c_str()));
}



void drawFA( TCanvas* c, vector<MnvH1D*> &hists,int iteration, string tag="" )
{
  gStyle->SetEndErrorSize(10);
  gStyle->SetErrorX(0);
  TFile f_genie("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/GENIEXSecExtract/cmt/GENIEXSECEXTRACT_CCQENuNeutron_Q2QE_minervame6D.root");
  MnvH1D *h_ds_dq2_genie = (MnvH1D*) f_genie.Get("ds_dq2_xsec");

  MnvPlotter plotter;
  gStyle->SetHistLineWidth(4);
  gStyle->SetErrorX(0);
  c->cd();
  c->SetLogx(false);
  c->SetLogy(false);
  c->SetLogz(false);

  TLegend* leg = new TLegend(.6,.6,.85,.90);
  leg->SetBorderSize(0);

  MnvH1D* data = (MnvH1D*) hists[6]->Clone("data");
  //data->Scale(1,"width");
  //vector<double> new_bins({0.025,0.0375,0.05,0.2,0.3,0.4,0.6,0.8,1.0,1.2,2.0,4.0,6.0});
  vector<double> new_bins({0.0375,0.05,0.2,0.3,0.4,0.6,0.8,1.0,1.2,2.0,4.0,6.0});
  MnvH1D* data_rebinned = (MnvH1D*) data->Rebin( new_bins.size()-1, "h_xs_rebinned", &new_bins[0] );

  MnvH1D* Fa_BBA05 = SolveFA( data, true );
  MnvH1D* Fa_BBBA07 = SolveFA( data,false);

  Fa_BBA05->Scale(-1);
  Fa_BBBA07->Scale(-1);
  Fa_BBA05->SetMarkerColor( kBlue ); Fa_BBA05->SetLineColor( kBlue ); Fa_BBA05->SetMarkerStyle(kFullCircle);
  Fa_BBBA07->SetMarkerColor( kRed ); Fa_BBBA07->SetLineColor( kRed );Fa_BBBA07->SetMarkerStyle(kFullCircle);


  Fa_BBA05->SetLineWidth(2);
  Fa_BBBA07->SetLineWidth(2);
  Fa_BBA05->SetMarkerSize(0.75);
  Fa_BBBA07->SetMarkerSize(0.75);

  //cout<<"Print Fa_BBA05: "<<endl;
  //for( int i = 0; i< 20;i++) cout<<Fa_BBA05->GetBinContent(i)<<endl;

  //Fa_BBA05->GetXaxis()->SetRangeUser(0.025, 6);
  //Fa_BBA05->Draw("pe1");
  //Fa_BBBA07->Draw("pe1same");
  //vector<TH1D*> v05=Fa_BBA05->GetVertErrorBand("Error")->GetHists();
  //vector<TH1D*> v07=Fa_BBBA07->GetVertErrorBand("Error")->GetHists();
  //for(int i = 0; i<2; i++)
  //{
  //  v05[i]->Draw("histlsame");
  //  //v07[i]->Draw("histlsame");
  //}
  TGraphErrors Fa_BBA05G, Fa_BBBA07G;
  for( int i = 1; i< Fa_BBA05->GetNbinsX()+1; i++ )
  {
    double Q2 = Fa_BBA05->GetBinCenter(i);
    Fa_BBA05G.SetPoint(i,GetWeightedQ2Center(Q2), Fa_BBA05->GetBinContent(i) );
    Fa_BBA05G.SetPointError(i, 0, Fa_BBA05->GetBinError(i) );

    Fa_BBBA07G.SetPoint(i, GetWeightedQ2Center(Q2), Fa_BBBA07->GetBinContent(i) );
    Fa_BBBA07G.SetPointError(i, 0, Fa_BBBA07->GetBinError(i) );
  }

  Fa_BBA05G.SetMarkerColor(kRed);Fa_BBA05G.SetLineColor(kRed);
  Fa_BBA05G.SetMarkerStyle(25);Fa_BBA05G.SetLineWidth(1.5);
  Fa_BBBA07G.SetMarkerColor(kBlue);Fa_BBBA07G.SetLineColor(kBlue);
  Fa_BBBA07G.SetMarkerStyle(32);Fa_BBBA07G.SetLineWidth(1.5);
  
  


  TGraph bnl,anl,fnal;
  TGraph dipole;
  TGraphErrors meyer, tejin, joint;
  bnl.SetName("BNL");
  anl.SetName("ANL");
  fnal.SetName("FINAL");

  Zexp::SetPar("BNL");
  for(int i = 1; i<200;i++) { double Q2 = i/50.; bnl.SetPoint(i,Q2, Zexp::FA(Q2)*-1); }
  Zexp::SetPar("ANL");
  for(int i = 1; i<200;i++) { double Q2 = i/50.; anl.SetPoint(i,Q2, Zexp::FA(Q2)*-1); }
  Zexp::SetPar("FINAL");
  for(int i = 1; i<200;i++) { double Q2 = i/50.; fnal.SetPoint(i,Q2, Zexp::FA(Q2)*-1); }
  for(int i = 0; i<200;i++) { 
    double Q2 = i/50.; 
    meyer.SetPoint(i, Q2, Zexp::MeyerFitCV(Q2) );
    meyer.SetPointError(i, 0, Zexp::MeyerFitError(Q2) );


    double tejin_par[4]={2.60463, -1.32693, -12.7091, 35.689};
    tejin.SetPoint(i, Q2, Zexp::Zexp(Q2, tejin_par) );
    tejin.SetPointError(i, 0, Zexp::TejinFitError(Q2) );


    double joint_par[4]={2.27384, 0.643142, -6.4269, 5.12256};
    joint.SetPoint(i, Q2, Zexp::Zexp(Q2, joint_par) );
    joint.SetPointError(i, 0, Zexp::JointFitError(Q2) );
  }

  cout<<"Printing Dipole"<<endl;
  for(int i = 0; i<200;i++) { double Q2 = i/50.; dipole.SetPoint(i,Q2, BBA05::FA(Q2,0.99)*-1); cout<<Q2<<", "<<BBA05::FA(Q2,0.99)*-1<<endl; }
  

  //bnl.SetLineColor(kGray); bnl.SetFillStyle(0);
  //anl.SetLineColor(kGray+1); anl.SetFillStyle(0);
  meyer.SetLineColor(kOrange); meyer.SetFillStyle(1001); meyer.SetLineWidth(2);
  meyer.SetFillColorAlpha(kOrange,0.3);
  dipole.SetLineColor(kGray); dipole.SetFillStyle(0); dipole.SetLineWidth(2);



  //bnl.Draw("c");
  //anl.Draw("c");
  meyer.GetXaxis()->SetRangeUser(xs_xmin,xs_xmax);
  meyer.GetYaxis()->SetRangeUser(0,1.7);
  meyer.GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
  meyer.GetYaxis()->SetTitle("-F_{A}");
  meyer.Draw("al3");
  //meyer.Draw("c");
  dipole.Draw("l");
  //meyer.GetYaxis()->SetRangeUser(0.01,2);
  Fa_BBA05G.Draw("p");
  Fa_BBBA07G.Draw("p");
  c->Update();


  //Fa_BBA05->Draw("pe1same");
  //Fa_BBBA07->Draw("pe1same");


  leg->AddEntry( &Fa_BBA05G, "F_{A} from BBBA05", "pl" );
  leg->AddEntry( &Fa_BBBA07G, "F_{A} from BBBA07", "pl" );
  //leg->AddEntry( &bnl, "Z-exp BNL", "l");
  //leg->AddEntry( &anl, "Z-exp ANL", "l");
  leg->AddEntry( &meyer, "Z-exp Meyer Fit", "l");
  leg->AddEntry( &dipole, "Dipole Ma=.99", "l");
  leg->Draw();
  
  c->SetLogx();
  c->SetGrid();
  c->Print(Form("plots/xs-FA-it_%d-%s.pdf", iteration, tag.c_str()));

  bnl.Draw("Ac");
  anl.Draw("c");
  fnal.Draw("c");
  c->Print(Form("plots/xs-Zexp-it_%d-%s.pdf", iteration, tag.c_str()));



  TFile f(Form("plots/FaFitHisto_it%d.root",iteration),"recreate");
 
  f.cd();


  Fa_BBA05G.Write("Fa_BBA05G");
  Fa_BBBA07G.Write("Fa_BBBA07G");

  Fa_BBA05->SetName("Fa_BBA05");
  Fa_BBBA07->SetName("Fa_BBBA07");

  Fa_BBA05->Write();
  Fa_BBBA07->Write();
  data_rebinned->Write();
  data->Write("h_xs_hydrogen");

  //Function to find FA for a given xsec
  FAFuncBBA05->Write("FAFuncBBA05");
  FAFuncBBBA07->Write("FAFuncBBBA07");

  FAJacobianBBA05->Write("FAJacobianFuncBBA05");
  FAJacobianBBBA07->Write("FAJacobianFuncBBBA07");
  f.Close();
}

void drawXs2D( TCanvas* c, vector<MnvH2D*> &hists, bool doRatio, int iteration )
{

  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  gStyle->SetEndErrorSize(10);
  gStyle->SetErrorX(0);
  //TFile f_genie("/minerva/data/users/tejinc/Documents/Neutrons/figures/cross_sections/GENIEXSECEXTRACT_CCQENuNeutron_Q2.root");
  //MnvH1D *h_ds_dq2_genie = (MnvH1D*) f_genie.Get("ds_dq2_xsec");

  MnvPlotter plotter;
  gStyle->SetHistLineWidth(4);
  gStyle->SetErrorX(0);
  c->cd();
  c->SetLogx(false);
  c->SetLogy(false);
  c->SetLogz(false);

  TLegend* leg = new TLegend(.8,.8,.95,.95);

  MnvH2D* data = (MnvH2D*) hists[6]->Clone("data");
  MnvH2D* mc = (MnvH2D*) hists[7]->Clone("mc");

  data->SetLineWidth(2);
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(1);
  data->SetMarkerColor( kBlack );
  
  
  data->GetXaxis()->SetTitle("Q^{2} (GeV^{2}) vs E_{#nu} (GeV) ");

  mc->SetLineWidth(2);



  if(false){
  MnvH1D* enu_data = data->ProjectionX("enu_data_clone");
  MnvH1D* enu_mc =   mc->ProjectionX("enu_mc_clone");
  enu_data->Scale(1e42,"width");
  enu_mc->Scale(1e42,"width");
  enu_data->GetYaxis()->SetTitle("Arb. Unit");
  enu_data->SetMarkerStyle(32);
  enu_data->SetLineWidth(1);
  enu_data->SetLineColor(kBlack);
  enu_data->Draw("ep");
  enu_mc->Draw("histlsame");
  string fname1=Form("plots/effcor-2d-enu-it_%d.pdf", iteration );
  c->Print(fname1.c_str());
  }

  //xs_data = data;//not sure why it's here..
  //mc->Scale( GetScale( mc ) );

  data->Scale(1,"width");
  mc->Scale(1,"width");
  if(doRatio)
  {
    data->Divide(data,mc);
    mc->DivideSingle(mc,mc);
  }

  TH2* h_sys =    (TH2*)     data->GetCVHistoWithError().Clone("dataErr");
  TH2* h_stat =   (TH2*)    data->GetCVHistoWithStatError().Clone("dataStat");
  TH2* h_mc =     (TH2*)      mc->GetCVHistoWithError().Clone("mc");
  TH2* h_mc_err = (TH2*)  mc->GetCVHistoWithError().Clone("mcErr");
  h_mc_err->SetFillStyle(1001);
  h_mc_err->SetFillColorAlpha(kRed,0.35);
  h_mc_err->SetMarkerSize(0);

  vector<pair<TH2*, const char*> > 
    histAndOpts({
      {h_mc,"histl"},
      {h_mc_err,"e2"},
      {h_sys,"e1"},
      {h_stat,"e1"}});

  double *multipliers = NULL;

  string xtitle="E_{#nu} (GeV)";
  string ytitle="Q^2 (GeV^{2})";

  GridCanvas* gc = plotXAxis1D( histAndOpts, xtitle,ytitle, multipliers);
  gc->SetGridx(true);
  gc->SetLogx();
  //gc->SetLogx(logx);
  if ( doRatio ) gc->SetYLimits(0,2);
  gc->Modified();
  leg->AddEntry(data,"Efficiency Corrected Data");
  leg->AddEntry(mc,"Efficiency Corrected MC");
  string fname=Form("plots/xs-2d-enu-q2-axis_x-ratio_%d-it_%d.pdf", doRatio, iteration );
  gc->Print(fname.c_str());

  GridCanvas* gc2 = plotYAxis1D( histAndOpts, ytitle,xtitle, multipliers);
  gc2->SetGridx(true);
  //gc2->SetLogx(logx);
  if ( doRatio ) gc2->SetYLimits(0,2);
  gc2->SetLogx();
  gc2->Modified();
  leg->AddEntry(data,"Efficiency Corrected Data");
  leg->AddEntry(mc,"Efficiency Corrected MC");
  fname=Form("plots/xs-2d-enu-q2-axis_y-ratio_%d-it_%d.pdf", doRatio, iteration );
  gc2->Print(fname.c_str());
}

void drawSystematics( CCQENuPlotUtils*utils, MnvH1D* h,  double pot_data, double pot_mc, string name, bool asfrac, int iteration)
{
  h->GetXaxis()->SetRangeUser(xs_xmin, xs_xmax);
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

  plotter->axis_maximum = 0.9;
  //if( !asfrac )plotter->axis_maximum = 1000;
  plotter->headroom = 1.0;
  plotter->legend_n_columns = 2;
  plotter->legend_text_size = 0.025;
  plotter->AddPlotLabel( Form("%s #bullet%s", "#Q^{2}" ,  "Errors"), 0.33, 0.94, 0.025 );

  vector<string> sysGroupNames({"Muon Reconstruction", "Recoil Reconstruction", "Low Recoil Fit", "XSection Models", "FSI Models", "Others", "Flux"});
  vector<string> sysGroupPrintNames({"Muon_Reconstruction", "Recoil_Reconstruction", "Low_Recoil_Fit", "XSection_Models", "FSI_Models", "Others", "Flux"});

  h->GetXaxis()->SetRangeUser(0.0001,10);
  h->GetXaxis()->SetTitle("Q^{2} (GeV)");
  for( int i = 0; i < sysGroupPrintNames.size(); i++ )
  {
    //plotter->axis_maximum = 1000;
    cout<<sysGroupPrintNames[i]<<endl;
    utils->writeNorm( false, pot_data, pot_mc, true ); 
    plotter->DrawErrorSummary(h, "TL", false, true, 0.00001, false, sysGroupNames[i],asfrac);
    c->Print(Form("plots/sys_%s_%s_it-%d.pdf",name.c_str(), sysGroupPrintNames[i].c_str(), iteration ) );
  }

    //plotter->axis_maximum = 1000;
    plotter->AddPlotLabel( Form("%s #bullet%s", "#Q^{2}" ,  "Errors"), 0.33, 0.94, 0.025 );
    utils->writeNorm( false, pot_data, pot_mc, true ); 
    plotter->DrawErrorSummary(h, "TL", true, true, 0.00001, false, "", asfrac);
    c->Print( Form("plots/sys_%s_all_it-%d.pdf",name.c_str(), iteration) );
    delete plotter;
    delete c;
}



void DrawXs( TFile *f, TFile *f_fk,TFile *f2,  int iteration = 4 )
{
  cout<<"Entering Draw Xs"<<endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  CCQENuPlotUtils* putils = new CCQENuPlotUtils( fluxHistoExists );
  
  //setRedPalette();
  //setBlackbodyPalette();
  //setCorrelationPalette(0.5);
  //setRainbowToWhitePalette();
  //
  //
  TVector2 *pot = (TVector2*)f->Get("pot");

  double pot_data = pot->X();
  double pot_mc = pot->Y();
  double pot_norm = pot_data/pot_mc;


  
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetEndErrorSize(2);
  gStyle->SetErrorX(2);
  gStyle->SetTitleOffset(1.2,"Y");

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.15);

  TCanvas* c = new TCanvas();
  vector<MnvH1D*> AllHistos = GetAllHistos( f ) ;
  vector<MnvH1D*> AllHistos_FullKine = GetAllHistos( f_fk ) ;
  vector<MnvH2D*> AllHistos2D = GetAllHistos2D(f2);
  FormatHistos( AllHistos );


  drawFA(c, AllHistos, iteration);

  drawUnfold(c, AllHistos, iteration );
  drawCorrelation(c, AllHistos, iteration);
  drawEnuQ2(c);


  drawXs(c,AllHistos, iteration);
  drawXsZexp(c,AllHistos, iteration);
  drawXsZexp(c,AllHistos, iteration, true);
  //drawXs(c,AllHistos_FullKine, iteration, "fullKine");

  //drawXs2D(c,AllHistos2D, true, iteration);
  //drawXs2D(c,AllHistos2D, false, iteration);
  drawSystematics( putils, AllHistos[6], pot_data, pot_mc,"xsec_data", true, iteration );
  drawSystematics( putils, AllHistos[7], pot_data, pot_mc,"xsec_mc", true, iteration );
  //drawSystematics( putils, AllHistos_FullKine[6], pot_data, pot_mc,"xsec_dataFK", true, iteration );
  //drawSystematics( putils, AllHistos_FullKine[7], pot_data, pot_mc,"xsec_mcFK", true, iteration );
}




void DrawFATests()
{
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->SetGrid();
  double x[400],y1[400],y2[400],y3[400];
  double xs1_antinu[400],xs1_nu[400];
  double xs2_antinu[400],xs2_nu[400];
  int n = 400;
  double xs0;

  double E = 6;
  double MA = 0.99;
  for( int i = 0; i<n;i++ )
  {
    double Q2 = i*.1;
    x[i]=Q2;
    double gmp = BBA05::GMp(Q2);
    double gep = BBA05::GEp(Q2);
    double gmn = BBA05::GMn(Q2);
    double gen = BBA05::GEn(Q2);
    y1[i]= pow(gen/gmn, 2);
    y2[i]= pow(gmn/gmp, 2);
    y3[i]= pow(gep/gmp, 2);


    xs1_antinu[i]=BBA05::XsLS(E, Q2, MA, false)*scale;
    xs1_nu[i]=BBA05::XsLS(E, Q2, MA, true)*scale;

    xs2_antinu[i]=BBBA07::XsLS(E, Q2, MA, false)*scale;
    xs2_nu[i]=BBBA07::XsLS(E, Q2, MA, true);
    if(i==0) 
    {
      xs0 = BBA05::XsLS(E, Q2, MA, false);
      cout<<"Testing at E=6, Q2=3.125e-03: "<<BBA05::XsLS(6, 3.125*pow(10,-3), MA, false)<<endl;
      cout<<"Testing at E=3, Q2=3.125e-03: "<<BBA05::XsLS(3, 3.125*pow(10,-3), MA, false)<<endl;
      cout<<"Testing at E=1, Q2=3.125e-03: "<<BBA05::XsLS(1, 3.125*pow(10,-3), MA, false)<<endl;
      cout<<"Testing at E=6, Q2=9.375e-03: "<<BBA05::XsLS(6, 9.375*pow(10,-3), MA, false)<<endl;
      cout<<"Testing at E=3, Q2=9.375e-03: "<<BBA05::XsLS(3, 9.375*pow(10,-3), MA, false)<<endl;
      cout<<"Testing at E=1, Q2=9.375e-03: "<<BBA05::XsLS(1, 9.375*pow(10,-3), MA, false)<<endl;
      cout<<"Testing at E=6, Q2=0: "<<BBA05::XsLS(6, 0., MA, false)<<endl;
      cout<<"Testing at E=3, Q2=0: "<<BBA05::XsLS(3, 0., MA, false)<<endl;
      cout<<"Testing at E=1, Q2=0: "<<BBA05::XsLS(1, 0., MA, false)<<endl;

      cout<<"BBBA05 nu: "<<xs1_nu[0]*pow(10,38)<<endl;
      cout<<"BBBA05 antinu: "<<xs1_antinu[0]*pow(10,38)<<endl;
      cout<<"BBBA07 nu: "<<xs2_nu[0]*pow(10,38)<<endl;
      cout<<"BBBA07 antinu: "<<xs2_nu[0]*pow(10,38)<<endl;
    }
  }
  TGraph* gr = new TGraph(n,x,y2);
  gr->SetTitle("(Gmn/Gmp)^{2}");
  gr->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
  gr->Draw("AC");
  c1->Print("testFA_gmngmp.pdf");

  TMultiGraph* mg = new TMultiGraph();
  gr = new TGraph(n,x,y1);
  gr->SetTitle("(Gen/Gmn)^{2}");
  gr->SetLineColor(kBlue);
  gr->GetYaxis()->SetRangeUser(0,0.1);
  gr->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
  TGraph *gr1 =new TGraph(n,x,y3);
  gr1->SetLineColor(kBlack);
  gr1->SetTitle("(Gep/Gmp)^{2}");
  gr1->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");

  mg->Add(gr,"C");
  mg->Add(gr1,"C");
  mg->SetTitle("(Gen/Gmn)^{2} and (Gep/Gmp)^{2}");

  mg->Draw("AC");
  mg->GetXaxis()->SetRangeUser(0,30);
  mg->GetYaxis()->SetRangeUser(0,0.1);
  mg->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
  //mg->GetXaxis()->SetTitleSize(0.05);
  c1->Update();

  TLegend* leg = new TLegend(.5,0.7,.85,0.85);
  leg->SetHeader("Legend");

  leg->AddEntry( gr, "(Gep/Gmp)^{2}", "l" );
  leg->AddEntry( gr1, "(Gen/Gmn)^{2}", "l" );
  leg->Draw();

  c1->Print("testFA_gengmn_gepgmp.pdf");

  int n2=100;
  auto gr_xs1_antinu =  new TGraph(n2,x,xs1_antinu);
  auto gr_xs1_nu =  new TGraph(n2,x,xs1_nu);
  auto gr_xs2_antinu =  new TGraph(n2,x,xs2_antinu);
  auto gr_xs2_nu =  new TGraph(n2,x,xs2_nu);

  gr_xs1_antinu->SetLineColor(kBlack);
  gr_xs1_nu->SetLineColor(kBlack);
  gr_xs1_nu->SetLineStyle(2);

  gr_xs2_antinu->SetLineColor(kBlue);
  gr_xs2_nu->SetLineColor(kBlue);
  gr_xs2_nu->SetLineStyle(2);


  gr_xs1_antinu->SetMarkerStyle(0);
  gr_xs1_nu->SetMarkerStyle(0);
  gr_xs2_antinu->SetMarkerStyle(0);
  gr_xs2_nu->SetMarkerStyle(0);

  TLegend *leg2 = new TLegend(.5,.5,.9,.9 );
  leg2->SetHeader("Legend");
  leg2->AddEntry(gr_xs1_antinu, "BBA05 antinu", "l");
  leg2->AddEntry(gr_xs1_nu, "BBA05 nu", "l");
  leg2->AddEntry(gr_xs2_antinu, "BBBA07 antinu", "l");
  leg2->AddEntry(gr_xs2_nu, "BBBA07 nu", "l");




  gr_xs1_antinu->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
  xs0*= pow(10,38);
  gr_xs1_antinu->SetTitle(Form( "d#sigma/dQ^{2}(0) = %.2f #times 10^{-38}", xs0));
  gr_xs1_antinu->Draw("AC");
  //gr_xs1_antinu->GetXaxis()->SetRangeUser(0,5);
  c1->Update();
  gr_xs2_antinu->Draw("C");
  gr_xs1_nu->Draw("C");
  gr_xs2_nu->Draw("C");
  leg2->Draw();
  //c1->SetLogx();
  c1->Print("testFA_xs.pdf");

  //Test FA

}

void GenerateAMatrix(  MnvH1D* input_data, string name, int iteration  )
{
  cout<<"Generating AMatrix"<<endl;
  //MnvH1D* data = (MnvH1D*) hists[6]->Clone("data");

  //double new_bins[10]={0.05,0.2,0.3,0.4,0.6,0.8,1.0,1.2,2.0,6.0};
  //MnvH1D* data_rebinned = data->Rebin(9, "data_rebinned", new_bins )


  MnvH1D* data = (MnvH1D*) input_data->Clone("this_data");
  cout<<"Solving FA"<<endl;
  MnvH1D* Fa_BBA05 = SolveFA( data, true );
  MnvH1D* Fa_BBBA07 = SolveFA( data,false);
  cout<<"Solved FA"<<endl;
  TGraphErrors Fa_BBA05G, Fa_BBBA07G;

  for( int i = 0; i< Fa_BBA05->GetNbinsX()+2; i++ )
  {
    double Q2 = Fa_BBA05->GetBinCenter(i);
    double Q2w = GetWeightedQ2Center( Q2 );
    
    Fa_BBA05G.SetPoint(i, Q2w , Fa_BBA05->GetBinContent(i) );
    Fa_BBA05G.SetPointError(i, 0, Fa_BBA05->GetBinError(i) );
    cout<<"Q2: "<<Q2<<", "<< Fa_BBA05->GetBinContent(i) <<endl;

    Fa_BBBA07G.SetPoint(i, Q2w , Fa_BBBA07->GetBinContent(i) );
    Fa_BBBA07G.SetPointError(i, 0, Fa_BBBA07->GetBinError(i) );
  }

  Fa_BBA05G.SetMarkerColor(kRed);Fa_BBA05G.SetLineColor(kRed);
  Fa_BBA05G.SetMarkerStyle(25);Fa_BBA05G.SetLineWidth(1.5);
  Fa_BBBA07G.SetMarkerColor(kBlue);Fa_BBBA07G.SetLineColor(kBlue);
  Fa_BBBA07G.SetMarkerStyle(32);Fa_BBBA07G.SetLineWidth(1.5);
  
 

  data->Scale(1,"width");
  //data_rebinned->Scale(1,"width");
  TMatrixD cov=data->GetTotalErrorMatrix( true,   false, false);
  TMatrixD matrix=data->GetTotalErrorMatrix( true,   false, false);


  for( int i = 0; i< data->GetNbinsX()+2; i++ )
  {
    for( int j = 0; j< data->GetNbinsX()+2; j++ ) matrix[i][j]=0;
  }
  TMatrixD jac_BBA05 = matrix;
  TMatrixD jac_BBBA07 = matrix;

  TMatrixD cov_BBA05 = matrix;
  TMatrixD cov_BBBA07 = matrix;

  TH1D h = data->GetCVHistoWithError();
  h.Reset();

  TFile fout(Form("plots/%s_it%d.root",name.c_str(), iteration),"recreate");
  cout<<"created fout"<<endl;
  fout.cd();
  bool UseKevin = false;
  for( int i = 1; i< h.GetNbinsX()+1; i++ )
  {
    cout<<i<<", "<<endl;
    double Q2 = GetWeightedQ2Center( h.GetBinCenter( i ) );
    bool isNeutrino = false;
    double dxs1=BBA05::dXsIntegral(Q2, Fa_BBA05->GetBinContent(i),isNeutrino, UseKevin);
    double dxs2=BBBA07::dXsIntegral(Q2,Fa_BBBA07->GetBinContent(i),isNeutrino, UseKevin);


    jac_BBA05[i][i]=1/dxs1;
    jac_BBBA07[i][i]=1/dxs2;
    //if( Q2>=0.05 )
    //{
    //  jac_BBA05[i][i]=1/dxs1;
    //  jac_BBBA07[i][i]=1/dxs2;
    //}

    TGraph g1,g2;
    g1.SetName(Form("xs_fa_05_Q2_%i",i));
    g2.SetName(Form("xs_fa_07_Q2_%i",i));
    g1.GetXaxis()->SetTitle("-F_{A}");
    g1.GetYaxis()->SetTitle("Xs");
    g2.GetXaxis()->SetTitle("-F_{A}");
    g2.GetYaxis()->SetTitle("Xs");

    for( int j = 0; j<140;j++) 
    {
      cout<<j<<", ";
      double fa = -(j-20)/99.;
      double xsj1 = BBA05::XsIntegral(Q2, fa, isNeutrino, UseKevin);
      double xsj2 = BBBA07::XsIntegral(Q2,fa, isNeutrino, UseKevin);
      g1.SetPoint(j, -fa, xsj1 );
      g2.SetPoint(j, -fa, xsj2 );
    }
    cout<<endl;

    g1.Write();
    g2.Write();
  }

  TH1D hdata = data->GetCVHistoWithError();
  TMatrixD errormatrix = data->GetTotalErrorMatrix( true,   false, false);
  TMatrixD errormatrix_shape = data->GetTotalErrorMatrix( true,   false, true);

  data->Write("data_mnvh1d");
  hdata.Write("data");
  errormatrix.Write("Data_Covariance");
  errormatrix_shape.Write("Data_Covariance_Shape");



  jac_BBA05.Write("Jacobian_BBA05");
  cov_BBA05 = jac_BBA05*cov*jac_BBA05;
  cov_BBA05.Write("Cov_BBA05");

  Fa_BBA05->ClearAllErrorBands();
  Fa_BBA05->PushCovMatrix("fa_cov",cov_BBA05);
  for(int i = 0; i< Fa_BBA05->GetNbinsX()+2;i++) Fa_BBA05->SetBinError(i, TMath::Sqrt( cov_BBA05[i][i] ) );
  Fa_BBA05->GetCVHistoWithStatError().Write("Fa_BBA05");



  jac_BBBA07.Write("Jacobian_BBBA07");
  cov_BBBA07 = jac_BBBA07*cov*jac_BBBA07;
  cov_BBBA07.Write("Cov_BBBA07");
  Fa_BBBA07->ClearAllErrorBands();
  Fa_BBBA07->PushCovMatrix("fa_cov",cov_BBBA07);
  for(int i = 0; i< Fa_BBBA07->GetNbinsX()+2;i++) Fa_BBBA07->SetBinError(i, TMath::Sqrt( cov_BBBA07[i][i] ) );
  Fa_BBBA07->GetCVHistoWithStatError().Write("Fa_BBBA07");


  TGraphErrors bodek_fit;
  int i=0;
  double Q2=0, dQ2 = 0.001,Q2Max=10;
  while(Q2<Q2Max)
  {
    bodek_fit.SetPoint(i, Q2, BBBA07::G("FA",Q2) );
    Q2+=dQ2; i++;
  }

  bodek_fit.Write("Bodek_Fit");

  Fa_BBA05G.Write("Fa_BBA05G");
  Fa_BBBA07G.Write("Fa_BBBA07G");
  fout.Close();
  //Test Plotting ds/dFA:
  cout<<endl;

}

int main(int argc, char* argv[])
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
 
  vector<string> par(7);
  string f_eff = argv[1];
  string f_eff_fk = argv[2];
  string f_mig = argv[3];
  string f_xs = argv[4];
  string f_xs_fk = argv[5];
  string f_xs_2D = argv[6];
  string s_it = argv[7];

  par[0] = f_eff;
  par[1] = f_eff_fk;
  par[2] = f_mig;
  par[3] = f_xs;
  par[4] = f_xs_fk;
  par[5] = f_xs_2D;
  par[6] = s_it;

  int iteration= stoi(s_it );

  for( int i = 0; i<par.size(); i++ ) cout<<"Parameter "<<i<<": "<<par[i]<<endl;
  //int iteration = atoi( argv[4] );
  //cout<<iteration<<endl;

  string f_flux_suffix = (neutrinoMode)? 
    "data/flux/flux-gen2thin-pdg14-minervame1D1M1NWeightedAve.root": 
    "data/flux/flux-gen2thin-pdg-14-minervame6A.root";
  string f_flux = Form("%s/%s", std::getenv("PLOTUTILSROOT"),f_flux_suffix.c_str() ) ; 
 
  TFile* file_eff = new TFile( f_eff.c_str(),"read");
  TFile* file_eff_fk = new TFile( f_eff_fk.c_str(),"read");
  TFile* file_mig = new TFile( f_mig.c_str(),"read");
  TFile* file_xs = new TFile( f_xs.c_str(),"read");
  TFile* file_xs_fk = new TFile( f_xs_fk.c_str(),"read");
  TFile* file_xs_2D = new TFile( f_xs_2D.c_str(),"read");
  TFile* file_flux = new TFile( f_flux.c_str(), "read");

 
  PlotUtils::FluxReweighter* frw = new PlotUtils::FluxReweighter(-14,false,FluxReweighter::minervame6A,FluxReweighter::gen2thin, FluxReweighter::g4numiv6,100);
  MnvH1D* h_flux = frw->GetFluxReweighted(-14);
  BBA05::SetPars(file_xs, file_xs_2D, h_flux);
  BBBA07::SetPars(file_xs, file_xs_2D, h_flux);


  vector<MnvH1D*> hists = GetAllHistos( file_xs ) ;
  MnvH1D* data = (MnvH1D*) hists[6]->Clone("test_data");

  weightedQ2 = (TH1D*) data->Clone("weightedQ2");
  weightedQ2->Reset();
  for( auto it : WeightedQ2 )
    weightedQ2->Fill( it.first, it.second );

  DrawXs( file_xs, file_xs_fk, file_xs_2D, iteration );

  return 0;

  GenerateAMatrix( data , "Fa_OriginalBinning", iteration );

  //vector<double> new_bins({0.05,0.2,0.3,0.4,0.6,0.8,1.0,1.2,2.0,4.0,10});
  vector<double> new_bins({0.0375,0.05,0.2,0.3,0.4,0.6,0.8,1.0,1.2,2.0,4.0,6.0});
  MnvH1D* data_rebinned = (MnvH1D*) data->Rebin(new_bins.size()-1, "data_rebinned", &new_bins[0] );
  GenerateAMatrix( data_rebinned , "Fa_Rebinned", iteration );



  //DrawFATests();

  DrawMig( file_mig );
  DrawEff( file_eff );
  //DrawEff( file_eff_fk, "fullKine" );
  return 0;
}
