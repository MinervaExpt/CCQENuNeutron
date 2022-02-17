//#include "myPlotStyle.h"
#include "TParameter.h"

#include "include/GeneralIncludes.h"


using namespace PlotUtils;
//forward declaration
//gROOT->SetBatch();

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



void DrawEnergy( TFile *f_orig )
{
  cout<<"DrawEnergy"<<endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  CCQENuPlotUtils *putils = new CCQENuPlotUtils( fluxHistoExists );
  CCQENuUtils *utils = new CCQENuUtils( false, fluxHistoExists );


  TFile *f1 = f_orig;
  //TFile *f2 = f_signal_weighted;

  vector<MnvH2D*> h_enuTrue_q2qe, h_enuFit_q2qe, h_enuCCQE_q2qe;
  vector<MnvH2D*> h_denuTrue_q2qe, h_denuFit_q2qe, h_denuCCQE_q2qe;
  vector<MnvH2D*> h_denuFracTrue_q2qe, h_denuFracFit_q2qe, h_denuFracCCQE_q2qe;

  utils->bookHistos(f1, &h_enuTrue_q2qe[0], "h_enuTrue_q2qe"); utils->bookHistos(f1, &h_enuFit_q2qe[0], "h_enuFit_q2qe"); utils->bookHistos(f1, &h_enuCCQE_q2qe[0], "h_enuCCQE_q2qe");
  utils->bookHistos(f1, &h_denuTrue_q2qe[0], "h_denuTrue_q2qe"); utils->bookHistos(f1, &h_denuFit_q2qe[0], "h_denuFit_q2qe"); utils->bookHistos(f1, &h_denuCCQE_q2qe[0], "h_denuCCQE_q2qe");
  utils->bookHistos(f1, &h_denuFracTrue_q2qe[0], "h_denuFracTrue_q2qe"); utils->bookHistos(f1, &h_denuFracFit_q2qe[0], "h_denuFracFit_q2qe"); utils->bookHistos(f1, &h_denuFracCCQE_q2qe[0], "h_denuFracCCQE_q2qe");

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
