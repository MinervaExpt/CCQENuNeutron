#include "include/CCQENuPlotUtils.h"
using namespace std;
int main()
{
  PlotUtils::MnvPlotter *plotter = new PlotUtils::MnvPlotter();
  //gStyle->SetErrorX(0.0001);
  //TFile* f_eff_new = new TFile( "/minerva/data/users/tejinc/hists/v10r8p12/effPurity/_Cut1_2p2hVariation0_plastic_minos_match_genie_var_atleast_one_proton_450MeV_ana-michel_LE-nu_minervaAllQERW_v10r8p12.root","read");
  TFile* f_eff_new = new TFile( "/minerva/data/users/tejinc/hists/v10r8p12/effPurity/transverseEffPlots_WithCuttransverseEb0_20180926_minerva17913C_Cut1_2p2hVariation0_plastic_tuned-bkgd_minos_match_genie_var_atleast_one_proton_450MeV_ana-michel_LE-nu_minerva17913C_v10r8p12.root");
  
  std::vector<TString> Variables, AxisTitles; 
  Variables.push_back("neutronmomentum");
  Variables.push_back("dpTx");
  Variables.push_back("dpTy");
  Variables.push_back("dpt");
  Variables.push_back("dalphat");
  Variables.push_back("dphit");
  Variables.push_back("muonmomentum");
  Variables.push_back("muontheta");
  Variables.push_back("protonmomentum");
  Variables.push_back("protontheta");

  AxisTitles.push_back("#it{p}_{n} (GeV/#it{c})");
  AxisTitles.push_back("#delta#it{p}_{Tx} (GeV/#it{c})");
  AxisTitles.push_back("#delta#it{p}_{Ty} (GeV/#it{c})");
  AxisTitles.push_back("#delta#it{p}_{T} (GeV/#it{c})");
  AxisTitles.push_back("#delta#it{#alpha}_{T} (degree)");
  AxisTitles.push_back("#delta#it{#phi}_{T} (degree)");
  AxisTitles.push_back("#it{p}_{#mu} (GeV/#it{c})");
  AxisTitles.push_back("#it{#theta}_{#mu} (degree)");
  AxisTitles.push_back("#it{p}_{p} (GeV/#it{c})");
  AxisTitles.push_back("#it{#theta}_{p} (degree)");



  TCanvas c("c","c", 800,800);
 
  for( int i = 0; i< Variables.size(); i++ )
  {
    const char *var = Variables[i].Data(), *xtitle=AxisTitles[i].Data();
    cout<<Form("h_%s_eff_efficiency_plastic_num", var )<<endl;
    PlotUtils::MnvH1D *h0_num = (PlotUtils::MnvH1D*) f_eff_new->Get( Form("h_%s_eff_efficiency_plastic_num", var ) );
    PlotUtils::MnvH1D *h0_den = (PlotUtils::MnvH1D*) f_eff_new->Get( Form("h_%s_eff_efficiency_plastic_den", var ) );
    //PlotUtils::MnvH1D *h0_eff = (PlotUtils::MnvH1D*) h0_num->Clone( Form("eff_%s_new", var ) );
    //h0_eff->Divide( h0_eff, h0_den );

    //TH1D h0_eff_all_errors = (TH1D) h0_eff->GetCVHistoWithError();
    //TH1D h0_eff_sys_errors = (TH1D) h0_eff->GetCVHistoWithError(false);

    //h0_eff->GetXaxis()->SetTitle(xtitle);
    //h0_eff->GetYaxis()->SetTitle("Efficiency");

    ////plotter.DrawMCWithErrorBand(h0_eff);
    //plotter->DrawMCWithErrorBand( &h0_eff_sys_errors);
    //c.Print( Form("be_test/eff_%s.eps", var) );
   
    

  }
  c.Close();

  cout<<"end"<<endl;
  

  return 0;
}

