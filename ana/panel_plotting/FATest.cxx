//#include "myPlotStyle.h"
#include "TParameter.h"
#include "include/GeneralIncludes.h"

#include "PlotUtils/FluxReweighter.h"

#include "TMinuit.h"
#include "TFitter.h"

#include "FA.h"
#include "FA05.h"

#include "TMultiGraph.h"


using namespace PlotUtils;








double MA=0.99;


void PrintValuesBBA05( vector<double>& Q2List )
{
  cout<<"Printing Values of BBA05"<<endl;
  for( auto Q2 : Q2List )
  {

    double FA = BBA05::FA( Q2, MA );
    double GVE= BBA05::GEp(Q2) - BBA05::GEn(Q2);
    double GVM= BBA05::GMp(Q2) - BBA05::GMn(Q2);
    //cout<<"{ "<<Q2<<", {FA->"<<FA<<", GVE->"<<GVE<<", GVM->"<<GVM<<"} },"<<endl;
    cout<<"{ FA->"<<FA<<", GVE->"<<GVE<<", GVM->"<<GVM<<"},"<<endl;
  }
}

void PrintValuesBBBA07( vector<double>& Q2List )
{
  cout<<"Printing Values of BBBA07"<<endl;
  for( auto Q2 : Q2List )
  {

    double FA = BBBA07::FA( Q2, MA );
    double GVE= BBBA07::GEp(Q2) - BBBA07::GEn(Q2);
    double GVM= BBBA07::GMp(Q2) - BBBA07::GMn(Q2);
    //cout<<"{ "<<Q2<<", {FA->"<<FA<<", GVE->"<<GVE<<", GVM->"<<GVM<<"} },"<<endl;
    cout<<"{ FA->"<<FA<<", GVE->"<<GVE<<", GVM->"<<GVM<<"},"<<endl;
  }
}







int main(int argc, char* argv[])
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
 

  string f_flux_suffix = (neutrinoMode)? 
    "data/flux/flux-gen2thin-pdg14-minervame1D1M1NWeightedAve.root": 
    "data/flux/flux-gen2thin-pdg-14-minervame6A.root";
  string f_flux = Form("%s/%s", std::getenv("PLOTUTILSROOT"),f_flux_suffix.c_str() ) ; 

  TFile f_genie("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/GENIEXSecExtract/cmt/GENIEXSECEXTRACT_CCQENuNeutron_Q2QE_minervame6D.root");
  MnvH1D *h_ds_dq2_genie = (MnvH1D*) f_genie.Get("ds_dq2_xsec");

  vector<double> Q2List({0.003125,0.009375,0.01875,0.03125,0.04375,0.075,0.125,0.175,0.25,0.35,.5,.7,.9,1.1,1.6,3,5,8});


  PrintValuesBBA05( Q2List );
  BBBA07::DeclareFittingPars();
  PrintValuesBBBA07( Q2List );

  return 0;
}
