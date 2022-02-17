#include "TFile.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"

using namespace PlotUtils;

vector<double> optimizeBinning( vector<double> &binEdges, vector<double> &statErrs )
{

  vector<double> newbinEdges;
  
  return newbinEdges;
}

void OptimizeBinning()
{
  TFile *f = new TFile("/minerva/app/users/tejinc/cmtuser/Minerva_v21r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/rootfiles/minervame1D/bkg_sub_hists_ptmu.root");

  MnvH2D* dthetaR_ptmu = (MnvH2D*) f->Get("h_dthetaR_ptmu_mc_2track_weighted_nobck_signal");
  MnvH2D* dthetaP_ptmu = (MnvH2D*) f->Get("h_dthetaP_ptmu_mc_2track_weighted_nobck_signal");

  TH1D


  


}
