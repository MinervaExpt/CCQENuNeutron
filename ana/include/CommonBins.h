#include "include/CCQENuUtils.h"
#include "PlotUtils/HyperDimLinearizer.h"

#include <vector>

using namespace CCQENU_ANA;
AnaBinning *binner = new AnaBinning();
CCQENuBinning *minmodbinner = new CCQENuBinning();
NeutronBlobBinning *neutbinner = new NeutronBlobBinning();

axis_binning GetBins( int nBins, double size )
{
  axis_binning bin;
  vector<double> bins;
  for( int i =0; i<nBins; i++ ) bins.push_back(i*size);
  bin.bin_edges = bins;
  bin.nbins = bins.size()-1;
  bin.min = bins.front();
  bin.max = bins.back();
  return bin;
};

axis_binning GetBins( double min, double max, double nBins )
{
  axis_binning bin;
  vector<double> bins;
  double size = (max-min)/nBins;
  for( int i =0; i<nBins+1; i++ ) bins.push_back(i*size+min);
  bin.bin_edges = bins;
  bin.nbins = nBins;
  bin.min = min;
  bin.max = max;
  return bin;
};


axis_binning nubins           = minmodbinner->GetNuBinsGeV(); 
axis_binning thetabins      = minmodbinner->GetMuonThetaBins();
axis_binning costhetabins   = minmodbinner->GetMuonCosThetaBinsMiniBoone();
axis_binning muonTbins      = binner->GetMuonEnergyBinsGeV();
axis_binning muonPbins      = binner->GetMuonEnergyBinsGeV();
//axis_binning muonPbins      = binner->GetMuonEnergyUniformBinsGeV();
//axis_binning muonPtbins      = minmodbinner->GetMuonPtUniformBinsGeV();
axis_binning muonPtbins     = minmodbinner->GetMuonPtBinsGeV();
axis_binning muonPzbins     = minmodbinner->GetMuonPzBinsGeV();
axis_binning Q2bins         = minmodbinner->GetQ2BinsGeV();
axis_binning Enubins        = minmodbinner->GetMuonPzBinsGeV();
axis_binning blobEnergyBinsMeV= minmodbinner->GetBlobEnergyBinsMeV();// up to 500 MeV: Events / 10 MeV
axis_binning dalphatbins    = minmodbinner->GetDalphaT();
axis_binning dphitbins    = minmodbinner->GetDphiT();
axis_binning dphitbins_sign = GetBins(-180,180,72);
axis_binning pnbins    = minmodbinner->GetPn();
axis_binning dptbins    = minmodbinner->GetDeltaPt();
axis_binning dptxbins    = minmodbinner->GetDeltaPtx();
axis_binning dptybins    = minmodbinner->GetDeltaPty();
axis_binning protonThetabins = minmodbinner->GetProtonThetaBins();
axis_binning protonKineticbins = minmodbinner->GetProtonEnergyBinsGeV();
axis_binning signbins          = minmodbinner->GetSignedBins();
axis_binning signeddalphatbins = minmodbinner->GetSignedDeltaAlphaT();
axis_binning signeddphitbins = minmodbinner->GetSignedDeltaPhiT();


axis_binning dthetaPerpbins = neutbinner->GetThetaPerpBinsDegreeOptimized();
axis_binning dthetaPerpbinsPositive = neutbinner->GetThetaPerpBinsDegreeOptimizedPositive();

axis_binning dthetaReactbins = neutbinner->GetThetaReactBinsDegreeOptimized();
axis_binning dthetaReactbinsPositive = neutbinner->GetThetaReactBinsDegreeOptimizedPositive();

axis_binning dthetaFineBins = neutbinner->GetThetaBinsDegreeFine();
axis_binning dthetaFineBinsUniform = GetBins(-180,180,180);

axis_binning dalphatbins_TKI    = minmodbinner->GetDalphaT_TKI();
axis_binning dphitbins_TKI    = minmodbinner->GetDphiT_TKI();
axis_binning pnbins_TKI    = minmodbinner->GetPn_TKI();
axis_binning dptbins_TKI    = minmodbinner->GetDeltaPt_TKI();
axis_binning dptxbins_TKI    = minmodbinner->GetDeltaPtx_TKI();
axis_binning dptybins_TKI    = minmodbinner->GetDeltaPty_TKI();
axis_binning protonThetabins_TKI = minmodbinner->GetProtonThetaBins_TKI();
axis_binning protonKineticbins_TKI = minmodbinner->GetProtonEnergyBinsGeV_TKI();
axis_binning signeddalphatbins_TKI = minmodbinner->GetSignedDeltaAlphaT_TKI(); 
axis_binning signeddphitbins_TKI = minmodbinner->GetSignedDeltaPhiT_TKI();

axis_binning ptbins3D = minmodbinner->Get3DPtBins();
axis_binning pzbins3D = minmodbinner->Get3DPzBins();
axis_binning recoilbins3D = minmodbinner->Get3DRecoilBins();



axis_binning leptonmomentumbins = neutbinner->GetLeptonPBinsGeV();

//custom bins : nblob/ntracks bins
axis_binning nDiscreteBins = GetBins(20,1);
axis_binning Uniform100BinsMeV = GetBins(100,10.);
axis_binning Uniform100BinsGeV = GetBins(100,10./1000.);
axis_binning Uniform50BinsMeV = GetBins(100,20.);
axis_binning Uniform50BinsGeV = GetBins(100,20./1000.);
axis_binning Uniform25BinsMeV = GetBins(100,40.);
axis_binning Uniform25BinsGeV = GetBins(100,40./1000.);

axis_binning UniformDThetaBins = GetBins(-100,100,50);
axis_binning UniformMuonPBins = GetBins(0,20,20);
axis_binning UniformVisibleEBins = GetBins(0,500,50); //MeV
axis_binning UniformVtxEBins = GetBins(0,100,50); //MeV

axis_binning UniformQ0Q3Bins = GetBins(0,0.8,8);//GeV

axis_binning BlobDistBins = GetBins(0,5000,50 );
axis_binning BlobEnergyBins = GetBins(0, 200, 20 );
axis_binning BlobEnergyBinsLow = GetBins(0, 10, 10 );
axis_binning BlobEnergyBinsFine = GetBins(0, 200, 200 );
axis_binning BlobNumBins = GetBins( 0, 10, 10 );
axis_binning MuonThetaBins = GetBins(0,20,10);
axis_binning MuonPhiBins = GetBins(-180,180, 18 );
axis_binning MuonEnergyBins = GetBins(0, 20, 10 );

axis_binning VtxEnergyBins = GetBins(0,50,50);
axis_binning Q2binsLow = GetBins(0,0.03,10);

axis_binning nHitsBins = GetBins(0,100,100);
axis_binning hitsEBins = GetBins(0,1000,100);
axis_binning nClusBins = GetBins(0,50,50);
axis_binning PartTypeBins = GetBins(0,10,10); //0: photon, 1: electron, 2:pi0, 3:charged pi, 4:proton, 5: neutron, 6 muo

axis_binning kineticBins = GetBins(0,2,50);
axis_binning kineticResBins = GetBins(-1,1,50);
axis_binning EnuResBins = GetBins(-1,1,50);
