#include<vector>
  

using namespace PlotUtils;
void makePlot()
{
  string file_loc = "/minerva/app/users/tejinc/cmtuser/Minerva_v21r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/rootfiles/minervameAntiNu/cuts_muonPCut-0_maxClusECut-0_vtxCut-0_blobMaxECut-0_n2DCut-0_nBlobCut-1_newBlobCut-1_ElasticQE/CV2ElasticNuwroSF/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal-Angles_CombinedPlaylists.root";
  string file_loc2 = "/minerva/app/users/tejinc/cmtuser/Minerva_v21r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/rootfiles/minervameAntiNu/cuts_muonPCut-0_maxClusECut-0_vtxCut-0_blobMaxECut-0_n2DCut-0_nBlobCut-1_newBlobCut-1_ElasticQE/CV2ElasticNuwroSF/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal_CombinedPlaylists.root";
  //string file_loc2 = "/minerva/app/users/tejinc/cmtuser/Minerva_v21r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/make_hists/test_muon_ela_SF_neutronW.root";
  TFile* f = new TFile(file_loc.c_str(), "READ");
  TFile* f2 = new TFile(file_loc2.c_str(), "READ");
  //MnvH1D* h = (MnvH1D*) f->Get("h_q2qe_region_00_qelike_qe_oth");
  MnvH1D* h = (MnvH1D*) f->Get("h_q2qe_qelike_qe_oth");
  MnvH2D* h2D = (MnvH2D*) f2->Get("h_dthetaPdthetaR_q2qe_qelike_qe_oth");
  MnvH1D* h2 = (MnvH1D*) h2D->ProjectionX();

  vector<string> vNames = h->GetVertErrorBandNames();
  for( vector<string>::iterator it = vNames.begin(); it!=vNames.end();it++) cout<<*it<<endl;
  cout<<"================="<<endl;
  vector<string> lNames = h->GetLatErrorBandNames();
  for( vector<string>::iterator it = lNames.begin(); it!=lNames.end();it++) cout<<*it<<endl;

  //only plot lat
  
  for( vector<string>::iterator it = lNames.begin(); it!=lNames.end();it++)
  {
    MnvLatErrorBand *eb = h2->GetLatErrorBand( *it );

    TCanvas c1;
    c1.cd();
    eb->SetLineWidth(3);
    eb->Draw("hist");
    TH1D* hlow = eb->GetHist(0);
    hlow->Draw("histlsame");
    TH1D* hhigh = eb->GetHist(1);
    hhigh->Draw("histsame");
    c1.SetLogx();
    c1.Print(Form("lat_wrong_%s.pdf",it->c_str() ) );


  }



}

















