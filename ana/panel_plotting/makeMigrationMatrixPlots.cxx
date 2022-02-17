
{
  using namespace PlotUtils;
  vector<string> tunes;
  tunes.push_back("CV2");
  tunes.push_back("CV2Elastic");
  tunes.push_back("CV2ElasticNuwroSF");
  tunes.push_back("CV2ElasticNuwroLFG");

  vector<string> names;
  names.push_back("mnvGENIEv2");
  names.push_back("Elastic Fix");
  names.push_back("Nuwro SF");
  names.push_back("Nuwro LFG");


    TCanvas *c = new TCanvas("c","c");
    c->cd();
    c->SetMargin(.10,.10,.10,.10);
    TFile *f = new TFile(Form("../rootfiles/minervameAntiNu/CV2/Migration_CombinedPlaylists.root"));
    MnvH2D* mc  = (MnvH2D*) f->Get("h_q2qe_angle_10_migration_qelike_qe_h");
    //TH2D migration("migration","migration",20,0,19,20,0,19);
    TH2D migration = mc->GetCVHistoWithStatError();
    migration.Reset();
    migration.GetXaxis()->SetTitle("reco");
    migration.GetYaxis()->SetTitle("truth");
    migration.GetZaxis()->SetRangeUser(0,1);
    for( int iy = 0; iy < mc->GetNbinsY(); iy++)
    {
      double sum = mc->Integral(0,-1,iy,iy);
      if(sum==0) sum=1;
      for( int ix = 0; ix < mc->GetNbinsX(); ix++ ){
        double value = mc->GetBinContent(ix,iy);
        value/=sum;
        migration.SetBinContent(ix,iy, value);
      }
    }
    migration.Draw("colz");
    c->SetLogx();
    c->SetLogy();

    gStyle->SetPalette(53);
    c->Print("plots/minervameAntiNu/Migration.pdf");


}
