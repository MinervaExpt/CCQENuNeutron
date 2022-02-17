{
#include <vector>
#include <map>
  using namespace PlotUtils;

  vector<string>tags;
  tags.push_back("0_20");
  //tags.push_back("full");
  tags.push_back("0_40");
  tags.push_back("40_70");
  tags.push_back("70_100");

  for( int it = 0; it<tags.size(); it++)
  {
    string tag = tags[it];
    //TFile* f = new TFile(Form("../root_files/inelastic_graph_%s_CombinedPlaylists.root",tag.c_str()),"read");
    TFile* f = new TFile(Form("../root_files/inelastic_graph_%s_minervame5A.root",tag.c_str()),"read");
    string type = "mc";

    TVector2 pot = (TVector2) f->Get("pot");
    double scale = pot->X()/pot->Y();
    //scale = 1;

    double rws[] = {0.1,0.3,.5,.6,.7,.8,.9,1,1.1,1.2,1.3,1.4,1.5,1.7,2};
    int n = sizeof(rws)/sizeof(rws[0]) ;
    //vector<double> reweights( rws, rws+ n );
    map<double, TH2D*> h_blobE, h_blobDist;
    map<double, TH2D*> h_blobE_data, h_blobDist_data;
    map<double, TH1D*> h_blobE1D, h_blobDist1D;
    map<double, TH1D*> h_blobE1D_data, h_blobDist1D_data;
    double nmc = 0;
    for( int i = 0; i< n; i++ )
    {
      cout<<i<<endl;
      double w = rws[i];
      h_blobE[w] = (TH2D*) f->Get( Form("h_blobE_q2qe_%0.1f_%s", type.c_str(), w ) );
      h_blobDist[w] = (TH2D*) f->Get( Form("h_blobDist_q2qe_%0.1f_%s", type.c_str(), w ) );
      if(w==1) nmc = h_blobE[w]->Integral();

      h_blobE_data[w] = (TH2D*) f->Get( Form("h_blobE_q2qe_%0.1f_data",  w ) );
      h_blobDist_data[w] = (TH2D*) f->Get( Form("h_blobDist_q2qe_%0.1f_data",  w ) );
      scale = h_blobE_data[w]->Integral()/h_blobE[w]->Integral();

      h_blobE[w]->Scale(scale);
      h_blobDist[w]->Scale(scale);

      h_blobE1D_data[w] = (TH1D*) h_blobE_data[w]->ProjectionX();
      h_blobDist1D_data[w] = (TH1D*) h_blobDist_data[w]->ProjectionX();

      h_blobDist1D_data[w]->Rebin(2);

      h_blobE1D[w] = (TH1D*) h_blobE[w]->ProjectionX();
      h_blobDist1D[w] = (TH1D*) h_blobDist[w]->ProjectionX();

      h_blobDist1D[w]->Rebin(2);

      h_blobE1D[w]->SetLineColor(i+1);
      h_blobDist1D[w]->SetLineColor(i+1);
      cout<<h_blobDist1D[w]<<endl;
    }

    TCanvas* c = new TCanvas("c","c");
    c->cd();

    string opt = "hist";
    TLegend *l = new TLegend( .7, .6, 1,1 );
    int imin = 0;
    TH1D* blobE = (TH1D*) h_blobE1D[1]->Clone("blobE");
    for( int i = 0; i< n; i++ )
    {
      double w = rws[i];
      //calc chi2
      double chi2 = 0;
      for( int j = 1; j <= h_blobE1D[w]->GetNbinsX(); j++ )
      {
        double vmc = h_blobE1D[w]->GetBinContent(j);
        double vdata = h_blobE1D_data[1]->GetBinContent(j);
        chi2+= ( vmc - vdata )*( vmc - vdata )/(  vdata ) ;
      }
      h_blobE1D[w]->Divide( blobE  );
      h_blobE1D[w]->GetYaxis()->SetRangeUser(.5,1.5);
      h_blobE1D[w]->GetXaxis()->SetRangeUser(0,100);
      h_blobE1D[w]->GetXaxis()->SetTitle("Blob E (MeV)");
      h_blobE1D[w]->Draw(opt.c_str() );

      l->AddEntry(h_blobE1D[w], Form("%.1f, #chi^{2}=%.1f",w,chi2 ) );
      if ( i==imin ) opt+="same";
    }
    if(type=="mc")
    {
      h_blobE1D_data[1]->Divide( blobE  );
      h_blobE1D_data[1]->Draw("esame");
    }
    l->Draw();
    //h_blobE[0]->Draw("colz");
    c->SetGrid();
    c->Print(Form("BlobE_%s_%s.pdf",type.c_str(), tag.c_str() ) );

    opt = "hist";
    int imin = 0;
    TH1D* blobDist = (TH1D*) h_blobDist1D[1]->Clone("blobDist");
    TLegend *l2 = new TLegend( .7, .6, 1,1 );
    for( int i = imin; i< n; i++ )
    {
      double w = rws[i];
      double chi2=0;
      for( int j = 1; j <= h_blobDist1D[w]->GetNbinsX(); j++ )
      {
        double vmc = h_blobDist1D[w]->GetBinContent(j);
        double vdata = h_blobDist1D_data[1]->GetBinContent(j);
        if( h_blobDist1D_data[1]->FindBin( 1500 ) == j ) break;
        chi2+= ( vmc - vdata )*( vmc - vdata )/(  vdata ) ;
        //cout<<vmc<<", "<<vdata<<endl;
      }

      h_blobDist1D[w]->Divide( blobDist );
      h_blobDist1D[w]->GetYaxis()->SetRangeUser(0.5,1.5);
      h_blobDist1D[w]->GetXaxis()->SetRangeUser(0,2500);
      h_blobDist1D[w]->GetXaxis()->SetTitle("Blob Dist (mm)");
      h_blobDist1D[w]->Draw(opt.c_str() );
      if ( i==imin ) opt+="same";
      l2->AddEntry(h_blobE1D[w], Form("%.1f, #chi^{2}=%.2f",w,chi2 ) );
    }
    if(type=="mc")
    {
      h_blobDist1D_data[1]->Divide( blobDist );
      h_blobDist1D_data[1]->Draw("esame");
    }
    l2->Draw();
    //h_blobDist[0]->Draw("colz");
    c->SetGrid();
    c->Print(Form("BlobDist_%s_%s.pdf",type.c_str(), tag.c_str() ));


  }

}
