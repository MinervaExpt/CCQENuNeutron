{
#include <vector>
#include <map>
  using namespace PlotUtils;
  vector<string>tags;
  tags.push_back("pre020");
  //tags.push_back("full");
  tags.push_back("0_40");
  tags.push_back("40_70");
  tags.push_back("70_100");

  for( int it = 0; it<tags.size(); it++)
  {
    string tag = tags[it];
    TFile* f = new TFile(Form("../root_files/inelastic_%s_CombinedPlaylists.root",tag.c_str()),"read");
    //TFile* f = new TFile("../root_files/inelastic_pre020_CombinedPlaylists.root","read");
    string type = "mc";

   // double rws[] = {0, .5,.6,.7,.8,.9,1,1.1,1.2,1.3,1.4,1.5,2,3};
    double rws[] = {0.1,0.3,.5,.6,.7,.8,.9,1,1.1,1.2,1.3,1.4,1.5,1.7,2};
    int n = sizeof(rws)/sizeof(rws[0]) ;
    //vector<double> reweights( rws, rws+ n );
    map<double, TH2D*> h_blobE, h_blobDist;
    map<double, TH2D*> h_blobE_data, h_blobDist_data;
    map<double, TH1D*> h_blobE1D, h_blobDist1D;
    TVector2 pot = (TVector2) f->Get("pot");
    int rebinQ2 = 3;
    double scale = pot->X()/pot->Y();
    for( int i = 0 ; i< n; i++ )
    {
      double w = rws[i];
      cout<<i<<", "<<w<<endl;

      h_blobE_data[w] = (TH2D*) f->Get( Form("h_blobE_q2qe_%0.1f_data", w ) );
      h_blobE_data[w]->GetXaxis()->SetTitle("Energy (MeV)");
      h_blobE_data[w]->GetYaxis()->SetTitle("Events/bin");
      h_blobE_data[w]->Rebin2D(1, rebinQ2);
      h_blobDist_data[w] = (TH2D*) f->Get( Form("h_blobDist_q2qe_%0.1f_data",  w ) );
      h_blobDist_data[w]->GetXaxis()->SetTitle("Dist (mm)");
      h_blobDist_data[w]->GetYaxis()->SetTitle("Events/bin");
      h_blobDist_data[w]->Rebin2D(5,rebinQ2);





      h_blobE[w] = (TH2D*) f->Get( Form("h_blobE_q2qe_%0.1f_%s", w, type.c_str() ) );
      h_blobE[w]->GetXaxis()->SetTitle("Energy (MeV)");
      h_blobE[w]->GetYaxis()->SetTitle("Events/bin");
      h_blobE[w]->Rebin2D(1, rebinQ2);
      h_blobDist[w] = (TH2D*) f->Get( Form("h_blobDist_q2qe_%0.1f_%s", w, type.c_str() ));
      h_blobDist[w]->GetXaxis()->SetTitle("Dist (mm)");
      h_blobDist[w]->GetYaxis()->SetTitle("Events/bin");
      h_blobDist[w]->Rebin2D(5,rebinQ2);

      scale = h_blobE_data[w]->Integral()/h_blobE[w]->Integral();
      h_blobE[w]->Scale(scale);
      h_blobDist[w]->Scale(scale);

      h_blobE1D[w] = (TH1D*) h_blobE[w]->ProjectionX();
      h_blobDist1D[w] = (TH1D*) h_blobDist[w]->ProjectionX();
      h_blobDist1D[w]->Rebin(2);
      h_blobE1D[w]->SetLineColor(i+1);
      h_blobDist1D[w]->SetLineColor(i+1);
      cout<<h_blobDist1D[w]<<endl;
    }



    TCanvas* c = new TCanvas("c","c");
    c->cd();

    string opt = "histl";
    TLegend *l = new TLegend( 0, 0, 1,1 );

    TLatex* text = new TLatex();
    text->SetTextColor(kRed);
    text->SetTextFont(82);
    text->SetTextSize(0.1);
    text->SetTextAlign(22);

    TLine * line  = new TLine(0,1,2500,1);
    line->SetLineColor(kBlack);
    line->SetLineWidth(1);
    line->SetLineStyle(2);

    int imin = 0;
    int NPlots = h_blobE[1]->GetNbinsY();
    int nX = int(TMath::Sqrt(NPlots))+1;
    cout<<"Dividing canvas: "<<nX<<", "<<nX-1<<endl;
    c->Divide(nX,nX-1,0,0);

    TH2D* h_blobE_base = (TH2D*) h_blobE[1.0]->Clone("h_blobE_base");
    for( int i = 0; i< n; i++ )
    {
      double w = rws[i];
      h_blobE[w]->Divide( h_blobE_base);
      for( int j = 0; j< NPlots; j++ )
      {
        TH1D* h = (TH1D*) h_blobE[w]->ProjectionX(Form("%.1f_%d",w,j+1), j+1,j+1);
        h->SetLineColor(i+1);
        h->GetYaxis()->SetRangeUser(0.5,1.5);
        c->cd(j+1);
        cout<<j+1<<endl;
        h->Draw(opt.c_str() );
      }

      l->AddEntry(h_blobE1D[w], Form("%.1f",w ) );
      if ( i==imin ) opt+="same";
    }

    h_blobE_data[1.0]->Divide( h_blobE_base );
    for( int j = 0; j< NPlots; j++ )
    {
      c->cd(j+1);
      if(type == "mc")
      {
        TH1D* hdata = (TH1D*) h_blobE_data[1.0]->ProjectionX(Form("data_%.1f_%d",w,j+1), j+1,j+1);
        hdata->SetMarkerStyle(0);
        hdata->Draw("esame");
        //hdata->Draw("esame");
      }

      double q2Low = h_blobE[1]->GetYaxis()->GetBinLowEdge(j+1);
      double q2High = h_blobE[1]->GetYaxis()->GetBinUpEdge(j+1);
      text->DrawLatexNDC(0.5,0.9, Form("(%.3f, %.3f)",q2Low,q2High));
      //gPad->SetGrid();
      line->Draw();
    }

    c->cd(NPlots+1);
    l->Draw();
    c->cd(NPlots+2);
    text->DrawLatexNDC(.5,.5,type.c_str() );
    //h_blobE[0]->Draw("colz");
    //c->SetGrid();
    c->Print(Form( "BlobE2D_%s_%s.pdf", type.c_str(), tag.c_str() ) );


    opt="histl";
    TH2D* h_blobDist_base = (TH2D*) h_blobDist[1.0]->Clone( "h_blobDist_base" );
    for( int i = imin; i< n; i++ )
    {
      double w = rws[i];
      h_blobDist[w]->Divide( h_blobDist_base );
      for( int j = 0; j< NPlots; j++ )
      {
        TH1D* h = (TH1D*) h_blobDist[w]->ProjectionX(Form("%.1f_%d",w,j+1), j+1,j+1);
        h->GetXaxis()->SetRangeUser(0,2500);
        h->SetLineColor(i+1);
        h->GetYaxis()->SetRangeUser(0.5,1.5);
        c->cd(j+1);
        h->Draw(opt.c_str() );
      }

      if ( i==imin ) opt+="same";
    }


    h_blobDist_data[1]->Divide( h_blobDist_base);
    for( int j = 0; j< NPlots; j++ )
    {
      c->cd(j+1);
      if(type == "mc")
      {
        TH1D* hdata = (TH1D*) h_blobDist_data[1]->ProjectionX(Form("data_%.1f_%d",w,j+1), j+1,j+1);
        hdata->SetMarkerStyle(0);
        hdata->Draw("esame");
        hdata->Draw("esame");
      }

      double q2Low = h_blobE[1]->GetYaxis()->GetBinLowEdge(j+1);
      double q2High = h_blobE[1]->GetYaxis()->GetBinUpEdge(j+1);
      text->DrawLatexNDC(0.5,0.9, Form("(%.3f, %.3f)",q2Low,q2High));
      //gPad->SetGrid();
      line->Draw();
    }

    c->cd(NPlots+1);
    l->Draw();
    c->cd(NPlots+2);
    text->DrawLatexNDC(.5,.5,type.c_str() );
    //h_blobDist[0]->Draw("colz");
    //c->SetGrid();
    c->Print(Form( "BlobDist2D_%s_%s.pdf", type.c_str(), tag.c_str()) );


  }





}
